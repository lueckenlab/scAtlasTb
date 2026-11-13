import re
import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
from pprint import pformat

from utils.io import read_anndata, write_zarr_linked


def check_same_categories(df, columns):
    first_categories = df[columns[0]].cat.categories
    if not all(df[col].cat.categories.equals(first_categories) for col in columns):
        logging.info(
            "Columns %s have different category sets; categories will be aligned before majority voting.",
            columns,
        )


def check_category_overlap(df, columns):
    from functools import reduce

    shared_categories = reduce(
        pd.Index.intersection,
        [pd.Index(df[col].dropna().unique()) for col in columns],
    )
    if len(shared_categories) == 0:
        categories_per_col = {
            col: pd.Index(df[col].dropna().unique()).tolist()
            for col in columns
        }
        raise ValueError(
            "No overlapping observed labels found across majority-vote columns. "
            f"Columns: {columns}. Observed labels per column: {pformat(categories_per_col)}"
        )


def mode_counts_per_row(arr, dtype=np.uint8):
    from tqdm import tqdm

    n_rows, n_cols = arr.shape
    
    if n_cols > 255:
        dtype = np.int16

    modes = np.empty(n_rows, dtype=arr.dtype)
    mode_counts = np.empty(n_rows, dtype=dtype)

    miniters = max(1, int(n_rows) / 10)
    for i in tqdm(range(n_rows), desc='Compute majority vote per row', miniters=miniters):
        # Filter out NaN codes (-1) to avoid selecting them as mode
        valid = arr[i][arr[i] != -1]
        
        if len(valid) == 0:
            # All values are NaN, set mode to -1
            modes[i] = -1
            mode_counts[i] = 0
        else:
            # Use bincount for efficiency: O(n_cols) instead of O(n_cols * log(n_cols))
            # This works because categorical codes are non-negative integers
            counts = np.bincount(valid.astype(int))
            max_idx = np.argmax(counts)
            modes[i] = max_idx
            mode_counts[i] = counts[max_idx]

    return modes, mode_counts


def get_majority_consensus(df, columns, new_key='majority_consensus', threshold=0.5):
    from functools import reduce

    agreement_col = f'{new_key}_agreement'
    low_agreement_col = f'{new_key}_low_agreement'
    n_vote_columns = len(columns)

    # Convert columns to categorical if not yet the case
    for col in columns:
        if not pd.api.types.is_categorical_dtype(df[col]):
            df[col] = df[col].astype('category')
    
    # Warn if categories aren't the same across columns
    check_same_categories(df, columns)
    check_category_overlap(df, columns)

    # Create a unified categorical dtype
    categories = reduce(pd.Index.union, [df[col].cat.categories for col in columns])
    cat_dtype = pd.CategoricalDtype(categories=categories)

    # Rebuild category codes from label values to avoid code/category order mismatches
    # when source categorical dtypes differ across columns.
    df_cat = pd.DataFrame(
        {
            col: pd.Categorical(df[col].astype('object'), dtype=cat_dtype).codes
            for col in columns
        },
        index=df.index,
    )
    
    majority_votes, mode_count = mode_counts_per_row(df_cat.to_numpy())
    
    logging.info('Assign majority votes to new column...')
    df[new_key] = pd.Categorical.from_codes(majority_votes, categories=categories)
    df[agreement_col] = mode_count / n_vote_columns
    
    min_agreement = 1 / n_vote_columns
    df.loc[df[agreement_col] <= min_agreement, new_key] = float('nan')
    # set to 0, because no agreement when all assignments different
    df.loc[df[agreement_col] <= min_agreement, agreement_col] = 0
    df[low_agreement_col] = df[agreement_col] <= threshold
    
    return df[[new_key, agreement_col, low_agreement_col]]


input_file = snakemake.input.zarr
output_file = snakemake.output.zarr
output_plots = Path(snakemake.output.plots)
output_plots.mkdir(parents=True, exist_ok=True)
column_patterns = snakemake.params.get('columns', [])
threshold = snakemake.params.get('threshold', 0.5)
new_key = 'majority_consensus'

logging.info('Read adata...')
adata = read_anndata(input_file, obs='obs')

logging.info(f'Parse column pattern: {column_patterns}')
columns = [
    col for col in adata.obs.columns for pattern in column_patterns
    if re.fullmatch(pattern, col)
]
columns = list(dict.fromkeys(columns))
assert len(columns) >= 2, (
    f'Insufficient columns found matching patterns: {column_patterns}. '
    f'Found {len(columns)} column(s), need at least 2.'
)
logging.info(f'parsed columns: {columns}')

logging.info('Majority voting...')
maj_df = get_majority_consensus(adata.obs, columns=columns, new_key=new_key, threshold=threshold)

logging.info('Compute consensus agreement stats...')
counts = maj_df.value_counts(
    subset=[new_key, f'{new_key}_agreement'],
    sort=False,
    dropna=False
).reset_index(name='count')
counts[f'total_cells_{new_key}_label'] = counts.groupby(new_key, observed=True)['count'].transform('sum')
counts[f'frac_cells_in_{new_key}_label'] = counts['count'] / counts[f'total_cells_{new_key}_label']
print(counts, flush=True)
counts.to_csv(output_plots / 'majority_consensus_agreement.tsv', sep='\t', index=False)

logging.info('Plot majority consensus stats...')
plot_df = pd.crosstab(
    maj_df[new_key],
    maj_df[f'{new_key}_low_agreement'],
    dropna=False,
    normalize='index'
)

if True not in plot_df.columns:
    plot_df[True] = 0
elif False not in plot_df.columns:
    plot_df[False] = 0

ax = plot_df.sort_values(True).plot(
    kind='barh',
    stacked=True,
    color=sns.color_palette("colorblind").as_hex(),
    figsize=(5, 5 + 0.1 * plot_df.shape[0]),
)
ax.legend(title='Low agreement', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title(f'Fraction of cells with low agreement in {new_key} labels')
plt.grid(False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.savefig(output_plots / 'majority_consensus_frac.png', bbox_inches='tight')

# add to adata
adata.obs[maj_df.columns] = maj_df
logging.info(f'Columns after adding majority consensus: {adata.obs.columns}')
assert maj_df.columns.isin(adata.obs.columns).all(), 'Not all majority consensus columns were added to adata.obs'

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata=adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs']
)