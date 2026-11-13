import pandas as pd
import numpy as np
import scib
import yaml
from tqdm import tqdm
from joblib import Parallel, delayed

import logging
logging.basicConfig(level=logging.INFO)
import warnings

# Workaround: suppress deprecation warning raised by some multiprocessing backends
# calling a now-deprecated `_destroy` helper (appears on Python 3.12+).
warnings.filterwarnings(
    "ignore",
    message=r"Call to deprecated function \(or staticmethod\) _destroy",
    category=DeprecationWarning,
)

from utils.io import read_anndata


input_file = snakemake.input[0]
setup_file = snakemake.input.setup
output_file = snakemake.output.tsv
covariate = snakemake.wildcards.covariate
sample_key = 'group'
n_threads = max(snakemake.threads, 1)

logging.info(f'Read setup file "{setup_file}"...')
with open(setup_file, 'r') as f:
    setup = yaml.safe_load(f)
n_permute = setup['n_permute']
is_numeric = setup['is_numeric']
na_strings = setup.get('na_strings', [])


logging.info('Read anndata file...')
adata = read_anndata(
    input_file,
    X='obsm/X_pca',
    obs='obs',
    uns='uns',
)

logging.info('Prepare PCA representation...')
if 'pca' not in adata.uns:
    raise KeyError('.uns["pca"] is missing, please make sure PCA is computed and PC loadings are saved in adata.uns as provided by scanpy')
pca_info = adata.uns['pca']
if 'variance' not in pca_info:
    raise KeyError('.uns["pca"]["variance"] is missing, please make sure PCA variance is available in AnnData uns metadata')

logging.info('Filter adata for non-NA covariate...')
adata = adata[
    ~adata.obs[covariate].isin(na_strings) &
    adata.obs[covariate].notna()
].copy()

pca_var = np.asarray(pca_info['variance'], dtype=np.float32)
if pca_var.ndim != 1:
    raise ValueError(f'Expected .uns["pca"]["variance"] to be 1D, got shape {pca_var.shape}')

n_pcs = adata.X.shape[1]
if pca_var.shape[0] != n_pcs:
    raise ValueError(
        f'PCA variance length mismatch: .uns["pca"]["variance"] has length {pca_var.shape[0]}, '
        f'but PCA matrix has {n_pcs} components'
    )

if sample_key not in adata.obs.columns:
    raise KeyError(f'Missing required "{sample_key}" column in obs.')
obs = adata.obs[list({sample_key, covariate})].copy()
n_covariate = obs[covariate].nunique()
logging.info(f'Using sample_key="{sample_key}" with {obs[sample_key].nunique()} unique samples')

# make sure the PCA embedding is an array
if not isinstance(adata.X, np.ndarray):
    adata.X = adata.X.toarray()
X_pca = adata.X.astype(np.float32)
del adata

if is_numeric:
    logging.info(f'Covariate "{covariate}" is numeric with {n_covariate} unique values')
    obs[covariate] = obs[covariate].astype(np.float32)
else:
    logging.info(f'Covariate "{covariate}" is categorical with {n_covariate} unique values')
    obs[covariate] = obs[covariate].astype('category')
    value_counts = obs[[sample_key, covariate]].drop_duplicates().value_counts(covariate)
    logging.info(f'Value counts for covariate "{covariate}":\n{value_counts[value_counts > 1]}')

    if value_counts.max() == 1:
        logging.info('Sample key is the same as covariate, skipping permutation...')
        n_permute = 0

observed_covariate = obs[covariate].to_numpy(copy=False)

if not is_numeric:
    sample_codes, uniques = pd.factorize(obs[sample_key], sort=False)
    sample_codes = sample_codes.astype(np.int32)
    cov_per_sample = (
        obs.groupby(sample_key, observed=True)[covariate]
        .first()
        .reindex(uniques)
    )
    logging.info(f'Aggregated covariate per sample:\n{cov_per_sample}')


def run_pcr(i, seed, **kwargs):
    if i == 0:
        covariate_array = observed_covariate
    else:
        rng = np.random.default_rng([i, seed])

        if is_numeric:
            # Numeric covariates are permuted across all cells, not at sample level.
            covariate_array = rng.permutation(observed_covariate)
        else:
            sample_map = cov_per_sample.to_numpy(copy=False)
            # Ensure permutation is different from original order
            permuted_cov = rng.permutation(sample_map)
            max_tries = 1000
            for _ in range(max_tries):
                if not np.array_equal(permuted_cov, sample_map):
                    break
                permuted_cov = rng.permutation(sample_map)
            if np.array_equal(permuted_cov, sample_map):
                return f"{covariate}-{i}", True, np.nan
            covariate_array = permuted_cov[sample_codes]

    return f"{covariate}-{i}", i > 0, scib.me.pc_regression(
        X_pca,
        pca_var=pca_var,
        covariate=covariate_array,
        **kwargs,
    )


logging.info(f'Principal component regression for {covariate=}')

total_jobs = n_permute + 1
pcr_scores = []

with tqdm(
    total=total_jobs,
    desc=f'Computing {n_permute} permutations with {n_threads=}',
    unit='perms',
    mininterval=0.5,
) as pbar:
    # chunks for job submission
    chunk_size = max(1, min(10 * n_threads, total_jobs))
    
    with Parallel(n_jobs=n_threads, backend='loky') as parallel:
        for i in range(0, total_jobs, chunk_size):
            chunk_jobs = [
                delayed(run_pcr)(
                    j,
                    n_threads=1,
                    seed=42,
                    verbose=False,
                    linreg_method='numpy'
                )
                for j in range(i, min(i + chunk_size, total_jobs))
            ]
            chunk_results = parallel(chunk_jobs)
            pcr_scores.extend(chunk_results)
            pbar.update(len(chunk_results))

# Set permuted score when covariate is the same as the group variable
if n_permute == 0:
    # permutations wouldn't change values in this case, impute same value as covariate score
    perm_score = (f'{covariate}-1', True, pcr_scores[0][2])
    pcr_scores.append(perm_score)

df = pd.DataFrame.from_records(
    pcr_scores,
    columns=['covariate', 'permuted', 'pcr'],
)
df = df.dropna(subset=['pcr']).reset_index(drop=True)

n_permute = int(df['permuted'].sum())
logging.info(f'Effective permutations used: {n_permute}')

# calculate summary stats
df['covariate'] = df['covariate'].str.split('-', expand=True)[0].astype('category')
df['covariate_type'] = 'numeric' if is_numeric else 'categorical'
df['n_covariates'] = n_covariate
df['perm_mean'] = df.loc[df['permuted'], 'pcr'].mean()
df['perm_std'] = df.loc[df['permuted'], 'pcr'].std()
df['z_score'] = (df['pcr'] - df['perm_mean']) / df['perm_std']
if n_permute == 0:
    df['p-val'] = np.nan
    df['signif'] = ''
elif n_permute < 100:
    df['signif'] = df['z_score'].apply(lambda x: '**' if x > 3 else '*' if x > 1.5 else '')
    df['p-val'] = df.query('permuted & z_score > 1.5').shape[0] / n_permute
else:
    stat = df['pcr'] - df['perm_mean']
    null = stat.loc[df['permuted']].values
    observed = stat.loc[~df['permuted']].values
    if observed.size != 1:
        raise ValueError('Expected exactly one non-permuted observation for p-value calculation')

    df['p-val'] = (np.sum(null >= observed[0]) + 1) / (n_permute + 1)
    df['signif'] = df['p-val'].apply(lambda x: '**' if x <= 0.01 else '*' if x <= 0.05 else '')

logging.info(f'PCR summary table:\n{df}')
df.to_csv(output_file, sep='\t', index=False)
