import logging
import numpy as np
import pandas as pd
import scanpy as sc
from pprint import pformat

from utils.io import read_anndata, write_zarr_linked
from utils.misc import dask_compute
from utils.accessors import match_genes

logging.basicConfig(level=logging.INFO)


Y_GENES = [
    "FAM197Y6",
    "FAM197Y7",
    "FAM41AY2",
    "LINC00279",
    "SRY",
    "TTTY1B",
]
X_GENES = ["XIST"]

Y_NONPAR_GENES = ["https://raw.githubusercontent.com/prabhakarlab/AIDA_Phase1/master/01_QualityControl/list_chrY_nonPAR_genes.txt"]
Y_PAR_GENES = ["https://raw.githubusercontent.com/prabhakarlab/AIDA_Phase1/master/01_QualityControl/list_chrY_PAR_genes.txt"]

input_file = snakemake.input[0]
output_file = snakemake.output[0]

layer = snakemake.params.get("layer")
is_normalized = snakemake.params.get("is_normalized", True)
params = snakemake.params.get("params", {})
donor_key = params.get("donor_key")
donors = params.get("donors")

x_genes = params.get("x_genes", X_GENES)
y_genes = params.get("y_genes", Y_GENES)
x_threshold = params.get("x_threshold", 0)
y_threshold = params.get("y_threshold", 4)
imbalance_frac = float(params.get("imbalance_frac", 0.1))

y_nonpar_genes = params.get("y_nonpar_genes", Y_NONPAR_GENES)
y_par_genes = params.get("y_par_genes", Y_PAR_GENES)

predict_key = params.get("predict_column", "sex")
reference_key = params.get("reference_key", predict_key)


def parse_genes(gene_list, adata):
    if isinstance(gene_list, str):
        gene_list = [gene_list]

    genes_use = match_genes(adata.var, gene_list, column='feature_name', return_index=False, as_list=True)
    assert len(genes_use) > 0, f"No genes found in adata.var['feature_name']:\n{pformat(gene_list)}"
    logging.info(f'Using {len(genes_use)} genes from list:\n{pformat(gene_list)}')
    return genes_use


def predict_sex_by_XY(
    adata,
    donor_key,
    predict_key,
    x_genes,
    y_genes,
    feature_column=None,
    default=np.nan,
    x_threshold=0,
    y_threshold=4,
    frac=0.1,
):
    """
    Predict donor sex based on aggregated expression of X- and Y-linked marker genes.

    This function computes per-cell mean expression of the specified X-linked genes
    and Y-linked genes, aggregates these values per donor, and then classifies each
    donor into sex categories using expression thresholds and X/Y imbalance ratios.

    The algorithm proceeds as follows:

    1. For each cell, compute the mean expression of the X marker genes (`x_genes`)
       and the Y marker genes (`y_genes`) using the expression matrix in ``adata.X``.
    2. Build a per-cell table containing the donor identifier (from ``adata.obs``)
       and the per-cell X and Y mean expression values.
    3. Aggregate this table by ``donor_key`` (summing the per-cell means) to obtain
       donor-level X and Y expression summaries.
    4. For each donor, derive the ratios

           y_x_frac = Y_expression / X_expression
           x_y_frac = X_expression / Y_expression

       and assign a sex label according to:

       * ``"female"`` if X expression is above ``x_threshold`` and Y expression is
         at or below ``y_threshold`` **or** if ``y_x_frac <= frac``.
       * ``"male"`` if X expression is at or below ``x_threshold`` and Y expression
         is above ``y_threshold`` **or** if ``x_y_frac < frac``.
       * ``"mix"`` if both X and Y expression are above their thresholds **or** if
         both ratios indicate similar contributions (``y_x_frac > frac`` and
         ``x_y_frac > frac``).
       * ``default`` for donors that do not match any of the above conditions.

    Parameters
    ----------
    adata
        AnnData-like object containing single-cell expression data. Uses
        ``adata.X`` for expression values, ``adata.var`` or ``adata.var_names`` to
        match marker genes, and ``adata.obs[donor_key]`` for donor identifiers.
    donor_key : str
        Column name in ``adata.obs`` indicating the donor ID used to group cells.
    predict_key : str
        Base name for the prediction; also used to construct the output column
        names for aggregated X and Y expression (``f"{predict_key}_x_exp"`` and
        ``f"{predict_key}_y_exp"``) and the final sex label column.
    x_genes : sequence of hashable
        Identifiers of X-linked marker genes. If ``feature_column`` is provided,
        these are matched against ``adata.var[feature_column]``; otherwise they are
        matched against ``adata.var_names``.
    y_genes : sequence of hashable
        Identifiers of Y-linked marker genes, matched in the same way as
        ``x_genes``.
    feature_column : str, optional
        Name of the column in ``adata.var`` used to match ``x_genes`` and
        ``y_genes``. If ``None``, gene identifiers are matched against
        ``adata.var_names``.
    default : str or None, optional
        Label assigned to donors that do not satisfy any of the sex classification
        rules. If ``None``, such donors will have a missing value in the
        prediction column.
    x_threshold : float, optional
        Minimum aggregated X-marker expression required to consider X expression
        "present" for a donor. Used in the classification rules described above.
    y_threshold : float, optional
        Minimum aggregated Y-marker expression required to consider Y expression
        "present" for a donor. Used in the classification rules described above.
    frac : float, optional
        Threshold on the imbalance between X and Y expression. A small value (e.g.
        0.1) requires that one chromosome's markers be much lower than the other's
        to classify a donor as unambiguously male or female based on expression
        ratios.

    Returns
    -------
    pandas.DataFrame
        A DataFrame indexed by donor (unique values of ``donor_key``) containing
        the aggregated X and Y expression columns (named
        ``f"{predict_key}_x_exp"`` and ``f"{predict_key}_y_exp"``) and a
        categorical column ``predict_key`` holding the predicted sex label for
        each donor.
    """
    from fast_array_utils import stats

    x_exp = f"{predict_key}_x_exp"
    y_exp = f"{predict_key}_y_exp"

    if not feature_column:
        feature_column = 'feature_name'
        adata.var[feature_column] = adata.var_names.astype(str)

    gene_mask = adata.var[feature_column].isin(x_genes+y_genes)
    print(f'Subset to {sum(gene_mask)} genes...', flush=True)
    adata = dask_compute(adata[:, gene_mask].copy(), layers='X')

    print('Scale data...', flush=True)
    sc.pp.scale(adata, zero_center=False)

    # Per-cell expression (sparse-safe)
    x_exp_cell = stats.mean(adata[:, adata.var[feature_column].isin(x_genes)].X, axis=1)
    y_exp_cell = stats.mean(adata[:, adata.var[feature_column].isin(y_genes)].X, axis=1)

    # Build per-cell table
    df = pd.DataFrame(
        {
            donor_key: adata.obs[donor_key].to_numpy(),
            x_exp: x_exp_cell,
            y_exp: y_exp_cell,
        }
    )

    # Aggregate per donor
    donor_exp = df.groupby(donor_key, sort=False, observed=True).sum()

    # Vectorized classification with safe division
    x_vals = donor_exp[x_exp].to_numpy(dtype=float)
    y_vals = donor_exp[y_exp].to_numpy(dtype=float)
    y_x_frac = np.divide(
        y_vals,
        x_vals,
        out=np.full_like(y_vals, np.inf),
        where=x_vals != 0,
    )
    x_y_frac = np.divide(
        x_vals,
        y_vals,
        out=np.full_like(x_vals, np.inf),
        where=y_vals != 0,
    )
    donor_exp[predict_key] = np.select(
        [
            ((donor_exp[x_exp] > x_threshold) & (donor_exp[y_exp] <= y_threshold)) | (y_x_frac <= frac),
            ((donor_exp[x_exp] <= x_threshold) & (donor_exp[y_exp] > y_threshold)) | (x_y_frac < frac),
            ((donor_exp[x_exp] > x_threshold) & (donor_exp[y_exp] > y_threshold)) | ((y_x_frac > frac) & (x_y_frac > frac)),
        ],
        ["female", "male", "mix"],
        default=default,
    )

    return donor_exp


def predict_sex_by_chrY_ratio(
    adata,
    donor_key,
    y_nonpar_genes,
    y_par_genes,
    feature_column=None,
    threshold=0.5,
    default=np.nan,
    predict_prefix="sex_chrY_",
    nonpar_col="nonPAR_UMIs",
    par_col="PAR_UMIs",
    ratio_col="nonPAR_to_PAR",
):
    """
    Predict donor sex using the ratio of chrY non-PAR to PAR UMIs.

    This function mirrors the provided R logic:

    1. Compute per-cell UMI sums for chrY non-PAR genes and chrY PAR genes.
    2. Aggregate these sums per donor (``donor_key``).
    3. Compute the ratio ``nonPAR_to_PAR = nonPAR_UMIs / PAR_UMIs`` per donor.
    4. Classify sex as ``"male"`` if the ratio is >= ``threshold`` (default 0.5),
       otherwise ``"female"``. If the ratio is undefined due to zero PAR counts,
       it is treated as ``inf`` and thus classified as ``"male"``.

    Parameters
    ----------
    adata
        AnnData-like object with expression in ``adata.X`` and donor IDs in
        ``adata.obs[donor_key]``.
    donor_key : str
        Column name in ``adata.obs`` indicating the donor ID used to group cells.
    y_nonpar_genes : sequence
        Gene identifiers for chrY non-PAR genes.
    y_par_genes : sequence
        Gene identifiers for chrY PAR genes.
    feature_column : str, optional
        Column in ``adata.var`` to match genes. If ``None``, matches against
        ``adata.var_names``.
    threshold : float, optional
        Cutoff for ``nonPAR_to_PAR`` above which donors are classified as male.
    default : str or None, optional
        Label for donors with missing ratios if desired.
    predict_key : str, optional
        Name of the output sex prediction column.
    nonpar_col, par_col : str, optional
        Column names for aggregated non-PAR and PAR UMIs in the output.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by donor containing ``nonpar_col``, ``par_col``,
        ``nonPAR_to_PAR``, and a ``predict_key`` column with values
        ``"male"``/``"female"``.
    """
    from fast_array_utils import stats

    if predict_prefix:
        nonpar_col = f"{predict_prefix}{nonpar_col}"
        par_col = f"{predict_prefix}{par_col}"
        ratio_col = f"{predict_prefix}{ratio_col}"

    if feature_column:
        nonpar_mask = adata.var[feature_column].isin(y_nonpar_genes)
        par_mask = adata.var[feature_column].isin(y_par_genes)
    else:
        nonpar_mask = adata.var_names.isin(y_nonpar_genes)
        par_mask = adata.var_names.isin(y_par_genes)

    nonpar_cell = stats.sum(adata[:, nonpar_mask].X, axis=1)
    par_cell = stats.sum(adata[:, par_mask].X, axis=1)

    df = pd.DataFrame(
        {
            donor_key: adata.obs[donor_key].to_numpy(),
            nonpar_col: nonpar_cell,
            par_col: par_cell,
        }
    )

    donor_chrY = df.groupby(donor_key, sort=False, observed=True).sum()

    nonpar_vals = donor_chrY[nonpar_col].to_numpy(dtype=float)
    par_vals = donor_chrY[par_col].to_numpy(dtype=float)
    ratio = np.divide(
        nonpar_vals,
        par_vals,
        out=np.full_like(nonpar_vals, np.inf),
        where=par_vals != 0,
    )
    donor_chrY[ratio_col] = ratio

    donor_chrY[predict_prefix] = np.where(
        donor_chrY[ratio_col] >= threshold,
        "male",
        "female",
    )

    if default is not None:
        # Only apply the default when both PAR and nonPAR counts are zero
        mask_nan = (
            ~np.isfinite(donor_chrY[ratio_col])
            & (nonpar_vals == 0)
            & (par_vals == 0)
        )
        donor_chrY.loc[mask_nan, predict_prefix] = default

    return donor_chrY


logging.info(f"Read file: {input_file}...")
adata = read_anndata(
    input_file,
    dask=True,
    backed=True,
    X=layer,
    obs="obs",
    var="var",
    uns="uns",
)

# check donors
assert donor_key in adata.obs.columns, f'"{donor_key}" not in adata.obs.columns'
columns = [donor_key] + ([reference_key] if reference_key in adata.obs.columns else [])
adata.obs = adata.obs[columns]

if not donors:
    donors = adata.obs[donor_key].unique().tolist()
else:
    adata = adata[adata.obs[donor_key].isin(donors)].copy()

# check all donors are in adata
missing_donors = set(donors) - set(adata.obs[donor_key].unique().tolist())
assert len(missing_donors) == 0, f"Donors not found in adata: {missing_donors}"

# parse genes
if 'feature_name' not in adata.var.columns:
    adata.var['feature_name'] = adata.var_names.astype(str)

x_genes = parse_genes(x_genes, adata)
y_genes = parse_genes(y_genes, adata)
y_nonpar_genes = parse_genes(y_nonpar_genes, adata)
y_par_genes = parse_genes(y_par_genes, adata)
all_genes = set(x_genes) | set(y_genes) | set(y_nonpar_genes) | set(y_par_genes)

if not is_normalized:
    logging.info('Normalizing and log-transforming data...')
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

logging.info(f'Predict sex by X and Y gene expression for {len(donors)} donors...')
donor_exp = predict_sex_by_XY(
    adata,
    donor_key=donor_key,
    predict_key=predict_key,
    x_genes=x_genes,
    y_genes=y_genes,
    default=None,
    feature_column='feature_name',
    x_threshold=x_threshold,
    y_threshold=y_threshold,
    frac=imbalance_frac,
)
adata.obs = adata.obs.join(donor_exp[predict_key], on=donor_key)

logging.info(f'Predict sex by chrY non-PAR to PAR ratio for {len(donors)} donors...')
predict_chry = f"{predict_key}_chrY"
donor_chrY = predict_sex_by_chrY_ratio(
    adata,
    donor_key=donor_key,
    y_nonpar_genes=y_nonpar_genes,
    y_par_genes=y_par_genes,
    feature_column='feature_name',
    threshold=params.get("chrY_threshold", 0.5),
    predict_prefix=predict_chry,
)
adata.obs = adata.obs.join(donor_chrY, on=donor_key)

adata.uns['predict_sex'] = {
    "X_Y_expression": {
        "predictions": donor_exp,
        "x_genes": x_genes,
        "y_genes": y_genes,
        "x_threshold": x_threshold,
        "y_threshold": y_threshold,
        "imbalance_frac": imbalance_frac,
    },
    "chrY_ratio": {
        "y_nonpar_genes": y_nonpar_genes,
        "y_par_genes": y_par_genes,
        "predictions": donor_chrY,
    },
}

logging.info(f'Predicted sex for {len(donors)} donors.')
df = adata.obs
if reference_key in df.columns:
    # subset to reference with non-missing values
    ref_mask = df[reference_key].isin(["male", "female"])

    # expression-based prediction
    exp_mismatch_donors = df[(df[reference_key] != df[predict_key])][donor_key].dropna().unique()
    donor_exp = donor_exp.loc[exp_mismatch_donors]
    accuracy = (df.loc[ref_mask, reference_key] == df.loc[ref_mask, predict_key]).mean()
    logging.info(f'X_Y_expression accuracy: {accuracy:.2%}\n{donor_exp.to_string()}')

    adata.uns['predict_sex']['X_Y_expression']['accuracy'] = accuracy
    adata.uns['predict_sex']['X_Y_expression']['mismatched'] = donor_exp

    # chrY-based prediction
    chrY_mismatch_donors = df[(df[reference_key] != df[f"{predict_key}_chrY"])][donor_key].dropna().unique()
    donor_chrY = donor_chrY.loc[chrY_mismatch_donors]
    accuracy = (df.loc[ref_mask, reference_key] == df.loc[ref_mask, f"{predict_key}_chrY"]).mean()
    logging.info(f'chrY_ratio accuracy: {accuracy:.2%}\n{donor_chrY.to_string()}')

    adata.uns['predict_sex']['chrY_ratio']['accuracy'] = accuracy
    adata.uns['predict_sex']['chrY_ratio']['mismatched'] = donor_chrY

    # subset to union of mismatches
    mismatched_donors = set(exp_mismatch_donors) | set(chrY_mismatch_donors)
    df = df[df[donor_key].isin(mismatched_donors)].copy()
else:
    logging.info(f'X_Y_expression:\n{donor_exp.to_string()}')
    logging.info(f'chrY_ratio:\n{donor_chrY.to_string()}')

value_counts = df.value_counts(dropna=False).reset_index().set_index(donor_key)
adata.uns['predict_sex']['prediction_counts'] = value_counts
logging.info(f'Predictions:\n{value_counts}')

logging.info(f'Write file: {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obs', 'uns'],
)
