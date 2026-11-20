from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc
import scanorama

from integration_utils import add_metadata, remove_slots, clean_categorical_column
from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg


# TODO: use scanpy/anndata directly
def merge_adata(adata_list, **kwargs):
    """Merge adatas from list while remove duplicated ``obs`` and ``var`` columns

    :param adata_list: ``anndata`` objects to be concatenated
    :param kwargs: arguments to be passed to ``anndata.AnnData.concatenate``
    """
    import anndata
    
    if len(adata_list) == 1:
        return adata_list[0]

    # Make sure that adatas do not contain duplicate columns
    for _adata in adata_list:
        for attr in ("obs", "var"):
            df = getattr(_adata, attr)
            dup_mask = df.columns.duplicated()
            if dup_mask.any():
                print(
                    f"Deleting duplicated keys `{list(df.columns[dup_mask].unique())}` from `adata.{attr}`."
                )
                setattr(_adata, attr, df.loc[:, ~dup_mask])

    return anndata.concat(adata_list, **kwargs)


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
batch_key = wildcards.batch
params = snakemake.params

hyperparams = params.get('hyperparams', {})
hyperparams = {} if hyperparams is None else hyperparams
hyperparams = {'seed': params.get('seed', 0)} | hyperparams

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='layers/normcounts',
    obs='obs',
    var='var',
    uns='uns',
    dask=True,
    backed=True,
)

clean_categorical_column(adata, batch_key)

# subset features
adata, _ = subset_hvg(adata, var_column='integration_features')
var = adata.var.copy()

batch_categories = adata.obs[batch_key].unique().tolist()
adatas = [
    adata[adata.obs[batch_key] == batch].copy()
    for batch in batch_categories
]

# run method
logging.info(f'Run Scanorama with parameters {pformat(hyperparams)}...')
corrected = scanorama.correct_scanpy(
    adatas,
    return_dimred=True,
    **hyperparams
)
adata = merge_adata(
    corrected,
    label=batch_key,
    keys=batch_categories,
    index_unique=None
)
# anndata.concat drops var columns, restore them
adata.var = var.loc[adata.var_names]

# save full feature output
# adata.obsm["X_full"] = adata.X

# save embedding output
adata.obsm["X_emb"] = adata.obsm["X_scanorama"]
del adata.obsm["X_scanorama"]

# prepare output adata
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    # files_to_keep=['obsm', 'uns'],
    files_to_keep=['X', 'obsm', 'var', 'uns'],
)