from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc

from integration_utils import clean_categorical_column, add_metadata, \
    remove_slots, get_hyperparams, PCA_PARAMS
from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg
from utils.misc import dask_compute


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch

hyperparams = params.get('hyperparams', {})
hyperparams = {} if hyperparams is None else hyperparams

pca_kwargs, hyperparams = get_hyperparams(
    hyperparams=hyperparams,
    model_params=PCA_PARAMS,
)
hyperparams = {'pynndescent_random_state': params.get('seed', 0)} | hyperparams

files_to_keep = ['obsm', 'obsp', 'uns']

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
adata, _ = subset_hvg(
    adata,
    var_column='integration_features',
    compute_dask=False,
)

# quickfix: remove batches with fewer than 3 cells
neighbors_within_batch = hyperparams.get('neighbors_within_batch', 3)
cells_per_batch = adata.obs.value_counts(batch_key)
min_batches = cells_per_batch[cells_per_batch < neighbors_within_batch].index.tolist()

if len(min_batches) > 0:
    files_to_keep.extend(['obs'])
    adata = adata[~adata.obs[batch_key].isin(min_batches)].copy()
    # Example alternative: adata = adata[adata.obs[batch_key].isin(min_batches)]

# recompute PCA according to user-defined hyperparameters
logging.info(f'Compute PCA with parameters {pformat(pca_kwargs)}...')
use_rep = 'X_pca'
sc.pp.pca(adata, **pca_kwargs)
dask_compute(adata, layers=use_rep)
del adata.X

# run method
logging.info(f'Run BBKNN with parameters {pformat(hyperparams)}...')
sc.external.pp.bbknn(
    adata,
    batch_key=batch_key,
    use_rep=use_rep,
    copy=False,
    **hyperparams,
)

# prepare output adata
X_pca = adata.obsm[use_rep]
adata = remove_slots(adata=adata, output_type=params['output_type'])
# keep the representation the graph was built on, as recorded in .uns['neighbors'],
# so that downstream steps (e.g. UMAP) don't recompute it from the full feature matrix
adata.obsm[use_rep] = X_pca
add_metadata(adata, wildcards, params)

logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=files_to_keep,
)