from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
from scanpy.pp import pca

import dask
dask.config.set(scheduler='threads')

from integration_utils import add_metadata, remove_slots, get_hyperparams, PCA_PARAMS, \
    clean_categorical_column
from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg
from utils.misc import dask_compute
from utils.processing import sc, USE_GPU

if USE_GPU:
    harmony_integrate = sc.pp.harmony_integrate
else:
    harmony_integrate = sc.external.pp.harmony_integrate


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch

dask.config.set({'num_workers': snakemake.threads})

hyperparams = params.get('hyperparams', {})
hyperparams = {} if hyperparams is None else hyperparams

pca_kwargs, hyperparams = get_hyperparams(
    hyperparams=hyperparams,
    model_params=PCA_PARAMS,
)
hyperparams = {'random_state': params.get('seed', 0)} | hyperparams

# set harmony var_use
keys = hyperparams.get('key', [])
if keys is None:
    keys = {batch_key}
elif isinstance(keys, str):
    keys = {batch_key, keys}
elif isinstance(keys, list):
    keys = set(keys).union({batch_key})
hyperparams['key'] = list(keys)

if USE_GPU:
    hyperparams['dtype'] = 'float32'

logging.info
(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='layers/normcounts',
    obs='obs',
    var='var',
    uns='uns',
    dask=True,
    backed=True,
)

for column in keys:
    clean_categorical_column(adata, column)

# subset features
adata, _ = subset_hvg(
    adata,
    var_column='integration_features',
    compute_dask=True
)

# recompute PCA according to user-defined hyperparameters
logging.info(f'Compute PCA with parameters {pformat(pca_kwargs)}...')
use_rep = 'X_pca'
# adata.X = adata.X.map_blocks(lambda x: x.toarray(), dtype=adata.X.dtype)
pca(adata, **pca_kwargs)
del adata.X
# dask_compute(adata, layers=use_rep)

# run method
logging.info(f'Run harmonypy with parameters {pformat(hyperparams)}...')
harmony_integrate(
    adata,
    basis=use_rep,
    adjusted_basis='X_emb',
    **hyperparams
)

# prepare output adata
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns'],
)