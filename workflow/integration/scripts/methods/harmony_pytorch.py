import os
print(os.environ['PYTORCH_CUDA_ALLOC_CONF'], flush=True)
import logging
logging.basicConfig(level=logging.INFO)
from pprint import pformat
from pathlib import Path
import torch
from harmony import harmonize

import dask
dask.config.set(scheduler='threads')

from integration_utils import clean_categorical_column, add_metadata, \
    remove_slots, get_hyperparams, PCA_PARAMS
from mcv import mcv_optimal_pcs_scanpy, plot_mcv_pca
from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg
from utils.misc import dask_compute
from utils.processing import sc, USE_GPU


input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_plot_dir = snakemake.output.plots
Path(output_plot_dir).mkdir(parents=True, exist_ok=True)

wildcards = snakemake.wildcards
params = dict(snakemake.params)
batch_key = wildcards.batch

dask.config.set({'num_workers': snakemake.threads})

hyperparams = params.get('hyperparams', {})
hyperparams = {} if hyperparams is None else hyperparams
params['hyperparams'] = hyperparams.copy()

scale = hyperparams.pop('scale', False)
pca_kwargs, hyperparams = get_hyperparams(
    hyperparams=hyperparams,
    model_params=PCA_PARAMS,
)
hyperparams = {'random_state': params.get('seed', 0)} | hyperparams


# set harmony var_use
keys = hyperparams.get('batch_key', [])
if keys is None:
    keys = {batch_key}
elif isinstance(keys, str):
    keys = {batch_key, keys}
elif isinstance(keys, list):
    keys = set(keys).union({batch_key})
hyperparams['batch_key'] = list(keys)

# check GPU
use_gpu = torch.cuda.is_available()
logging.info(f'GPU available: {use_gpu}')

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

if pca_kwargs.get('n_comps') == 'mcv':
    adata.layers['raw'] = read_anndata(
        input_file,
        X='layers/counts',
        dask=True,
        backed=True,
    ).X

for column in keys:
    clean_categorical_column(adata, column)

# subset features
adata, _ = subset_hvg(
    adata,
    var_column='integration_features',
    compute_dask=False
)
logging.info(f'Subset features: {adata.shape}')

# Param optimisation
if pca_kwargs.get('n_comps') == 'mcv':
    logging.info('Find optimal n PCs via molecular cross validation (MCV)...')
    optimal_k, mcv_summary = mcv_optimal_pcs_scanpy(
        adata,
        raw_name='raw',
        max_pcs=100,
        scale_for_pca=scale,
    )
    logging.info(f'Optimal n_comps for PCA: {optimal_k}')
    pca_kwargs['n_comps'] = optimal_k
    params['hyperparams']['n_comps'] = optimal_k

    plot_mcv_pca(mcv_summary, figdir=output_plot_dir)

if USE_GPU:
    sc.get.anndata_to_GPU(adata)

if scale:
    logging.info('Scale counts...')
    sc.pp.scale(adata, zero_center=False)

# recompute PCA according to user-defined hyperparameters
logging.info(f'Compute PCA with parameters {pformat(pca_kwargs)}...')
pca_kwargs['svd_solver'] = pca_kwargs.get('svd_solver', 'covariance_eigh')
sc.pp.pca(adata, **pca_kwargs)
adata = dask_compute(adata, layers='X_pca')
del adata.X

if USE_GPU:
    import cupy as cp
    adata.obsm['X_pca'] = cp.asnumpy(adata.obsm['X_pca'])

# run method
logging.info(f'Run Harmony pytorch with parameters {pformat(hyperparams)}...')
try:
    adata.obsm['X_emb'] = harmonize(
        adata.obsm['X_pca'],
        adata.obs,
        use_gpu=use_gpu,
        n_jobs=snakemake.threads,
        **hyperparams
    )
except torch.OutOfMemoryError as e:
    logging.warning(f'GPU OOM with error: {e}')
    logging.info('Fallback to CPU...')
    adata.obsm['X_emb'] = harmonize(
        adata.obsm['X_pca'],
        adata.obs,
        use_gpu=False,
        n_jobs=snakemake.threads,
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
