from pathlib import Path
import logging
from pprint import pformat
import numpy as np
import anndata as ad
import scanpy as sc

from utils.accessors import subset_hvg
from utils.io import read_anndata, write_zarr_linked, write_zarr
from utils.processing import compute_neighbors, _filter_genes
from utils.misc import dask_compute, apply_layers


def compute_pca(adata, matrix):
    X_pca, _, variance_ratio, variance = sc.tl.pca(matrix, return_info=True)
    adata.obsm['X_pca'] = X_pca
    adata.uns['pca'] = {
        'variance_ratio': variance_ratio,
        'variance': variance,
    }


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params
neighbor_args = params.get('neighbor_args', {})
unintegrated_layer = params.get('unintegrated_layer', 'X')
corrected_layer = params.get('corrected_layer', 'X')
var_key = params.get('var_mask', 'highly_variable')
output_type = params.get('output_type', 'embed')
recompute_neighbors = params.get('recompute_neighbors', False)

PERSIST_MATRIX_THRESHOLD = params.get('PERSIST_MATRIX_THRESHOLD', 5e5)

files_to_keep = ['uns', 'var']
slot_map = {}

# determine output types
output_type = read_anndata(
    input_file,
    uns='uns',
    verbose=False
).uns.get('output_type', output_type)

logging.info(f'Output type: {output_type}')
kwargs = dict(
    obs='obs',
    obsm='obsm',
    obsp='obsp',
    var='var',
    uns='uns',
)

is_h5ad = input_file.endswith('.h5ad')
if is_h5ad:
    kwargs |= dict(
        X=corrected_layer,
        layers='layers',
        dask=True,
        backed=True
    )
    files_to_keep.extend(['obs', 'X', 'layers'])

if output_type == 'full':
    kwargs |= dict(X=corrected_layer, dask=True, backed=True)
    slot_map |= dict(X=corrected_layer)

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, **kwargs)
adata.uns['output_type'] = output_type # ensure that output type is set for run.py

if is_h5ad:
    adata.layers['unintegrated'] = read_anndata(
        input_file,
        X=unintegrated_layer,
        dask=True,
        backed=True,
    ).X

all_obs_names = adata.obs_names.copy()
all_var_names = adata.var_names.copy()

# remove cells without labels
n_obs = adata.n_obs
# logging.info('Filtering out cells without labels')
# TODO: only for metrics that require labels?
# logging.info(f'Before: {adata.n_obs} cells')
# adata = adata[(adata.obs[label_key].notna() | adata.obs[label_key] != 'NA') ]
# logging.info(f'After: {adata.n_obs} cells')
# if adata.is_view:
#     adata = adata.copy()
force_neighbors = (
    n_obs > adata.n_obs
    or not {'connectivities', 'distances'}.issubset(adata.obsp)
    or recompute_neighbors
)

# set HVGs
if var_key is None:
    var_key = "highly_variable"
    logging.info(
        f'var_key set to None, assuming no feature selection desired and setting new "{var_key}" to True'
        f'\n{pformat(adata.var.columns)}'
    )
    adata.var[var_key] = True
else:
    if var_key not in adata.var.columns:
        if adata.var.shape[1] == 0:
            # Integration output has empty .var (features already subsetted by
            # var_mask during integration_prepare); treat all remaining genes
            # as selected, matching the var_key=None branch above.
            logging.info(
                f'adata.var has no columns; assuming integration already applied '
                f'"{var_key}" mask, setting all True'
            )
            adata.var[var_key] = True
        else:
            raise AssertionError(
                f'"{var_key}" not in adata var columns: {adata.var.columns.tolist()}'
            )

logging.info(f'Set "{var_key}" in adata.var: {adata.var[var_key].sum()} HVGs')

logging.info('Compute PCA...')
if output_type == 'full':
    sc.pp.pca(
        adata,
        mask_var=var_key,
        svd_solver='covariance_eigh',
    )
    adata = dask_compute(adata, layers='X_pca')
    adata.obsm['X_emb'] = adata.obsm['X_pca']
    force_neighbors = True
    files_to_keep.append('obsm/X_pca')
elif output_type == 'embed':
    logging.info('Run PCA on embedding...')
    compute_pca(adata, matrix=adata.obsm['X_emb'])
    files_to_keep.append('obsm/X_pca')

if force_neighbors:
    files_to_keep.append('obsp')

logging.info(f'Computing neighbors for output type {output_type} force_neighbors={force_neighbors}...')
compute_neighbors(
    adata,
    output_type,
    force=force_neighbors,
    check_n_neighbors=False,
    **neighbor_args
)

# remove any obsm keys except for X_pca and X_emb to save space
for key in list(adata.obsm.keys()):
    if key in ['X_pca', 'X_emb']:
        continue
    del adata.obsm[key]

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
if is_h5ad:
    write_zarr(adata, output_file)
else:
    write_zarr_linked(
        adata,
        in_dir=input_file,
        out_dir=output_file,
        files_to_keep=files_to_keep,
        slot_map=slot_map,
    )
