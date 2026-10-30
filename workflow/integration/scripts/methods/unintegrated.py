from pathlib import Path
import scanpy as sc
import logging
logging.basicConfig(level=logging.INFO)

import dask
dask.config.set(scheduler='threads')

from integration_utils import add_metadata, remove_slots, clean_categorical_column
from utils.io import read_anndata, write_zarr_linked
from utils.processing import assert_neighbors
from utils.accessors import subset_hvg
from utils.misc import dask_compute

input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

dask.config.set({'num_workers': snakemake.threads})

adata = read_anndata(
    input_file,
    X='layers/normcounts',
    obs='obs',
    var='var',
    obsp='obsp',
    obsm='obsm',
    varm='varm',
    uns='uns',
    dask=True,
    backed=True,
)

clean_categorical_column(adata, wildcards.batch)

# subset features
adata, _ = subset_hvg(
    adata,
    var_column='integration_features',
    compute_dask=False,
)

# prepare output adata
files_to_keep = ['obsm', 'uns']

if 'X_pca' not in adata.obsm:
    logging.info('Compute PCA...')
    sc.pp.pca(adata, svd_solver='covariance_eigh')
adata.obsm['X_emb'] = adata.obsm['X_pca']
del adata.obsm['X_pca']
dask_compute(adata, layers='X_emb')

logging.info(adata.__str__())
logging.info(adata.uns.keys())
try:
    assert_neighbors(adata)
    print('Neighbors params:', adata.uns['neighbors'], flush=True)
except AssertionError as e:
    logging.info(f'Recompute neighbors due to {e}...')
    sc.pp.neighbors(adata, use_rep='X_emb')
    print('Neighbors params:', adata.uns['neighbors'], flush=True)
    files_to_keep.extend(['obsp', 'uns'])

del adata.X
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

# write file
logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=files_to_keep,
    slot_map={
        'X': 'layers/normcounts',
        # 'obsm/X_full': 'layers/normcounts',
    },
)