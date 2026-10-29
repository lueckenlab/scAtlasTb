"""
Highly variable gene selection with Scanpy
"""
import logging
logging.basicConfig(level=logging.INFO)
import math
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)

from utils.io import read_anndata, write_zarr_linked
from utils.processing import _filter_genes

input_file = snakemake.input[0]
output_file = snakemake.output[0]
dask = snakemake.params.get('dask', True)
args = snakemake.params.get('args') or {}

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='X',
    var='var',
    backed=dask,
    dask=dask,
    stride=100_000,
)
logging.info('%s', adata)
var = adata.var.copy()

if 'min_cells' not in args and 'min_counts' not in args:
    # Use a stronger adaptive default pre-filter to reduce downstream HVG cost.
    args['min_cells'] = max(10, min(200, int(math.ceil(0.001 * adata.n_obs))))

logging.info('Filter gene args: %s', args)
logging.info('Determine nonzero genes...')
adata.var['nonzero_genes'] = _filter_genes(adata, **args)

write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['var/nonzero_genes'],
)