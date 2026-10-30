from pathlib import Path
import scanpy as sc
import warnings
warnings.simplefilter("ignore", UserWarning)
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr, write_zarr_linked, ALL_SLOTS
from utils.misc import dask_compute
import scatlastb_utils as atl


sc.set_figure_params(dpi=100, frameon=False)
input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr
output_bulk_zarr = snakemake.output.bulks

sample_key = snakemake.params.get('sample_key')
aggregate = snakemake.params.get('aggregate', 'sum')
layer = snakemake.params.get('layer', 'X')

logging.info(f'Read "{input_zarr}"...')
n_obs = read_anndata(input_zarr, X=layer, dask=True, backed=True, verbose=False).n_obs
dask = snakemake.params.get('dask', n_obs > 5e5)

adata = read_anndata(
    input_zarr,
    X=layer,
    obs='obs',
    var='var',
    backed=dask,
    dask=dask,
    stride=min(500_000, int(n_obs / 5)),
)
adata.obs_names_make_unique()

logging.info(f'Pseudobulk by "{sample_key}"...')
if sample_key is None or not str(sample_key).strip() or sample_key == 'None':
    logging.info('No sample_key provided. Skipping pseudobulk and keeping full cell resolution.')
    adata.obs['group'] = adata.obs_names
    adata_bulk = adata
else:
    sample_key = [x.strip() for x in str(sample_key).split(',') if x.strip()]
    missing_keys = [key for key in sample_key if key not in adata.obs.columns]
    if missing_keys:
        raise KeyError(
            f'sample_key contains columns not found in adata.obs: {missing_keys}. '
            f'Available columns: {list(adata.obs.columns)}'
        )
    adata.obs['group'] = adata.obs[sample_key].astype(str).agg('-'.join, axis=1)
    adata_bulk = atl.pp.pseudobulk(
        adata,
        group_key='group',
        agg=aggregate,
        use_legacy=True,
    )

# preprocess pseudobulks
sc.pp.normalize_total(adata_bulk)
sc.pp.log1p(adata_bulk)
logging.info(adata_bulk.__str__())

logging.info(f'Write "{output_bulk_zarr}"...')
write_zarr(adata_bulk, output_bulk_zarr, compute=True)

logging.info(f'Write "{output_zarr}"...')
if Path(input_zarr).suffix == '.h5ad':
    obs = adata.obs
    del adata
    logging.info('Read full h5ad file...')
    adata = read_anndata(input_zarr, dask=True, backed=True, verbose=False)
    adata.obs = obs

write_zarr_linked(
    adata,
    in_dir=input_zarr,
    out_dir=output_zarr,
    files_to_keep=['obs'],
)