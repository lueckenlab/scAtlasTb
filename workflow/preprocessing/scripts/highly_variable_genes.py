"""
Highly variable gene selection with Scanpy
"""
import logging
logging.basicConfig(level=logging.INFO)
import time
import warnings
warnings.filterwarnings("ignore", message="The frame.append method is deprecated and will be removed from pandas in a future version.")
import anndata as ad
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)
from pprint import pformat

from utils.io import read_anndata, write_zarr_linked
from utils.accessors import _filter_batch
from utils.processing import sc, USE_GPU

input_file = snakemake.input[0]
output_file = snakemake.output[0]
args = snakemake.params.get('args') or {}

if isinstance(args, dict):
    args.pop('subset', None) # don't support subsetting

run_hvg = args is not False
hvg_args = dict(args) if isinstance(args, dict) else {}


def _log_timing(start_ts, step):
    elapsed = time.perf_counter() - start_ts
    logging.info(f'{step} finished in {elapsed:.2f}s')

logging.info(f'HVG args:\n{pformat(args)}')

dask = snakemake.params.get('dask', True) # get global dask flag
dask = hvg_args.pop('dask', dask) # overwrite with hvg-specific dask flag

logging.info('args: %s', args)
logging.info('dask: %s', dask)

logging.info('Read %s...', input_file)
layer = 'layers/normcounts'
if hvg_args.get('flavor') == 'seurat_v3':
    layer = 'layers/counts'
adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var',
    uns='uns',
    backed=True,
    dask=True,
)
logging.info('%s', adata)
var = adata.var.copy()

# add metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}
adata.uns['preprocessing']['highly_variable_genes'] = args

# Remove any previous highly_variable_genes columns
for col in list(var.columns):
    if col.startswith('highly_variable'):
        var.drop(columns=col, inplace=True)

# Determine HVG column name based on args (use parameter hash in filename as suffix)
hvg_column_name = 'highly_variable'
for key in sorted(hvg_args.keys()):
    hvg_column_name += f'-{key}={hvg_args[key]}'


def _write_empty_output(adata_obj, message, linked):
    logging.info(message)
    target_var = var if linked else adata_obj.var
    target_var[hvg_column_name] = True
    target_var['highly_variable'] = True

    if linked:
        empty_adata = ad.AnnData(var=var, uns=adata_obj.uns)
        write_zarr_linked(
            empty_adata,
            in_dir=input_file,
            out_dir=output_file,
            files_to_keep=['uns', 'var']
        )
    else:
        adata_obj.write_zarr(output_file)

    raise SystemExit(0)


if adata.n_obs == 0:
    _write_empty_output(adata, 'No data, write empty file...', linked=False)

if not run_hvg:
    logging.info('No highly variable gene parameters provided, including all genes...')
    var['highly_variable'] = True
else:
    batch_key = hvg_args.get('batch_key')
    min_cells = hvg_args.pop('min_cells', max(hvg_args.get('n_bins', 20), 200))
    large_batch_threshold = int(hvg_args.pop('large_batch_threshold', 20))
    persist_dask = hvg_args.pop('persist_dask', dask)

    # filter genes and cells that would break HVG function
    t0 = time.perf_counter()
    batch_mask = _filter_batch(
        adata,
        batch_key=batch_key,
        min_cells=min_cells,
    )
    adata = adata[batch_mask, adata.var['nonzero_genes']].copy()
    _log_timing(t0, 'Batch/cell filtering')
    
    # Handle case where filtering resulted in empty data
    if adata.n_obs == 0:
        _write_empty_output(adata, 'No data after filtering, write empty file...', linked=True)
    
    # Gather batch metadata for heuristics/logging.
    n_batches = adata.obs[batch_key].nunique() if batch_key in adata.obs.columns else 0
    logging.info(f'Filtered data shape: n_obs={adata.n_obs}, n_vars={adata.n_vars}, n_batches={n_batches}')

    # For many batches, dask task overhead commonly dominates HVG runtime.
    if dask and n_batches >= large_batch_threshold:
        logging.info(
            'Detected many batches (n_batches=%s >= %s); switching to in-memory compute for HVG.',
            n_batches,
            large_batch_threshold,
        )
        dask = False
    
    if dask:
        logging.info('Dask HVG path: persist_dask=%s', persist_dask)
        t0 = time.perf_counter()
        if persist_dask:
            adata.X = adata.X.persist()
        _log_timing(t0, 'Dask preparation')
    else:
        logging.info('Compute X before HVG computation...')
        t0 = time.perf_counter()
        adata.X = adata.X.compute()
        _log_timing(t0, 'Compute X')
    
    # make sure data is on GPU for rapids_singlecell
    if USE_GPU:
        logging.info('Move data to GPU for RAPIDS SingleCell...')
        t0 = time.perf_counter()
        sc.get.anndata_to_GPU(adata)
        _log_timing(t0, 'Move data to GPU')
    
    logging.info('Select features with arguments: %s...', hvg_args)
    t0 = time.perf_counter()
    sc.pp.highly_variable_genes(adata, **hvg_args)
    _log_timing(t0, 'HVG computation')

    # add HVG info back to adata
    hvg_column_map = {
        'highly_variable': False,
        'means': 0,
        'dispersions': 0,
        'dispersions_norm': 0,
        'highly_variable_nbatches': 0,
        'highly_variable_intersection': False,
    }
    for column, default_value in hvg_column_map.items():
        if column not in adata.var.columns:
            continue
        dtype = adata.var[column].dtype
        var[column] = default_value
        var[column] = var[column].astype(dtype)
        var.loc[adata.var_names, column] = adata.var[column]

# rename highly_variable column
var[hvg_column_name] = var['highly_variable']

logging.info('Write to %s...', output_file)
adata = ad.AnnData(var=var, uns=adata.uns)
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['uns', 'var']
)