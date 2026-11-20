"""
Highly variable gene selection with Scanpy
"""
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
import warnings
warnings.filterwarnings("ignore", message="The frame.append method is deprecated and will be removed from pandas in a future version.")
import anndata as ad
from dask import array as da
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)
from pprint import pformat

from utils.io import read_anndata, write_zarr_linked
from utils.accessors import _filter_batch
from utils.processing import sc, USE_GPU

input_file = snakemake.input[0]
output_file = snakemake.output[0]
args = snakemake.params.get('args', {})

if args is None:
    args = {}
elif isinstance(args, dict):
    args.pop('subset', None) # don't support subsetting

logging.info(f'HVG args:\n{pformat(args)}')

dask = snakemake.params.get('dask', True) # get global dask flag
if isinstance(args, dict):
    dask = args.pop('dask', dask) # overwrite with hvg-specific dask flag

logging.info(f'args: {args}')
logging.info(f'dask: {dask}')

logging.info(f'Read {input_file}...')
layer = 'layers/normcounts'
if args and args.get('flavor') in ('seurat_v3', 'seurat_v3_paper'):
    # both seurat_v3 flavors expect raw counts
    layer = 'layers/counts'
adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var',
    uns='uns',
    backed=dask,
    dask=dask,
)
logging.info(adata.__str__())
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
if isinstance(args, dict):
    for key in sorted(args.keys()):
        hvg_column_name += f'-{key}={args[key]}'

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.var[hvg_column_name] = True
    adata.var['highly_variable'] = True
    adata.write_zarr(output_file)
    exit(0)

if args == False:
    logging.info('No highly variable gene parameters provided, including all genes...')
    var['highly_variable'] = True
else:
    min_cells = args.pop('min_cells', max(args.get('n_bins', 20), 200))

    # filter genes and cells that would break HVG function
    batch_mask = _filter_batch(
        adata,
        batch_key=args.get('batch_key'),
        min_cells=min_cells,
    )
    adata = adata[batch_mask, adata.var['nonzero_genes']].copy()
    
    # Handle case where filtering resulted in empty data
    if adata.n_obs == 0:
        logging.info('No data after filtering, write empty file...')
        var[hvg_column_name] = True
        var['highly_variable'] = True
        adata = ad.AnnData(var=var, uns=adata.uns)
        write_zarr_linked(
            adata,
            in_dir=input_file,
            out_dir=output_file,
            files_to_keep=['uns', 'var']
        )
        exit(0)
    
    # Keep batches contiguous for count-based operations (after filtering for efficiency).
    batch_key = args.get('batch_key') if isinstance(args, dict) else None
    if batch_key in adata.obs.columns:
        adata.obs_names_make_unique()
        adata = adata[adata.obs.sort_values(batch_key, kind='stable').index].copy()
    
    # Persist dask arrays to reset task graph before expensive HVG computation.
    if dask:
        logging.info('Persist dask arrays before HVG computation...')
        adata.X = adata.X.rechunk((200_000, -1)).persist()
    
    # make sure data is on GPU for rapids_singlecell
    if USE_GPU:
        sc.get.anndata_to_GPU(adata)
    
    logging.info(f'Select features with arguments: {args}...')
    sc.pp.highly_variable_genes(adata, **args)

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

logging.info(f'Write to {output_file}...')
files_to_keep = ['uns', 'var']
adata = ad.AnnData(var=var, uns=adata.uns)
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=files_to_keep
)