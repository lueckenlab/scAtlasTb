"""
Highly variable gene selection
- HVG by group -> take union of HVGs from each group
- allow including user-specified genes
"""
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")
from tqdm.dask import TqdmCallback
from dask import config as da_config
da_config.set(num_workers=snakemake.threads)
import numpy as np
import anndata as ad

from utils.io import read_anndata, write_zarr_linked
from utils.accessors import _filter_batch, match_genes
from utils.misc import dask_compute
from utils.processing import sc, USE_GPU
rsc = sc


input_file = snakemake.input[0]
output_file = snakemake.output[0]
args = snakemake.params.get('args', {})
extra_hvg_args = snakemake.params.get('extra_hvgs', {})
overwrite_args = snakemake.params.get('overwrite_args', {})
union_over = extra_hvg_args.get('union_over')
extra_genes = extra_hvg_args.get('extra_genes', [])
remove_genes = extra_hvg_args.get('remove_genes', [])
min_cells = extra_hvg_args.pop('min_cells', max(args.get('n_bins', 20), 200))

logging.info(f'Extra HVG args:\n{pformat(extra_hvg_args)}')
logging.info(f'Overwrite args:\n{pformat(overwrite_args)}')

dask = snakemake.params.get('dask', True) # get global dask flag
dask = args.pop('dask', dask) # overwrite with highly_variable-specific dask flag
dask = extra_hvg_args.pop('dask', dask) # overwrite with extra_hvgs-specific dask flag
logging.info(f'Dask enabled: {dask}')

hvg_column_name = 'extra_hvgs'
use_gpu = USE_GPU

if not args:
    args = {}
elif isinstance(args, dict):
    args.pop('subset', None) # don't support subsetting
if isinstance(overwrite_args, dict):
    args |= overwrite_args
    for key in sorted(overwrite_args.keys()):
        hvg_column_name += f'-{key}={overwrite_args[key]}'

logging.info(f'args: {args}')
logging.info(f'dask: {dask}')

logging.info(f'Read {input_file}...')
# both seurat_v3 flavors expect raw counts
layer = 'layers/counts' if args.get('flavor') in ('seurat_v3', 'seurat_v3_paper') else 'layers/normcounts'
adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var',
    uns='uns',
    backed=True,
    dask=True,
)
logging.info(adata.__str__())
for col in adata.var.columns:
    if col.startswith('extra_hvgs'):
        del adata.var[col]  # remove old HVG columns

# if adata.n_obs > 2e6:
#     use_gpu = False
#     sc = scanpy

# add metadata
if 'preprocessing' not in adata.uns:
    adata.uns['preprocessing'] = {}
adata.uns['preprocessing'][hvg_column_name] = args | extra_hvg_args

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.var[hvg_column_name] = True
    adata.write_zarr(output_file)
    exit(0)

# make a copy of var for later use since we'll be subsetting adata for HVG selection
var = adata.var.copy()

# filter genes and cells that would break HVG function
batch_mask = _filter_batch(
    adata,
    batch_key=args.get('batch_key'),
    min_cells=min_cells,
)
logging.info(f'Before filtering for HVG selection: {adata.shape}')
adata = adata[batch_mask, adata.var['nonzero_genes']].copy()
logging.info(f'After filtering for HVG selection: {adata.shape}')

# Handle case where filtering resulted in empty data
if adata.n_obs == 0:
    logging.info('No data after filtering, write empty file...')
    var[hvg_column_name] = True
    adata = ad.AnnData(var=var, uns=adata.uns)
    write_zarr_linked(
        adata,
        in_dir=input_file,
        out_dir=output_file,
        files_to_keep=['uns', 'var']
    )
    exit(0)

# Keep batches contiguous for count-based operations.
batch_key = args.get('batch_key') if isinstance(args, dict) else None
if batch_key in adata.obs.columns:
    adata.obs_names_make_unique()
    adata = adata[adata.obs.sort_values(batch_key, kind='stable').index].copy()

# workaround for CxG datasets
feature_col = 'feature_name' if 'feature_name' in var.columns else None

# remove user-specified genes
if remove_genes:
    remove_genes, remove_mask = match_genes(
        adata.var,
        remove_genes,
        column=feature_col,
        return_mask=True
    )
    logging.info(f'Remove {len(remove_genes)} genes (subset data)...')
    adata = adata[:, ~remove_mask].copy()

adata.var[hvg_column_name] = False

if union_over is not None:
    logging.info(f'Compute highly variable genes per {union_over} with args={args}...')
    if isinstance(union_over, str):
        union_over = [union_over]
    union_values = adata.obs[union_over].astype(str) \
        .replace(['nan', 'unknown'], np.nan).dropna() \
        .apply(lambda x: '--'.join(x), axis=1)
    
    # remove groups with fewer than min_cells
    value_counts = union_values.value_counts()
    remove_groups = value_counts[value_counts < min_cells]
    union_values = union_values[~union_values.isin(remove_groups.index)]
    
    if union_values.empty:
        logging.info('No valid union_over groups found; running HVG selection on all cells...')
        union_over = None
    else:
        # set union_over values in adata.obs
        adata.obs['union_over'] = union_values
        adata.obs['union_over'] = adata.obs['union_over'].astype('category')
        
        logging.info(adata.obs['union_over'].value_counts(dropna=False))
        if len(remove_groups) > 0:
            logging.info(f'Removed groups with fewer than {min_cells} cells:\n{remove_groups}')
        
        for group in tqdm(
            adata.obs['union_over'].dropna().unique(),
            miniters=1,
            desc='Computing HVGs per group',
        ):
            _ad = adata[adata.obs['union_over'] == group].copy()
            
            if not dask:
                logging.info('Compute matrix...')
                _ad = dask_compute(_ad)
            
            if use_gpu:
                rsc.get.anndata_to_GPU(_ad)

            sc.pp.highly_variable_genes(_ad, **args)
            
            # get union of gene sets
            adata.var[hvg_column_name] = adata.var[hvg_column_name] | _ad.var['highly_variable']
            del _ad
        
        logging.info(f'Computed {adata.var[hvg_column_name].sum()} highly variable genes.')

if union_over is None:
    # default gene selection
    logging.info(f'Select features for all cells with arguments: {args}...')
    if use_gpu:
        sc.get.anndata_to_GPU(adata)
    
    if not dask:
        logging.info('Compute matrix...')
        adata = dask_compute(adata)
    else: # persist
        adata.X = adata.X.rechunk((200_000, -1)).persist()

    with TqdmCallback(desc=f'Select features with arguments: {args}...', miniters=1):
        sc.pp.highly_variable_genes(adata, **args)
    adata.var[hvg_column_name] = adata.var['highly_variable']

# set extra_hvgs in full dataset
var[hvg_column_name] = False
var.loc[adata.var_names, hvg_column_name] = adata.var[hvg_column_name]

# add user-provided genes
if extra_genes:
    n_genes = len(extra_genes)
    extra_genes = match_genes(var, extra_genes, column=feature_col)
    
    if len(extra_genes) < n_genes:
        logging.warning(f'Only {len(extra_genes)} of {n_genes} user-provided genes found in data...')
    if len(extra_genes) == 0:
        logging.info('No extra user genes added...')
    else:
        logging.info(f'Add {len(extra_genes)} user-provided genes...')
        var.loc[extra_genes, hvg_column_name] = True

# recreate AnnData object for full feature space
adata = ad.AnnData(var=var, uns=adata.uns)

logging.info(f'Write to {output_file}...')
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['uns', 'var']
)
