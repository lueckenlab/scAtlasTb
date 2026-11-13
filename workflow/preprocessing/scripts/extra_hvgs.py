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

def _write_empty_result(adata, output_file, hvg_column_name):
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

def _compute_group_mapping(obs, union_keys):
    keys = [union_keys] if isinstance(union_keys, str) else union_keys
    return (
        obs[keys]
        .astype(str)
        .replace(['nan', 'unknown'], np.nan)
        .dropna()
        .agg('--'.join, axis=1)
    )

def _build_ordered_index(adata, valid_mapping, batch_key, min_cells):
    """Filter and sort cells for HVG selection.

    Orders by group label and/or `batch_key`.
    """
    logging.info(f'Before filtering for HVG selection: {adata.shape}')
    batch_mask = _filter_batch(adata, batch_key=batch_key, min_cells=min_cells)

    if batch_key:
        assert batch_key in adata.obs.columns, f"Column '{batch_key}' not found in adata.obs"

    if valid_mapping is not None:
        sort_cols = ['_group_label', batch_key] if batch_key else ['_group_label']
        adata.obs['_group_label'] = valid_mapping
        return adata.obs.loc[valid_mapping.index].sort_values(sort_cols, kind='stable').index.tolist()

    filtered_idx = adata.obs_names[batch_mask].tolist()
    if batch_key:
        return adata.obs.loc[filtered_idx].sort_values(batch_key, kind='stable').index.tolist()
    return filtered_idx


input_file = snakemake.input[0]
output_file = snakemake.output[0]
args = snakemake.params.get('args', {})
batch_key = args.get('batch_key') if isinstance(args, dict) else None
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
layer = 'layers/counts' if args.get('flavor') == 'seurat_v3' else 'layers/normcounts'
adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var',
    uns='uns',
    backed=True,
    dask=True,
)
adata.obs_names_make_unique()
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

    # Normalize union_over to a list and map each cell to its group label
    union_keys = [union_over] if isinstance(union_over, str) else union_over
    group_mapping = _compute_group_mapping(adata.obs, union_keys)

    # Filter out groups below the minimum cell threshold (operate on filtered obs)
    group_counts = group_mapping.value_counts()
    valid_groups = group_counts[group_counts >= min_cells]
    logging.info(f'Union over {union_keys}:\n{group_counts}')

    n_removed = group_mapping.size - valid_groups.sum()
    if n_removed > 0:
        logging.info(f'Removed {n_removed} cells from groups with < {min_cells} cells')

    if valid_groups.empty:
        logging.info('No valid union_over groups found; running HVG selection on all cells...')
        union_over = None
    else:
        valid_mapping = group_mapping[group_mapping.isin(valid_groups.index)]

        logging.info(f'Before filtering for HVG selection: {adata.shape}')
        # build ordered index and single-subset adata (rows + nonzero genes)
        new_index = _build_ordered_index(adata, valid_mapping, batch_key, min_cells)
        adata = adata[new_index, adata.var['nonzero_genes']].copy()
        logging.info(f'After filtering for HVG selection: {adata.shape}')

        if adata.n_obs == 0:
            _write_empty_result(adata, output_file, hvg_column_name)

        # Compute HVGs per group and accumulate results via union
        adata.X = adata.X.rechunk((200_000, -1))
        for group, cell_idx in tqdm(valid_mapping.groupby(valid_mapping), miniters=1, desc='Computing HVGs per group'):
            _ad = adata[cell_idx.index].copy()

            if not dask:
                _ad = dask_compute(_ad, layers='X')
            if use_gpu:
                rsc.get.anndata_to_GPU(_ad)

            sc.pp.highly_variable_genes(_ad, **args)
            adata.var[hvg_column_name] |= _ad.var['highly_variable']
            del _ad

        logging.info(f'Total highly variable genes across all groups: {adata.var[hvg_column_name].sum()}')

else:
    logging.info(f'Select features for all cells with arguments: {args}...')

    logging.info(f'Before filtering for HVG selection: {adata.shape}')
    # build ordered index and subset using helper (handles batch mask internally)
    new_index = _build_ordered_index(adata, None, batch_key, min_cells)
    adata = adata[new_index, adata.var['nonzero_genes']].copy()
    logging.info(f'After filtering for HVG selection: {adata.shape}')
    
    if adata.n_obs == 0:
        _write_empty_result(adata, output_file, hvg_column_name)

    if use_gpu:
        sc.get.anndata_to_GPU(adata)

    adata.X = adata.X.rechunk((200_000, -1))
    if not dask:
        logging.info('Compute matrix...')
        adata = dask_compute(adata, layers='X')

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
