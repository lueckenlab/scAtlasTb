import logging
from pathlib import Path
logging.basicConfig(level=logging.INFO)
import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
import gc
import faulthandler
faulthandler.enable()
from scipy.sparse import issparse
from tqdm import tqdm
import pandas as pd
import yaml
import scanpy as sc
from anndata.experimental import AnnCollection
from anndata import AnnData
from dask import array as da
from dask import config as da_config
# from dask.diagnostics import ProgressBar

da_config.set(
    **{
        'num_workers': snakemake.threads,
        'array.slicing.split_large_chunks': False,
    }
)
logging.info(f"Dask using {da_config.get('num_workers')} workers")

from utils.io import read_anndata, link_zarr, write_zarr
from utils.misc import apply_layers, ensure_sparse


def read_adata(
    file,
    file_id=None,
    backed=False,
    dask=False,
    stride=100_000,
    chunks=(-1, -1),
    log=True,
    **slots
):
    if file_id is None:
        file_id = file
    if log:
        logging.info(f'Read {file}...')
    try:
        adata = read_anndata(
            file,
            backed=backed,
            dask=dask,
            stride=stride,
            chunks=chunks,
            **slots,
            verbose=False,
        )
        adata.obs['file_id'] = file_id
        adata = ensure_sparse(adata)
    except Exception as e:
        raise Exception(f'Error reading {file}') from e
    
    return adata


def get_common_obs_columns(adatas):
    if not adatas:
        return []

    common = set(adatas[0].obs.columns)
    for adata in adatas[1:]:
        common &= set(adata.obs.columns)

    return [column for column in adatas[0].obs.columns if column in common]


def keep_only_obs_columns(adata, columns):
    if columns:
        adata.obs = adata.obs.loc[:, columns].copy()
    else:
        adata.obs = pd.DataFrame(index=adata.obs_names)
    return adata


def remove_slots(adata):
    for slot in ['X', 'layers', 'obsm', 'obsp', 'varm', 'varp']:
        if hasattr(adata, slot):
            delattr(adata, slot)


def persist_batch_adata(adata):
    """Persist dask arrays in adata to reset task graph."""
    dask_layers = ['X', 'raw'] + list(adata.layers.keys())
    return apply_layers(
        adata,
        func=lambda x: x if not isinstance(x, da.Array) else x.persist(),
        layers=dask_layers,
    )


def calculate_batch_boundaries(adatas, batch_cell_size):
    """
    Calculate batch boundaries based on target cell count.
    
    Returns list of (start_idx, end_idx, total_cells) tuples.
    """
    batches = []
    batch_start = 0
    cells_in_batch = 0
    
    for idx, adata in enumerate(adatas):
        cells_in_batch += adata.n_obs
        is_last_file = (idx == len(adatas) - 1)
        
        if cells_in_batch >= batch_cell_size or is_last_file:
            batches.append((batch_start, idx + 1, cells_in_batch))
            batch_start = idx + 1
            cells_in_batch = 0
    
    return batches


def reduce_pairwise_with_persist(adatas, merge_strategy):
    """
    Reduce a list of AnnData objects in pairwise rounds.

    Algorithm
    ---------
    1. Start with all batch-level merged objects.
    2. In each round, merge neighbors in pairs: (0,1), (2,3), ...
    3. Persist each newly merged pair immediately to reset Dask graph depth.
    4. If a round has an odd item, carry it unchanged to the next round.
    5. Repeat until one AnnData remains.
    6. Persist once more before returning to ensure a compact final graph.

    Why this helps
    --------------
    Pairwise tree reduction avoids building one long concatenation chain,
    which keeps graph depth and scheduler overhead lower for many inputs.
    """
    level = adatas
    round_idx = 0

    while len(level) > 1:
        round_idx += 1
        next_level = []
        n_pairs = (len(level) + 1) // 2
        logging.info(f'Reduction round {round_idx}: {len(level)} nodes -> {n_pairs}')

        for i in tqdm(range(0, len(level), 2), unit='pair'):
            pair = level[i:i + 2]

            # Odd-node carry: keep as-is for the next round.
            if len(pair) == 1:
                next_level.append(pair[0])
                continue

            merged = sc.concat(pair, join=merge_strategy, axis=0)
            merged = persist_batch_adata(merged)
            next_level.append(merged)

        level = next_level

    return persist_batch_adata(level[0])


def merge_adatas_with_batching(adatas, merge_strategy, batch_cell_size):
    """
    Merge adatas in batches by cell count and persist after each batch to reduce dask layer depth.
    
    Parameters
    ----------
    adatas : list of AnnData
        List of AnnData objects to merge
    merge_strategy : str
        Join strategy ('inner', 'outer')
    batch_cell_size : int
        Target number of cells per batch
        
    Returns
    -------
    AnnData
        Merged AnnData object with persisted dask arrays
    """
    if len(adatas) <= 1:
        return adatas[0] if adatas else None
    
    # Calculate batch boundaries
    batch_boundaries = calculate_batch_boundaries(adatas, batch_cell_size)
    total_cells = sum(adata.n_obs for adata in adatas)
    
    logging.info(
        f'Merging {len(adatas)} files ({total_cells:,} cells) '
        f'in {len(batch_boundaries)} batches of ~{batch_cell_size:,} cells'
    )
    
    # Phase 1: merge file-level inputs into cell-count-based batches.
    # Each batch merge is persisted immediately to keep each batch graph shallow.
    merged_batches = []

    for batch_num, (start, end, n_cells) in enumerate(
        tqdm(batch_boundaries, desc='Merging batches', unit='batch'), start=1
    ):
        batch_slice = adatas[start:end]
        merged_batch = sc.concat(batch_slice, join=merge_strategy, axis=0)
        merged_batch = persist_batch_adata(merged_batch)
        merged_batches.append(merged_batch)
        
        logging.info(
            f'[Batch {batch_num}/{len(batch_boundaries)}] '
            f'{len(batch_slice)} files, {n_cells:,} cells -> {merged_batch.shape}'
        )

    # Phase 2: reduce persisted batches with pairwise rounds.
    # This keeps final-merge graph growth controlled for many batches.
    return reduce_pairwise_with_persist(merged_batches, merge_strategy)


dataset = snakemake.wildcards.dataset
files = snakemake.input
out_file = snakemake.output.zarr

merge_strategy = snakemake.params.get('merge_strategy', 'inner')
keep_all_columns = snakemake.params.get('keep_all_columns', False)
allow_duplicate_obs = snakemake.params.get('allow_duplicate_obs', False)
allow_duplicate_vars = snakemake.params.get('allow_duplicate_vars', False)
new_indices = snakemake.params.get('new_indices', False)
backed = snakemake.params.get('backed', False)
dask = snakemake.params.get('dask', False)
persist = snakemake.params.get('persist', False)
stride = snakemake.params.get('stride', 100_000)
slots = snakemake.params.get('slots', {})
if slots is None:
    slots = {}
    
logging.info(f'Read slots: {slots}')

if len(files) == 1:
    link_zarr(in_dir=files[0], out_dir=out_file)
    exit(0)

# subset to non-empty datasets
def get_shape(file):
    if file.endswith('.zarr'):
        zattr_path = Path(file) / slots.get('X', 'X') / '.zattrs'
        if zattr_path.exists():
            shape = read_anndata(file, X='X', dask=True, backed=True, verbose=False).shape
    # # doesn't work when using file symlinking
    #         with open(zattr_path, 'r') as f:
    #             zattrs = yaml.safe_load(f) or {}
    #         if 'shape' in zattrs:
    #             return list(zattrs['shape'])
    else:
        shape = read_anndata(file, obs='obs', var='var', verbose=False).shape
    return list(shape)  # tuples fail when saving zarr

original_shapes = {
    file_id: get_shape(file)
    for file_id, file
    in zip(files.keys(), files)
}

files = {
    file_id: file
    for file_id, file
    in zip(files.keys(), files)
    if original_shapes[file_id][0] > 0  # check that n_obs > 0
}

if len(files) == 0:
    logging.info('All adatas are empty, skip concatenation...')
    AnnData().write_zarr(out_file)
    exit(0)


if dask:
    logging.info('Loading all files as backed Dask AnnDatas...')
    adatas = [
        read_adata(
            file_path,
            file_id=file_id,
            backed=True,
            dask=True,
            stride=stride,
            **slots,
            log=False,
        )
        for file_id, file_path in tqdm(files.items(), desc='Read files', miniters=1)
    ]
    obs_col_check = list(zip(files.keys(), (adata.obs.shape[1] for adata in adatas)))
    insufficient_files = [fid for fid, n_cols in obs_col_check if n_cols <= 20]
    if insufficient_files:
        logging.warning(f'Files with insufficient obs columns (<=20): {insufficient_files}')
        for fid, n_cols in obs_col_check:
            logging.info(f'  {fid}: {n_cols} obs columns')
    assert all(n_cols > 20 for _, n_cols in obs_col_check), 'not all adatas have more than 20 obs columns'

    common_obs_columns = [] if keep_all_columns else get_common_obs_columns(adatas)
    print(common_obs_columns)
    assert len(common_obs_columns) > 20, 'not enough common obs columns found across adatas, cannot merge with current settings'
    for adata_item in adatas:
        keep_only_obs_columns(adata_item, common_obs_columns)
    
    # Merge adatas with batching to reduce dask layer depth
    if persist and len(adatas) > 2:
        adata = merge_adatas_with_batching(adatas, merge_strategy, batch_cell_size=stride)
    else:
        logging.info('Concatenating adatas with lazy counts...')
        adata = sc.concat(adatas, join=merge_strategy, axis=0)
        if persist:
            logging.info('Persisting concatenated adata to reset Dask graph...')
            adata = persist_batch_adata(adata)
    
    logging.info(adata.__str__())
    logging.info(f'raw present: {adata.raw is not None}, layers: {adata.layers.keys()}')
    
    if persist and len(adatas) == 1:
        logging.info('Persisting single file (batching not applicable)...')
        adata = persist_batch_adata(adata)

    
elif backed:
    logging.info('Read all files in backed mode...')
    adatas = [
        read_adata(file_path, file_id=file_id, backed=backed, dask=dask)
        for file_id, file_path in files.items()
    ]
    common_obs_columns = [] if keep_all_columns else get_common_obs_columns(adatas)
    for adata_item in adatas:
        keep_only_obs_columns(adata_item, common_obs_columns)
    dc = AnnCollection(
        adatas,
        join_obs='outer',
        join_obsm=None,
        join_vars=merge_strategy,
        indices_strict=not backed,
    )
    logging.info(dc.__str__())

    logging.info('Subset AnnDataCollection, returning a View...')
    adata = dc[:].to_adata()
    assert adata.X is not None
    
    for _ad in adatas:
        remove_slots(_ad)
        gc.collect()

else:

    adatas = []
    adata = None
    for file_id, file_path in tqdm(files.items()):
        logging.info(f'Read {file_path}...')
        _ad = read_adata(
            file_path,
            file_id=file_id,
            **slots,
            backed=backed,
            dask=dask
        )
        logging.info(f'{file_id} shape: {_ad.shape}')
        
        if adata is None:
            adata = _ad
            if keep_all_columns:
                keep_only_obs_columns(adata, [])
            continue

        common_obs_columns = [] if keep_all_columns else [
            column for column in adata.obs.columns if column in _ad.obs.columns
        ]
        keep_only_obs_columns(adata, common_obs_columns)
        keep_only_obs_columns(_ad, common_obs_columns)
        
        # merge adata
        adata = sc.concat([adata, _ad], join=merge_strategy)
        logging.info(f'after merge:\n{adata}')

        remove_slots(_ad)
        adatas.append(_ad)
        gc.collect()

if not allow_duplicate_obs:
    # Remove duplicate obs_names, keeping the first occurrence
    duplicates = adata.obs_names.duplicated(keep='first')
    if duplicates.any():
        logging.info(f'Removing {duplicates.sum()} duplicate obs_names...')
        adata = adata[~duplicates].copy()


# Assert no duplicate var_names in the merged object
if not allow_duplicate_vars:
    assert not adata.var_names.duplicated().any(), "Duplicate var_names found in merged AnnData object"

# merge lost annotations
logging.info('Add gene info...')
var_dfs = (
    read_anndata(file, var='var', verbose=False).var
    for file in files.values()
)
for var in var_dfs:
    # get intersection of var_names
    var_names = list(set(adata.var_names).intersection(set(var.index)))

    # conver categorical columns to str
    categorical_columns = var.select_dtypes(include='category').columns
    var[categorical_columns] = var[categorical_columns].astype(str)

    # add new columns to adata.var
    adata.var.loc[var_names, var.columns] = var.loc[var_names, :]

# fix dtypes
adata.var = adata.var.infer_objects()
logging.info(adata.var)

if keep_all_columns:
    logging.info('Merging observation columns from all datasets, all columns will be kept...')
    obs_dfs = (
        read_anndata(file, obs='obs', verbose=False).obs
        for file in files.values()
    )
    merged_obs = pd.concat(obs_dfs, axis=0, join='outer', ignore_index=False)
    
    if not allow_duplicate_obs:
        # remove duplicates the same way as above
        merged_obs = merged_obs[~merged_obs.index.duplicated(keep='first')]
        
        if not merged_obs.index.equals(adata.obs_names):
            # Ensure indices match
            extra = merged_obs.index.difference(adata.obs_names)
            if len(extra):
                logging.info(f'Dropping {len(extra)} extra observation indices not in adata.obs')
                merged_obs = merged_obs.drop(extra)
            
            missing = adata.obs_names.difference(merged_obs.index)
            if len(missing):
                raise ValueError(f'Merged observation table is missing {len(missing)} indices from adata.obs')
            
            # Reorder to match adata.obs_names
            merged_obs = merged_obs.loc[adata.obs_names]
    
    adata.obs = adata.obs.combine_first(merged_obs)

# Store information about the files that went into the merge
adata.uns['merge'] = {
    'files': list(files.values()),
    'file_ids': list(files.keys()),
    'n_files': len(files),
    'merge_strategy': merge_strategy,
    'allow_duplicate_obs': allow_duplicate_obs,
    'allow_duplicate_vars': allow_duplicate_vars,
    'reindexed': new_indices,
    'original_shapes': original_shapes,
}

# set new indices
if new_indices:
    adata.obs[f'obs_names_before_{dataset}'] = adata.obs_names
    adata.obs_names = dataset + '-' + adata.obs.reset_index(drop=True).index.astype(str)

# add uns data
logging.info(adata.__str__())

logging.info(f'Write to {out_file}...')
if isinstance(adata.X, da.Array):
    graph = adata.X.__dask_graph__()
    n_tasks = len(graph)
    if hasattr(graph, "layers"):
        n_layers = len(graph.layers)
        logging.info(f'Dask graph: {n_layers} layers, {n_tasks} tasks')
    else:
        logging.info(f'Dask graph: {n_tasks} tasks')
write_zarr(adata, out_file)
