import pandas as pd
import numpy as np
from tqdm import tqdm
from pprint import pformat

LARGE_SLOTS = ['layers', 'obsm', 'obsp']


def check_same_index(dfs: dict):
    iterator = iter(dfs.items())
    key1, df1 = next(iterator)
    
    for key, df in iterator:
        if key == key1:
            continue
        diff = df.index.difference(df1.index)
        assert len(diff) == 0, \
            f'Index must be the same\n {len(diff)} differing indices \
            when comparing {key} against {key1}\n' 


def check_columns_equal(dfs: dict, col: str):
    def _check_columns_equal(s1, s2, col):
        assert s1.index.equals(s2.index), \
            f'Indices do not match for column comparison: {col}'
        # Fast path: check if both Series are the same object or have the same underlying data
        if s1 is s2:
            return True
        if s1._values is s2._values:
            return True
        # Use numpy for fast comparison if dtypes are compatible
        if s1.dtype == s2.dtype and pd.api.types.is_numeric_dtype(s1.dtype):
            return np.array_equal(s1.values, s2.values, equal_nan=True)
        # Fallback to pandas equals
        return s1.equals(s2)

    iterator = iter(dfs.values())
    df1 = next(iterator)
    if all(_check_columns_equal(df[col], df1[col], col) for df in iterator):
        return []
    return [col]


def get_same_columns(dfs: dict, n_threads=1):
    check_same_index(dfs)

    # collect columns that exist in all files
    iterator = iter(dfs.values())
    df1 = next(iterator)
    same_columns = set(
        df1.columns.tolist()
    ).intersection(
        *[df.columns.tolist() for df in iterator]
    )
    
    # determine which columns do not have the same values across all files
    cols_to_remove = []
    for col in tqdm(same_columns, desc='Check for equal columns'):
        cols_to_remove.extend(check_columns_equal(dfs, col))
    
    print(f'Shared columns that are not the same across datasets:\n{pformat(cols_to_remove)}', flush=True)
    return [col for col in same_columns if col not in cols_to_remove]


def analyze_slot_columns(dfs_dict, n_threads=1):
    same = set(get_same_columns(dfs_dict, n_threads=n_threads))
    multicol = {
        col
        for col, count in (
            pd.Series(
                [col for df in dfs_dict.values() for col in df.columns]
            )
            .loc[lambda s: ~s.isin(same)]
            .value_counts()
            .items()
        )
        if count > 1
    }
    return {'same': same, 'multicol': multicol}


def merge_df(
    df_current,
    file_id,
    df_previous,
    column_suffix_map=None,
    sep='_',
):
    """Merge dataframes, adding suffix only to columns that appear in multiple files.
    
    Parameters
    ----------
    column_suffix_map : dict, optional
        'same': set of columns identical across all files (excluded)
        'multicol': set of columns needing suffix
    """
    if column_suffix_map is None:
        column_suffix_map = {'same': set(), 'multicol': set()}
    
    same_cols = column_suffix_map.get('same', set())
    multicol_cols = column_suffix_map.get('multicol', set())
    
    if df_previous is None:
        return df_current.rename(
            columns={
                col: f'{col}{sep}{file_id}'
                for col in df_current.columns.tolist()
                if col not in same_cols and col in multicol_cols
            }
        )
    
    idx_diff = df_current.index.difference(df_previous.index).sort_values()
    n_diff = len(idx_diff)
    assert n_diff == 0, \
        f'Index must be the same\n {n_diff} differing indices: {idx_diff}\n' \
        f'current index: {df_current.index.sort_values()} \nprevious index: {df_previous.index.sort_values()}'
    
    unique_columns = [col for col in df_current.columns if col not in same_cols]
    df_current = df_current[unique_columns].rename(
        columns={
            col: f'{col}{sep}{file_id}' 
            for col in unique_columns 
            if col in multicol_cols
        }
    )
    df = pd.concat([df_previous, df_current], axis=1)
    return df


def _reload_large_slots(adata, file_path, large_merge_slots):
    """Reload large slots from disk to avoid view-based chained operations."""
    from utils.io import read_anndata
    adata_tmp = read_anndata(
        file_path,
        **{slot: slot for slot in large_merge_slots},
        dask=True,
        backed=True,
        verbose=False
    )
    for slot in large_merge_slots:
        setattr(adata, f'_{slot}', getattr(adata_tmp, f'_{slot}'))


def align_file_indices(adatas, files, merge_slots):
    """
    Align AnnData objects to match reference: same obs_names, intersection of var_names.
    
    Reorders/subsets all adatas to match the first one's obs_names and the intersection
    of all var_names. Reloads large slots (layers, obsm, obsp) from disk before subsetting
    to avoid expensive view chains. Converts resulting views to copies.
    
    Parameters
    ----------
    adatas : dict[str, AnnData]
        Mapping from file_id to AnnData. First item is the reference.
    files : dict[str, str]
        Mapping from file_id to file path on disk.
    merge_slots : Iterable[str]
        Slot names to merge (e.g. {'X', 'layers', 'obsm'}).
    
    Returns
    -------
    dict[str, list[str]]
        Maps large slot names to file_ids where data was reloaded from disk.
    
    Raises
    ------
    ValueError
        If obs_names sets differ across files (permutation OK, different content not OK).
    """
    import logging
    from utils.io import read_anndata
    
    adata_iter = iter(adatas.items())
    first_key, first_adata = next(adata_iter)
    
    first_obs_names = first_adata.obs_names
    first_var_names = first_adata.var_names
    large_merge_slots = [slot for slot in LARGE_SLOTS if slot in merge_slots]
    file_copy_map = {slot: [] for slot in large_merge_slots}
    
    # Compute intersection of var_names
    var_intersection = set(first_adata.var_names)
    remaining_adatas = []
    for file_id, adata in adata_iter:
        var_intersection &= set(adata.var_names)
        remaining_adatas.append((file_id, adata))
    
    # Subset first adata if var_names differ
    if var_intersection != set(first_var_names):
        logging.info(f"var_names differ. Keeping {len(var_intersection)} common vars.")
        first_var_names = first_adata.var_names[first_adata.var_names.isin(var_intersection)]
        if large_merge_slots:
            _reload_large_slots(first_adata, files[first_key], large_merge_slots)
            for slot in large_merge_slots:
                file_copy_map[slot].append(first_key)
        adatas[first_key] = first_adata[:, first_var_names]
    
    # Align remaining adatas
    for file_id, adata in remaining_adatas:
        file_path = files[file_id]
        modified = False
        
        # Check/reorder obs_names
        if not adata.obs_names.equals(first_obs_names):
            if set(adata.obs_names) != set(first_obs_names):
                missing = set(first_obs_names) - set(adata.obs_names)
                extra = set(adata.obs_names) - set(first_obs_names)
                n_missing_examples = min(10, len(missing))
                n_extra_examples = min(10, len(extra))
                raise ValueError(
                    f"obs_names mismatch in '{file_id}':\n"
                    f"  Missing: {len(missing)} (e.g. {sorted(missing)[:n_missing_examples]})\n"
                    f"  Extra: {len(extra)} (e.g. {sorted(extra)[:n_extra_examples]})\n"
                )
            if large_merge_slots:
                _reload_large_slots(adata, file_path, large_merge_slots)
                for slot in large_merge_slots:
                    file_copy_map[slot].append(file_id)
            logging.info(f"Reordering obs_names: {file_id}")
            adata = adata[first_obs_names, :]
            modified = True
        
        # Subset var_names if needed
        if not adata.var_names.equals(first_var_names):
            if not modified and large_merge_slots:
                _reload_large_slots(adata, file_path, large_merge_slots)
                for slot in large_merge_slots:
                    if file_id not in file_copy_map[slot]:
                        file_copy_map[slot].append(file_id)
            logging.info(f"Subsetting var_names: {file_id}")
            adata = adata[:, first_var_names]
            modified = True
        
        adatas[file_id] = adata if not modified or not adata.is_view else adata.copy()
    
    return file_copy_map

