import numpy as np
import pandas as pd
import typing
import hashlib
import warnings


def get_use_gpu(config):
    use_gpu = bool(config.get('use_gpu', False))
    if isinstance(use_gpu, str):
        use_gpu = use_gpu.lower() == 'true'
    return use_gpu


def remove_outliers(adata, extrema='max', factor=10, rep='X_umap'):
    if factor == 0:
        return adata
    umap = adata.obsm[rep]
    if extrema == 'max':
        abs_values = np.abs(umap.max(axis=1))
    elif extrema == 'min':
        abs_values = np.abs(umap.min(axis=1))
    outlier_mask = abs_values < factor * abs_values.mean()
    return adata[outlier_mask]


def trim_umap(embed, old_key='X_umap', new_key='X_umap', threshold=1e-5):
    """Code provided by @moinfar, modified by @mumichae"""
    x_min, x_max = np.quantile(embed.obsm[old_key][:, 0], (threshold, 1-threshold))
    x_min, x_max = float(x_min), float(x_max)
    y_min, y_max = np.quantile(embed.obsm[old_key][:, 1], (threshold, 1-threshold))
    y_min, y_max = float(y_min), float(y_max)
    embed.obsm[new_key] = np.vstack([embed.obsm[new_key][:, 0].clip(x_min, x_max), embed.obsm[old_key][:, 1].clip(y_min, y_max)]).T


def all_but(_list, is_not):
    return [x for x in _list if x != is_not]


def unique_dataframe(df):
    if df.empty:
        return df
    # hashable_columns = [
    #     col for col in df.columns
    #     if all(isinstance(df[col].iloc[i], typing.Hashable) for i in range(df.shape[0]))
    # ]
    # duplicated = df[hashable_columns].duplicated()
    duplicated = df.astype(str).duplicated()
    return df[~duplicated].reset_index(drop=True)


def expand_dict(_dict):
    """
    Create a cross-product on a dictionary with literals and lists
    :param _dict: dictionary with lists and literals as values
    :return: list of dictionaries with literals as values
    """
    df = pd.DataFrame({k: [v] if isinstance(v, list) else [[v]] for k, v in _dict.items()})
    for col in df.columns:
        df = df.explode(col)
    dict_list = df.apply(lambda row: dict(zip(df.columns, row)), axis=1)

    def remove_chars(s, chars="{} ',"):
        for c in chars:
            s = s.replace(c, '')
        return s

    wildcards = df.apply(
        lambda row: '-'.join(
            [
                remove_chars(f'{col[0]}:{x}') for col, x in zip(df.columns, row)
            ]
        ),
        axis=1
    )
    return zip(wildcards, dict_list)


def expand_dict_and_serialize(_dict: dict, do_not_expand: list = None, hash_exclude: list = None):
    """
    Create a cross-product on a dictionary with literals and lists
    :param _dict: dictionary with lists and literals as values
    :return: list of dictionaries with literals as values
    """
    import jsonpickle
    import hashlib
    
    if do_not_expand is None:
        do_not_expand = []
    
    if hash_exclude is None:
        hash_exclude = []

    df = pd.DataFrame(
        {
            k: [v]
            if isinstance(v, list) and k not in do_not_expand
            else [[v]] for k, v in _dict.items()
        }
    )
    for col in df.columns:
        df = df.explode(col)
    dict_list = df.apply(lambda row: dict(zip(df.columns, row)), axis=1)

    wildcards = []
    for d in dict_list:
        hash_dict = {k: v for k, v in d.items() if k not in hash_exclude}
        wildcards.append(
            hashlib.blake2b(
                jsonpickle.encode(hash_dict).encode('utf-8'),
                digest_size=5,
            ).hexdigest()
        )

    return zip(wildcards, dict_list)


def unlist_dict(dictionary):
    return {
        k: v[0] if isinstance(v, list) and len(v) == 1 else v
        for k, v in dictionary.items()
    }

def unpack_dict_in_df(df, column):
    """
    Given a column in a pandas dataframe containing dictionies, extract these to top level
    
    :param df: pandas dataframe
    :param column: column name containing dictionaries 
    """
    return df.drop(columns=column).assign(**df[column].dropna().apply(pd.Series, dtype=object))


def ifelse(statement, _if, _else):
    if statement:
        return _if
    else:
        return _else


def check_sparse(matrix, sparse_type=None):
    import types
    from anndata.experimental import CSRDataset, CSCDataset
    from scipy.sparse import csr_matrix, csc_matrix
    from sparse import COO
    from dask import array as da
    
    if sparse_type is None:
        sparse_type = (csr_matrix, csc_matrix, CSRDataset, CSCDataset, COO)
    elif not isinstance(sparse_type, tuple):
        sparse_type = (sparse_type,)
    
    # convert to type for functions
    sparse_type = [type(x(0)) if isinstance(x, types.FunctionType) else x for x in sparse_type]
    sparse_type = tuple(sparse_type)
    
    if isinstance(matrix, da.Array):
        return isinstance(matrix._meta, sparse_type)
    return isinstance(matrix, sparse_type)


def check_sparse_equal(a, b):
    from scipy.sparse import csr_matrix
    a = a if check_sparse(a) else csr_matrix(a)
    b = b if check_sparse(b) else csr_matrix(b)
    if a.shape != b.shape:
        warnings.warn(f'Shape mismatch: {a.shape} != {b.shape}')
    return a.shape == b.shape and (a != b).nnz == 0


def ensure_sparse(adata, layers: [str, list] = None, sparse_type=None, **kwargs):
    
    def to_sparse(matrix, sparse_type=None):
        from scipy.sparse import csr_matrix
        from dask import array as da
        
        if matrix is None:
            return matrix
        
        if sparse_type is None:
            sparse_type = csr_matrix

        def convert_block(block):
            if isinstance(sparse_type, type) and isinstance(block, sparse_type):
                return block
            if sparse_type is csr_matrix and hasattr(block, 'tocsr'):
                return block.tocsr()
            return sparse_type(block)

        if check_sparse(matrix, sparse_type):
            return matrix
        elif isinstance(matrix, da.Array):
            return matrix.map_blocks(convert_block, dtype=matrix.dtype)
        return convert_block(matrix)

    return apply_layers(adata, func=to_sparse, layers=layers, sparse_type=sparse_type, **kwargs)


def ensure_dense(adata, layers: [str, list] = None, **kwargs):
    
    def to_dense(matrix):
        from dask import array as da
        
        if isinstance(matrix, da.Array):
            return matrix.map_blocks(np.array)
        if check_sparse(matrix):
            return matrix.toarray()
        return matrix
    
    return apply_layers(adata, func=to_dense, layers=layers, **kwargs)


def dask_compute(adata, layers: [str, list] = None, verbose: bool = True, **kwargs):
    from dask import array as da
    
    def compute_layer(x, persist=False):
        from tqdm.dask import TqdmCallback
        from contextlib import nullcontext
        
        if not isinstance(x, da.Array):
            return x

        if any(dim == 0 for dim in x.shape):
            return np.empty(x.shape, dtype=x.dtype)

        context = TqdmCallback(desc='Dask compute', miniters=10, mininterval=5) if verbose else nullcontext()
        with context:
            if persist:
                x = x.persist()
            x = x.compute()
        return x
    
    return apply_layers(
        adata,
        func=compute_layer,
        layers=layers,
        verbose=verbose,
        **kwargs,
    )


def apply_layers(adata, func, layers:[str, list, bool] = None, verbose: bool = False, **kwargs):
    if layers is None or layers is True:
        layers = ['X', 'raw'] + list(adata.layers.keys())
    elif isinstance(layers, str):
        layers = [layers]
    elif layers is False:
        return adata
    
    for layer in layers:
        if verbose:
            print(f'Apply function {func.__name__} to {layer}...', flush=True)
        if layer == 'X':
            adata.X = func(adata.X, **kwargs)
        elif layer in adata.layers:
            adata.layers[layer] = func(adata.layers[layer], **kwargs)
        elif layer in adata.obsm:
            adata.obsm[layer] = func(adata.obsm[layer], **kwargs)
        elif layer == 'raw':
            if adata.raw is None:
                continue
            raw_x = func(adata.raw.X, **kwargs)
            if raw_x is adata.raw.X:
                continue
            adata_raw = adata.raw.to_adata()
            adata_raw.X = raw_x
            adata.raw = adata_raw
        elif verbose:
            print(f'Layer {layer} not found, skipping...', flush=True)
    return adata


def merge(dfs: list, verbose: bool = True, **kwargs):
    """
    Merge list of dataframes
    :param dfs: list of dataframes
    :param kwargs: arguments passed to pd.merge
    :return: merged dataframe
    """
    from functools import reduce

    merged_df = reduce(
        lambda x, y: pd.merge(x, y, **kwargs),
        dfs
    )
    if verbose:
        print(merged_df)
    return merged_df


def create_hash(string: str, digest_size: int = 5):
    string = string.encode('utf-8')
    return hashlib.blake2b(string, digest_size=digest_size).hexdigest()
