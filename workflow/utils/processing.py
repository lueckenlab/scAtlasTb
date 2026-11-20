import traceback
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse
from tqdm import tqdm
try:
    import subprocess
    if subprocess.run('nvidia-smi', shell=True).returncode != 0:
        logging.info('No GPU found...')
        raise ImportError()
    import rapids_singlecell as sc
    import cupy as cp
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    rmm.reinitialize(
        managed_memory=True,
        pool_allocator=False,
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)
    logging.info('Using rapids_singlecell...')
    USE_GPU = True
except ImportError as e:
    import scanpy as sc
    logging.info('Importing rapids failed, using scanpy...')
    USE_GPU = False

from .assertions import assert_neighbors
from .io import to_memory
from .misc import dask_compute
from .accessors import _filter_genes


def compute_neighbors(adata, output_type=None, force=False, check_n_neighbors=False, **kwargs):
    """ Compute kNN graph based on output type.

    :param adata: integrated anndata object
    :param output_type: string of output type
    :param force: force re-computation of kNN graph
    :param kwargs: additional arguments for sc.pp.neighbors
    :return: anndata with kNN graph based on output type
    """
    
    if not output_type:
        output_type = 'full'

    try:
        logging.info(adata.__str__())
        assert not force, 'force neighbors computation'
        assert_neighbors(adata, check_n_neighbors=check_n_neighbors)
        logging.info(f'kNN graph already computed for {output_type}. Using pre-computed {output_type} kNN graph')
        return
    except AssertionError as e:
        logging.info(e.__str__())
        logging.info(f'Compute kNN graph for {output_type}...')

    if output_type == 'knn':
        assert_neighbors(adata, check_params=False)
        if 'params' not in adata.uns['neighbors']:
            adata.uns['neighbors']['params'] = adata.uns['neighbors'].get('params', {})
        adata.uns['neighbors']['params'] |= dict(use_rep='X')
    
    elif output_type == 'embed':
        assert 'X_emb' in adata.obsm, 'Embedding key "X_emb" not found'
        kwargs |= dict(use_rep='X_emb')
        sc.pp.neighbors(adata, **kwargs)
    
    elif output_type == 'full':
        if 'X_emb' not in adata.obsm:
            logging.info('Load corrected counts...')
            adata_tmp = adata
            if 'highly_variable' in adata_tmp.var.columns:
                adata_tmp = adata_tmp[:, adata_tmp.var['highly_variable']]
            adata_tmp = dask_compute(adata_tmp.copy(), layers='X')
            
            logging.info('Compute PCA on corrected counts...')
            sc.pp.pca(adata_tmp)
            adata.obsm['X_emb'] = adata_tmp.obsm['X_pca']
            del adata_tmp
        kwargs |= dict(use_rep='X_emb')
        logging.info(f'Neighbors with {kwargs}...')
        sc.pp.neighbors(adata, **kwargs)
    
    else:
        raise ValueError(f'Invalid output type {output_type}')
    
    assert_neighbors(adata, check_n_neighbors=False)


def filter_genes(adata, batch_key=None, compute=True, **kwargs):
    """
    Filter anndata based on .X matrix
    """
    from dask import array as da
    import dask.config
    
    # filter cells from batches with too few cells
    if batch_key is not None:
        cells_per_batch = adata.obs[batch_key].value_counts()
        if cells_per_batch.min() < 2:
            adata = adata[adata.obs[batch_key].isin(cells_per_batch[cells_per_batch > 1].index)]
        
    # apply filter
    gene_mask = _filter_genes(adata, **kwargs)
    if any(gene_mask == False):
        logging.info(f'Subset to {sum(gene_mask)}/{adata.n_vars} filtered genes...')
        adata = adata[:, gene_mask]
    
    if adata.is_view:
        adata = adata.copy()    
    
    # convert dask array to csr matrix
    if compute and isinstance(adata.X, da.Array):
        with dask.config.set(**{'array.slicing.split_large_chunks': True}):
            adata = dask_compute(adata, layers='X')

    return adata


def get_pseudobulks(adata, group_key, agg='sum', dtype='float32', sep='--', group_cols=None):
    from dask import array as da
    
    def aggregate(x, agg):
        if agg == 'sum':
            return x.sum(0)
        elif agg == 'mean':
            return x.mean(0)
        else:
            raise ValueError(f'invalid aggregation method "{agg}"')
    
    def _get_pseudobulk_matrix(adata, group_key, agg, dtype=None):
        dtype = dtype or np.float32
        X = adata.X
        value_counts = adata.obs[group_key].value_counts()
        value_counts = value_counts[value_counts >= 2]
        groups = value_counts.index

        print(f'Aggregate {len(groups)} pseudobulks...', flush=True)

        if isinstance(X, da.Array):
            from dask import config as dask_config
            from tqdm.dask import TqdmCallback

            df = adata.obs[[group_key]].reset_index(drop=True).query(f'`{group_key}` in @groups')
            df[group_key] = pd.Categorical(df[group_key], categories=groups, ordered=True)

            sorted_df = df.sort_values(by=group_key)

            # Derive chunk sizes from the sorted order — NOT from value_counts() which is
            # sorted by count descending and would silently assign cells to wrong groups
            sorted_group_sizes = (
                sorted_df[group_key]
                .value_counts()
                .reindex(sorted_df[group_key].cat.categories)
            )

            print(f'Sort and rechunk dask array by "{group_key}"...', flush=True)
            with dask_config.set(**{'array.slicing.split_large_chunks': False}):
                X = X[sorted_df.index.values]
                X = X.rechunk((tuple(sorted_group_sizes.values), -1))
                pseudobulks = X.map_blocks(lambda x: aggregate(x, agg), dtype=dtype)

            print(f'Compute aggregation...', flush=True)
            with TqdmCallback(
                bar_format='{percentage:3.0f}% |{bar}| {elapsed}<{remaining}\n',
                mininterval=10,
                delay=10,
            ):
                pseudobulks = pseudobulks.compute()

        elif isinstance(X, (scipy.sparse.spmatrix, np.ndarray)):
            import scanpy as sc
            pbulk_adata = sc.get.aggregate(adata, by=group_key, func=[agg])
            pbulk_adata = pbulk_adata[pbulk_adata.obs[group_key].isin(groups)].copy()
            pseudobulks = pbulk_adata.layers[agg]
            groups = pbulk_adata.obs_names

        else:
            raise ValueError(f'invalid type "{type(X)}"')

        pseudobulks = scipy.sparse.csr_matrix(pseudobulks, dtype='float32')
        return pseudobulks, groups

    if isinstance(group_key, list):
        group_keys = group_key
        group_key = sep.join(group_keys)
        adata.obs[group_key] = adata.obs[group_keys].astype(str).apply(
            lambda x: sep.join(x),
            axis=1
        )

    # call pseudobulk function
    pbulks, groups = _get_pseudobulk_matrix(adata, group_key=group_key, agg=agg, dtype=dtype)

    # aggregate metadata
    obs = aggregate_obs(adata, group_key, groups, include_cols=group_cols)
    print(obs, flush=True)

    return ad.AnnData(pbulks, obs=obs, var=adata.var)


def aggregate_obs(
    adata: ad.AnnData,
    group_key: str,
    groups: list,
    include_cols: list = None,
):

    print(f'Aggregate metadata by {group_key}...', flush=True)

    def smart_categorical_mode(x):
        n_vals = len(x)
        n_unique = x.nunique()
        
        # If all values are unique, just return the first one
        if n_unique == n_vals:
            return x.iloc[0]
        
        if n_unique / n_vals < 0.8:  # threshold can be tuned
            codes = x.cat.codes.values
            # Handle missing values (-1 codes)
            codes = codes[codes >= 0]
            if len(codes) == 0:
                return None
            counts = np.bincount(codes)
            mode_code = np.argmax(counts)
            return x.cat.categories[mode_code]
        
        # If very high uniqueness ratio, use value_counts (avoid bincount)
        return x.value_counts().index[0]
    

    df = adata.obs
    if include_cols is not None:
        include_cols = list(set(include_cols+[group_key]))
        df = df[include_cols].copy()
    bool_columns = df.select_dtypes('bool').columns.tolist()
    num_columns = df.select_dtypes('number').columns.tolist()
    cat_columns = df.select_dtypes(exclude=['number', 'bool']).columns
    cat_columns = [col for col in cat_columns if col != group_key]
    df[cat_columns] = df[cat_columns].astype('category')

    if bool_columns or num_columns or cat_columns:
        g = df.groupby(group_key)
        df = pd.concat([
            g[bool_columns+num_columns].mean(),
            g[cat_columns].agg(smart_categorical_mode),
        ], axis=1)
    else:
        df = df.groupby(group_key).first()
    df[bool_columns] = df[bool_columns] > 0.5 # mode for bool
    df['n_agg'] = adata.obs.groupby(group_key).size()

    print(f'Set aggregated metadata...', flush=True)
    obs = pd.DataFrame(groups, columns=[group_key]).merge(df, on=group_key, how='left')
    obs.index = obs[group_key].values

    print(f'Convert metadata to categorical...', flush=True)
    # convert category columns to string
    for col in obs.columns:
        if obs[col].dtype.name == 'category':
            obs[col] = obs[col].astype(str).astype('category')

    return obs
