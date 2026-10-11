import warnings
import numpy as np
import anndata as ad
from dask import array as da
import logging
import re
logging.basicConfig(level=logging.INFO)

from .io import to_memory
from .misc import dask_compute


# deprecated
def select_layer(adata, layer, force_dense=False, force_sparse=False, dtype='float32'):
    from scipy.sparse import csr_matrix, issparse
    from dask.array import Array as DaskArray
    
    # select matrix
    if layer == 'X' or layer is None:
        matrix = adata.X
    elif layer in adata.layers:
        matrix = adata.layers[layer]
    elif layer in ['raw', 'counts']:
        try:
            assert not isinstance(adata.raw, type(None))
            matrix = adata.raw[:, adata.var_names].X
        except AssertionError as e:
            raise ValueError(f'Cannot find layer "{layer}" and no counts in adata.raw') from e
    else:
        raise ValueError(f'Invalid layer {layer}')

    if isinstance(matrix, DaskArray):
        return matrix

    if force_dense and force_sparse:
        raise ValueError('force_dense and force_sparse cannot both be True')

    if force_dense:
        matrix = np.asarray(matrix.todense()) if issparse(matrix) else matrix
        return matrix.astype(dtype)

    if force_sparse:
        return matrix if issparse(matrix) else csr_matrix(matrix, dtype=dtype)

    return matrix


def select_neighbors(adata, output_type):
    neighbors_key = f'neighbors_{output_type}'
    adata.uns['neighbors'] = adata.uns[neighbors_key]
    adata.obsp['connectivities'] = adata.obsp[adata.uns[neighbors_key]['connectivities_key']]
    adata.obsp['distances'] = adata.obsp[adata.uns[neighbors_key]['distances_key']]
    return adata


def _filter_batch(adata, batch_key=None, min_cells=10):
    """
    Filter cells from batches with too few cells
    """
    mask = np.ones(adata.n_obs, dtype=bool)
    if batch_key is not None:
        cells_per_batch = adata.obs[batch_key].value_counts()
        if cells_per_batch.min() < min_cells:
            batches = cells_per_batch[cells_per_batch >= min_cells].index
            print(f'Keeping {len(batches)} out of {len(cells_per_batch)} batches with at least {min_cells} cells')
            mask = adata.obs[batch_key].isin(batches)
    return mask


def _filter_genes(adata, verbose=True, return_varnames=False, **kwargs):
    """
    :return: boolean mask of genes that pass the filter
    """
    import scanpy as sc
    from dask import array as da
    from tqdm.dask import TqdmCallback
    from contextlib import nullcontext
    
    # Handle empty adata
    if adata.n_obs == 0:
        gene_subset = np.zeros(adata.n_vars, dtype=bool)
        if return_varnames:
            gene_subset = adata.var_names[gene_subset]
        return gene_subset
    
    X = adata.X
    context = TqdmCallback(desc='Filter genes', miniters=1) if verbose and isinstance(X, da.Array) else nullcontext()
    
    with context:
        gene_subset, _ = sc.pp.filter_genes(X, **kwargs)
    
    if return_varnames:
        gene_subset = adata.var_names[gene_subset]
    
    return gene_subset


def subset_hvg(
    adata: ad.AnnData,
    to_memory: [str, list, bool] = 'X',
    var_column: str = 'highly_variable',
    compute_dask: bool = True,
    add_column: str = 'highly_variable',
    min_cells: int = 0,
) -> (ad.AnnData, bool):
    """
    Subset to highly variable genes
    :param adata: anndata object
    :param to_memory: layers to convert to memory
    :return: subsetted anndata object, bool indicating whether subset was performed
    """
    import sparse
    from tqdm.dask import TqdmCallback
    
    if var_column is None or var_column == 'None':
        var_column = 'mask'
        adata.var[var_column] = True
    elif var_column not in adata.var.columns and adata.var.shape[1] == 0:
        # Integration has already subsetted features; no .var columns survived.
        # Treat as if the mask was applied (all remaining genes selected),
        # mirroring the var_column=None branch above.
        adata.var[var_column] = True
    assert var_column in adata.var.columns, f'Column {var_column} not found in adata.var'
    assert adata.var[var_column].dtype == bool, f'Column {var_column} is not boolean'
    
    if add_column is not None and add_column != var_column:
        adata.var[add_column] = adata.var[var_column]
    
    if min_cells > 0:
        print(f'Filter to genes that are in at least {min_cells} cells...', flush=True)
        if isinstance(adata.X, da.Array):
            low_count_mask = adata.X.map_blocks(sparse.COO).sum(axis=0) < min_cells
            with TqdmCallback(desc=f"Determine genes with < {min_cells} cells"):
                low_count_mask = low_count_mask.compute().todense()
        else:
            filtered_mask = _filter_genes(adata, min_cells=min_cells)
        adata.var[var_column] &= filtered_mask
    
    if adata.var[var_column].sum() == adata.var.shape[0]:
        warnings.warn('All genes are highly variable, not subsetting')
        subsetted = False
    else:
        subsetted = True
        print(
            f'Subsetting to {adata.var[var_column].sum()} features from {var_column}...',
            flush=True
        )
        adata._inplace_subset_var(adata.var_names[adata.var[var_column]])
    
    if compute_dask:
        adata = dask_compute(adata, layers=to_memory)
    else:
        adata = adata_to_memory(
            adata,
            layers=to_memory,
            verbose=False,
        )
    
    return adata, subsetted


def adata_to_memory(
    adata: ad.AnnData,
    layers: [str, list, bool] = None,
    verbose: bool = False,
) -> ad.AnnData:
    if layers is None or layers is True:
        layers = ['X', 'raw'] + list(adata.layers.keys())
    elif isinstance(layers, str):
        layers = [layers]
    elif layers is False:
        return adata
    
    for layer in layers:
        if verbose:
            print(f'Convert {layer} to memory...', flush=True)
        if layer in adata.layers:
            adata.layers[layer] = to_memory(adata.layers[layer])
        elif layer == 'X':
            adata.X = to_memory(adata.X)
        elif layer == 'raw':
            if adata.raw is None:
                continue
            adata_raw = adata.raw.to_adata()
            adata_raw.X = to_memory(adata.raw.X)
            adata.raw = adata_raw
        elif verbose:
            print(f'Layer {layer} not found, skipping...', flush=True)
    return adata


def parse_gene_names(adata, gene_list):
    var_names = adata.var_names.astype(str)
    genes_not_in_var = [str(g) for g in gene_list if g not in var_names]
    gene_list = [g for g in gene_list if g in var_names]

    if genes_not_in_var:
        mask = var_names.str.contains(pat='|'.join(re.escape(g) for g in genes_not_in_var))
        gene_list += var_names[mask].tolist()

    return gene_list


def match_genes(var_df, gene_list, column=None, return_index=True, as_list=False, return_mask=False):
    """
    Match genes from a provided list against a variable (gene) annotation table.

    This function supports three kinds of gene specifications:
      1. Direct gene identifiers contained in `gene_list`.
      2. File paths pointing to plain-text files, one gene identifier per line.
      3. HTTP(S) URLs pointing to plain-text resources, one gene identifier per line.

    Any entry in `gene_list` that resolves to an existing local file path or a valid URL
    is read and replaced by the list of gene identifiers extracted from that resource.
    All collected gene identifiers (original and loaded) are then used to construct a
    regular-expression pattern which is matched against the genes in `var_df`.

    Parameters
    ----------
    var_df
        A pandas-like DataFrame or Series containing gene annotations, typically
        ``adata.var``. Matching is performed either on its index or on a specified
        column (or columns).
    gene_list
        Iterable of strings representing gene identifiers to match. Entries may
        also be:
          * Paths to local text files (one gene per line), or
          * HTTP(S) URLs to text resources (one gene per line).
        For such entries, the file/URL is read and the contained genes are added
        to the list of genes to match.
    column
        Optional name(s) of column(s) in ``var_df`` to use for matching. Can be:
          * A single column name (str), or
          * A list/tuple of column names to match against.
        If ``None``, the index of ``var_df`` is used. Matching is performed
        against any of the specified columns (OR operation).
    return_index
        If ``True``, return the index labels corresponding to the matched genes.
        If ``False``, return the matched values themselves (e.g. a Series or
        column slice from ``var_df``).
    as_list
        If ``True``, convert the final result (index or values) to a Python
        ``list`` before returning. If ``False``, return it in its native type
        (e.g. pandas Index or Series).
    return_mask
        If ``True``, return a tuple of (genes, mask) where mask is a boolean
        array indicating which rows in ``var_df`` matched the gene list.
        If ``False``, return only the matched genes.

    Returns
    -------
    genes
        The set of matched genes. Depending on the combination of flags:
          * If ``return_index=True`` and ``as_list=False``: a pandas Index of
            matched gene indices.
          * If ``return_index=True`` and ``as_list=True``: a ``list`` of matched
            gene index labels.
          * If ``return_index=False`` and ``as_list=False``: a pandas Series or
            similar object of matched gene values.
          * If ``return_index=False`` and ``as_list=True``: a ``list`` of matched
            gene values.
    mask : ndarray, optional
        If ``return_mask=True``, returns a boolean array of shape (n_genes,)
        indicating which rows matched the gene list. Otherwise, not returned.

    Raises
    ------
    Exception
        Propagates any exception raised while reading from URLs or if regex-based
        matching against ``var_df`` fails. Additional context is logged using the
        module logger.
    """
    import urllib
    from pathlib import Path

    def is_url(url):
        try:
            result = urllib.parse.urlparse(url)
            return all([result.scheme, result.netloc])
        except ValueError:
            return False

    genes_from_path = dict()
    for gene in gene_list:
        if Path(gene).exists():
            with open(gene, 'r') as f:
                genes_from_path[gene] = f.read().splitlines()
        elif is_url(gene):
            try:
                with urllib.request.urlopen(gene) as f:
                    genes_from_path[gene] = f.read().decode('utf-8').splitlines()
            except Exception as e:
                logging.error(f'Error reading gene list from URL {gene}...')
                raise e

    for path, genes in genes_from_path.items():
        logging.info(f'Gene list from {path}: {len(genes)} genes')
        gene_list.extend(genes)
        gene_list.remove(path)

    # Normalize column parameter to list
    if column is None:
        columns = [None]
    elif isinstance(column, str):
        columns = [column]
    else:
        columns = list(column)
    
    # Get the genes series from the dataframe or specified columns
    # Create a mask for matching across all specified columns (OR operation)
    mask = np.zeros(len(var_df), dtype=bool)
    
    # Handle empty gene list - return empty result instead of matching all genes
    if not gene_list:
        logging.warning('Gene list is empty after processing, returning empty result')
        genes = var_df.index.to_series()[:0]  # Return empty Series
    else:
        try:
            pattern = '|'.join(re.escape(g) for g in gene_list)
            
            # Match against each column and combine with OR
            for col in columns:
                if col is None:
                    genes_series = var_df.index.to_series()
                else:
                    genes_series = var_df[col]
                
                col_mask = genes_series.astype(str).str.contains(pattern, regex=True)
                mask |= col_mask.values  # OR operation
            
            # Get results from first column for the genes output
            genes_series = var_df.index.to_series() if columns[0] is None else var_df[columns[0]]
            genes = genes_series[mask].drop_duplicates()
            
        except Exception as e:
            logging.error(f'Error: {e}')
            logging.error(f'Gene list: {gene_list}')
            logging.error(f'Pattern: {pattern}')
            logging.error(f'Gene names: {var_df.index}')
            raise e

    if return_index:
        genes = genes.index

    if as_list:
        genes = genes.tolist()
    
    if return_mask:
        return genes, mask
    
    return genes
