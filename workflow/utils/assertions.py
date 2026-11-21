import logging
import anndata as ad
import numpy as np


class EmptyDataError(Exception):
    """All observations were filtered out and the module is configured to fail on empty output"""


def _yaml_block(keys: list, value=None) -> str:
    """Render nested config keys as an indented YAML block, optionally ending with a value"""
    lines = [f'{"  " * (i + 1)}{key}:' for i, key in enumerate(keys)]
    if value is not None:
        lines[-1] += f' {value}'
    return '\n'.join(lines)


def check_empty(
    n_input: int,
    n_remaining: int,
    fail_on_empty: bool = True,
    reason: str = '',
    param: str = None,
    value=None,
    module: str = None,
    dataset: str = None,
) -> bool:
    """
    Decide how to proceed after filtering, before any expensive computation.

    An input that is already empty never fails, since there is nothing to filter.
    If the filters remove everything from a non-empty input, raise EmptyDataError
    (fail_on_empty=True) or log a warning so that the caller writes an empty output.

    :param n_input: number of observations before filtering
    :param n_remaining: number of observations that pass the filters
    :param fail_on_empty: raise EmptyDataError if all observations are filtered out
    :param reason: which filter removed everything, including the observed values
    :param param: module parameter that caused the filtering, as dot-separated keys
        relative to the module configuration, e.g. 'rank_genes_groups.min_cells'
    :param value: current value of param
    :param module: module name, to point to the module configuration in the config file
    :param dataset: dataset name, to point to the dataset entry of the module configuration
    :return: True if the caller should write an empty output, False to continue
    """
    if n_input == 0:
        logging.warning('Input data is empty, writing empty output.')
        return True
    if n_remaining > 0:
        return False

    message = f'All {n_input} observations were filtered out: {reason}'
    if not fail_on_empty:
        logging.warning(f'{message}\nWriting empty output, since fail_on_empty is false.')
        return True

    module_keys = [module] if module is not None else []
    options = []
    if param is not None:
        options.append(
            f'change `{param.split(".")[-1]}` (currently {value}):\n\n'
            + _yaml_block(module_keys + param.split('.'), value)
        )
    options.append(
        'write an empty output instead:\n\n'
        + _yaml_block(module_keys + ['fail_on_empty'], 'false')
    )
    fix = '\n\nor\n\n'.join(f'- {option}' for option in options)

    location = ''
    if module is not None:
        # show where the module configuration is located in the config file
        location = (
            f'The {module} module configuration is located either under\n\n'
            + _yaml_block(['DATASETS', dataset or '<dataset>', module])
            + '\n\nfor this dataset, or under\n\n'
            + _yaml_block(['defaults', module])
            + '\n\nfor all datasets.\n\n'
        )

    raise EmptyDataError(
        f'\n\n{message}\n\n'
        'If you did not expect this, first check that the input data is what you expect\n'
        '(e.g. the values in the columns used for filtering).\n\n'
        f'If the data is correct, update the module configuration in your config file and either\n\n{fix}\n\n'
        f'{location}'
    )


def assert_pca(adata: ad.AnnData, check_varm=True):
    assert 'X_pca' in adata.obsm, f'PCA layer missing from obsm:\n{adata}'
    if check_varm:
        assert 'PCs' in adata.varm, f'PCS missing from varm:\n{adata}'
    assert 'pca' in adata.uns, f'pca missing from uns:\n{adata}'
    assert 'variance' in adata.uns['pca'], f'variance missing from uns[pca]:\n{adata.uns}'
    assert 'variance_ratio' in adata.uns['pca'], f'variance_ratio missing from uns[pca]:\n{adata.uns}'


def assert_neighbors(
    adata,
    neighbors_key='neighbors',
    conn_key='connectivities',
    dist_key='distances',
    check_params=True,
    check_n_neighbors=False,
):
    assert neighbors_key in adata.uns, f'neighbors key "{neighbors_key}" not on .uns'
    assert 'connectivities_key' in adata.uns[neighbors_key], f'"connectivities_key" not in .uns["{neighbors_key}"]'
    assert 'distances_key' in adata.uns[neighbors_key], f'"distances_key" not in .uns["{neighbors_key}"]'
    assert adata.uns[neighbors_key]['connectivities_key'] == conn_key, f'"{conn_key} is not saved as conectivities_key for "{neighbors_key}": {adata.uns[neighbors_key]}'
    assert adata.uns[neighbors_key]['distances_key'] == dist_key, f'"{dist_key} is not saved as distances_key for "{neighbors_key}": {adata.uns[neighbors_key]}'
    assert conn_key in adata.obsp, f'"{conn_key}" not in .obsp {adata.obsp.keys()}'
    assert dist_key in adata.obsp, f'"{dist_key}" not in .obsp {adata.obsp.keys()}'
    if check_params:
        assert 'params' in adata.uns[neighbors_key], f'"params" missing from .uns["{neighbors_key}"]'
        assert 'use_rep' in adata.uns[neighbors_key]['params'], f'"use_rep" missing from .uns["{neighbors_key}"]["params"]'
    
    # check that all cells have the same number of neighbors
    n_neighbors = np.unique(adata.obsp['distances'].nonzero()[0], return_counts=True)[1]
    n_neighbors = np.unique(n_neighbors)
    try:
        assert len(n_neighbors) == 1, f'Cells do not have the same number of neighbors {n_neighbors}'
    except AssertionError as e:
        if check_n_neighbors:
            raise e
        print(f'Warning: {e}')