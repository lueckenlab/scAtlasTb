import numpy as np
import scanpy as sc
import pandas as pd


def write_metrics(filename, metric_names, scores, output_types, **kwargs):
    """
    Write metrics for output type specific scores
    :param filename: file to write to
    :param metric_names: list of metric names
    :param scores: list of scores, order must match that of metric_names
    :param output_types: list of output types, order must match that of metric_names
    :param kwargs: additional information to add to output
    """
    meta_names = list(kwargs.keys())
    meta_values = [kwargs[col] for col in meta_names]

    records = [
        (*meta_values, metric, output_type, score)
        for metric, score, output_type
        in zip(metric_names, scores, output_types)
    ]
    df = pd.DataFrame.from_records(
        records,
        columns=meta_names + ['metric_name', 'output_type', 'score']
    )
    df = df.explode('score')
    df.to_csv(filename, sep='\t', index=False)


def rename_categories(adata, obs_col):
    s = adata.obs[obs_col].astype('category')
    s = s.cat.rename_categories({i for i, _ in enumerate(s.cat.categories)})
    return s.to_numpy()


def select_neighbors(adata, output_type):
    # neighbors_key = f'neighbors_{output_type}'
    neighbors_key = 'neighbors'
    adata.uns['neighbors'] = adata.uns[neighbors_key]
    
    connectivities_key = adata.uns[neighbors_key]['connectivities_key']
    assert connectivities_key in adata.obsp, f'Connectivities key "{connectivities_key}" missing from adata.obsp {adata}'
    adata.obsp['connectivities'] = adata.obsp[connectivities_key]
    
    distances_key = adata.uns[neighbors_key]['distances_key']
    assert distances_key in adata.obsp, f'Distances key "{distances_key}" missing from adata.obsp {adata}'
    adata.obsp['distances'] = adata.obsp[distances_key]
    
    return adata


def scanpy_to_neighborsresults(adata):
    from scib_metrics.nearest_neighbors import NeighborsResults

    n_neighbors = adata.uns["neighbors"]["params"]["n_neighbors"]
    n_obs = adata.n_obs
    dist_matrix = adata.obsp["distances"].tocsr()

    # scib_metrics expects each cell as neighbor 0 (distance 0), n_neighbors columns total.
    # scanpy stores n_neighbors-1 off-diagonal neighbors (self excluded); some backends
    # (e.g. rapids) may include self. Force the convention regardless of backend.
    indices = np.full((n_obs, n_neighbors), -1, dtype=int)
    distances = np.full((n_obs, n_neighbors), np.inf, dtype=float)
    indices[:, 0] = np.arange(n_obs)
    distances[:, 0] = 0.0

    indptr = dist_matrix.indptr
    for i in range(n_obs):
        sl = slice(indptr[i], indptr[i + 1])
        row = dist_matrix.indices[sl]
        d = dist_matrix.data[sl]
        keep = row != i                            # drop self if the backend stored it
        row, d = row[keep], d[keep]
        order = np.argsort(d)[: n_neighbors - 1]   # K-1 nearest non-self
        k = len(order)
        indices[i, 1 : 1 + k] = row[order]
        distances[i, 1 : 1 + k] = d[order]

    return NeighborsResults(indices=indices, distances=distances)
