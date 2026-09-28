import numpy as np
from .utils import select_neighbors, rename_categories, scanpy_to_neighborsresults


def clisi(adata, output_type, batch_key, label_key, **kwargs):
    import scib

    adata = adata[adata.obs[label_key].notna()].copy()
    score = scib.me.clisi_graph(
        adata,
        batch_key=batch_key,
        label_key=label_key,
        type_='knn',
        subsample=0.5 * 100,
        scale=True,
    )
    return (score, 'cLISI graph')


def clisi_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    if output_type == 'knn':
        return np.nan
    
    adata = select_neighbors(adata, output_type)
    adata = adata[adata.obs[label_key].notna()].copy()
    labels = rename_categories(adata, label_key)
    
    return scib_metrics.clisi_knn(
        X=scanpy_to_neighborsresults(adata),
        labels=labels
    ).mean()


def ilisi(adata, output_type, batch_key, label_key, **kwargs):
    import scib

    score = scib.me.ilisi_graph(
        adata,
        batch_key=batch_key,
        type_='knn',
        subsample=0.5 * 100,
        scale=True,
    )
    return (score, 'iLISI graph')


def ilisi_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics
    
    if output_type == 'knn':
        return np.nan
    
    adata = select_neighbors(adata, output_type)
    batches = rename_categories(adata, batch_key)

    return scib_metrics.ilisi_knn(
        X=scanpy_to_neighborsresults(adata),
        batches=batches
    ).mean()
