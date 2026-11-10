import numpy as np
import scanpy as sc
from .utils import select_neighbors, rename_categories, scanpy_to_neighborsresults


def nmi(adata, output_type, batch_key, label_key, cluster_keys, **kwargs):
    import scib

    adata = select_neighbors(adata, output_type)
    adata = adata[adata.obs[label_key].notna()].copy()
    # scib.cl.cluster_optimal_resolution(
    #     adata=adata,
    #     label_key=label_key,
    #     cluster_key=cluster_key,
    #     metric=scib.me.nmi,
    #     use_rep=None,
    # )
    cluster_key = cluster_keys[0]
    nmi_max = 0
    for col in cluster_keys:
        nmi = scib.me.nmi(adata, label_key, col)
        if nmi > nmi_max:
            nmi_max = nmi
            cluster_key = col

    score = scib.me.nmi(
        adata=adata,
        label_key=label_key,
        cluster_key=cluster_key,
    )
    return (score, 'NMI')

def nmi_leiden_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    adata = select_neighbors(adata, output_type)
    labels = rename_categories(adata, label_key)

    scores = scib_metrics.nmi_ari_cluster_labels_leiden(
        X=scanpy_to_neighborsresults(adata),
        labels=labels,
        optimize_resolution=True,
    )
    return scores['nmi']


def nmi_kmeans_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    if output_type == 'knn':
        return np.nan

    labels = rename_categories(adata, label_key)
    X = adata.obsm['X_emb'] if output_type == 'embed' else adata.obsm['X_pca']
    X = X if isinstance(X, np.ndarray) else np.asarray(X.todense())

    scores = scib_metrics.nmi_ari_cluster_labels_kmeans(
        X=X,
        labels=labels,
    )
    return scores['nmi']
