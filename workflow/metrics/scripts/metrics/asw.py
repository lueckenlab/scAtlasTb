import numpy as np

from .utils import rename_categories


def asw_batch(adata, output_type, batch_key, label_key, **kwargs):
    import scib

    if output_type == 'knn':
        return np.nan

    return scib.me.silhouette_batch(
        adata,
        batch_key=batch_key,
        label_key=label_key,
        embed='X_emb' if output_type == 'embed' else 'X_pca',
    )


def asw_batch_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    if output_type == 'knn':
        return np.nan

    X = adata.obsm['X_emb'] if output_type == 'embed' else adata.obsm['X_pca']
    X = X if isinstance(X, np.ndarray) else X.toarray()
    labels = rename_categories(adata, label_key)
    batches = rename_categories(adata, batch_key)

    return scib_metrics.silhouette_batch(
        X=X,
        batch=batches,
        labels=labels,
    )


def asw_label(adata, output_type, batch_key, label_key, **kwargs):
    import scib

    if output_type == 'knn':
        return np.nan

    adata = adata[adata.obs[label_key].notna()].copy()
    return scib.me.silhouette(
        adata,
        label_key=label_key,
        embed='X_emb' if output_type == 'embed' else 'X_pca',
    )


def asw_label_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    if output_type == 'knn':
        return np.nan

    X = adata.obsm['X_emb'] if output_type == 'embed' else adata.obsm['X_pca']
    X = X if isinstance(X, np.ndarray) else X.toarray()
    labels = rename_categories(adata, label_key)

    return scib_metrics.silhouette_label(
        X=X,
        labels=labels,
    )
