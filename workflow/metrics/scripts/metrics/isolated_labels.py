import numpy as np
from .utils import select_neighbors


def replace_last(source_string, replace_what, replace_with, cluster_key='leiden'):
    """
    Helper function for parsing clustering columns
    adapted from
    https://stackoverflow.com/questions/3675318/how-to-replace-some-characters-from-the-end-of-a-string/3675423#3675423
    """
    if not source_string.startswith(cluster_key):
        return source_string
    if not source_string.endswith(replace_what):
        return source_string
    head, _sep, tail = source_string.rpartition(replace_what)
    return head + replace_with + tail


def isolated_label_f1(adata, output_type, batch_key, label_key, cluster_keys, **kwargs):
    import scib

    adata = select_neighbors(adata, output_type)
    
    # rename cluster columns to expected format
    cluster_rename_map = {col: replace_last(col, '_1', '') for col in cluster_keys}
    adata.obs.rename(columns=cluster_rename_map, inplace=True)
    
    if any('leiden' in key for key in cluster_keys):
        cluster_key = 'leiden'
    elif any('louvain' in key for key in cluster_keys):
        cluster_key = 'louvain'
    else:
        raise NotImplementedError(f'No leiden or louvain cluster key found in {cluster_keys}')

    return scib.me.isolated_labels_f1(
        adata[adata.obs[label_key].notna()].copy(),
        label_key=label_key,
        batch_key=batch_key,
        cluster_key=cluster_key,
        embed=None,
    )


def isolated_label_asw(adata, output_type, batch_key, label_key, **kwargs):
    import scib

    if output_type == 'knn':
        return np.nan

    adata = adata[adata.obs[label_key].notna()].copy()
    return scib.me.isolated_labels_asw(
        adata,
        label_key=label_key,
        batch_key=batch_key,
        embed='X_emb' if output_type == 'embed' else 'X_pca',
    )


def isolated_label_asw_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    if output_type == 'knn':
        return np.nan

    X = adata.obsm['X_emb'] if output_type == 'embed' else adata.obsm['X_pca']
    X = X if isinstance(X, np.ndarray) else X.toarray()

    return scib_metrics.isolated_labels(
        X=X,
        labels=adata.obs[label_key],
        batch=adata.obs[batch_key],
    )
