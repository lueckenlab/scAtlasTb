import numpy as np
from .utils import rename_categories



def bras_batch(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    if output_type == 'knn':
        return np.nan
    
    X = adata.obsm['X_emb'] if output_type == 'embed' else adata.obsm['X_pca']
    X = X if isinstance(X, np.ndarray) else X.toarray()
    labels = rename_categories(adata, label_key)
    batches = rename_categories(adata, batch_key)

    return scib_metrics.bras(X, labels, batches, chunk_size=256, metric='cosine', between_cluster_distances='mean_other')




