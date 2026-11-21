from .utils import select_neighbors, rename_categories, scanpy_to_neighborsresults


def graph_connectivity(adata, output_type, batch_key, label_key, **kwargs):
    import scib

    adata = select_neighbors(adata, output_type)
    adata = adata[adata.obs[label_key].notna()].copy()
    return scib.me.graph_connectivity(
        adata,
        label_key=label_key
    )


def graph_connectivity_y(adata, output_type, batch_key, label_key, **kwargs):
    import scib_metrics

    adata = select_neighbors(adata, output_type)
    adata = adata[adata.obs[label_key].notna()].copy()
    labels = rename_categories(adata, label_key)

    return scib_metrics.graph_connectivity(
        X=scanpy_to_neighborsresults(adata),
        labels=labels
    )
