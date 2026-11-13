######### UMAP and embedding plots #########

use rule umap from preprocessing as integration_compute_umap with:
    input:
        anndata=rules.integration_postprocess.output.zarr,
        rep=rules.integration_postprocess.output.zarr,
    output:
        zarr=directory(out_dir / f'{paramspace.wildcard_pattern}.zarr'),
    threads:
        lambda wildcards: max(
            1,
            mcfg.get_from_parameters(wildcards, 'threads', exclude=['output_type']),
        )
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt, attempt_to_cpu=2),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt, attempt_to_cpu=2),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt, attempt_to_cpu=2),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt, attempt_to_cpu=2),


def get_colors(wildcards):
    dataset = wildcards.dataset
    labels = mcfg.get_from_parameters(wildcards, 'label')
    labels = labels if isinstance(labels, list) else [labels]
    batch = mcfg.get_from_parameters(wildcards, 'batch')
    batch = batch if isinstance(batch, list) else [batch]
    
    umap_colors = mcfg.get_from_parameters(wildcards, 'umap_colors', default=[])
    umap_colors += mcfg.get_from_parameters(wildcards, 'plots', default={}).get('colors', [])
    return [*labels, *batch, *umap_colors]


use rule plots from preprocessing as integration_plot_umap with:
    input:
        anndata=rules.integration_compute_umap.output.zarr,
    output:
        plots=directory(image_dir / 'umap' / paramspace.wildcard_pattern.replace('--output_type', '/output_type')),
    params:
        color=get_colors,
        basis='X_umap',
        min_cells_per_category=0.0001,
        plot_centroids=lambda w: mcfg.get_from_parameters(w, 'plots', default={}).get('plot_centroids', []),
        gene_chunk_size=lambda w: mcfg.get_from_parameters(w, 'plots', default={}).get('plot_gene_chunk_size', 12),
        add_outline=lambda w: mcfg.get_from_parameters(w, 'plots', default={}).get('add_outline', False),
        max_label_length=lambda w: mcfg.get_from_parameters(w, 'plots', default={}).get('max_label_length', 10),
    threads:
        lambda wildcards: max(
            1,
            mcfg.get_from_parameters(wildcards, 'threads', exclude=['output_type']),
        )
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
