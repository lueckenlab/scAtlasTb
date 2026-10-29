def get_cluster_keys(wildcards):
    return mcfg.get_output_files(
        '{algorithm}_{resolution}_{level}',
        subset_dict=dict(wildcards),
        all_params=True
    )


use rule plots from preprocessing as clustering_plot_umap with:
    input:
        zarr=rules.clustering_merge.output.zarr,
    output:
        plots=directory(mcfg.image_dir / 'dataset~{dataset}' / 'file_id~{file_id}' / 'umap')
    params:
        basis='X_umap',
        color=lambda wildcards: mcfg.get_from_parameters(wildcards, 'umap_colors', default=[]) +
            mcfg.get_from_parameters(wildcards, 'plots', default={}).get('colors', []),
        plot_centroids=lambda wildcards: get_cluster_keys(wildcards) +
            mcfg.get_from_parameters(wildcards, 'plots', default={}).get('plot_centroids', []),
    threads:
        lambda wildcards: max(
            1,
            min(
                mcfg.get_from_parameters(wildcards, 'threads', default=5),
                len(mcfg.get_from_parameters(wildcards, 'umap_colors', default=[]))
            )
        )
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb',attempt=attempt, factor=1),
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),


def get_mem_mb(attempt, profile):
    mem_mb = mcfg.get_resource(profile=profile, resource_key='mem_mb', attempt=attempt)
    try:
        mem_mb = int(mem_mb)
    except ValueError:
        return mem_mb
    return min(100_000, int(mem_mb // 5))


use rule plot_evaluation from clustering as clustering_plot_evaluation with:
    input:
        zarr=rules.clustering_merge.output.zarr,
    output:
        plots=directory(mcfg.image_dir / 'dataset~{dataset}' / 'file_id~{file_id}' / 'evaluation')
    params:
        cluster_keys=get_cluster_keys,
        covariates=lambda wildcards: mcfg.get_from_parameters(wildcards, 'umap_colors', default=[]) +
            mcfg.get_from_parameters(wildcards, 'plots', default={}).get('colors', []),
        max_categories=lambda wildcards: mcfg.get_from_parameters(wildcards, 'plots', default={}).get('max_categories', 100),
    conda:
        get_env(config, 'scanpy')
    retries: 1
    threads:
        lambda wildcards: max(
            1,
            min(
                mcfg.get_from_parameters(wildcards, 'threads', default=5),
                len(mcfg.get_from_parameters(wildcards, 'umap_colors', default=[]))
            )
        )
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
        mem_mb=lambda w, attempt: get_mem_mb(attempt=attempt, profile='cpu'),


rule plots_all:
    input:
        mcfg.get_output_files(rules.clustering_plot_umap.output),
    localrule: True

rule plot_evaluation_all:
    input:
        mcfg.get_output_files(rules.clustering_plot_evaluation.output)
    localrule: True