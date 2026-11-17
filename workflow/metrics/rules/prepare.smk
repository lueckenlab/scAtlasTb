attempt_to_cpu = 2

rule prepare:
    input:
        lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'prepare.zarr'),
    params:
        # labels=lambda wildcards: mcfg.get_from_parameters(wildcards, 'label', single_value=False),
        neighbor_args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors', default={}),
        recompute_neighbors=lambda wildcards: mcfg.get_from_parameters(wildcards, 'recompute_neighbors', default=False),
        unintegrated_layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'unintegrated', default='X'),
        corrected_layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'corrected', default='X'),
        var_mask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'var_mask', default='highly_variable'),
        output_type=lambda wildcards: mcfg.get_from_parameters(wildcards, 'output_type', default='embed'),
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=not config.get("use_gpu", False))
    resources:
        partition=lambda w, attempt: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        qos=lambda w, attempt: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu, factor=3),
        gpu=lambda w, attempt: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
    script:
        '../scripts/prepare.py'


use rule pca from preprocessing as metrics_pca with:
    input:
        lambda wildcards: mcfg.get_input_file(**wildcards),
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'pca.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(
            wildcards,
            'pca',
            default=dict(
                mask_var=mcfg.get_from_parameters(wildcards, 'var_mask', default='highly_variable')
            )
        ),
        use_gpu=config.get("use_gpu", False),
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'unintegrated', default='X'),
        subset=True,
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=not config.get("use_gpu", False))
    resources:
        partition=lambda w, attempt: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        qos=lambda w, attempt: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu, factor=3),
        gpu=lambda w, attempt: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),


rule pca_all:
    input:
        pca=[
            mcfg.get_output_files(rules.metrics_pca.output, subset_dict=dict(dataset=dataset))
            for dataset in mcfg.get_datasets()
            if sum(mcfg.get_from_parameters(dict(dataset=dataset), 'comparison', single_value=False)) > 0
        ]
    localrule: True

# Clustering for cluster-based metrics

def scale_mem_mb(wildcards, attempt, factor=3, min_mb=1_000, profile=None):
    # TODO: add to utils
    if profile is None:
        profile = mcfg.get_profile(wildcards)
    mem_mb = mcfg.get_resource(profile=profile, resource_key='mem_mb', attempt=attempt)

    try:
        mem_mb = int(mem_mb)
    except ValueError:
        return mem_mb
    
    mem_mb = int(mem_mb // factor)
    return max(min_mb, mem_mb)


use rule cluster from clustering as metrics_cluster with:
    input:
        zarr=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'cluster_resolutions' / '{algorithm}--{resolution}--{level}.zarr'),
    params:
        clustering_args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'clustering', default={}).get('kwargs', {}),
        overwrite=lambda wildcards: mcfg.get_from_parameters(wildcards, 'clustering', default={}).get('overwrite', True),
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=not config.get("use_gpu", False))
    resources:
        partition=lambda w, attempt: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        qos=lambda w, attempt: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu, factor=3),
        gpu=lambda w, attempt: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),


use rule merge from clustering as metrics_cluster_collect with:
    input:
        zarr=rules.prepare.output.zarr,
        cluster_anno=expand(
            rules.metrics_cluster.output.zarr,
            resolution=[2 * (x + 1) / 10 for x in range(10)],
            algorithm='leiden',
            level=1,
            allow_missing=True,
        ),
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'clustered.zarr'),


rule cluster_all:
    input:
        clusters=[
            mcfg.get_output_files(rules.metrics_cluster_collect.output, subset_dict=dict(dataset=dataset))
            for dataset in mcfg.get_datasets()
            if sum(mcfg.get_from_parameters(dict(dataset=dataset), 'needs_clustering', single_value=False)) > 0
            and not mcfg.get_from_parameters(dict(dataset=dataset), 'clustering', default={}).get('precomputed_key')
        ]
    localrule: True


# Gene scoring for gene-based metrics

rule score_genes:
    input:
        zarr=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / 'prepare' / paramspace.wildcard_pattern / 'gene_scores.zarr'),
    params:
        unintegrated_layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'unintegrated', default='X'),
        var_mask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'var_mask', default='highly_variable'),
        gene_sets=lambda wildcards: mcfg.get_gene_sets(wildcards.dataset),
        n_permutations=lambda wildcards: mcfg.get_from_parameters(wildcards, 'gene_score', default={}).get('n_permutations', 20),
        n_quantiles=lambda wildcards: mcfg.get_from_parameters(wildcards, 'gene_score', default={}).get('n_quantiles', 5),
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=not config.get("use_gpu", False))
    resources:
        partition=lambda w, attempt: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        qos=lambda w, attempt: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu, factor=3),
        gpu=lambda w, attempt: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
    script:
        '../scripts/score_genes.py'


rule score_genes_all:
    input:
        gene_scores=[
            mcfg.get_output_files(rules.score_genes.output, subset_dict=dict(dataset=dataset))
            for dataset in mcfg.get_datasets()
            if sum(mcfg.get_from_parameters(dict(dataset=dataset), 'use_gene_set', single_value=False)) > 0
        ]
    localrule: True


rule prepare_all:
    input:
        prepare=mcfg.get_output_files(rules.prepare.output),
        pca=rules.pca_all.input,
        gene_scores=rules.score_genes_all.input,
        clusters=rules.cluster_all.input,
    localrule: True
