class MetricNotDefinedError(RuntimeError):
    
    def __init__(self, wildcards):
        super().__init__(
            f'Metric {wildcards.metric} is not defined for dataset {wildcards.dataset} and file_id {wildcards.file_id}'
        )


def get_metric_input(wildcards):
    files = dict(zarr=rules.prepare.output.zarr)
    if mcfg.get_from_parameters(wildcards, 'needs_clustering', default=False) and \
       not mcfg.get_from_parameters(wildcards, 'clustering', default={}).get('precomputed_key'):
        files['zarr'] = rules.metrics_cluster_collect.output.zarr
    if mcfg.get_from_parameters(wildcards, 'use_gene_set', default=False):
        files['zarr'] = rules.score_genes.output.zarr
    if mcfg.get_from_parameters(wildcards, 'comparison', default=False):
        files['raw'] = rules.metrics_pca.output.zarr
    return files


def get_mem_mb(wildcards, attempt):
    profile = mcfg.get_profile(wildcards)
    mem_mb = mcfg.get_resource(profile=profile, resource_key='mem_mb', attempt=attempt)
    input_type = mcfg.get_from_parameters(wildcards, 'input_type')
    comparison = mcfg.get_from_parameters(wildcards, 'comparison', default=False)

    try:
        mem_mb = int(mem_mb)
    except ValueError:
        return mem_mb
    
    if input_type == 'full' or comparison:
        return mem_mb
    return max(1_000, int(mem_mb // 3))


rule run:
    message:
       """
       Metrics: Evaluate metric={wildcards.metric} on dataset={wildcards.dataset} and file_id={wildcards.file_id}
       input: {input}
       output: {output}
       wildcards: {wildcards}
       resources: gpu={resources.gpu} mem_mb={resources.mem_mb}
       """
    input:
        unpack(get_metric_input)
    output:
        metric=mcfg.out_dir / paramspace.wildcard_pattern / 'label={label}--batch={batch}' / 'metric={metric}.tsv'
    benchmark:
        mcfg.out_dir / paramspace.wildcard_pattern / '.benchmark' / 'label={label}--batch={batch}' / 'metric={metric}.tsv'
    params:
        metric_type=lambda wildcards: mcfg.get_from_parameters(wildcards, 'metric_type', default=MetricNotDefinedError(wildcards)),
        allowed_output_types=lambda wildcards: mcfg.get_from_parameters(wildcards, 'allowed_output_types', default=MetricNotDefinedError(wildcards)),
        input_type=lambda wildcards: mcfg.get_from_parameters(wildcards, 'input_type', default=MetricNotDefinedError(wildcards)),
        clustering=lambda wildcards: mcfg.get_from_parameters(wildcards, 'clustering', default={}),
        use_covariate=lambda wildcards: mcfg.get_from_parameters(wildcards, 'use_covariate', default=False),
        use_gene_set=lambda wildcards: mcfg.get_from_parameters(wildcards, 'use_gene_set', default=False),
        covariates=lambda wildcards: mcfg.get_from_parameters(wildcards, 'covariate', default=[]),
        gene_sets=lambda wildcards:
            mcfg.get_gene_sets(wildcards.dataset)
            if mcfg.get_from_parameters(wildcards, 'use_gene_set', default=False)
            else {},
        env=lambda wildcards: mcfg.get_from_parameters(wildcards, 'env', check_null=True, default=MetricNotDefinedError(wildcards)),
    conda:
        lambda wildcards, params: get_env(config, params.env)
    threads:
        lambda wildcards: mcfg.get_from_parameters(wildcards, 'threads', default=1, as_type=int),
    resources:
        partition=lambda w: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w)),
        qos=lambda w: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w)),
        mem_mb=lambda w, attempt: get_mem_mb(w, attempt),
        gpu=lambda w: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w)),
    script:
        '../scripts/run.py'


rule run_all:
    input:
        rules.prepare_all.input,
        mcfg.get_output_files(rules.run.output, all_params=True),
    localrule: True
