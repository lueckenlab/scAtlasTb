attempt_to_cpu = 2

def get_threads(wildcards, attempt):
    use_gpu = mcfg.get_profile(wildcards) == 'gpu'
    if use_gpu and attempt <= attempt_to_cpu:
        return 1
    return max(
        1,
        mcfg.get_from_parameters(
            wildcards,
            'threads',
            exclude=['output_type'],
            default=1
        )
    )


rule prepare:
    message:
        """
        Prepare: Prepare dataset={wildcards.dataset} with file_id={wildcards.file_id} and var_mask={wildcards.var_mask}
        input: {input}
        output: {output}
        params: norm_counts={params.norm_counts} raw_counts={params.raw_counts} save_subset={params.save_subset}
        threads: {threads}
        """
    input:
        anndata=lambda wildcards: mcfg.get_input_file(wildcards.dataset, wildcards.file_id)
    output:
        zarr=directory(out_dir / 'prepare' / 'dataset~{dataset}' / 'file_id~{file_id}' / 'var_mask~{var_mask}.zarr'),
    params:
        norm_counts=lambda wildcards: mcfg.get_from_parameters(wildcards, 'norm_counts', exclude=['output_type']),
        raw_counts=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_counts', exclude=['output_type']),
        save_subset=lambda wildcards: mcfg.get_from_parameters(wildcards, 'save_subset', exclude=['output_type']),
    conda:
        get_env(config, 'scanpy')
    threads:
        lambda wildcards: mcfg.get_for_dataset(
            dataset=wildcards.dataset,
            query=[mcfg.module_name, 'threads'],
            default=1
        )
    resources:
        partition=lambda w, attempt: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        qos=lambda w, attempt: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu, factor=3),
        gpu=lambda w, attempt: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
    script:
        '../scripts/prepare.py'


rule prepare_all:
    input:
        mcfg.get_output_files(rules.prepare.output),
    localrule: True


integration_run_pattern = 'run_method/' + paramspace.wildcard_pattern.replace('--output_type~{output_type}', '')

use rule run_method from integration as integration_run_method with:
    message:
       """
       Integration: Run integration method {wildcards.method} on dataset={wildcards.dataset} file_id={wildcards.file_id}
       input: {input}
       output: {output}
       wildcards: {wildcards}
       resources: gpu={resources.gpu} mem_mb={resources.mem_mb} partition={resources.partition} qos={resources.qos}
       threads: {threads}
       """
    input:
        zarr=rules.prepare.output.zarr,
    output:
        zarr=directory(out_dir / integration_run_pattern / 'adata.zarr'),
        model=touch(directory(out_dir / integration_run_pattern / 'model')),
        plots=touch(directory(out_dir / integration_run_pattern / 'plots')),
    benchmark:
        out_dir / integration_run_pattern / 'benchmark.tsv'
    params:
        norm_counts=lambda wildcards: mcfg.get_from_parameters(wildcards, 'norm_counts', exclude=['output_type']),
        raw_counts=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_counts', exclude=['output_type']),
        output_type=lambda wildcards: mcfg.get_from_parameters(wildcards, 'output_types', exclude=['output_type']),
        hyperparams=lambda wildcards: mcfg.get_from_parameters(wildcards, 'hyperparams_dict', exclude=['output_type']),
        seed=lambda wildcards: mcfg.get_from_parameters(wildcards, 'seed', exclude=['output_type'], default=0),
        env=lambda wildcards: mcfg.get_from_parameters(wildcards, 'env', exclude=['output_type']),
    threads:
        lambda wildcards, attempt: get_threads(wildcards, attempt)
    resources:
        partition=lambda w, attempt: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        qos=lambda w, attempt: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu, factor=3),
        gpu=lambda w, attempt: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        time="2-00:00:00",


def update_neighbors_args(wildcards):
    args = mcfg.get_for_dataset(
        dataset=wildcards.dataset,
        query=['integration', 'neighbors'],
        default={}
    ).copy()
    output_type = wildcards.output_type
    if output_type == 'full':
        args |= {'use_rep': 'X_pca'}
    elif output_type == 'embed':
        args |= {'use_rep': 'X_emb'}
    elif output_type == 'knn':
        args = False
    else:
        raise ValueError(f'Unknown output_type {output_type}')
    return args


use rule neighbors from preprocessing as integration_postprocess with:
    input:
        zarr=rules.integration_run_method.output.zarr,
        done=rules.integration_run_method.output.model,
    output:
        zarr=directory(out_dir / 'postprocess' / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        args=update_neighbors_args,
        extra_uns=lambda wildcards: {'output_type': wildcards.output_type},
    retries: 2
    resources:
        partition=lambda w, attempt: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        qos=lambda w, attempt: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu, factor=3),
        gpu=lambda w, attempt: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=attempt_to_cpu),


rule run_all:
    input:
        mcfg.get_output_files(rules.integration_run_method.output),
        mcfg.get_output_files(rules.integration_postprocess.output),
    localrule: True
