checkpoint split_batches:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards)
    output:
        batches=directory(mcfg.out_dir / 'scatter' /  params.wildcard_pattern / '_batches'),
    params:
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch_key'),
        chunk_size=lambda wildcards: mcfg.get_from_parameters(wildcards, 'chunk_size', default=100_000),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/split_batches.py'


def get_checkpoint_output(wildcards):
    return f'{checkpoints.split_batches.get(**wildcards).output[0]}/{{batch}}.txt'


def get_from_checkpoint(wildcards, pattern=None):
    checkpoint_output = get_checkpoint_output(wildcards)
    if pattern is None:
        pattern = checkpoint_output
    return expand(
        pattern,
        batch=glob_wildcards(checkpoint_output).batch,
        allow_missing=True
    )


def get_mem_mb(attempt, profile, factor=2):
    mem_mb = mcfg.get_resource(profile=profile, resource_key='mem_mb', attempt=attempt)
    try:
        mem_mb = int(mem_mb)
    except ValueError:
        return mem_mb
    return max(1_000, int(mem_mb // factor))


def get_profile_from_use_gpu(wildcards):
    return 'gpu' if mcfg.get_from_parameters(wildcards, 'use_gpu', default=False) else 'cpu'


rule scrublet:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        batch=get_checkpoint_output,
    output:
        tsv=mcfg.out_dir / 'scatter' / params.wildcard_pattern / 'scrublet' / '{batch}.tsv',
    params:
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch_key', check_query_keys=False),
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'counts', default='X'),
    conda:
        lambda wildcards: get_env(
            config,
            'qc',
            gpu_env='rapids_singlecell',
            no_gpu=not mcfg.get_from_parameters(wildcards, 'use_gpu', default=False)
        )
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile=get_profile_from_use_gpu(w), resource_key='partition', attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile=get_profile_from_use_gpu(w), resource_key='qos', attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile=get_profile_from_use_gpu(w), resource_key='gpu', attempt=attempt),
        mem_mb=lambda w, attempt: get_mem_mb(attempt=attempt, profile=get_profile_from_use_gpu(w)),
    script:
        '../scripts/scrublet.py'


rule doubletdetection:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        batch=get_checkpoint_output
    output:
        tsv=mcfg.out_dir / 'scatter' / params.wildcard_pattern / 'doubletdetection' / '{batch}.tsv',
    params:
        batch_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'batch_key', check_query_keys=False),
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'counts', default='X'),
    conda:
        get_env(config, 'qc')
    threads: 3
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
        mem_mb=lambda w, attempt: get_mem_mb(attempt, profile='cpu'),
    # shadow: "minimal"
    script:
        '../scripts/doubletdetection.py'


def collect_results(wildcards):
    methods = mcfg.get_from_parameters(wildcards, 'methods', default=['scrublet'])
    files = {'zarr': mcfg.get_input_file(**wildcards)}
    if 'scrublet' in methods:
        files['scrublet'] = get_from_checkpoint(wildcards, rules.scrublet.output.tsv)
    if 'doubletdetection' in methods:
        files['doubletdetection'] = get_from_checkpoint(wildcards, rules.doubletdetection.output.tsv)
    return files


rule collect:
    input:
        unpack(collect_results)
        # zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        # scrublet=lambda wildcards: get_from_checkpoint(wildcards, rules.scrublet.output.tsv),
        # doubletdetection=lambda wildcards: get_from_checkpoint(wildcards, rules.doubletdetection.output.tsv),
    output:
        zarr=directory(mcfg.out_dir / f'{params.wildcard_pattern}.zarr'),
    params:
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'counts', default='X'),
    localrule: True
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb')
    script:
        '../scripts/collect.py'
