
rule scimilarity:
    input:
        unpack(lambda wildcards: get_input(wildcards, model_name='scimilarity'))
    output:
        zarr=directory(mcfg.out_dir / 'scimilarity' / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'scimilarity', default={}).get('layer', 'X'),
        var_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'scimilarity', default={}).get('var_key'),
        model_params=lambda wildcards: mcfg.get_from_parameters(wildcards, 'scimilarity', default={}).get('model_params', {}),
    threads:
        lambda wildcards: mcfg.get_from_parameters(wildcards, 'threads', default=10, as_type=int)
    conda:
        get_env(config, 'scimilarity', env_dir='envs')
    resources:
        partition=mcfg.get_resource(profile='gpu',resource_key='partition'),
        gpu=mcfg.get_resource(profile='gpu',resource_key='gpu'),
        qos=mcfg.get_resource(profile='gpu',resource_key='qos'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt,factor=1),
    script:
        '../scripts/scimilarity.py'


rule scimilarity_all:
    input:
        [
            mcfg.get_output_files(rules.scimilarity.output, subset_dict=dict(dataset=dataset))
            for dataset in mcfg.get_datasets()
            if mcfg.get_from_parameters(dict(dataset=dataset), 'scimilarity', default={}).get('model') is not None
        ]
    localrule: True