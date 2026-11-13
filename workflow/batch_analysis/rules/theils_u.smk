rule theils_u:
    input:
        zarr=lambda wildcards: get_file(wildcards, 'data'),
    output:
        tsv=mcfg.out_dir / paramspace.wildcard_pattern / 'theils_u.tsv',
        plot=mcfg.image_dir / paramspace.wildcard_pattern / 'theils_u_heatmap.svg',
    params:
        covariates=lambda wildcards: mcfg.get_from_parameters(wildcards, 'covariates', default=[]),
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample', check_query_keys=False),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/theils_u.py'


rule theils_u_all:
    input:
        mcfg.get_output_files(rules.theils_u.output),
    localrule: True