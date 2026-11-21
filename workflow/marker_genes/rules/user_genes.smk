def get_plot_user_input(wildcards):
    files = dict(zarr=mcfg.get_input_file(**wildcards))
    # add group_order file as input, so that changes to the file trigger reruns
    group_order = mcfg.get_from_parameters(wildcards, 'plot', default={}).get('group_order')
    if isinstance(group_order, dict):
        group_order = group_order.get(wildcards.group)
    if isinstance(group_order, str) and Path(group_order).is_file():
        files['group_order'] = group_order
    return files


rule plot_user:
    input:
        unpack(get_plot_user_input)
    output:
        plots=directory(mcfg.image_dir / paramspace.wildcard_pattern / 'group={group}' / 'user_markers'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'plot', default={}),
        layer=lambda wildcards: mcfg.get_from_parameters(wildcards, 'layer', default='X'),
        markers=lambda wildcards: get_marker_gene_set(mcfg, wildcards),
    conda:
        get_env(config, 'scanpy')
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt),
    script:
        '../scripts/plot_user.py'

