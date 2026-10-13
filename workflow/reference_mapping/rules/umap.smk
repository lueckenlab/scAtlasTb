def update_neighbors_args(wildcards):
    args = mcfg.get_for_dataset(
        dataset=wildcards.dataset,
        query=['reference_mapping', 'neighbors'],
        default={}
    ).copy()
    args['use_rep'] = 'X_emb'
    return args


def get_reference_mapped(wildcards):
    if mcfg.get_from_parameters(wildcards, 'scarches', default={}).get('model') is not None:
        return rules.scarches.output.zarr
    elif mcfg.get_from_parameters(wildcards, 'scimilarity', default={}).get('model') is not None:
        return rules.scimilarity.output.zarr
    else:
        raise ValueError(f"No reference mapping method specified for dataset {wildcards.dataset}")


use rule neighbors from preprocessing as reference_mapping_neighbors with:
    input:
        zarr=get_reference_mapped,
    output:
        zarr=directory(mcfg.out_dir / 'postprocess' / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        args=update_neighbors_args,
        extra_uns=lambda wildcards: {'output_type': 'embed'},
    retries: 2
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt),



use rule umap from preprocessing as reference_mapping_compute_umap with:
    input:
        anndata=rules.reference_mapping_neighbors.output.zarr,
        rep=rules.reference_mapping_neighbors.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / f'{paramspace.wildcard_pattern}.zarr'),
    threads: 10
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt, attempt_to_cpu=2),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt, attempt_to_cpu=2),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt, attempt_to_cpu=2),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt, attempt_to_cpu=2),


use rule plots from preprocessing as reference_mapping_plot_umap with:
    input:
        anndata=rules.reference_mapping_compute_umap.output.zarr,
    output:
        plots=directory(mcfg.image_dir / 'umap' / paramspace.wildcard_pattern),
    params:
        color=lambda wildcards: mcfg.get_for_dataset(wildcards.dataset, query=[mcfg.module_name, 'umap_colors'], default=[]),
        basis='X_umap',
        min_cells_per_category=0.0001,
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
