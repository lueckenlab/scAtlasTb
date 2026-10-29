"""
Assemble anndata with different Preprocessing outputs
"""

use rule normalize from preprocessing as preprocessing_normalize with:
    input:
        lambda wildcards: mcfg.get_input_file(wildcards.dataset, wildcards.file_id)
    output:
        zarr=directory(mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / 'normalized.zarr'),
    params:
        raw_counts=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_counts'),
        gene_id_column=lambda wildcards: mcfg.get_from_parameters(wildcards, 'gene_id_column'),
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'normalize', default={}),
        dask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'dask', default=True),
    conda:
        lambda w: get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=mcfg.get_profile(w) == 'cpu')
    threads:
        lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_threads', default=4)
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w),resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w),resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w),resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w),resource_key='mem_mb',attempt=attempt, factor=1),
    benchmark:
        mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / '.benchmark_normalized.tsv'


use rule filter_genes from preprocessing as preprocessing_filter_genes with:
    input:
        zarr=rules.preprocessing_normalize.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / 'filtered_genes.zarr'),
    params:
        dask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'dask', default=True),
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'filter', default=dict(min_cells=1)),
    threads:
        lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_threads', default=4)
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb',attempt=attempt, factor=1),
    benchmark:
        mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / '.benchmark_filtered_genes.tsv'

use rule highly_variable_genes from preprocessing as preprocessing_highly_variable_genes with:
    input:
        zarr=rules.filter_genes.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / 'highly_variable_genes--{hvg_args}.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'hvg_args_dict', default={}),
        dask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'dask', default=False),
    conda:
        lambda w: get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=mcfg.get_profile(w) == 'cpu')
    threads:
        lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_threads', default=4)
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w),resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w),resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w),resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w),resource_key='mem_mb',attempt=attempt, factor=1),
    benchmark:
        mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / '.benchmark_hvg--{hvg_args}.tsv'


use rule extra_hvgs from preprocessing as preprocessing_extra_hvgs with:
    input:
        zarr=rules.filter_genes.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / 'extra_hvgs--{overwrite_args}.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'highly_variable_genes', default={}),
        extra_hvgs=lambda wildcards: mcfg.get_from_parameters(wildcards, 'extra_hvgs', default={}),
        overwrite_args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'overwrite_args_dict', default={}),
        dask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'dask', default=False),
    conda:
        lambda w: get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=mcfg.get_profile(w) == 'cpu')
    threads:
        lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_threads', default=4)
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w),resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w),resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w),resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w),resource_key='mem_mb',attempt=attempt, factor=1),
    benchmark:
        mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / '.benchmark_extra_hvgs--{overwrite_args}.tsv'


use rule pca from preprocessing as preprocessing_pca with:
    input:
        zarr=lambda wildcards: mcfg.get_output_files(
            rules.preprocessing_highly_variable_genes.output.zarr,
            subset_dict=wildcards,
            wildcard_names=['dataset', 'file_id', 'hvg_args'],
        )[0],
    output:
        zarr=directory(mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / 'pca.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'pca', default={}),
        scale=lambda wildcards: mcfg.get_from_parameters(wildcards, 'scale', default=False),
        dask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'dask', default=True),
    conda:
        lambda w: get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=mcfg.get_profile(w, default='gpu') == 'cpu')
    threads:
        lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_threads', default=4)
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w, default='gpu'),resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w, default='gpu'),resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w, default='gpu'),resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w, default='gpu'),resource_key='mem_mb',attempt=attempt, factor=1),
    benchmark:
        mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / '.benchmark_pca.tsv'


use rule neighbors from preprocessing as preprocessing_neighbors with:
    input:
        zarr=rules.preprocessing_pca.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / 'neighbors.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'neighbors', default={}),
    threads:
        lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_threads', default=4)
    conda:
        lambda w: get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=mcfg.get_profile(w, default='gpu') == 'cpu')
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w, default='gpu'),resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w, default='gpu'),resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w, default='gpu'),resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w, default='gpu'),resource_key='mem_mb',attempt=attempt),
    benchmark:
        mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / '.benchmark_neighbors.tsv'


use rule umap from preprocessing as preprocessing_umap with:
    input:
        zarr=rules.preprocessing_neighbors.output.zarr,
        rep=rules.preprocessing_pca.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / 'umap.zarr'),
    params:
        args=lambda wildcards: mcfg.get_from_parameters(wildcards, 'umap', default={}),
    conda:
        lambda w: get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=mcfg.get_profile(w, default='gpu') == 'cpu')
    threads:
        lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_threads', default=4)
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w, default='gpu'),resource_key='partition',attempt=attempt),
        qos=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w, default='gpu'),resource_key='qos',attempt=attempt),
        gpu=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w, default='gpu'),resource_key='gpu',attempt=attempt),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile=mcfg.get_profile(w, default='gpu'),resource_key='mem_mb',attempt=attempt),
    benchmark:
        mcfg.out_dir / 'preprocessed' / paramspace.wildcard_pattern / '.benchmark_umap.tsv'


def collect_files(wildcards):
    # Configuration: step_name -> (rule_name, wildcard_names)
    base_wildcards = ['dataset', 'file_id']
    step_config = {
        'normalize': (rules.preprocessing_normalize, base_wildcards),
        'highly_variable_genes': (rules.preprocessing_highly_variable_genes, base_wildcards + ['hvg_args']),
        'extra_hvgs': (rules.preprocessing_extra_hvgs, base_wildcards + ['overwrite_args']),
        'pca': (rules.preprocessing_pca, base_wildcards),
        'neighbors': (rules.preprocessing_neighbors, base_wildcards),
        'umap': (rules.preprocessing_umap, base_wildcards),
    }
    
    # Build file dictionary
    file_dict = {
        step: mcfg.get_output_files(
            _rule.output.zarr,
            subset_dict=wildcards,
            wildcard_names=step_wildcards,
        )
        for step, (_rule, step_wildcards) in step_config.items()
    }
    
    # Filter by assembly configuration if present
    assembly_config = mcfg.get_from_parameters(wildcards, 'assemble')
    if assembly_config is not None:
        file_dict = {k: v for k, v in file_dict.items() if k in assembly_config}
    
    return file_dict


use rule assemble from preprocessing as preprocessing_assemble with:
    input:
        unpack(collect_files)
    output:
        zarr=directory(mcfg.out_dir / f'{paramspace.wildcard_pattern}.zarr'),
    resources:
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
    retries: 0
    conda:
        get_env(config, 'scanpy')


rule assemble_all:
    message:
        """
        Collect all preprocessing output for all datasets and methods
        """
    input:
        assemble_input=[
            file for dataset in mcfg.get_datasets()
            for file_list in collect_files(dict(dataset=dataset)).values()
            for file in file_list
        ],
        assemble_output=mcfg.get_output_files(rules.preprocessing_assemble.output.zarr),
    localrule: True