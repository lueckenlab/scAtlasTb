metric_wildcards = mcfg.get_wildcard_names() + ['metric', 'parse_file_id']

rule merge_per_file:
    message:
        """
        Merge all metrics for all datasets and methods
        {params.wildcards_string}
        """
    input:
        metrics=lambda wildcards: mcfg.get_output_files(rules.run.output, subset_dict=wildcards, all_params=True),
        benchmark=lambda wildcards: mcfg.get_output_files(rules.run.benchmark, subset_dict=wildcards, all_params=True),
    output:
        tsv=mcfg.out_dir / 'results' / paramspace.wildcard_pattern / 'metrics.tsv',
        extra_columns=mcfg.out_dir / 'results' / paramspace.wildcard_pattern / 'extra_columns.txt',
    params:
        wildcards=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude=['dataset', 'file_id'], as_df=True, wildcard_names=metric_wildcards),
        wildcards_string=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude=['dataset', 'file_id'], as_df=True, wildcard_names=metric_wildcards).to_string(index=False),
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=1000,
    retries: 0
    # group:
    #     'metrics_postprocess'
    localrule: True
    script: '../scripts/merge.py'


use rule merge_per_file as merge with:
    message:
        """
        Merge all metrics for {wildcards}
        {params.wildcards_string}
        """
    input:
        metrics=lambda wildcards: mcfg.get_output_files(rules.run.output, subset_dict=wildcards, all_params=True),
        benchmark=lambda wildcards: mcfg.get_output_files(rules.run.benchmark, subset_dict=wildcards, all_params=True),
    output:
        tsv=mcfg.out_dir / 'results' / 'metrics.tsv',
        extra_columns=mcfg.out_dir / 'results' / 'extra_columns.txt',
    params:
        wildcards=mcfg.get_wildcards(as_df=True, wildcard_names=metric_wildcards),
        wildcards_string=mcfg.get_wildcards(as_df=True, wildcard_names=metric_wildcards).to_string(index=False),
    # group:
    #     'metrics_postprocess'
    localrule: True


use rule merge_per_file as merge_per_dataset with:
    message:
        """
        Merge all metrics for {wildcards}
        {params.wildcards_string}
        """
    input:
        metrics=lambda wildcards: mcfg.get_output_files(rules.run.output, subset_dict=dict(wildcards), all_params=True),
        benchmark=lambda wildcards: mcfg.get_output_files(rules.run.benchmark, subset_dict=dict(wildcards), all_params=True),
    output:
        tsv=mcfg.out_dir / 'results' / 'dataset~{dataset}' / 'metrics.tsv',
        extra_columns=mcfg.out_dir / 'results' / 'dataset~{dataset}' / 'extra_columns.txt',
    params:
        wildcards=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='dataset', as_df=True, wildcard_names=metric_wildcards),
        wildcards_string=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='dataset', as_df=True, wildcard_names=metric_wildcards).to_string(index=False),
    # group:
    #     'metrics_postprocess'
    localrule: True


use rule merge_per_file as merge_per_batch with:
    message:
        """
        Merge all metrics for {wildcards}
        {params.wildcards_string}
        """
    input:
        metrics=lambda wildcards: mcfg.get_output_files(rules.run.output, subset_dict=dict(wildcards), all_params=True),
        benchmark=lambda wildcards: mcfg.get_output_files(rules.run.benchmark, subset_dict=dict(wildcards), all_params=True),
    output:
        tsv=mcfg.out_dir / 'results' / 'dataset~{dataset}' / 'batch~{batch}' / 'metrics.tsv',
        extra_columns=mcfg.out_dir / 'results' / 'dataset~{dataset}' / 'batch~{batch}' / 'extra_columns.txt',
    params:
        wildcards=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='batch', as_df=True, wildcard_names=metric_wildcards),
        wildcards_string=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='batch', as_df=True, wildcard_names=metric_wildcards).to_string(index=False),
    # group:
    #     'metrics_postprocess'
    localrule: True


use rule merge_per_file as merge_per_label with:
    message:
        """
        Merge all metrics for {wildcards}
        {params.wildcards_string}
        """
    input:
        metrics=lambda wildcards: mcfg.get_output_files(rules.run.output, subset_dict=dict(wildcards), all_params=True),
        benchmark=lambda wildcards: mcfg.get_output_files(rules.run.benchmark, subset_dict=dict(wildcards), all_params=True),
    output:
        tsv=mcfg.out_dir / 'results' / 'dataset~{dataset}' / 'label~{label}' / 'metrics.tsv',
        extra_columns=mcfg.out_dir / 'results' / 'dataset~{dataset}' / 'label~{label}' / 'extra_columns.txt',
    params:
        wildcards=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='label', as_df=True, wildcard_names=metric_wildcards),
        wildcards_string=lambda wildcards: mcfg.get_wildcards(subset_dict=wildcards, exclude='label', as_df=True, wildcard_names=metric_wildcards).to_string(index=False),
    # group:
    #     'metrics_postprocess'
    localrule: True



def get_collect_input(wildcards):
    files = dict(
        zarr=mcfg.get_input_file(**wildcards),
        metrics=rules.merge_per_file.output.tsv,
    )

    if sum(mcfg.get_from_parameters(wildcards, 'use_gene_set', single_value=False)) > 0:
        gene_scores = mcfg.get_output_files(rules.score_genes.output, subset_dict=wildcards)
        files['gene_scores'] = gene_scores
    
    return files


rule collect:
    message:
        """
        Collect all metrics for {wildcards}
        input: {input}
        """
    input:
        unpack(get_collect_input)
    output:
        zarr=directory(mcfg.out_dir / f'{paramspace.wildcard_pattern}.zarr'),
    conda:
        get_env(config, 'scanpy')
    localrule: True
    script:
        '../scripts/collect.py'


rule collect_all:
    message:
        """
        Collect all metrics for all datasets and methods
        """
    input:
        mcfg.get_output_files(rules.collect.output),
    localrule: True


rule merge_all:
    input:
        # rules.merge.output,
        mcfg.get_output_files(rules.merge_per_dataset.output, wildcard_names=['dataset']),
        mcfg.get_output_files(rules.merge_per_batch.output, wildcard_names=['dataset', 'batch']),
        mcfg.get_output_files(rules.merge_per_label.output, wildcard_names=['dataset', 'label']),
        # mcfg.get_output_files(rules.merge_per_file.output, wildcard_names=['file_id']),
    localrule: True