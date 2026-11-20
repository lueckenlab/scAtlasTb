rule download:
    message:
        """
        Downloading file for {wildcards}
        """
    output:
        h5ad=out_dir / 'download' / '{dataset}.h5ad'
    params:
        dataset_df=lambda w: dataset_df.query(f'dataset == "{w.dataset}"')[DOWNLOAD_COLUMNS].reset_index(drop=True)
    conda:
        get_env(config, 'scanpy', env_dir='../../../envs')
    script:
        '../scripts/download.py'


rule download_all:
    input:
        expand(
            rules.download.output,
            dataset=dataset_df[dataset_df['url'].apply(lambda x: not Path(x).is_file())] \
                .to_dict(orient='list').get('dataset')
        )
    localrule: True


def get_files_for_metadata_harm(wildcards):
    meta = unlist_dict(
            get_wildcards(dataset_df, columns=all_but(dataset_df.columns,'subset'), wildcards=wildcards)
        )
    file_path = meta['url']
    anno_file = meta.get('annotation_file')
    schema_file = config.get('schema_file')

    # download file if not present
    if not Path(str(file_path)).exists():
        file_path = rules.download.output.h5ad

    file_map = {'h5ad': file_path}

    # include annotation file if specified
    if pd.notna(anno_file) and pd.notnull(anno_file):
        file_map['annotation_file'] = anno_file

    # include schema file if specified
    if pd.notna(schema_file) and pd.notnull(schema_file):
        file_map['schema'] = schema_file

    return file_map


rule harmonize_metadata:
    input:
        unpack(get_files_for_metadata_harm)
    params:
        meta=lambda wildcards: unlist_dict(
            get_wildcards(dataset_df, columns=all_but(dataset_df.columns,'subset'), wildcards=wildcards)
        ),
        dask=True,
    output:
        zarr=directory(out_dir / 'harmonize_metadata' / '{dataset}.zarr'),
        # plot=image_dir / 'harmonize_metadata' / 'counts_sanity--{dataset}.png',
    conda:
        get_env(config, 'scanpy', env_dir='../../../envs')
    threads: 5
    resources:
        mem_mb=lambda wildcards, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),
    # shadow: 'shallow'
    script:
        '../scripts/harmonize_metadata.py'


rule harmonize_metadata_all:
    input:
        expand(rules.harmonize_metadata.output,dataset=dataset_df['dataset'])
    localrule: True


use rule merge from load_data as load_data_merge_study with:
    input:
        lambda wildcards: expand(
            rules.harmonize_metadata.output.zarr,
            dataset=get_wildcards(dataset_df,'dataset',wildcards)['dataset']
        )
    output:
        zarr=directory(out_dir / 'merged' / 'study' / '{study}.zarr'),
    params:
        dataset=lambda wildcards: wildcards.study,
        merge_strategy='inner',
        keep_all_columns=True,
        backed=False,
        dask=False,
    threads: 5
    resources:
        mem_mb=lambda wildcards, attempt: get_resource(config,profile='cpu',resource_key='mem_mb', attempt=attempt),


rule merge_study_all:
    input: expand(rules.load_data_merge_study.output,**get_wildcards(dataset_df,['study']))
    localrule: True
