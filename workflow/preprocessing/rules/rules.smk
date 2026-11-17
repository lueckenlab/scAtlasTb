"""
Preprocessing steps
Only the unique outputs per step are saved for storage efficiency. For assembled zarr files, see `assemble.smk`.
"""
from utils.environments import get_env


rule normalize:
    input: '{dataset}.h5ad'
    output:
        zarr=directory('{dataset}_normalized.zarr'),
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=not config.get("use_gpu", False))
    script:
        '../scripts/normalize.py'


rule filter_genes:
    input: '{dataset}.h5ad'
    output:
        zarr=directory('{dataset}_normalized.zarr'),
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/filter_genes.py'


rule highly_variable_genes:
    input:
        zarr='{dataset}.h5ad',
    output:
        zarr=directory('{dataset}_highly_variable_genes.zarr')
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=not config.get("use_gpu", False))
    script:
        '../scripts/highly_variable_genes.py'


rule extra_hvgs:
    input:
        zarr='{dataset}.h5ad',
    output:
        zarr=directory('{dataset}_extra_hvgs.zarr')
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=not config.get("use_gpu", False))
    script:
        '../scripts/extra_hvgs.py'


rule pca:
    input:
        zarr='{dataset}.h5ad',
    output:
        zarr=directory('{dataset}_pca.zarr')
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=not config.get("use_gpu", False))
    script:
        '../scripts/pca.py'


rule neighbors:
    input:
        zarr='{dataset}.h5ad'
    output:
        zarr=directory('{dataset}_neighbors.zarr')
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=not config.get("use_gpu", False))
    script:
        '../scripts/neighbors.py'


rule umap:
    input:
        zarr='{dataset}.h5ad',
        rep='{dataset}.h5ad',
    output:
        zarr=directory('{dataset}_umap.zarr')
    conda:
        get_env(config, 'scanpy', gpu_env='rapids_singlecell', no_gpu=not config.get("use_gpu", False))
    script:
        '../scripts/umap.py'


### Assemble ###

rule assemble:
    input:
        unpack(
            dict(
                counts='{dataset}.h5ad',
                normalize='{dataset}_normalized.zarr',
                highly_variable_genes='{dataset}_highly_variable_genes.zarr',
                pca='{dataset}_pca.zarr',
                neighbors='{dataset}_neighbors.zarr',
                umap='{dataset}_umap.zarr',
            )
        )
    output:
        zarr=directory('{dataset}_preprocessed.zarr')
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/assemble.py'


### Plots ###

rule plots:
    input:
        anndata='{dataset}.h5ad',
    output:
        plots=directory('{dataset}_plot'),
    params:
        basis='X_pca'
    threads: 10
    conda:
        get_env(config, 'scanpy')
    script:
        '../scripts/plot.py'
