# Single Cell Atlasing Toolbox 🧰

[![Documentation][badge-docs]][documentation]
[![Tests][badge-tests]][tests]

[badge-tests]: https://img.shields.io/github/actions/workflow/status/HCA-integration/scAtlasTb/test.yaml?label=tests
[badge-docs]: https://readthedocs.org/projects/scatlastb/badge/?version=latest

**Toolbox of Snakemake pipelines for easy-to-use analyses and benchmarks for building integrated atlases**

This toolbox provides multiple modules that can be easily combined into custom workflows that leverage the file management of [Snakemake](https://snakemake.readthedocs.io/en/v7.31.1/).
This allows for an efficient and scalable way to run analyses on large datasets that can be easily configured by the user.

## Getting started

Please refer to the [documentation][] for a detailed guide. In short:

### System requirements

* Linux (preferred) or MacOS, with [conda](https://github.com/conda-forge/miniforge) and git installed
* All software dependencies are installed via the conda environments in `envs/`
* The demo runs on a standard desktop computer (CPU only). An NVIDIA GPU is optional and only needed for deep-learning-based and GPU-accelerated methods.

See [system requirements](https://scatlastb.readthedocs.io/en/latest/getting_started/installation.html#system-requirements) for details.

### Installation

```
git clone https://github.com/HCA-integration/scAtlasTb.git
cd scAtlasTb
for env in snakemake scanpy bbknn scanorama scib plots funkyheatmap; do
    bash envs/install_environment.sh -f envs/$env.yaml
done
conda run -n funkyheatmap Rscript -e "install.packages('funkyheatmap', repos='https://cloud.r-project.org')"
```

This installs the environments needed for the demo. Typical install time: ~10 minutes (~8 GB disk space).
On Apple Silicon, install the `funkyheatmap` environment with `CONDA_SUBDIR=osx-64` (see [troubleshooting](https://scatlastb.readthedocs.io/en/latest/principles/troubleshooting.html#working-on-apple-silicon)).
See the [installation guide](https://scatlastb.readthedocs.io/en/latest/getting_started/installation.html) for installing all environments.

### Demo

The repository ships a small demo dataset (`data/pbmc68k.h5ad`, 700 cells with 3 simulated batches).
Run preprocessing, data integration and integration benchmarking on it with:

```
conda activate snakemake
bash run_example.sh preprocessing_all integration_all metrics_all -c 4
```

Expected run time: ~15 minutes with 4 CPU cores (measured on a laptop with Apple M1 Pro and 16 GB RAM).
The [quickstart](https://scatlastb.readthedocs.io/en/latest/getting_started/quickstart.html) describes the expected output.

### Usage on your own data

Write a configuration file for your data following the [configuration guide](https://scatlastb.readthedocs.io/en/latest/getting_started/configure_workflow.html) and [call the pipeline](https://scatlastb.readthedocs.io/en/latest/getting_started/call_pipeline.html) with it.
Each module is described in detail in the [module documentation](https://scatlastb.readthedocs.io/en/latest/modules/index_data_preparation.html).


## 🧰 Which Modules does the Toolbox Support?

The modules are located under `workflow/` and can be run independently or combined into a more complex workflow.

Helper modules (`subset`, `split_data`, `relabel`, `collect`, `uncollect`) can be inserted at any point to reshape inputs and outputs.

<details>
<summary><b>Click to expand the full list of modules</b></summary>

| Module                 | Description                                                               |
|------------------------|---------------------------------------------------------------------------|
| `load_data`            | Loading datasets from URLs and converting them to AnnData objects         |
| `exploration`          | Exploration and quality control of datasets                               |
| `batch_analysis`       | Exploration and quality control of batches within datasets                |
| `qc`                   | Semi-automated quality control of datasets using [sctk AutoQC](https://teichlab.github.io/sctk/notebooks/automatic_qc.html) |
| `doublets`             | Identifying and handling doublets in datasets                             |
| `merge`                | Merging datasets                                                          |
| `filter`               | Filtering datasets based on specified criteria                            |
| `subset`               | Creating subsets of datasets                                              |
| `relabel`              | Relabeling data points in datasets                                        |
| `split_data`           | Splitting datasets into training and testing sets                         |
| `preprocessing`        | Preprocessing of datasets (normalization, feature selection, PCA, kNN graph, UMAP) |
| `integration`          | Running single cell batch correction methods on datasets                  |
| `metrics`              | Calculating [scIB metrics](https://scib.readthedocs.io/en/latest/), mainly for benchmarking of integration methods  |
| `clustering`           | Multi-resolution and hierarchical clustering of datasets                  |
| `label_harmonization`  | Providing alignment between unharmonized labels using [CellHint](https://cellhint.readthedocs.io/en/latest) |
| `label_transfer`       | Transfer annotations of annotated cells to unannotated cells              |
| `majority_voting`      | Consensus voting across multiple cell type assignments                    |
| `celltype_prediction`  | Predict cell types from reference model e.g. celltypist                   |
| `reference_mapping`    | Map query datasets to reference atlases                                   |
| `marker_genes`         | Identify marker genes for cell types                                      |
| `collect`              | Collect multiple input anndata objects into a single anndata object       |
| `uncollect`            | Distribute slots of an anndata object to multiple anndata objects         |
| `common`               | Common utilities and helper functions for workflows                       |

</details>

## 👀 TL;DR What does a full workflow look like?

The heart of the configuration is captured in a YAML (or JSON) configuration file.
Here is a simplified version of the demo configuration in `configs/quickstart.yaml` containing the `preprocessing`, `integration` and `metrics` modules:

```yaml
output_dir: data/out
images: images

use_gpu: false

DATASETS:

  my_task: # custom task/workflow name

    # input specification: map of module name to map of input file name to input file path
    input:
      preprocessing:
        file_1: data/pbmc68k.h5ad
        # file_2: ... # more files if required
      integration: preprocessing # all outputs of module will automatically be used as input
      metrics: integration

    # module configuration
    preprocessing:
      highly_variable_genes:
        n_top_genes: 2000
      pca:
        n_comps: 50
      assemble:
        - normalize
        - highly_variable_genes
        - pca

    # module configuration
    integration:
      raw_counts: layers/counts
      norm_counts: layers/normcounts
      batch: batch
      methods:
        unintegrated:
        harmonypy:
        scanorama:
          batch_size: 100

    # module configuration
    metrics:
      unintegrated: layers/normcounts
      batch: batch
      label: bulk_labels
      metrics:
        - nmi
        - graph_connectivity
```

Which allows you to call the pipeline as follows (after saving the config to `my_config.yaml`):

```
snakemake --configfile my_config.yaml --snakefile workflow/Snakefile --use-conda preprocessing_all integration_all metrics_all -nq
```

giving you the following dryrun output:

```
Job stats:
job                                    count
-----------------------------------  -------
integration_all                            1
integration_barplot_per_dataset            3
integration_benchmark_per_dataset          1
integration_compute_umap                   6
integration_plot_umap                      6
integration_postprocess                    6
integration_prepare                        1
integration_run_method                     3
metrics_all                                1
metrics_barplot                            3
metrics_barplot_per_dataset                3
metrics_cluster                           60
metrics_cluster_collect                    6
metrics_collect                            6
metrics_funkyheatmap                       1
metrics_funkyheatmap_per_dataset           1
metrics_merge                              1
metrics_merge_per_batch                    1
metrics_merge_per_dataset                  1
metrics_merge_per_file                     6
metrics_merge_per_label                    1
metrics_prepare                            6
metrics_run                               12
preprocessing_all                          1
preprocessing_assemble                     1
preprocessing_filter_genes                 1
preprocessing_highly_variable_genes        1
preprocessing_neighbors                    1
preprocessing_normalize                    1
preprocessing_pca                          1
preprocessing_plot_pca                     1
preprocessing_plot_umap                    1
preprocessing_umap                         1
total                                    146
```

💖 Beautiful, right? Check out the [documentation][] to learn how to set up your own workflow!

## Release notes

See the [changelog][].

## Contact

If you found a bug, please use the [issue tracker][].

## Citation

```bibtex
@article {Mueller2026.07.30.741695,
	author = {Mueller, Michaela F. and Cujba, Ana-Maria and Romanovskaia, Daria and Cohen, Carla J. and Bright, Chelsea A. and Lance, Christopher and Ram{\'\i}rez-Su{\'a}stegui, Ciro and Strobl, Daniel C. and Yuan, Hao and Hulsen, Janneke and Naas, Julia and Limbeck, Katharina and Kock, Kian Hong and Halle, Lennard and Knoll, Rainer and Kfuri-Rubens, Raphael and Aguilar-Fern{\'a}ndez, Sergio and Parikh, Shrey and Shitov, Vladimir A. and Said, Wamia and Kasper, Maria and Snelling, Sarah J. B. and Teichmann, Sarah A. and Reynolds, Gary and Prabhakar, Shyam and Villani, Alexandra-Chloe and Theis, Fabian J. and Luecken, Malte D.},
	title = {Building optimized single-cell reference atlases with scAtlasTb},
	year = {2026},
	doi = {10.64898/2026.07.30.741695},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2026/08/02/2026.07.30.741695},
	journal = {bioRxiv}
}
```

[issue tracker]: https://github.com/HCA-integration/scAtlasTb/issues
[tests]: https://github.com/HCA-integration/scAtlasTb/actions/workflows/test.yaml
[documentation]: https://scatlastb.readthedocs.io
[changelog]: https://scatlastb-utils.readthedocs.io/en/latest/changelog.html
