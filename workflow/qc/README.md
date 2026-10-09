# QC
Quality control of single cell datasets with [sctk autoqc](https://teichlab.github.io/sctk/notebooks/automatic_qc.html) thresholds.
This module provides workflows that quantify the quality of single cell datasets and provide automated thresholds for filtering out low quality cells.
QC metrics and thresholds are visualised in different plots and tables.

## Input

### AnnData file
AnnData file in h5ad or zarr format with the following:

+ raw count matrix under adata.X (or the slot given by `counts`, e.g. `layers/counts`)
+ any categorical or continuous metadata under adata.obs


### Configuration

```yaml
out_dir: test/out
images: test/images

DATASETS:
  dataset_name:
    input:
      qc:
        test: test/input/pbmc68k.h5ad
        test2: test/input/pbmc68k.h5ad
    qc:
      hue:  # colors for colouring joint plots and stratify removed cellls plots
        - phase
        - bulk_labels
      thresholds:  # mapping of file_id -> custom thresholds to overwrite autoqc thresholds
        test:
          n_counts_max: 500
          n_genes_min: 100
          percent_mito_max: 0.5
      alternative_thresholds: # stricter thresholds to flag ambiguous cells
        test:
          n_counts_max: 1000
          n_genes_min: 200
          percent_mito_max: 0.2
      thresholds_file: test/input/user_thresholds.tsv # path to a file with custom thresholds, will overwrite thresholds above
```

TODO: explain sctk autoqc threshold naming


### Defining custom thresholds

You can define custom thresholds on a per file_id basis that will overwrite the autoqc thresholds. There are two ways of defining custom thresholds. You can either make use of the `thresholds` section of the configuration file or create a TSV file that allow more customisation as well as alternative thresholds. When specifying thresholds directly in the yaml, you need to provide a mapping of file_id to a dictionary of thresholds (see example above).

If you wnat to define the thresholds via TSV file, these values will overwrite the user-specified in the config file. Providing a file is useful if you want to customise many file_ids. The TSV file should have the following columns:

* file_id (mandatory): the file_id of the file for which the thresholds are defined
* threshold_type: the type of user input. Internally, the table will be subset to 'user' and 'alternative', which will be passed on accordingly.
* n_counts_max: the maximum number of counts per cell
* n_counts_min: the minimum number of counts per cell
* n_genes_max: the maximum number of genes per cell
* n_genes_min: the minimum number of genes per cell
* percent_mito_max: the maximum percentage of mitochondrial genes per cell
* percent_mito_min: the minimum percentage of mitochondrial genes per cell


### Alternative thresholds

Alternative thresholds are used to define ambiguous cells. This is useful if you are not sure about which threshold to choose and allows you to flag cells that would have been removed by a stricter threshold. These annotations can be useful for downstream analysis.

You can specify alternative thresholds the same way as user-provided thresholds under `alternattive_thresholds` in the config or via the `thresholds_file` with `threshold_type` == 'alternative'.


## Output

The output of the QC workflow is a set of plots and tables that visualise the quality of the dataset and the thresholds used for filtering cells. The output plots are saved under the images directory from the configuration file.

### Outputs include:

* Joint QC plots under `<images_dir>/dataset~{dataset}/file_id~{file_id}/joint_plots/`
  * `main.svg`: scatter plots for all configured QC metric pairs without hue stratification (only created if no valid `hue` is configured)
  * `hue={hue}.svg`: scatter plots for all configured QC metric pairs colored by each configured hue and by `qc_status`
  * `density.svg`: one combined density figure covering all configured QC metric pairs, with regular-scale and log-scale panels side by side (only created with `plot_params: {plot_density: true}`)
* QC removal plots under `<images_dir>/dataset~{dataset}/file_id~{file_id}/removed/`
  * `by={hue}.svg`: removed-cell summaries stratified by each configured hue
  * `cells_passed_all.svg`: overall counts of passed, failed, and ambiguous cells
  * `per_metric_violin.svg`: per-metric QC distributions split by QC status
* Table of thresholds used for filtering cells at the dataset and file level
* Thresholds TSV containing the autoQC, user-provided, alternative and updated thresholds
  * Updated thresholds are based on the autoQC thresholds, overwritten by the user-provided thresholds (alternative thresholds are not considered in this case)
  * The file contains number and fraction of kept and removed cells as well as the corresponding QC metrics on a per file_id basis
  * The combined thresholds are stored under `<images_dir>/dataset~{dataset}/thresholds.tsv`
* QC status summary table
  * Table that contains the number of cells of the QC status categories 'passed', 'failed', 'ambiguous' on a per file_id basis
  * The combined summary is stored under `<images_dir>/dataset~{dataset}/qc_stats.tsv`
* An anndata file in zarr format containing the QC metrics and status of each cell
  * The metrics are stored under `adata.obs[['n_counts', 'n_genes', 'percent_mito']]`
  * The QC status is stored under `adata.obs['qc_status']`
  * The zarr file is under `<out_dir>/qc/dataset~{dataset}/file_id~{file_id}.zarr`