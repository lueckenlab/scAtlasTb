# Batch Analysis

This module provides an exploratory framework for understanding the technical effects of batch variables at different hierarchical levels in single-cell datasets. It quantifies the strength of technical covariates by measuring the linear variance they explain using principal component regression (PCR). This systematic evaluation helps to assess the impact of each covariate, making it easier to decide which variable should be treated as the batch for subsequent batch correction or data integration steps.

## Features
- Principal component regression (PCR) analysis for quantifying linear effects of technical covariates on data
- Theil's U analysis for quantifying the pairwise association between covariates at sample level
- Pseudobulk generation and PCA plotting for sample-level exploration of batch structure
- Optional preprocessing pipeline: normalization, gene filtering, PCA
- Parallelized computation of batch PCR for scalability
- Flexible configuration for covariates, permutations, and sample keys

## Environments

The following conda environments are used by different steps:
- `scanpy` (for preprocessing, PCA, plotting)
- `scib` (for PC regression and batch PCR)

## Configuration

Configure the module under your dataset key using the `batch_analysis` section. Common keys:

- `sample`: column(s) in `.obs` to use as the sample key (can be comma-separated for composite keys). It should represent the smallest common grouping of a technical effect and the covariate of interest. An error will be thrown if this is not the case.
- `covariates`: list of covariate columns in `.obs` to test for batch effects. These covariates are also used to color the pseudobulk PCA plots and to generate the Theil's U heatmap.
- `permute_covariates`: (optional) list of covariates to permute for computing a z-score. If not specified, all `covariates` will be used for permutations. The covariates will be permuted per sample to compute empirical null distributions 
- `n_permutations`: number of permutations for each covariate (default: 10)
- `na_strings`: (optional) values treated as missing in covariates (default: `['NA', 'NaN', 'nan', '', 'unknown']`)
- `raw_counts`: (optional) slot with counts that are aggregated into pseudobulks (default: `X`)
- `max_threads`: (optional) maximum number of threads for pseudobulk aggregation and permutations (default: 1)
- `pca_plot`: (optional) settings for the pseudobulk PCA plots: `plot_centroids`, `plot_gene_chunk_size`, `size` (default: 50)
- Step-specific overrides (e.g., `normalize`, `highly_variable_genes`, `pca`) for preprocessing

> **Note:** Preprocessing steps such as normalization, and PCA are optional. A step is only executed if its corresponding key (`normalize`, `highly_variable_genes`, `pca`) is defined in your configuration and the result is used by a later step: PCA is only computed if `pca` is defined, gene filtering and highly variable gene selection only if both `highly_variable_genes` and `pca` are defined, and normalization only if all three keys are defined.

### Example configuration

```yaml
DATASETS:
  BATCH_ANALYSIS_PREPROCESSED:
    input:
      batch_analysis:
        file_1: test/input/blood_pca.zarr
    batch_analysis:
      covariates:
        - sample
        - donor
        - assay
        - sex
        - disease
        - self_reported_ethnicity
      permute_covariates:
        - assay
        - sex
        - disease
        - self_reported_ethnicity
      n_permutations: 1000
      sample: sample,donor

  BATCH_ANALYSIS_UNPROCESSED:
    input:
      batch_analysis:
        file_1: test/input/pbmc68k.h5ad
    batch_analysis:
      sample: batch, bulk_labels
      covariates:
        - bulk_labels
        - batch
        - is_cd14_mono
      normalize:
      highly_variable_genes:
      pca:
```

The example configuration above demonstrates how to set up the `batch_analysis` module for two datasets:

- **`BATCH_ANALYSIS_PREPROCESSED`**: Uses a dataset with precomputed PCA in the input AnnData object. It specifies multiple covariates (e.g., `sample`, `donor`, `assay`, `sex`, `disease`, `self_reported_ethnicity`) to test for batch effects, and defines which covariates to permute for significance testing. The number of permutations is set to 1000, and a composite sample key (`sample,donor`) is used.

- **`BATCH_ANALYSIS_UNPROCESSED`**: Uses a dataset without the necessary PCA information and configures the workflow to perform normalization, highly variable gene selection, and PCA as preprocessing steps. It sets `batch` and `bulk_labels` as the composite sample key, and tests covariates such as `bulk_labels`, `batch`, and `is_cd14_mono` for batch effects.


## Workflow steps

The batch_analysis workflow consists of the following steps:

1. **Preprocessing** (optional)
  - Steps: normalize, filter genes, HVG selection, PCA
  - Each preprocessing step is optional and will only be executed if its corresponding key is defined in the configuration and its output is consumed downstream. For example, if only `pca` is defined, only PCA will be computed on the input; if `normalize`, `highly_variable_genes` and `pca` are defined, all steps are performed. This allows users to skip preprocessing if their input data is already preprocessed and contains the necessary PCA information for batch PCR analysis.
  - Uses rules from the preprocessing module, with dataset-specific overrides.
2. **Prepare data**:
  - Sets sample key for pseudobulk aggregation and PCR analysis based on the configured `sample` key. The sample key represents the smallest common grouping of a technical effect and the covariate of interest (e.g., `batch` or `donor`).
  - Aggregates cells into pseudobulk samples using the configured `sample` key.
  - Recomputes PCA on the pseudobulk data and colors the PCA plots by the configured `covariates`.
3. **Covariate setup** (`determine_covariates`):
  - Determines which covariates to test and sets up permutation schemes.
4. **Batch PC regression** (`batch_pcr`):
   - Runs principal component regression for each covariate and permutation, computes z-scores.
5. **Collect**:
   - Aggregates per-covariate results into a single table.
6. **Plot**:
  - Generates barplots and violin plots summarizing PCR and permutation results.
  - Generates a Theil's U heatmap to visualize the pairwise association between covariates.

## Output

### Pseudobulk
- Pseudobulk AnnData/Zarr output: `<output_dir>/batch_analysis/prepare/dataset~<dataset>/file_id~<file_id>/pseudobulks.zarr`
- Pseudobulk PCA plots: `<images>/batch_analysis/dataset~<dataset>/file_id~<file_id>/pca_plots`

### Principal regression analysis
- Per-covariate PCR results: `<output_dir>/batch_analysis/dataset~<dataset>/file_id~<file_id>/batch_pcr/{covariate}.tsv`
- Aggregated results: `<output_dir>/batch_analysis/dataset~<dataset>/file_id~<file_id>/batch_pcr.tsv`
- Plots: `<images>/batch_analysis/dataset~<dataset>/file_id~<file_id>/batch_pcr_bar.png`, `<images>/batch_analysis/dataset~<dataset>/file_id~<file_id>/batch_pcr_violin.png`

Each result file contains columns such as:
- `covariate`: tested covariate
- `pcr`: PCR score
- `permuted`: whether the score is from a permutation
- `n_covariates`: number of unique values in the covariate
- `z_score`: z-score of observed PCR vs. permutations
- `p-val`, `signif`: empirical p-value and significance stars


### Theil's U
- Results: `<output_dir>/batch_analysis/dataset~<dataset>/file_id~<file_id>/theils_u.tsv`
- Plots: `<images>/batch_analysis/dataset~<dataset>/file_id~<file_id>/theils_u_heatmap.png`

## Testing

Activate the snakemake environment and run the test workflow:

```bash
conda activate snakemake
bash test/run_test.sh -n
```

This test requires precomputed objects from the `preprocessing` tests.
See the `preprocessing` module for more details.