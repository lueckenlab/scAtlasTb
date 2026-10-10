# Sample representation

This module computes sample-level representations from single-cell data, e.g. pseudobulk expression, cell type composition or learned sample embeddings.
Each method returns one representation (embedding and/or sample-by-sample distance matrix) per sample, which is then used to build a sample kNN graph, UMAP and plots.
Most methods are run through [patpy](https://github.com/lueckenlab/patpy).

## Environments

The following environments are needed, depending on which methods are configured:

- [`scanpy`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scanpy.yaml): data preparation and pseudobulk aggregation (`prepare` rule), UMAP and plots
- [`sample_representation`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/sample_representation.yaml): `pseudobulk`, `composition`, `cell_type_pseudobulk`, distance histograms
- [`pilot`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/pilot.yaml): `pilot`
- [`scarches`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scarches.yaml): `scpoli`
- [`scvi-tools`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scvi-tools.yaml): `mrvi`
- [`gloscope`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/gloscope.yaml): `gloscope` (R)
- [`scitd`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scitd.yaml): `scitd` (R)

## Configuration

```yaml
DATASETS:
  test:
    input:
      sample_representation:
        pbmc: test/input/pbmc68k.h5ad
    sample_representation:
      sample_key: bulk_labels  # .obs column(s) defining samples, comma-separated for composite keys
      cell_type_key: batch  # .obs column with cell type labels
      raw_counts: layers/counts  # counts slot used by feature-based methods
      norm_counts: layers/normcounts  # normalized slot used for pseudobulk aggregation in the prepare step
      aggregate: sum  # aggregation function for pseudobulks: sum or mean
      use_rep:  # embeddings used by embedding-based methods
        - obsm/X_pca
      var_mask: highly_variable  # .var column used to subset genes for count-based input
      min_cells_per_sample: 50
      min_cells_per_cell_type: 10
      colors:  # .obs columns to color the sample UMAP and PCA plots by
        - bulk_labels
      dask: false  # read data lazily with dask in the prepare step
      methods:  # each method with optional hyperparameters
        pseudobulk:
        composition:
        cell_type_pseudobulk:
        pilot:
        scpoli:
          n_epochs: 100
        mrvi:
          max_epochs: 100
        gloscope:
        scitd:
```

### Parameters

- **`methods`** (required): methods to run, as keys with optional hyperparameters. Available methods are listed in `params.tsv`: `pseudobulk`, `composition`, `cell_type_pseudobulk`, `pilot`, `scpoli`, `mrvi`, `gloscope`, `scitd`.
- **`sample_key`**: `.obs` column(s) that define a sample. Multiple columns can be given comma-separated and are concatenated. If not set, every cell is treated as its own sample.
- **`cell_type_key`**: `.obs` column with cell type labels, required by `composition`, `cell_type_pseudobulk`, `pilot`, `scpoli`, `mrvi` and `scitd`.
- **`use_rep`**: one or more embeddings (e.g. `obsm/X_pca`) for methods that accept an embedding as input (`pseudobulk`, `cell_type_pseudobulk`, `pilot`, `gloscope`). Embedding-based runs are skipped if `use_rep` is not set.
- **`raw_counts`** (default: `X`): slot with the counts used by feature-based methods.
- **`norm_counts`** (default: `X`): slot that is aggregated into pseudobulks in the preparation step.
- **`var_mask`**: boolean `.var` column used to subset genes. Only applies to count-based input (`X`, `layers/...`, `raw/...`).
- **`aggregate`** (default: `sum`): pseudobulk aggregation, `sum` or `mean`.
- **`min_cells_per_sample`**, **`min_cells_per_cell_type`**: minimum number of cells. Currently passed to the preparation step but not applied by the code; only samples with fewer than 2 cells are removed.
- **`colors`**: `.obs` columns used to color the sample-level plots.
- **`dask`**: whether to read the input lazily with dask in the preparation step.

## Output

Wildcards and file naming follow the general pattern `dataset~<dataset>/file_id~<file_id>/input_type~<input_type>/method~<method>--hyperparams~<hyperparams>--var_mask~<var_mask>--use_rep~<use_rep>`.

- `<output_dir>/sample_representation/prepare/dataset~<dataset>/file_id~<file_id>/prepare.zarr`: input with the sample column `obs['group']`
- `<output_dir>/sample_representation/prepare/dataset~<dataset>/file_id~<file_id>/bulks.zarr`: normalized and log-transformed pseudobulks per sample
- `<output_dir>/sample_representation/run_method/<pattern>.zarr`: sample-level AnnData (one observation per sample) with the representation in `.obsm` (`X_emb`, `X_pca` and/or `distances`) and a sample kNN graph in `.obsp`
- `<output_dir>/sample_representation/<pattern>.zarr`: same object with `obsm['X_umap']`
- `<images>/sample_representation/<pattern>/`: `distances_histplot.png`, `umap/` and `embeddings/` (PCA) plots

## Testing

Activate the snakemake environment and run the test workflow:

```bash
conda activate snakemake
bash test/run_test.sh -n
```
