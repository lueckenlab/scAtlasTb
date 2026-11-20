# Label Harmonisation

This module harmonises author cell type labels across datasets with [CellHint](https://github.com/Teichlab/cellhint), assigning each cell a harmonised label and grouping related labels across datasets, and produces tree plots, heatmaps, UMAPs and marker gene dot plots per harmonised group.

## Testing

Activate the snakemake environment and call `test/run_test.sh` with run specific Snakemake parameters.

```
conda activate snakemake
bash test/run_test.sh -n  # dry run
bash test/run_test.sh -c2  # actual run with max 2 cores
```

## Input

AnnData file (`.h5ad` or `.zarr`) containing:

+ `.obs[dataset_key]`: dataset/study of origin of each cell
+ `.obs[author_label_key]`: author cell type labels to be harmonised across datasets
+ `.obsm[use_rep]` (default `X_pca`): low-dimensional representation used by CellHint. If `X_pca` is missing, PCA is computed on the highly variable genes, which requires `.X` and a boolean `.var['highly_variable']`
+ `.X` normalised expression (needed for `cellhint: use_pct: true` and for marker gene dot plots)
+ `.obsm['X_umap']` for UMAP plots, unless `recompute_umap` or `recompute_neighbors` is set

## Config file

```yaml
DATASETS:
  dataset_name:  # can replace with a representative name
    input:
      label_harmonization: anndata_file_path.h5ad  # can be output from another module
    label_harmonization:
      dataset_key: study  # .obs column with dataset of origin
      author_label_key: cell_type  # .obs column with author labels to harmonise
      cellhint:  # parameters passed to cellhint.harmonize
        use_rep: X_pca
        use_pct: false
      subsample: 0.8  # optional fraction of cells to subsample before harmonisation
      force_scale: false  # scale .X before PCT even if data is marked as scaled
      recompute_neighbors: false  # recompute kNN graph (and UMAP) for plotting
      recompute_umap: false  # recompute UMAP for plotting
      plot_colors:  # additional .obs columns to colour UMAPs by
        - study
      marker_genes: blood  # comma-separated keys of MARKER_GENES for dot plots

MARKER_GENES:
  blood:
    - CD14
    - ITGB2
```

A complete example is in `workflow/label_harmonization/test/config.yaml`.
