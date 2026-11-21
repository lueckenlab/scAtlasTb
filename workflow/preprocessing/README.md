# Preprocessing

This module runs a configurable single-cell preprocessing workflow and assembles selected results into one output zarr file.

Implemented steps:

1. Normalize (`normalize.py`)
2. Mark non-zero genes for HVG-safe filtering (`filter_genes.py`)
3. Highly variable genes (`highly_variable_genes.py`)
4. Optional extra HVGs (`extra_hvgs.py`)
5. PCA (`pca.py`)
6. Neighbors graph (`neighbors.py`)
7. UMAP (`umap.py`)
8. PCA/UMAP plots (`plot.py`)
9. Assembly (`assemble.py`)

Rules are declared in `rules/rules.smk`, parameterized in `rules/assemble.smk`, and plotting rules are in `rules/plots.smk`.

## Testing

You can run the module on a small test dataset using:

```bash
snakemake --configfile test/config_no_gpu.yaml --use-conda --cores 8

# or, for the GPU-enabled variant
snakemake --configfile test/config.yaml --use-conda --cores 8
```

## Quickstart configuration

Global settings/defaults plus the `no_hvg` scenario from `test/config.yaml`:

```{eval-rst}
.. literalinclude:: ../../workflow/preprocessing/test/config.yaml
   :language: yaml
   :lines: 1-19,96-108
```

## Full configuration example

```yaml
DATASETS:
  dataset_name:
    input:
      preprocessing: adata.h5ad
    preprocessing:
      raw_counts: X
      dask: true
      n_threads: 10
      highly_variable_genes:
        n_top_genes: 2000
        batch_key: batch
      extra_hvgs:
        union_over: [lineage]
        extra_genes: [CCR7, PTPRC]
        min_cells: 200         # minimum cells per batch, counted within each lineage
        union_min_cells: 500   # minimum cells per lineage (default: min_cells)
      pca:
        n_comps: 50
      neighbors:
        n_neighbors: 15
      umap:
        min_dist: 0.5
      colors: [batch, cell_type, CCR7]
      plot_centroids: [cell_type]
      plot_gene_chunk_size: 12
      assemble:
        - normalize
        - highly_variable_genes
        - extra_hvgs
        - pca
        - neighbors
        - umap
```

## Cell and gene thresholds

`min_cells` means different things depending on the step:

| Key | Applies to | Default |
|---|---|---|
| `filter.min_cells` | Genes: minimum number of cells a gene must be expressed in to be considered for HVG selection (the feature space is not reduced) | `max(10, min(200, ceil(0.001 * n_obs)))`, unless `filter` sets `min_cells` or `min_counts` |
| `highly_variable_genes.min_cells` | Batches of `batch_key`: cells in smaller batches, or with a missing batch, are excluded from HVG selection | `max(n_bins, 200)` |
| `extra_hvgs.min_cells` | Batches of the effective `batch_key` (after `overwrite_args`). With `union_over`, batches are counted within each group, so a batch that is small in one group is excluded only from that group | `max(n_bins, 200)` |
| `extra_hvgs.union_min_cells` | `union_over` groups: groups with fewer cells (after batch filtering) are skipped; if none remain, HVGs are computed on all cells | `extra_hvgs.min_cells` |

Batch filtering exists because batch-aware HVG selection fails on very small batches: `cell_ranger` and `seurat_v3` raise errors (or, for `seurat_v3` with a single cell, crash the process) when a batch has only a handful of cells.

All step-specific keys (`normalize`, `filter`, `highly_variable_genes`, `extra_hvgs`, `pca`, `neighbors`, `umap`) as well as the global knobs `dask`, `n_threads` and `resources` can be set under `defaults: preprocessing:` and overridden per dataset.
The per-step behaviour and all parameters are described in the *Functional description* section of the module documentation page.
