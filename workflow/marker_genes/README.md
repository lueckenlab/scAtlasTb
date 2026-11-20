# Marker Genes

This module computes marker genes for user-defined groupings (e.g. cell type labels, batches, cell cycle phases) using the Wilcoxon rank-sum test or other methods available in [scanpy](https://scanpy.scverse.org/en/stable/api/generated/scanpy.tl.rank_genes_groups.html).
It also supports plotting user-provided gene sets alongside the statistically derived marker genes.

## Environments

- [`scanpy`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scanpy.yaml)

## Configuration

Example configurations for the marker genes module can be found in `workflow/marker_genes/test/config.yaml`.

### Example: multiple groupings with rank genes groups

```yaml
DATASETS:
  test:
    input:
      marker_genes:
        test_1: test/input/pbmc68k.h5ad
        test_2: test/input/pbmc68k.h5ad
    marker_genes:
      sample: louvain
      marker_genes: default
      rank_genes_groups:
        reference: rest
        n_genes: 100
        method: wilcoxon
      plot:
        n_genes: 20
        min_logfoldchange: 3
      groups:
        - bulk_labels
        - phase
        - batch
```

### Example: user-provided marker gene sets with a layer

```yaml
DATASETS:
  test2:
    input:
      marker_genes: test/input/pbmc68k.h5ad
    marker_genes:
      marker_genes: T_cells, Myeloid, B_cells, Other
      groups:
        - bulk_labels
      plot:
        n_genes: 10
        n_groups_per_split: 5
      layer: layers/counts
```

### Example: defining custom marker gene sets

```yaml
MARKER_GENES:
  default:
    "Pan-immune": [ "CD44" ]
    "T": [ "CD3D", "CD4", "CD8A" ]
    "NK": [ "KLRB1", "NCR1", "NCAM1", "GNLY" ]
  T_cells:
    "T": [ "CD3D", "CD4", "CD8A" ]
    "Naive T": [ "CCR7", "LEF1", "CD27", "SELL" ]
```

## Parameters

- **`groups`**: One or more `.obs` columns by which marker genes are computed. Each entry produces a separate analysis run.

- **`rank_genes_groups`**: Arguments passed directly to [`sc.tl.rank_genes_groups`](https://scanpy.scverse.org/en/stable/api/generated/scanpy.tl.rank_genes_groups.html).
  - `method`: Statistical test to use (e.g. `wilcoxon`, `logreg`, `t-test`). Default: `wilcoxon`.
  - `reference`: Reference group for comparison (e.g. `rest`). Default: `rest`.
  - `n_genes`: Number of top genes to compute per group. Default: all genes.

- **`marker_genes`**: User-provided marker gene sets to visualise.
  - String (comma-separated): names of gene sets defined under `MARKER_GENES` in the config (e.g. `default` or `T_cells, Myeloid`).
  - Dictionary: inline gene set mapping `{group_name: [gene1, gene2, ...]}`.

- **`plot`**: Parameters controlling visualisation of ranked gene results.
  - `n_genes`: Number of top genes per group to include in plots. Default: `10`.
  - `min_logfoldchange`: Minimum log-fold change to filter genes before plotting. Default: none.
  - `n_groups_per_split`: Number of cell groups to include per plot panel. Default: `100 / n_genes`.

- **`sample`**: `.obs` column used for pseudo-bulk aggregation before ranking. If set, counts are summed per `(group, sample)` pair before the statistical test. Default: none (single-cell level).

- **`layer`**: AnnData layer to use as the expression matrix. Default: `X`.

## Output

- **`<out_dir>/marker_genes/dataset~<dataset>/file_id~<file_id>.zarr`**: AnnData file containing the ranked gene results stored in `adata.uns` under the key `marker_genes_group=<group>` for each configured grouping.

- **`<out_dir>/marker_genes/groups/<params>/group=<group>--marker_genes.tsv`**: TSV file with differential expression statistics per cluster, including columns `gene`, `z-score`, `logfoldchange`, `pval`, `pval_adj`, `pct_within`, `pct_outside`, and `-log10 pvalue`.

- **`<image_dir>/marker_genes/<params>/group=<group>/rank_plot.png`**: Rank plot showing the top marker genes per group.

- **`<image_dir>/marker_genes/<params>/group=<group>/dotplot/`**: Dot plots of top marker gene expression and log-fold changes, split into panels for readability.

- **`<image_dir>/marker_genes/<params>/group=<group>/matrixplot/`**: Matrix plots of top marker gene expression and log-fold changes.

- **`<image_dir>/marker_genes/<params>/group=<group>/user_markers/`**: Dot plots of user-provided marker gene sets.
