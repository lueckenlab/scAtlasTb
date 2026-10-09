# Marker Genes

This module computes differentially expressed marker genes for one or more cell groupings (e.g. cell type or cluster columns in `.obs`) with `scanpy.tl.rank_genes_groups`, optionally on per-sample pseudobulks, and visualises both the computed markers and user-defined marker gene sets.

## Environments

- [`scanpy`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scanpy.yaml)

## Configuration

```yaml
DATASETS:
  test:
    input:
      marker_genes:
        test_1: test/input/pbmc68k.h5ad
    marker_genes:
      groups:              # .obs columns to compute markers for (one job per column)
        - bulk_labels
        - phase
      sample: louvain      # optional: aggregate pseudobulks per group x sample
      layer: X             # matrix to use, e.g. X or layers/counts
      marker_genes: default  # optional: comma-separated keys of MARKER_GENES
      rank_genes_groups:   # passed to scanpy.tl.rank_genes_groups
        reference: rest
        n_genes: 100
        method: wilcoxon
      plot:                # passed to the scanpy plotting functions
        n_genes: 20
        min_logfoldchange: 3
        n_groups_per_split: 5

MARKER_GENES:
  default:
    "T": [ "CD3D", "CD4", "CD8A" ]
    "B": [ "CD19", "CD79A" ]
```

A full example is in `workflow/marker_genes/test/config.yaml`.

### Parameters

- **`groups`**: `.obs` column(s) to group cells by. Each column is processed separately.
- **`sample`** (optional): `.obs` column with sample/donor IDs. If set, counts are summed into pseudobulks per group × sample before testing.
- **`layer`** (default: `X`): expression matrix to use for testing and user marker plots.
- **`rank_genes_groups`**: keyword arguments for [`scanpy.tl.rank_genes_groups`](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.tl.rank_genes_groups.html) (e.g. `method`, `reference`, `n_genes`).
- **`plot`**: keyword arguments for the scanpy `rank_genes_groups` plotting functions. Special keys: `n_genes` (genes per group, default 10), `min_logfoldchange` (filter for plotted genes), `n_groups_per_split` (number of groups per dot/matrix plot).
- **`marker_genes`**: either a comma-separated string of keys into the top-level `MARKER_GENES` dictionary, or a dictionary of gene sets. Gene sets may be a list of genes, a dictionary `{label: [genes]}`, or a nested dictionary of such dictionaries.

## Output

- `<output_dir>/marker_genes/dataset~<dataset>/file_id~<file_id>.zarr`: AnnData with `.uns['marker_genes_group=<group>']` (the `rank_genes_groups` result) for every configured group.
- `<output_dir>/marker_genes/groups/dataset~<dataset>/file_id~<file_id>/group=<group>--marker_genes.tsv`: table of marker genes per group with scores, log fold changes, p-values and expression fractions.
- `<images>/marker_genes/dataset~<dataset>/file_id~<file_id>/group=<group>/`: rank plot, dot plots and matrix plots of computed markers, and dot plots of user-defined marker genes (`user_markers/`).
