# Marker Genes

```{include} ../../workflow/marker_genes/README.md
:heading-offset: 1
```

## Functional description

### Inputs

* **File formats:** `.h5ad` or `.zarr` (AnnData), configured under `input: marker_genes:` (see {ref}`architecture`).
* **Expression matrix:** the slot named by `layer` (default `X`; e.g. `layers/counts`). The matrix is used as stored — no normalisation is applied, so normalised data should be selected for cell-level tests, raw counts for pseudobulk tests.
* **`.obs`:** one grouping column per entry of `groups` (the parameter space is exploded so each column is its own `group` wildcard); optionally the `sample` column.
* **`.var`:** if a `feature_name` column exists it replaces `var_names` (gene symbols in outputs and for marker matching).
* **Marker gene sets:** resolved by `utils/marker_genes.py::get_marker_gene_set` from config `marker_genes` (string of comma-separated keys into top-level `MARKER_GENES`, or a dictionary).

### Processing steps

1. **`marker_genes_rank_genes_groups`** (script `scripts/rank_genes_groups.py`, environment `scanpy`). One job per file × `group` column.
   ```text
   read obs[[group, sample?]]
   use dask/backed reading if sample is set, or n_obs > 1e6, or n_unique(group[,sample]) > 1e5
   read X=<layer>, obs, var, uns
   if sample is set:
       pseudo_group = f"{group}_{sample}"
       adata = get_pseudobulks(adata, group_key=pseudo_group, agg="sum")   # pseudobulks with ≥2 cells
       rank_genes_groups.n_genes = n_vars                                  # report all genes
   groups = categories of <group> with > 1 observation (cells or pseudobulks)
   scanpy.tl.rank_genes_groups(adata, groupby=<group>, groups=groups, pts=True,
                               use_raw=False, key_added=f"marker_genes_group={group}",
                               **rank_genes_groups)
   ```
   Defaults of `scanpy.tl.rank_genes_groups` apply for unspecified arguments (e.g. `method='t-test'`, `reference='rest'`). For `.h5ad` input the original (non-pseudobulked) data is reloaded so that only `.uns` changes; the result is written zarr-linked with `X` mapped to the configured `layer`. A long-format table is derived from the result: per group, `gene`, `z-score`, `logfoldchange`, `pval`, `pval_adj`, `pct_within` (`pts`), `pct_outside` (`pts_rest`; for duplicated gene names the maximum is kept) and `-log10 pvalue`, indexed by `cluster`.

2. **`marker_genes_collect`** (script `scripts/collect.py`, environment `scanpy`). Merges `.uns` of the input file and of all per-group outputs (later keys overwrite earlier ones) into one object and writes `uns` only.

3. **`marker_genes_plot`** (script `scripts/plot.py`, environment `scanpy`). Reads the per-group output lazily, extracts results with `scanpy.get.rank_genes_groups_df(log2fc_min=plot.min_logfoldchange)`, takes the top `plot.n_genes` (default 10) genes per group, subsets cells of tested groups and these genes into memory and computes PCA (`scanpy.pp.pca`). Groups are partitioned into `max(1, n_groups // n_groups_per_split)` splits (default `n_groups_per_split = 100 // n_genes`). Produces `scanpy.pl.rank_genes_groups` (rank plot; an empty figure is saved on error) and, per split, `rank_genes_groups_dotplot` and `rank_genes_groups_matrixplot` coloured by expression and by log fold change (`cmap='bwr'`, `vmin=-4`, `vmax=4`), with group totals added.

4. **`marker_genes_plot_user`** (script `scripts/plot_user.py`, environment `scanpy`). Reads `layer`, subsets to the union of user marker genes present in `var_names` and loads them densely. List-valued gene sets are combined into one set `markers`; each dictionary-valued superset becomes its own gene set. For each gene set, genes not found are dropped (empty labels removed, empty sets skipped) and `scanpy.pl.dotplot(groupby=<group>, standard_scale='var', use_raw=False)` is drawn with totals. Exits without plots if no genes match. `plot.n_groups_per_split` is accepted but not used for splitting here.

### Outputs

Relative to `<output_dir>/marker_genes/`:

* `dataset~<dataset>/file_id~<file_id>.zarr` — `.uns['marker_genes_group=<group>']` for all groups (scanpy `rank_genes_groups` structure: `names`, `scores`, `logfoldchanges`, `pvals`, `pvals_adj`, `pts`, `pts_rest`, `params`).
* `groups/dataset~<dataset>/file_id~<file_id>/group=<group>.zarr` — per-group result (in `.uns`), linked to the input.
* `groups/dataset~<dataset>/file_id~<file_id>/group=<group>--marker_genes.tsv` — marker table (see step 1).

Relative to `<images>/marker_genes/dataset~<dataset>/file_id~<file_id>/group=<group>/`:

* `rank_plot.png`
* `dotplot/split=<i>_expression.png`, `dotplot/split=<i>_logfoldchanges.png`
* `matrixplot/split=<i>_expression.png`, `matrixplot/split=<i>_logfoldchanges.png`
* `user_markers/<gene_set>.png`

### Environments

* `scanpy` — all rules (no GPU path).
