# Exploration

```{include} ../../workflow/exploration/README.md
:heading-offset: 1
```

## Functional description

### Inputs

- One or more AnnData files per dataset (`input: exploration: {file_id: path}`). These are typically harmonized outputs of `load_data`. The scripts rely on the harmonized schema:
  - `.obs`: the columns named by `sample` and `donor`, all `summary_columns`, plus the fixed columns `sample` (in `summary_stats` and `barcode_matching`), `author_annotation` and `cell_type` (in `marker_genes`),
  - `.var['feature_name']`: gene symbols, used for marker matching and mitochondrial genes,
  - `.X`: raw counts, used for QC metrics and marker dot plots.
- `marker_genes`: comma-separated keys into the top-level `MARKER_GENES` mapping (resolved by `utils.marker_genes.get_marker_gene_set`).
- `marker_genes` and `barcode_matching` read the file with `anndata.read_zarr` / `zarr.open`, so they require `.zarr` input. `summary_stats` uses `read_anndata` (`.zarr` or `.h5ad`).

### Processing steps

All rules run per `dataset~<dataset>/file_id~<file_id>` (see {ref}`architecture`), except `summary_stats_all`, which runs per dataset.

1. **`exploration_summary_stats`** (script `scripts/summary_stats.py`, environment `scanpy`):

   ```text
   adata := read(X, obs, var)
   bar plot (barh) of obs[sample].value_counts()  -> per_sample.png
   bar plot (barh) of obs[donor].value_counts()   -> per_donor.png
   var['mito'] := feature_name starts with "MT-"
   sc.pp.calculate_qc_metrics(adata, qc_vars=['mito'])
   stats := study (= file_id), n_cells, n_genes, n_samples, n_donors,
            median/mean cells per sample, median cells per donor,
            median total_counts per cell, median n_genes_by_counts per cell
   for c in summary_columns: assert c in obs; stats[c] := ','.join(unique values)
   ```

   Empty files produce empty plots.
2. **`exploration_summary_stats_all`** (script `scripts/summary_plot.py`, environment `scanpy`, local rule): concatenates the per-file `stats.tsv` of a dataset, indexed by `study`. It draws one horizontal bar subplot per statistic (studies sorted by `n_cells`, two columns) and writes the sums of `n_cells`, `n_samples` and `n_donors` over all files.
3. **`exploration_marker_genes`** (script `scripts/marker_genes.py`, environment `scanpy`): sets `var_names := var['feature_name']`, keeps only marker genes present in the data, and draws two `sc.pl.dotplot(use_raw=False, standard_scale='var')` panels, grouped by `obs['author_annotation']` and by `obs['cell_type']`. It fails if `author_annotation` has no values. Empty files produce an empty figure. Note: the rule calls `get_marker_gene_set` without `flatten=True`. When `marker_genes` names a key of `MARKER_GENES`, the script therefore receives the nested mapping `{key: {group: [genes]}}`, and `isin` matches the group names instead of the gene lists. For marker genes to be matched, the per-group gene lists must be at the top level.
4. **`exploration_barcode_matching`** (script `scripts/barcode_matching.py`, environment `plots`): checks whether cell barcodes are shared between samples (e.g. to detect duplicated or demultiplexing artefacts).

   ```text
   obs := read obs from zarr
   samples := up to 100 random samples of obs['sample'] (np.random.seed(42))
   barcodes[s] := obs_names of cells in sample s
   membership := for each barcode, the set of samples containing it
   counts := value counts of membership patterns
             keep patterns with <= 5 samples and size >= min(50, max count)
   upsetplot.UpSet(counts).plot()
   ```

   The commented-out rules `get_barcodes` and `barcode_matching` (R) in `rules/barcode_matching.smk` are a legacy alternative and are not used by default. In that version, `scripts/get_barcodes.py` writes `(barcode, sample)` pairs with the barcode suffix after the first `-` removed, and `scripts/barcode_matching.R` casts them to a barcode × sample membership table and draws a `ComplexUpset::upset` plot (`min_size` parameter, saved with `ggsave`, 300 dpi).

### Outputs

All outputs are images and tables under `<images>/exploration/`:

| Path | Content |
|---|---|
| `summary/dataset~…/file_id~…/stats.tsv` | per-file statistics (one row) |
| `summary/dataset~…/file_id~…/per_sample.png`, `per_donor.png` | cells per sample / donor |
| `summary/dataset~…/per_dataset_stats.tsv` | stats of all files of a dataset |
| `summary/dataset~…/per_dataset_stats_aggregated.tsv` | summed `n_cells`, `n_samples`, `n_donors` |
| `summary/dataset~…/per_dataset_stats.png` | per-study bar charts |
| `marker_genes/dataset~…/file_id~….png` | marker gene dot plots |
| `barcode_matching/dataset~…/file_id~….png` | UpSet plot of shared barcodes |

No AnnData objects are written.

### Environments

- `scanpy`: summary statistics, summary plot, marker genes.
- `plots`: barcode matching (`upsetplot`; `r-complexupset` for the legacy R script).
