# Preprocessing

```{include} ../../workflow/preprocessing/README.md
:heading-offset: 1
```

## Functional description

### Inputs

* Input files are configured under `DATASETS.<dataset>.input.preprocessing` (mapping `file_id` → path). Both `.h5ad` and `.zarr` are accepted; every rule writes `.zarr`.
* Count matrix: read from the slot given by `raw_counts` (default `X`, e.g. `layers/counts`). It must contain un-normalised counts.
* `gene_id_column` (optional): `.var` column that replaces `var_names` before processing (names are cast to `str` and made unique).
* `.obs` columns referenced by the configuration: `highly_variable_genes.batch_key`, `extra_hvgs.overwrite_args.batch_key`, `extra_hvgs.union_over`, and plot `colors`/`plot_centroids`. `.var['feature_name']` is used (if present) for gene matching in `extra_hvgs` and as gene labels in plots.
* `.uns` is carried through; an existing neighbour graph (`.obsp['connectivities']`, `.obsp['distances']`, `.uns['neighbors']`) is only read when `neighbors: false`.

Global knobs (overridable per step inside the step dictionary where noted):

| key | default | effect |
| --- | --- | --- |
| `dask` | `true` for normalize/filter/PCA, `false` for HVG/extra HVGs | read matrices lazily as dask arrays; a `dask` entry inside `normalize`, `highly_variable_genes`, `extra_hvgs` or `pca` overrides the global value for that step |
| `n_threads` | `10` | threads per job (also dask workers) |
| `resources` | `cpu` (normalize, HVG, extra HVGs), `gpu` (PCA, neighbors, UMAP) | resource profile; also selects the GPU environment, see Environments |
| `scale` | `false` | scale HVG matrix before PCA |
| `assemble` | all steps | steps linked into the final file |

Parameter sweeps: `highly_variable_genes` and `extra_hvgs.overwrite_args` may contain lists; `PreprocessingConfig` expands them into all combinations, each becoming a value of the `hvg_args` / `overwrite_args` wildcard (a 10-character blake2b hash of the combination; only the HVG keys `n_top_genes, min_disp, max_disp, min_mean, max_mean, span, n_bins, flavor, batch_key` enter the hash, other keys such as `dask` or `min_cells` do not). `highly_variable_genes: false` is kept as a distinct value meaning "all genes". Wildcards and file naming follow {ref}`architecture`.

### Processing steps

Rules in DAG order (rule names as in the full pipeline):

1. **`preprocessing_normalize`** (script `scripts/normalize.py`, environment `scanpy` / `rapids_singlecell`)
2. **`preprocessing_filter_genes`** (`scripts/filter_genes.py`, `scanpy`)
3. **`preprocessing_highly_variable_genes`** (`scripts/highly_variable_genes.py`, `scanpy` / `rapids_singlecell`), one job per `hvg_args`
4. **`preprocessing_extra_hvgs`** (`scripts/extra_hvgs.py`, `scanpy` / `rapids_singlecell`), one job per `overwrite_args`
5. **`preprocessing_pca`** (`scripts/pca.py`, `scanpy` / `rapids_singlecell`)
6. **`preprocessing_neighbors`** (`scripts/neighbors.py`, `scanpy` / `rapids_singlecell`)
7. **`preprocessing_umap`** (`scripts/umap.py`, `scanpy` / `rapids_singlecell`)
8. **`preprocessing_plot_pca`**, **`preprocessing_plot_umap`** (`scripts/plot.py`, `scanpy`)
9. **`preprocessing_assemble`** (`scripts/assemble.py`, `scanpy`)

Each step writes only the slots it creates and links all other slots to its input zarr (`write_zarr_linked`), so intermediate files are small. Every script handles empty inputs (`n_obs == 0`) by writing a valid empty zarr with the expected slots.

When the GPU environment is active and `nvidia-smi` succeeds, `utils.processing` imports `rapids_singlecell` as `sc` (with an RMM managed-memory allocator) and data are moved to the GPU (`sc.get.anndata_to_GPU`) before the heavy calls; otherwise `scanpy` is used with identical arguments.

#### Normalize

```text
X ← read(raw_counts, backed/dask if dask)
if gene_id_column: var_names ← var[gene_id_column]; make unique
ensure_sparse(X)
if input is .h5ad: layers['counts'] ← X; raw ← adata
scanpy.pp.normalize_total(adata, **normalize)    # e.g. target_sum
scanpy.pp.log1p(adata)
layers['normcounts'] ← X
uns['preprocessing'] |= {normalization: normalize, log-transformed: True}; uns['log1p'] = {base: None}
```

For zarr inputs the counts are not copied: `layers/counts`, `raw/X` are linked to the `raw_counts` slot of the input, `raw/var` to `var`, and `X` to the new `layers/normcounts`.

Example:

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      raw_counts: X
      gene_id_column: gene_id
      normalize:
        target_sum: 1e4
```

#### Filter genes

Reads `X` (normalised counts) and computes `var['nonzero_genes'] = scanpy.pp.filter_genes(X, **filter)[0]` (default `filter: {min_cells: 1}`). This is only a mask used by the HVG steps to avoid genes that break HVG selection; the feature space is not reduced.

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      filter:
        min_cells: 3
```

#### Highly variable genes

```text
args ← hvg_args_dict (one expanded combination); drop 'subset'; pop 'dask'
column ← 'highly_variable' + ''.join(f'-{k}={v}' for k in sorted(args))
X ← layers/counts if flavor == 'seurat_v3' else layers/normcounts
if args is False:  var['highly_variable'] ← True        # no HVG selection
else:
    min_cells ← args.pop('min_cells', max(args.n_bins or 20, 200))
    drop batches of batch_key with < min_cells cells (utils.accessors._filter_batch)
    subset genes to var['nonzero_genes']
    sort cells by batch_key (contiguous batches); persist dask array (chunks 200 000 × all genes)
    scanpy.pp.highly_variable_genes(adata, **args)
    map highly_variable, means, dispersions, dispersions_norm,
        highly_variable_nbatches, highly_variable_intersection back to full var
        (genes that were filtered out get False / 0)
var[column] ← var['highly_variable']
```

Previous `highly_variable*` columns of the input are removed. Arguments are stored in `uns['preprocessing']['highly_variable_genes']`. If no HVG config is given (`{}`), `scanpy` defaults are used and the column is just `highly_variable`. List-valued parameters produce one output per combination; only the first combination feeds PCA/neighbors/UMAP and becomes the canonical `highly_variable` column in the assembled file.

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      highly_variable_genes:
        n_top_genes: 2000
        flavor: seurat_v3
        batch_key: sample
```

#### Extra HVGs

Additional gene mask, independent of the PCA input.

```text
args ← highly_variable_genes config (drop subset) | overwrite_args
column ← 'extra_hvgs' + ''.join(f'-{k}={v}' for k in sorted(overwrite_args))
min_cells ← extra_hvgs.min_cells or max(args.n_bins or 20, 200)
read layers/counts (seurat_v3) or layers/normcounts; filter batches (< min_cells) and nonzero_genes
sort cells by batch_key
remove genes matching remove_genes (subset before HVG selection)
if union_over:
    group ← '--'.join(obs[union_over]) with 'nan'/'unknown' as missing; drop groups < min_cells
    for each group: scanpy.pp.highly_variable_genes(group, **args); column |= group HVGs
    (falls back to all cells if no valid group remains)
else: scanpy.pp.highly_variable_genes(adata, **args)
map back to full var (False for filtered genes); set column True for matched extra_genes
```

`extra_genes` and `remove_genes` entries are matched with `utils.accessors.match_genes` against `var['feature_name']` (if present) or `var_names`; entries may be gene names, regular expressions, paths to existing local text files or HTTP(S) URLs with one gene per line (entries that are neither an existing file nor a URL are treated as gene patterns). Metadata: `uns['preprocessing'][<column>] = args | extra_hvgs`.

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      extra_hvgs:
        union_over: [lineage]
        extra_genes: [CCR7, PTPRC]
        remove_genes: [MALAT1]
        min_cells: 200
        overwrite_args:
          n_top_genes: 3000
          flavor: seurat_v3
```

#### PCA

Input is the HVG output of the first `hvg_args` combination.

```text
X ← normalised counts; mask ← pca.mask_var (default 'highly_variable')
adata_pca ← subset_hvg(adata, mask)            # no subset if all genes are True
if scale: scanpy.pp.scale(adata_pca)
scanpy.pp.pca(adata_pca, svd_solver='covariance_eigh' (default), **pca)   # 'zero_center: None' string → None
obsm['X_pca'], uns['pca'] ← adata_pca; varm['PCs'] mapped to all genes (NaN for non-HVGs)
uns['preprocessing'] |= {pca: pca args, scaled: scale}
```

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      scale: true
      pca:
        n_comps: 50
        svd_solver: covariance_eigh
```

#### Neighbors

```text
if neighbors is False:                      # reuse existing graph from input
    read obsp; set uns['neighbors'] keys and use_rep ('X_pca' if present else 'X'); assert_neighbors
else:
    use_rep ← neighbors.use_rep or ('X_pca' with n_pcs = #PCs if present else 'X')
    if use_rep == 'X': densify X;  if use_rep == 'X_pca' missing: scanpy.pp.pca(**uns.preprocessing.pca)
    n_neighbors ← min(neighbors.n_neighbors or 15, n_obs)
    sc.pp.neighbors(adata, **neighbors)     # on GPU failure: retry with scanpy.pp.neighbors
```

Writes `obsp['connectivities']`, `obsp['distances']`, `uns['neighbors']` (plus `obsm`/`varm` if PCA had to be computed). The same script is reused by the integration module (`integration_postprocess`).

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      neighbors:
        n_neighbors: 15
        metric: cosine
        use_rep: X_pca
```

#### UMAP

`scanpy.tl.umap(adata, init_pos='random', **umap)` on the graph referenced by `umap.neighbors_key` (default `neighbors`). Missing `connectivities_key`/`distances_key`/`params` entries are filled with defaults, `n_neighbors` is inferred from the distance matrix if absent. If the representation named in `uns[neighbors_key]['params']['use_rep']` is not in the graph file, it is read from the `rep` input (PCA file), or PCA (50 components) is computed on `X`. On GPU (cuML), connectivities are replaced by distances before the call and the embedding is clipped to the [1e-5, 1−1e-5] quantiles per axis (`utils.misc.trim_umap`). Writes `obsm['X_umap']` and `uns`.

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      umap:
        min_dist: 0.3
        spread: 1.0
```

#### Plots

`preprocessing_plot_pca` (basis `X_pca`, input PCA file) and `preprocessing_plot_umap` (basis `X_umap`, input UMAP file) call `scripts/plot.py`:

```text
colors ← colors + plot_centroids (deduplicated)
obs colors: keep columns with > 1 unique value; categorical/string columns: 'NaN','None','','nan','unknown' → missing,
            drop categories with ≤ 1e-4 · n_cells cells
gene colors: remaining entries matched to var_names (exact, else substring/regex), sorted, chunked by plot_gene_chunk_size (default 12)
remove embedding outliers (UMAP only, outlier_factor=100: cells with |max| or |min| ≥ 100 × mean)
if n_cells > 1e6: plot a random 70 % of cells
point size ← clip(200 000 / n_obs, min 0.5 or 1, max 200)
for each color (parallel threads): scanpy.pl.embedding(basis, color, palette) → <color>.png
for each gene chunk: scanpy.pl.embedding(color=genes, ncols=4) → genes_group=<i>.png
```

Palettes: > 102 categories `turbo` (legend removed), > 20 categories `godsnot_102`, numeric `coolwarm` if negative values else `plasma`. Categorical legends show group sizes; for columns in `plot_centroids` (≤ 102 categories, UMAP only) category numbers are drawn at the per-category median coordinates and prefixed in the legend. Figure titles list the wildcards and the number of cells.

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      colors: [batch, cell_type, CCR7, PTPRC]
      plot_centroids: [cell_type]
      plot_gene_chunk_size: 12
```

#### Assemble

`collect_files` gathers the outputs of `normalize`, `highly_variable_genes` (all `hvg_args`), `extra_hvgs` (all `overwrite_args`), `pca`, `neighbors` and `umap`, restricted to the entries listed in `assemble` (unknown entries such as `counts` are ignored; counts are always available through the linked input). The first collected file provides `obs`/`var`/`uns`; slots of the other files are linked (zarr) or copied (h5ad) into the output:

| `assemble` entry | slots linked |
| --- | --- |
| `normalize` | `X`, `layers` (`counts`, `normcounts`), `raw`, `uns/log1p`, `uns/preprocessing/{normalization,log-transformed}` |
| `highly_variable_genes` | `var/highly_variable-<args>` for every combination; the first (default) one also as `var/highly_variable`; `uns/preprocessing/highly_variable_genes` metadata |
| `extra_hvgs` | `var/extra_hvgs-<args>` per combination; default one also as `var/extra_hvgs`; matching `uns/preprocessing` entries |
| `pca` | `obsm/X_pca`, `uns/pca`, `varm/PCs` |
| `neighbors` | `obsp/connectivities`, `obsp/distances`, `uns/neighbors` |
| `umap` | `obsm/X_umap` |

Wildcards are stored in `uns['wildcards']` (`preprocessing_<wildcard>`). An empty first input produces an empty output.

```yaml
DATASETS:
  dataset_name:
    preprocessing:
      assemble:
        - normalize
        - highly_variable_genes
        - extra_hvgs
        - pca
        - neighbors
        - umap
```

### Outputs

Under `<output_dir>/preprocessing/` (`<pattern>` = `dataset~<dataset>/file_id~<file_id>`, see {ref}`architecture`):

* `preprocessed/<pattern>/normalized.zarr`, `filtered_genes.zarr`, `highly_variable_genes--<hvg_args>.zarr`, `extra_hvgs--<overwrite_args>.zarr`, `pca.zarr`, `neighbors.zarr`, `umap.zarr`: per-step intermediates.
* `<pattern>.zarr`: assembled result (default target), containing depending on `assemble`:
  `X` (log-normalised), `layers['counts']`, `layers['normcounts']`, `raw`, `var['nonzero_genes']` (via linking), `var['highly_variable']`, `var['highly_variable-<args>']`, `var['means'|'dispersions'|'dispersions_norm'|'highly_variable_nbatches'|...]` (from the HVG file used as base), `var['extra_hvgs']`, `var['extra_hvgs-<args>']`, `obsm['X_pca']`, `varm['PCs']`, `obsm['X_umap']`, `obsp['connectivities']`, `obsp['distances']`, `uns['pca']`, `uns['neighbors']`, `uns['log1p']`, `uns['preprocessing']` (arguments of each step and `scaled`), `uns['wildcards']`.

Under `<images>/preprocessing/<pattern>/`:

* `pca/<color>.png`, `umap/<color>.png` for every valid `.obs` color, and `genes_group=<i>.png` for gene panels.

### Environments

* `scanpy`: filter genes, plots, assemble; and all other steps when `use_gpu: false` or the step's `resources` profile is `cpu`.
* `rapids_singlecell` (GPU only): normalize, highly variable genes, extra HVGs, PCA, neighbors and UMAP when `use_gpu: true` and the step's profile is not `cpu`. By default only PCA, neighbors and UMAP use the `gpu` profile; set `resources: gpu` to also run normalisation and HVG selection on the GPU. Inside this environment the scripts still fall back to `scanpy` if no GPU is detected at runtime.
