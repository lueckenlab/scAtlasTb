# Common: commonly used workflows

Module for workflows that are commonly used in other modules.

* dependency graph (rule and job graphs)
* adding annotations to obs (`rules/annotate.smk`, not included by default)
* plotting results from TSV file and from AnnData objects

`common` is not configured per task. Every module imports `workflow/common/Snakefile` (`use rule * from common as common_*`), which includes `rules/dependency_graph.smk` and `rules/plots.smk`; the plotting rules are additionally imported as Snakemake modules by `integration` and `metrics`, which override their inputs, outputs and params (`use rule barplot from plots as ...`). The inputs and params listed below are the defaults in the rule definitions, which mainly serve as templates and for the module's own tests (`test_barplot`, `test_umap`, ...).

## Plots

Rules in `rules/plots.smk`

Facetted plots for comparing computational and quality metrics.
Currently supported plots are:

* `barplot`
* `swarmplot`

### Input

TSV file containing the columns needed by facetting and metrics columns
Which columns are required for facetting the plot can be specified in the `params` directive of the rule.
The rule can be imported as a [Snakemake module](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules)
so that the inputs, outputs and parameters can be modified.

#### Example

The rules expect the benchmark output provided by the Snakemake directive.
The snakemake benchmark columns are as follows:

```tsv
method    s    h:m:s    max_rss    max_vms max_uss    max_pss io_in    io_out    mean_load    cpu_time
```

The entries are explained nicely
here: https://stackoverflow.com/questions/46813371/meaning-of-the-benchmark-variables-in-snakemake

An example configuration of the plot parameters:

```snakemake
params:
    metric='s',
    category='method',
    hue=None,
    facet_row=None,
    facet_col=None,
    title='Barplot',
    description=lambda wildcards: ' '.join([f'{key}={value}' for key, value in wildcards.items()]),
    dodge=True,
```

## Functional description

### Inputs

* **Dependency graphs:** no input file required, as the graphs are built on the given target: the resolved Snakemake `config` and a target rule name (wildcards `images`, `target`).
* **AnnData plots:** `{filename}.h5ad` (read with `scanpy.read` or `utils/io.py::read_anndata`); `.obsm[basis|use_rep]`, `.obs[color|groupby]`, `.X`.
* **Table plots:** a TSV with one row per observation (e.g. Snakemake benchmark files or metric tables) containing the `metric`, `category` and facet/hue columns (see [Plots](#plots)).

### Processing steps

1. **`common_save_config`** (local, inline Python): dumps the in-memory `config` dict to `.snakemake/{images}_{target}/config.json`.
2. **`common_rulegraph`** / **`common_dag`** (local, shell): run `snakemake {target} --configfile config.json --rulegraph` (resp. `--dag`), extract the `digraph` block with `sed` and render it with Graphviz `dot -Tpng` (`rankdir=TB` for the rule graph, `LR` for the job DAG). Requires `snakemake` and `dot` on the `PATH`.
3. **`common_dependency_graph`**: aggregate target requiring both PNGs; each module defines its own `dependency_graph` rule pointing to `<images>/<module>` with `target='all'`.
4. **`common_dotplot`** (script `scripts/dotplot.py`, env `scanpy`): `scanpy.pl.dotplot(adata, **params)` — all rule params are forwarded (template: `var_names`, `groupby`, `use_raw=False`, `standard_scale='var'`, `dendrogram=False`, `swap_axes=False`); saved at 200 dpi.
5. **`common_embedding`** (script `scripts/embedding.py`, env `scanpy`): requires `basis`; densifies `.obsm[basis]` if sparse; drops colour columns with more than 128 categories; `scanpy.pl.embedding(adata, **params)` with `frameon=False, vector_friendly=True, fontsize=9`.
6. **`common_umap`** (script `scripts/umap.py`, env `scanpy`, `scanpy_rapids` if `use_gpu: true`):
   ```text
   use_rep = params.use_rep if in .obsm else 'X' (warning)
   densify use_rep
   if uns['neighbors'] missing or computed on another representation:
       sc.pp.neighbors(use_rep, method='rapids'), fallback to default method
   sc.tl.umap(method='rapids'), fallback to default
   save obsm['X_umap'] to .npy
   remove outlier cells whose max (then min) UMAP coordinate exceeds 10x the mean (plot only)
   drop colours with > 128 categories; sc.pl.umap(**params)
   ```
7. **`common_barplot`** (script `scripts/barplot.py`, env `plots`, local): `seaborn.catplot(kind='bar', x=metric, y=category, row=facet_row, col=facet_col, hue=hue)`; categories sorted by their minimum metric value (descending); `hue` is dropped if it has more than 6 levels; title = `title` + `description` (default: wildcards as `key=value`); optional `xlim`/`ylim`.
8. **`common_swarmplot`** (script `scripts/swarmplot.py`, env `plots`, local): `seaborn.catplot(kind='swarm', x=category, y=metric, hue=hue, row/col facets, s=10)`, categories ordered as for the bar plot, x labels rotated by 90°, optional `ylim`.

Two further rule files are present but **not** included by `common/Snakefile` (they are only usable when included explicitly):

* `rules/annotate.smk` — **`add_obs`** (script `scripts/add_obs.py`): reads an h5ad and a TSV, indexes the TSV by `tsv_column` and replaces `.obs` by `pandas.merge(adata.obs, annotation)` (merge on shared column names; `h5ad_column` is not used), writes h5ad.
* `rules/convert.smk` — **`zarr_to_h5ad`** (script `scripts/convert_zarr_h5ad.py`): `anndata.read_zarr(file.zarr).write_h5ad(file.h5ad)`.

### Outputs

* `<images>/<module>/rule_graphs/all.png`, `<images>/<module>/job_graphs/all.png`, `.snakemake/<images>_<target>/config.json`.
* `{filename}_dotplot.png`, `{filename}_embedding.png`, `{filename}_umap.png` + `{filename}_coordinates.npy`.
* `{out_dir}/barplot_{metric}.png`, `{out_dir}/swarmplot_{metric}.png` (paths are overridden by the importing modules).

### Environments

* [`scanpy`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scanpy.yaml) — dotplot, embedding, UMAP (and `zarr_to_h5ad`).
* `scanpy_rapids` — UMAP when `use_gpu: true` (no `envs/scanpy_rapids.yaml` is shipped, so with `env_mode: from_yaml` this environment must be provided by the user; inside the script, RAPIDS failures fall back to the default implementations).
* [`plots`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/plots.yaml) — bar and swarm plots.
