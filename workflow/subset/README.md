# Subset data

Subsetting cells or samples to create smaller, representative datasets for benchmarking and testing.

This module supports common subsetting strategies such as sampling within samples (per-library) or selecting full samples/datasets.

```mermaid
flowchart TD
  in(["AnnData with sample column"])
  shuf["Shuffle samples<br/>(reproducible, seed)"]
  bs["Add whole samples until n_cells is exceeded<br/>(skip samples with fewer than 100 cells)"]
  ws["Draw k random cells per sample until n_cells<br/>k = per_sample or n_cells / n_samples<br/>(skip samples with fewer than k cells)"]
  sub["Subset to selected cells<br/>(all genes kept)"]
  out(["Subsetted AnnData"])
  in --> shuf
  shuf -->|by_sample| bs
  shuf -->|within_sample| ws
  bs --> sub
  ws --> sub
  sub --> out
```

*Conceptual overview of the main steps of the module. See the [functional description](#functional-description) below for details.*

## Configuration

Example dataset configuration (keys match the implementation and test config):

```yaml
DATASETS:
  by_sample_subset:
    input:
      subset: path/to/input.h5ad
    subset:
      sample_key: sample
      seed: 42
      n_cells: 1000
      strategy: by_sample
  
  within_sample_subset:
    input:
      subset: path/to/input.h5ad
    subset:
      sample_key: sample
      seed: 42
      n_cells: 1000
      per_sample: 100
      strategy: within_sample
```

Notes about parameter names and mapping:
- `sample_key`: the key in `adata.obs` that defines the samples/libraries to subset over.
- `per_sample`: number of cells to keep per sample when using `within_sample` strategy.
- `seed`: random seed for reproducible subsetting.
- `strategy`: must match one of the implemented strategies: `by_sample` or `within_sample`.

## Strategies

### `by_sample`
  - Selects entire samples (libraries) and keeps adding shuffled samples until the `n_cells` (n_cell_max) is reached.
  - Samples smaller than `min_cells_per_sample` (100 cells, not configurable) are skipped.
  - Useful for testing behavior across full-sample composition and batch effects.

### `within_sample`
  - Samples exactly `per_sample` (n_cells_per_sample) cells per sample; samples with fewer cells are skipped. If `per_sample` is not set, the code computes floor(n_cells / n_samples).
  - Samples are shuffled to avoid systematic bias; sampling is done without replacement.
  - If `n_cells` is not set, there is no upper limit on the total number of cells; `per_sample` must then be set, otherwise no cells are selected and the job fails.


## Functional description

### Inputs

* **File formats:** `.h5ad` or `.zarr`, configured under `input: subset:` (see [architecture](https://scatlastb.readthedocs.io/en/latest/principles/architecture.html)). The whole object is opened lazily (`read_anndata(backed=True, dask=True)`); only `.obs` is used for the decision.
* **`.obs[sample_key]`** (config `sample_key`, no default; the job fails with an assertion if the column is missing).
* **Config keys:** `strategy` (`by_sample` or `within_sample`), `n_cells` (if unset: no upper limit, `np.iinfo(int).max`), `per_sample` (only `within_sample`), `seed` (default `42`). `label` is accepted in the configuration but not used by the script.

### Processing steps

1. **`subset_subset`** (script `scripts/run.py`, functions in `scripts/subset_functions.py`, environment `scanpy`). One job per task × input file.
   Samples are first shuffled reproducibly: `obs[sample_key].value_counts().sample(frac=1, random_state=seed)`.

   *`by_sample`* — keeps complete samples:
   ```text
   selected = []; n = 0
   for sample, count in shuffled samples:
       if selected and n > n_cells: stop
       if count < 100: skip                       # min_cells_per_sample, not configurable
       selected.append(sample); n += count
   mask = obs[sample_key] in selected
   ```
   The last sample that crosses `n_cells` is still added, so the result can exceed `n_cells` by up to one sample.

   *`within_sample`* — draws a fixed number of cells per sample:
   ```text
   k = per_sample or int(n_cells / n_samples)
   for sample in shuffled samples:
       cells = obs_names of sample
       if len(cells) < k: skip sample             # small samples are dropped, not down-weighted
       draw k cells without replacement (pandas .sample, random_state=seed)
       if total drawn >= n_cells: stop
   obs['subset'] = True for drawn cells, False otherwise; mask = obs['subset']
   ```
   Because `run.py` replaces an unset `n_cells` by the maximum integer, setting neither `n_cells` nor `per_sample` makes `k` larger than any sample, so no cells are selected and the job fails.
   (A third entry, `scarf_TopACeDo`, exists in the strategy map but is an unimplemented stub.)

   The AnnData is subset to the mask and `.uns['subset']` is set. Only `.obs` and `.uns` are written; all other slots of a `.zarr` input are symlinked and the cell mask is stored alongside them so they are subset on read (`write_zarr_linked(subset_mask=(mask, all genes))`). For `.h5ad` input the full subset is written instead.

### Outputs

Relative to `<output_dir>/subset/`:

* `dataset~<dataset>/file_id~<file_id>.zarr` — AnnData restricted to the selected cells (all genes kept) with
  * `.uns['subset'] = {strategy, n_cells, sample_key, per_sample, seed}`,
  * `.obs['subset']` (all `True`) — only for `within_sample`.

For `.zarr` inputs, all slots other than `.obs` and `.uns` are linked to the input arrays instead of being copied (see step 1); `.h5ad` inputs are written as a full subset copy.

### Environments

* [`scanpy`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scanpy.yaml)

## Testing

A test config is included under `workflow/subset/test/config.yaml` showing usage patterns for both strategies. Example test commands used by the project:

```
conda activate snakemake
bash test/run_test.sh -n       # dry run
bash test/run_test.sh -c2     # actual run with max 2 cores
```
