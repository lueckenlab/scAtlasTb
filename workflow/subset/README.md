# Subset data

Subsetting cells or samples to create smaller, representative datasets for benchmarking and testing.

This module supports common subsetting strategies such as sampling within samples (per-library) or selecting full samples/datasets.

## Environments

Install only the environments you need. For this module you need:

- [`scanpy`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scanpy.yaml)

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
  - Samples smaller than `min_cells_per_sample` (default 100) are skipped.
  - Useful for testing behavior across full-sample composition and batch effects.

### `within_sample`
  - Samples exactly `per_sample` (n_cells_per_sample) cells per sample; samples with fewer cells are skipped. If `per_sample` is not set, the code computes floor(n_cells / n_samples).
  - Samples are shuffled to avoid systematic bias; sampling is done without replacement.
  - If `n_cells` is None the function will default to keeping all cells (no downsampling).


## Testing

A test config is included under `workflow/subset/test/config.yaml` showing usage patterns for both strategies. Example test commands used by the project:

```
conda activate snakemake
bash test/run_test.sh -n       # dry run
bash test/run_test.sh -c2     # actual run with max 2 cores
```

## Output

Updated AnnData and summary files:

- `<out_dir>/subset/dataset~<datasets>/file_id~<file_id>.zarr` — AnnData with subsetting metadata and reduced observation set (for zarr inputs a reference/linking strategy is attempted when possible).
- adata.obs['subset'] — boolean mask (True for retained cells), only set by the `within_sample` strategy.
- adata.uns['subset'] — dictionary containing the subsetting strategy and parameters used (strategy, n_cells, sample_key, per_sample, seed).

> Note: For zarr input stores the pipeline will attempt to reference or link original arrays to avoid copying large arrays where supported by the filesystem/backend. If linking is unsupported some arrays may be copied.
