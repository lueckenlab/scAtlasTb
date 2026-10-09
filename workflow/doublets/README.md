# Doublet detection

This module runs per-batch doublet detection using one or more callers (currently [`scrublet`](https://github.com/swolock/scrublet) and [`doubletdetection`](https://github.com/JonathanShor/DoubletDetection)). The pipeline runs doublet callers within each library/batch and writes scores and predictions back into the AnnData object.

## Environments

The following environments are useful for running the module. Install only the ones you need.

- [`scanpy`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scanpy.yaml): splitting batches and collecting results
- [`qc`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/qc.yaml): running `scrublet` and `doubletdetection`
- [`rapids_singlecell`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/rapids_singlecell.yaml) (optional): GPU-accelerated `scrublet` when `use_gpu: true` is set for the dataset

# Configuration

```yaml
DATASETS:
  Lee2020:
    input:
      doublets: 
        Lee2020: test/input/load_data/harmonize_metadata/Lee2020.zarr
    doublets:
      counts: X
      batch: donor
      chunk_size: 10_000

  test:
    input:
      doublets:
        test: test/input/pbmc68k.h5ad
        test2: test/input/pbmc68k.h5ad
    doublets:
      counts: layers/counts
      methods:
        - scrublet
        - doubletdetection

defaults:
  datasets:
    - test
    - Lee2020
```

* `counts`: Slot in anndata that contains raw (unnormalized) counts. Examples: `X`, `raw/X`, or `layers/<layer_name>`.
* `batch`: Column in `obs` that contains batch/library IDs. Doublet detection is executed separately per batch.
* `chunk_size`: Number of cells used to group batches for more efficient parallel processing. Default: `100_000`.
* `methods`: Optional list of doublet callers to run for this dataset (`scrublet`, `doubletdetection`). Default: `['scrublet']`.
* `use_gpu`: Whether to run `scrublet` with `rapids_singlecell` on GPU. Default: `false`.

> Note: `counts` is resolved relative to the input object (e.g., anndata.X or anndata.layers). Methods are run per-batch; choose `batch` to reflect library-level grouping so cross-library doublets are not considered.

## Output

Results are written into the AnnData `.obs`:

* `<out_dir>/doublets/dataset~<dataset>/file_id~<file_id>.zarr` — Updated AnnData with caller-specific score and prediction columns in `.obs`:
  - scrublet:
    - `scrublet_score`
    - `scrublet_prediction`
  - doubletdetection:
    - `doubletdetection_score`
    - `doubletdetection_prediction`
* `<out_dir>/doublets/scatter/dataset~<dataset>/file_id~<file_id>/<method>/<batch>.tsv` — Intermediate per-batch scores and predictions of each caller.

Batches with fewer than 100 cells are skipped and receive a score and prediction of `0`.