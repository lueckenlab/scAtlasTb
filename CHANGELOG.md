# Changelog

## 13.11.2026 Preprocessing: adaptive gene filtering and SVG plots

PR: https://github.com/lueckenlab/scAtlasTb/pull/396

- Plots are written as SVG and label the PCA axes with explained variance.
- Default threads lowered from 10 to 4; HVG computation skips dask when there are many batches (`large_batch_threshold`).
- HVG scripts share common helpers.

### BREAKING CHANGE
- Without `min_cells` or `min_counts`, genes are filtered with `min_cells = max(10, min(200, ceil(0.001 * n_obs)))` instead of `1`.
- `rule all` no longer builds plots; run `plots_all` for them.

## 13.11.2026 Utils: operator queries and simpler sparse dask reads

PR: https://github.com/lueckenlab/scAtlasTb/pull/394

- Wildcard filtering (`subset_by_query`) accepts operator dicts: `eq`, `ne`, `in`, `notin`.
- `_filter_batch` ignores NaN batches and returns an all-`False` mask when no batch is left.
- Sparse dask reads no longer cast index dtypes; sparse and forced-sparse slots (`X`, `layers`, `raw/X`) are read as CSR blocks through one code path.
- `write_zarr` uses `ensure_sparse`, which now also converts dask arrays of `sparse.COO` blocks.

### BREAKING CHANGE
- `write_zarr` converts dense `X`, `layers` and `raw.X` to CSR before writing. Previously, dense matrices were stored dense.

## 02.04.2026 Metrics prepare optimisation

PR: https://github.com/HCA-integration/scAtlasTb/pull/361

Metrics module preparation was optimized to reduce memory and I/O overhead and avoid unnecessary recomputation.

- `prepare_all` now aggregates dedicated subtargets (`pca_all`, `score_genes_all`, `cluster_all`) instead of forcing all heavy preparation steps for every dataset.
- Clustering is skipped automatically for clustering-based metrics when `clustering.precomputed_key` is set.
- Prepared AnnData now keeps only required embeddings in `.obsm` (`X_pca`, `X_emb`) to reduce output size.
- Gene-score random permutations are precomputed once and reused across gene sets in `score_genes`, reducing repeated work.

## 05.05.2025 Storage optimisation when subsetting data

PR: https://github.com/HCA-integration/scAtlasTb/pull/233

Subset dataset by only saving coordinates of the entries that are being subset to. This omits saving a new copy of the data for subsets.

**Performance boost:** Massive storage savings when e.g. data is split after a parameter grid search is performed on those splits. This also saves time spent writing the unnecessary copies, and loading a large matrix and subsetting it is sufficiently fast (since data is not fully loaded).

### BREAKING CHANGE
Zarr files that have been subset via masking, cannot be read by the `anndata.read_zarr`, but instead relies on `utils.read_anndata`, which internally subsets the slots before creating the correct anndata object.

## 14.12.2023 extend relabel model

Extend the model to allow for merging of existing columns.

The old way of writing the config:

```yaml
DATASET:
...
      relabel:
        mapping:
          file:  test/input/mapping_test.tsv
          order:
            - bulk_labels
            - lineage
```

becomes:

```yaml
...
      relabel:
        new_columns:
          file:  test/input/mapping_test.tsv
          order:
            - bulk_labels
            - lineage
```

And the extension is specified as follows:

```yaml
...
      relabel:
        new_columns:
          ...
        merge_columns:
          file:  test/input/merge_test.tsv
          sep: '-'
```


## 7.11.2023 Optimise batch PCR analysis

- parallelise permutations per covariate withing script
- use unfiltered files for merging DCP columns
- optimise Preprocessing io

## 19.10.2023 Handle multiple inputs

### New feature: allow multiple inputs to a module
Allow the user to define multiple inputs for a given task ("dataset").

```yaml
DATASETS:
  test_data:
    input:
      preprocessing: 'data.h5ad'
```

The config still work, but it the output files will now have an additional hash code for the input file.

Additionally, the user can now define a mapping of an input file name and a mapping

```yaml
DATASETS:
  test_data:
    input:
      preprocessing:
        file1: 'data.h5ad'
        file2: 'data2.h5ad'
```


This gets particularly handy, when channeling multiple outputs of a module to another module


```yaml
DATASETS:
  test_data:
    input:
      integration: 'data.h5ad'  # provides multiple output files
      metrics: integration  # will take all the output files of integration as input with human readable input file IDs
```


### Breaking changes

1. The output folder structure will look very different now. In case you are using a custom `outputs.yaml` file, please make sure to adapt it to the latest file patterns in `configs/outputs.yaml`

2. The module `integration_per_lineage` is not supported at the moment and will be replaced by a `split_dataset` module combined with `integration`.

### Under the hood

New classes (`workflow/utils/ModuleConfig.py`, `workflow/utils/WildcardParameters.py`, `workflow/utils/InputFiles.py`, `workflow/integration/IntegrationConfig.py`) are now created to handle config inputs and parse wildcard combinations for the workflow.
They make writing new modules easier and improve code reuse.
Most of the functions in `workflow/utils/` will become redundant (they're still included just in case for now).

For example, the integration `Snakefile` goes from this:

<details>
<summary>Old code</summary>

```python

from utils.misc import all_but, unique_dataframe
from utils.config import get_hyperparams, get_resource, get_params_from_config, set_defaults, get_datasets_for_module, get_for_dataset
from utils.wildcards import expand_per, get_params, get_wildcards, wildcards_to_str
from utils.environments import get_env

module_name = 'integration'
config = set_defaults(config,module_name)
out_dir = Path(config['output_dir']) / module_name
image_dir = Path(config['images']) / module_name

# ... 

parameters = pd.read_table(workflow.source_path('params.tsv'))
parameters['output_type'] = parameters['output_type'].str.split(',')
parameters = get_params_from_config(
    config=get_datasets_for_module(config, module_name),
    module_name=module_name,
    config_params=['methods', 'label', 'batch', 'norm_counts', 'raw_counts'],
    wildcard_names=['dataset', 'method', 'label', 'batch', 'norm_counts', 'raw_counts'],
    defaults=config['defaults'],
    explode_by=['method', 'batch'],
).merge(parameters,on='method')

# subset to datasets that have module defined
parameters = parameters[~parameters['method'].isnull()]

# TODO: remove redundant wildcards
# parameters['label'] = np.where(parameters['use_cell_type'], parameters['label'], 'None')
# parameters = unique_dataframe(parameters)

hyperparams_df = get_hyperparams(config,module_name=module_name)
parameters = parameters.merge(hyperparams_df,on=['dataset', 'method'],how='left')
wildcard_names = ['dataset', 'batch', 'label', 'method', 'hyperparams']

# write hyperparameter mapping
Path(out_dir).mkdir(parents=True, exist_ok=True)
unique_dataframe(
    hyperparams_df[['method', 'hyperparams', 'hyperparams_dict']]
).to_csv(out_dir / 'hyperparams.tsv', sep='\t', index=False)

paramspace = Paramspace(
    parameters[wildcard_names],
    filename_params=['method', 'hyperparams'],
    filename_sep='--',
)
```

</details>


to this:

<details>
<summary>New code</summary>

```python
from utils.environments import get_env
from IntegrationConfig import IntegrationConfig

mcfg = IntegrationConfig(
    module_name='integration',
    config=config,
    parameters=workflow.source_path('params.tsv'),
    config_params=['methods', 'batch', 'label', 'norm_counts', 'raw_counts'],
    wildcard_names=['method', 'batch', 'label'],
    rename_config_params={'methods': 'method'},
    explode_by=['method', 'batch'],
)

out_dir = mcfg.out_dir
image_dir = mcfg.image_dir
paramspace = mcfg.get_paramspace()
```
</details>