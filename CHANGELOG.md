# Changelog

<!--
Format:
- One bullet per PR under the year heading, newest first: `**DD.MM.YYYY** Short summary ([#PR](link))`.
  The date is the day the PR was merged.
- Add indented sub-bullets only when needed (e.g. notable details, `**Breaking:**` notes).
- Larger, manually written entries on breaking changes get their own `###` section below the year's list
  (with bold labels such as `**Breaking changes**` to structure it), and the list bullet links to that section.
-->

This changelog lists new features, major changes and breaking changes.
Bug fixes and minor changes are listed in the [merged pull requests](https://github.com/lueckenlab/scAtlasTb/pulls?q=is%3Apr+is%3Amerged).

## 2026

- **21.11.2026** Marker genes: `plot.group_order` and matrix plots for user marker genes, configurable `min_cells`/`min_groups_size` and `fail_on_empty` ([#395](https://github.com/lueckenlab/scAtlasTb/pull/395))
  - **Breaking:** `rank_genes_groups` fails by default when `min_cells` or `min_groups_size` filter out all cells. Set `fail_on_empty: false` to write an empty result instead.
- **13.11.2026** Utils: operator queries in `WildcardParameters`, NaN-safe batch filter, simpler sparse dask reads ([#394](https://github.com/lueckenlab/scAtlasTb/pull/394))
  - Wildcard filtering (`subset_by_query`) accepts operator dicts: `eq`, `ne`, `in`, `notin`.
  - `_filter_batch` ignores NaN batches and returns an all-`False` mask when no batch is left.
  - Sparse dask reads no longer cast index dtypes; sparse and forced-sparse slots (`X`, `layers`, `raw/X`) are read as CSR blocks through one code path.
  - `write_zarr` uses `ensure_sparse`, which now also converts dask arrays of `sparse.COO` blocks.
  - **Breaking:** `write_zarr` converts dense `X`, `layers` and `raw.X` to CSR before writing. Previously, dense matrices were stored dense.
- **01.06.2026** QC: consolidate plots into one SVG per QC metric pair, with joint plots and removed-cell summaries per `file_id` ([#372](https://github.com/lueckenlab/scAtlasTb/pull/372))
- **20.05.2026** Batch analysis: add pseudobulk-level batch analysis ([#367](https://github.com/lueckenlab/scAtlasTb/pull/367))
- **14.04.2026** Batch analysis: Theil's U heatmap for associations between covariates ([#364](https://github.com/lueckenlab/scAtlasTb/pull/364))
- **02.04.2026** Metrics: optimise preparation to reduce memory and I/O overhead and avoid unnecessary recomputation ([#361](https://github.com/lueckenlab/scAtlasTb/pull/361))
  - `prepare_all` now aggregates dedicated subtargets (`pca_all`, `score_genes_all`, `cluster_all`) instead of forcing all heavy preparation steps for every dataset.
  - Clustering is skipped automatically for clustering-based metrics when `clustering.precomputed_key` is set.
  - Prepared AnnData now keeps only required embeddings in `.obsm` (`X_pca`, `X_emb`) to reduce output size.
  - Gene-score random permutations are precomputed once and reused across gene sets in `score_genes`.
- **02.04.2026** Metrics: more configurable funkyheatmap (`funkyheatmap` block with `group_col`, `weight_batch`, `n_top`, `scale`, `dpi`) ([#362](https://github.com/lueckenlab/scAtlasTb/pull/362))
  - **Breaking:** `overwrite_file_id` is renamed to `parse_file_id`.
- **15.03.2026** Preprocessing: support additional HVG flavors ([#358](https://github.com/lueckenlab/scAtlasTb/pull/358))
- **15.03.2026** Utils: shorter output file names and `defaults: datasets` to restrict which datasets are run ([#354](https://github.com/lueckenlab/scAtlasTb/pull/354))
  - **Breaking:** File names longer than 150 characters are shortened by hashing the wildcards of each module, so output paths of long names change.
- **11.01.2026** Cell type prediction: predict donor sex from XIST and chrY gene expression, and from the ratio of non-pseudoautosomal to pseudoautosomal chrY genes ([#323](https://github.com/lueckenlab/scAtlasTb/pull/323), [#326](https://github.com/lueckenlab/scAtlasTb/pull/326))

## 2025

- **19.10.2025** Preprocessing: streamline output paths ([#299](https://github.com/lueckenlab/scAtlasTb/pull/299), [#304](https://github.com/lueckenlab/scAtlasTb/pull/304))
  - **Breaking:** The final output is written to `<wildcards>.zarr` instead of `<wildcards>/preprocessed.zarr`, and intermediate files are written under `preprocessed/`. The legacy path is no longer created.
- **19.10.2025** Subset: streamline output ([#298](https://github.com/lueckenlab/scAtlasTb/pull/298))
  - **Breaking:** The config key `sample` is renamed to `sample_key`, and the `strategy` wildcard is removed from output paths.
- **09.10.2025** Reference mapping: new module for mapping to a reference with scArches ([#292](https://github.com/lueckenlab/scAtlasTb/pull/292))
- **01.10.2025** Integration: molecular cross-validation (MCV) to choose the number of PCs for Harmony ([#289](https://github.com/lueckenlab/scAtlasTb/pull/289))
- **16.09.2025** Batch PCR: optional preprocessing pipeline (normalisation, gene filtering, PCA) before computing PCR ([#286](https://github.com/lueckenlab/scAtlasTb/pull/286))
- **12.09.2025** Metrics: add BRAS metric from `scib_metrics` ([#284](https://github.com/lueckenlab/scAtlasTb/pull/284))
- **10.09.2025** Metrics: configurable `output_type` for input files without `.uns['output_type']` ([#285](https://github.com/lueckenlab/scAtlasTb/pull/285))
  - The `methods` key under `metrics` is deprecated in favour of `metrics`.
- **27.08.2025** Preprocessing: dask-based PCA ([#279](https://github.com/lueckenlab/scAtlasTb/pull/279))
- **25.08.2025** Integration: add scVI-VitKL model ([#267](https://github.com/lueckenlab/scAtlasTb/pull/267))
- **05.05.2025** Storage optimisation when subsetting data: subsets only save the coordinates of the selected entries instead of a new copy of the data ([#233](https://github.com/lueckenlab/scAtlasTb/pull/233))
  - Massive storage savings, e.g. when data is split for a parameter grid search, and no time spent writing unnecessary copies. Loading a large matrix and subsetting it is sufficiently fast, since data is not fully loaded.
  - **Breaking:** Zarr files that have been subset via masking cannot be read by `anndata.read_zarr`. Use `utils.read_anndata` instead, which subsets the slots before creating the AnnData object.
- **23.02.2025** Metrics: gene set scoring, gene set-based metrics, custom covariates and a collected metrics AnnData as default output ([#221](https://github.com/lueckenlab/scAtlasTb/pull/221))
- **24.01.2025** Integration: add DRVI and sysVI ([#216](https://github.com/lueckenlab/scAtlasTb/pull/216))
- **14.01.2025** Filter: remove cells by pandas query with `remove_by_query` ([#213](https://github.com/lueckenlab/scAtlasTb/pull/213))

## 2024

- **13.12.2024** Metrics: Moran's I for covariates and gene scores, with baselines and bootstrapping ([#184](https://github.com/lueckenlab/scAtlasTb/pull/184))
  - **Breaking:** The metric `pcr` is renamed to `pcr_regression`.
- **16.10.2024** Sample representation: framework for multiple patient representation methods and hyperparameters, with plots ([#199](https://github.com/lueckenlab/scAtlasTb/pull/199))
  - **Breaking:** New output paths for `sample_representation`, and the `subset` module now supports multiple subset strategies.
- **25.09.2024** Uncollect: new module to split a file into multiple files, each with a dedicated subset of slots ([#196](https://github.com/lueckenlab/scAtlasTb/pull/196))
- **19.09.2024** Label transfer: majority voting ([#192](https://github.com/lueckenlab/scAtlasTb/pull/192))
- **19.09.2024** Integration: rename counts layers to match the preprocessing module ([#193](https://github.com/lueckenlab/scAtlasTb/pull/193))
  - **Breaking:** Integration outputs use `layers/counts` and `layers/normcounts` instead of `layers/raw_counts` and `layers/norm_counts`.
- **17.09.2024** Collect: new module to collect multiple input files into a single file, optimised for zarr by linking large slots ([#191](https://github.com/lueckenlab/scAtlasTb/pull/191))
- **14.05.2024** Marker genes: new module ([#158](https://github.com/lueckenlab/scAtlasTb/pull/158))
- **03.04.2024** Integration: feature masking via `var_mask` as an expandable wildcard ([#138](https://github.com/lueckenlab/scAtlasTb/pull/138))
- **11.03.2024** Merge: new module to merge multiple input files ([#121](https://github.com/lueckenlab/scAtlasTb/pull/121))
- **05.03.2024** Filter: new module ([#111](https://github.com/lueckenlab/scAtlasTb/pull/111))
- **19.02.2024** Preprocessing: read data with dask and use dask-ml for normalisation, filtering and PCA ([#98](https://github.com/lueckenlab/scAtlasTb/pull/98))
- **11.01.2024** QC: automatic thresholds with AutoQC ([#90](https://github.com/lueckenlab/scAtlasTb/pull/90))
- **02.01.2024** Integration and preprocessing: more efficient loading of count matrices without loading all matrices ([#87](https://github.com/lueckenlab/scAtlasTb/pull/87))

## 2023

- **14.12.2023** Extend relabel module to allow merging of existing columns, see [details](#14122023-extend-relabel-module) ([#78](https://github.com/lueckenlab/scAtlasTb/pull/78))
  - **Breaking:** `relabel.mapping` is renamed to `relabel.new_columns`.
- **11.12.2023** Utils: backed reading of `.X` and `.layers` for h5py and zarr files ([#75](https://github.com/lueckenlab/scAtlasTb/pull/75))
- **05.12.2023** Doublets: new module ([#73](https://github.com/lueckenlab/scAtlasTb/pull/73))
- **30.11.2023** Sample representation: new module, starting with pseudobulk ([#71](https://github.com/lueckenlab/scAtlasTb/pull/71))
- **23.11.2023** QC: new module with QC plots, configurable thresholds and doublet detection ([#70](https://github.com/lueckenlab/scAtlasTb/pull/70))
- **14.11.2023** Split data: new module ([#69](https://github.com/lueckenlab/scAtlasTb/pull/69))
  - **Breaking:** The `integration_per_lineage` module is removed. Use `split_data` combined with `integration` instead.
- **10.11.2023** Integration and metrics: refactor output types ([#68](https://github.com/lueckenlab/scAtlasTb/pull/68))
  - **Breaking:** `output_type` is a wildcard of the integration module, metrics use a single output type per file, and kNN graphs are saved under standard slots.
- **07.11.2023** Optimise batch PCR analysis ([#67](https://github.com/lueckenlab/scAtlasTb/pull/67))
  - Parallelise permutations per covariate within script.
  - Use unfiltered files for merging DCP columns.
  - Optimise preprocessing I/O.
- **31.10.2023** Preprocessing: GPU support with `rapids-singlecell` ([#63](https://github.com/lueckenlab/scAtlasTb/pull/63))
- **19.10.2023** Allow multiple inputs per module, see [details](#19102023-allow-multiple-inputs-to-a-module) ([#55](https://github.com/lueckenlab/scAtlasTb/pull/55))
  - **Breaking:** New output folder structure; `integration_per_lineage` is not supported at the moment.
- **03.09.2023** Integration: refactor into separate `clustering` and `integration_per_lineage` modules ([#50](https://github.com/lueckenlab/scAtlasTb/pull/50))
- **10.08.2023** Integration and metrics: store AnnData as zarr, and save the unintegrated matrix in a separate file ([#37](https://github.com/lueckenlab/scAtlasTb/pull/37))

### 14.12.2023 Extend relabel module

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

### 19.10.2023 Allow multiple inputs to a module

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


**Breaking changes**

1. The output folder structure will look very different now. In case you are using a custom `outputs.yaml` file, please make sure to adapt it to the latest file patterns in `configs/outputs.yaml`

2. The module `integration_per_lineage` is not supported at the moment and will be replaced by a `split_dataset` module combined with `integration`.

**Under the hood**

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
