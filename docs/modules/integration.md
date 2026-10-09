# Integration

```{include} ../../workflow/integration/README.md
:heading-offset: 1
```

## Functional description

### Inputs

* Input files are configured under `DATASETS.<dataset>.input.integration` (`file_id` → path, `.h5ad` or `.zarr`).
* `.X`/`.layers`: the slots named by `norm_counts` (log-normalised expression) and `raw_counts` (raw counts) are both required; the prepare step asserts that neither is `None`. If `raw_counts` points to `raw/X`, `raw/var` is used as feature table.
* `.var[<var_mask>]`: boolean feature mask, one run per configured value; its existence is asserted. If `var_mask` is not set, all genes are used.
* `.obs[<batch>]` (mandatory) and `.obs[<label>]` (only for methods with `use_cell_type: true`), plus any covariates named in method hyperparameters (e.g. `covariates` for ComBat, `categorical_covariate_keys` for scVI-family models, `system_key` for SysVI, `key`/`batch_key` for Harmony).
* `.obsm['X_pca']`, `.obsp`/`.uns['neighbors']` (optional): reused by `unintegrated`; all other methods recompute what they need.

Method metadata comes from `params.tsv` (columns `output_type`, `use_cell_type`, `no_scale`, `env`, `resources`, `cpu_env`, `do_not_expand`). `IntegrationConfig` builds the parameter space:

```text
for each method in params.tsv:
    if use_gpu is false: env ← cpu_env (if set); resources ← cpu
    explode comma-separated output_type into one row per output type
for each dataset × configured method: expand list-valued hyperparameters into all combinations
    (keys listed in do_not_expand stay lists); hyperparams wildcard ← hash of combination
    (method without hyperparameters → 'None'); mapping written to <output_dir>/integration/hyperparams.tsv
explode batch, label, var_mask
label ← 'None' for methods with use_cell_type == false   # avoids redundant runs
drop rows whose output_type is not in the dataset's output_types (default: all of the method's types)
```

The `no_scale` column of `params.tsv` is currently not read by any rule or script.

Wildcards: `dataset`, `file_id`, `batch`, `var_mask`, `method`, `hyperparams`, `label`, `output_type` (see {ref}`architecture`).

### Processing steps

#### Methods

| method | what it does | input | output type(s) | env (GPU / CPU) |
| --- | --- | --- | --- | --- |
| `unintegrated` | no correction; reuses input `X_pca` or computes `scanpy.pp.pca` on integration features as `X_emb`; reuses input kNN graph if valid, else `scanpy.pp.neighbors(use_rep='X_emb')` | normalized | full, embed, knn | `rapids_singlecell` / `scanpy` |
| `bbknn` | `scanpy.pp.pca(**PCA params)` then batch-balanced kNN `scanpy.external.pp.bbknn(batch_key, use_rep='X_pca', pynndescent_random_state=seed)` ([BBKNN](https://github.com/Teichlab/bbknn)) | normalized | knn | `bbknn` / `bbknn` |
| `combat` | `scanpy.pp.combat(key=batch, covariates=...)` linear empirical-Bayes correction (ComBat) | normalized | full | `scanpy` / `scanpy` |
| `harmony_pytorch` | optional scaling and PCA, then `harmony.harmonize(X_pca, obs, batch_key=[batch, ...], random_state=seed)` ([harmony-pytorch](https://github.com/lilab-bcb/harmony-pytorch)) | normalized (+ counts for MCV) | embed | `rapids_singlecell` / `scanpy` |
| `harmonypy` | PCA, then `harmony_integrate(basis='X_pca', adjusted_basis='X_emb', key=[batch, ...])` from `rapids_singlecell` on GPU or `scanpy.external` ([harmonypy](https://github.com/slowkow/harmonypy)) | normalized | embed | `rapids_singlecell` / `scanpy` |
| `scanorama` | split by batch, `scanorama.correct_scanpy(return_dimred=True, seed=seed)` ([Scanorama](https://github.com/brianhie/scanorama)), concatenate | normalized | embed, full | `scanorama` / `scanorama` |
| `scvi` | `scvi.model.SCVI` VAE with batch (and covariates) ([scvi-tools](https://scvi-tools.org)); latent mean as `X_emb` | counts | embed | `scvi-tools` / `scvi-tools` |
| `scanvi` | scVI pre-training, then `SCANVI.from_scvi_model(labels_key=label, unlabeled_category='nan')` semi-supervised training | counts | embed | `scvi-tools` / `scvi-tools` |
| `scvi_vitkl` | same script logic as `scvi`, run in the `scvi-vitkl` environment (fork of scvi-tools, accepts additional model parameters `library_n_hidden`, `use_additive_background`, ...) | counts | embed | `scvi-vitkl` / `scvi-vitkl` |
| `drvi` | `drvi.model.DRVI` disentangled VAE ([DRVI](https://github.com/theislab/drvi)), `is_count_data=True` | counts | embed | `scvi-tools` / `scvi-tools` |
| `sysvi` | `scvi.external.SysVI` cVAE for substantial batch effects ([sysVI](https://github.com/theislab/sysVI)); `system_key` (mandatory) is the model batch, `batch` is added as categorical covariate | normalized | embed | `scvi-tools` / `scvi-tools` |
| `scpoli` | `scarches.models.scPoli` prototype-based cVAE ([scArches](https://github.com/theislab/scarches)), condition = batch, cell type = label | counts | embed | `scarches` / `scarches` |
| `scgen` | `scarches.models.scgen` VAE, `batch_removal(batch_key, cell_label_key=label)` | normalized (dense) | embed, full | `scarches` / `scarches` |

Output types: **full** — a corrected feature matrix in `.X` (kNN graph later computed on a PCA of it); **embed** — a corrected low-dimensional embedding in `.obsm['X_emb']` (kNN graph later computed on it); **knn** — a corrected kNN graph in `.obsp['connectivities']`/`.obsp['distances']` and `.uns['neighbors']`, used as is. Methods with several output types are run once; the output types differ only in postprocessing.

#### 1. `integration_prepare`

Script `scripts/prepare.py`, environment `scanpy`, one job per `dataset × file_id × var_mask`. No normalisation, scaling or HVG selection is performed here; it only harmonises slots:

```text
if var_mask is set: assert var_mask in input .var
norm ← read(norm_counts, backed, dask)
var['integration_features'] ← var[var_mask] (all True if var_mask is None)
var['integration_features'] &= genes with ≥ 1 non-zero cell among masked genes   # _filter_genes(min_cells=1)
raw ← read(raw_counts)
assert identical var_names and n_obs for norm and raw
if save_subset: subset both matrices to integration_features and write them; else link them
output: X ≙ layers['normcounts'] ← norm_counts, layers['counts'] ← raw_counts, obs/var/obsm/obsp/uns from input
```

#### 2. `integration_run_method`

Script `scripts/methods/<method>.py` (selected by the `method` wildcard), environment from `params.tsv` (`env`, or `cpu_env` when `use_gpu: false`). Common steps in every script (helpers in `scripts/methods/integration_utils.py`):

```text
read layers/normcounts or layers/counts (see table) + obs/var/uns, lazily via dask
clean_categorical_column(batch [, label, covariates])        # cast to str → category
subset_hvg(var_column='integration_features')                 # feature subset
split hyperparams into model/PCA/train arguments (get_hyperparams with *_MODEL_PARAMS / PCA_PARAMS lists)
run method (seeded with `seed`)
remove_slots(): drop layers and X_pca; keep X only for full outputs (or deep-learning models, linked)
add_metadata(): uns['integration'] = {method, label_key, batch_key, output_type, hyperparams[, model_history]}
                uns['wildcards'] |= integration_<wildcard>
write only new slots, link the rest to the prepared file
```

Method-specific behaviour:

* **PCA-based** (`bbknn`, `harmonypy`, `harmony_pytorch`): hyperparameters named in `PCA_PARAMS` (`n_comps`, `svd_solver`, `mask_var`, `zero_center`, ...) go to `scanpy.pp.pca`, the rest to the method. `harmony_pytorch` defaults to `svd_solver='covariance_eigh'`, applies `scanpy.pp.scale(zero_center=False)` if hyperparameter `scale: true`, and with `n_comps: mcv` selects the number of PCs by molecular cross-validation (`scripts/methods/mcv.py`: assign each non-zero count entry with probability 0.5 to one of two matrices (seed 42), normalise both to 10 000 + log1p (+ scale if `scale`), PCA with 100 components on the first, choose k ∈ {2…9, 10…28 step 2, 30…100 step 5} minimising the MSE between the rank-k reconstruction and the second matrix; plot `plots/pca_mcv.png`). Harmony batch keys are the union of `batch` and the `key`/`batch_key` hyperparameter. `harmonypy` sets `dtype='float32'` on GPU. `bbknn` removes batches with fewer cells than `neighbors_within_batch` (default 3) before running, so its output may contain fewer cells.
* **scvi-tools models** (`scvi`, `scanvi`, `scvi_vitkl`, `drvi`, `sysvi`): `scvi.settings.seed = seed`, `torch.set_float32_matmul_precision('medium')`. Keys in `SCVI_MODEL_PARAMS` / `SYSVI_MODEL_PARAMS` / `DRVI_MODEL_PARAMS` go to the model constructor; `categorical_covariate_keys`/`continuous_covariate_keys` go to `setup_anndata(batch_key=batch, ...)`; all remaining keys are passed to `model.train()` with default `check_val_every_n_epoch=1`. For `scanvi`, training keys prefixed `scanvi_` are passed (without prefix) to the scANVI stage, which uses a checkpoint callback on `elbo_validation` and loads the best checkpoint on failure. `drvi` sets `plan_kwargs={'n_epochs_kl_warmup': max_epochs (default 400)}` and plots latent-dimension statistics. `sysvi` maps `early_stopping: true` to validation every epoch. Models are saved to `model/`, training curves to `plots/training_metrics_<group>.png`, and `model.history` to `uns['integration']['model_history']`.
* **scArches models** (`scpoli`, `scgen`): `torch.manual_seed(seed)`, threads via `torch.set_num_threads`. scPoli defaults: `condition_keys=[batch]`, `cell_type_keys=[label]`, `unknown_ct_names=['NA']`, `early_stopping_kwargs` (metric `val_prototype_loss`, patience 20, `reduce_lr`, `lr_patience` 13, `lr_factor` 0.1), `pretrain_epochs = int(0.8 · n_epochs)`; embedding = `model.get_latent(mean=True)`, condition embeddings in `uns['scpoli_<batch>_embeddings']`. scGen densifies `X`, trains, and stores `batch_removal` output as corrected `X` and `obsm['latent_corrected']` as `X_emb`.

Resources: GPU methods (`resources` = `gpu`) run with 1 thread on the GPU profile for the first 2 attempts and fall back to the CPU profile (and `threads` threads) on later retries; time limit 2 days. Snakemake records runtime/memory in `benchmark.tsv`.

#### 3. `integration_postprocess`

Re-uses the preprocessing `neighbors` rule (`scripts/neighbors.py` from the preprocessing module, environment `scanpy` / `rapids_singlecell`, 2 retries), once per output type:

```text
args ← integration.neighbors config (default {}; n_neighbors default 15, capped at n_obs)
full  → use_rep='X_pca'; PCA is computed on the corrected X (arguments from uns['preprocessing']['pca'] if present)
embed → use_rep='X_emb'
knn   → args=False: keep the method's graph, only validate it (assert_neighbors)
sc.pp.neighbors(**args); uns['output_type'] ← output_type
```

#### 4. `integration_compute_umap`, `integration_plot_umap`

`integration_compute_umap` re-uses the preprocessing `umap` rule (`scripts/umap.py`, `scanpy` / `rapids_singlecell`) on the postprocessed file (see the preprocessing module for details). `integration_plot_umap` re-uses the preprocessing `plots` rule (`scripts/plot.py`, `scanpy`) with basis `X_umap`, colours = label(s) + batch(es) + `umap_colors` + `plots.colors`, centroids from `plots.plot_centroids`, gene chunk size `plots.plot_gene_chunk_size` (default 12), and minimum category size 0.01 % of cells.

#### 5. `integration_benchmark`, `integration_benchmark_per_dataset`, `integration_barplot`, `integration_barplot_per_dataset`

The Snakemake benchmark files of all `integration_run_method` jobs are concatenated with their wildcards into one TSV over all datasets and one per dataset (local rules). `integration_barplot*` (common `barplot` rule, `scripts/barplot.py` in `plots`) plot the metrics `s` (wall-clock seconds), `max_uss` (peak unique memory, MB) and `mean_load` (CPU load) per method, faceted by dataset.

### Outputs

Under `<output_dir>/integration/` (`<pattern>` = `dataset~<dataset>/file_id~<file_id>/batch~<batch>/var_mask~<var_mask>/method~<method>--hyperparams~<hash>--label~<label>--output_type~<output_type>`):

* `prepare/dataset~<dataset>/file_id~<file_id>/var_mask~<var_mask>.zarr`: `X`, `layers['normcounts']`, `layers['counts']`, `var['integration_features']`.
* `run_method/<pattern without --output_type>/adata.zarr`: raw method output — `X` (full), `obsm['X_emb']` (embed), `obsp['connectivities']`, `obsp['distances']`, `uns['neighbors']` (knn), `uns['integration']`, `uns['wildcards']`; plus `model/` (saved model, deep-learning methods), `plots/` (training curves, MCV, DRVI latent statistics) and `benchmark.tsv`.
* `postprocess/<pattern>.zarr`: adds `obsp['connectivities']`, `obsp['distances']`, `uns['neighbors']`, `uns['output_type']` (and `obsm['X_pca']`, `varm['PCs']` for full outputs).
* `<pattern>.zarr`: final file with `obsm['X_umap']` (default target, input for downstream modules such as metrics and clustering).
* `hyperparams.tsv` (hash → hyperparameters), `input_files.tsv`, `output_files.tsv`, `integration.benchmark.tsv`, `dataset~<dataset>/integration.benchmark.tsv`.

Under `<images>/integration/`:

* `umap/<pattern with /output_type~...>/<color>.png` and `genes_group=<i>.png`.
* `benchmark/<metric>.png` and `benchmark/dataset~<dataset>/metric~<metric>.png` for `s`, `max_uss`, `mean_load`.

### Environments

* `scanpy`: prepare, UMAP plots; CPU environment of `unintegrated`, `harmony_pytorch`, `harmonypy`; `combat`; postprocess and UMAP when `use_gpu: false`.
* `rapids_singlecell` (only with `use_gpu: true`): `unintegrated`, `harmony_pytorch`, `harmonypy`, postprocess (neighbors) and UMAP computation.
* `bbknn`, `scanorama`: the respective methods (CPU and GPU mode).
* `scvi-tools`: `scvi`, `scanvi`, `drvi`, `sysvi`; `scvi-vitkl`: `scvi_vitkl`; `scarches`: `scpoli`, `scgen`. These use a GPU when available (`resources: gpu`), but have no CPU-specific environment.
* `plots`: benchmark barplots.
