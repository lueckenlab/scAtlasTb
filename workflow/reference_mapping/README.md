# Reference mapping

This module enables projection of new single-cell RNA-seq datasets onto pre-trained variational autoencoder (VAE) models or foundation models. This allows for:

- Transfer learning from large reference datasets to smaller query datasets
- Integration of new data with existing atlases
- Consistent cell type annotation across studies
- Leveraging pre-trained embeddings for downstream analysis

## Environments

The following environments are needed for reference mapping:

- [`scvi-tools`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scvi-tools.yaml): mapping with scArches (`scarches` rule)
- [`scanpy`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/scanpy.yaml) or [`rapids_singlecell`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/rapids_singlecell.yaml): neighbors, UMAP and plots (reused from the preprocessing module)

## Supported Models

The module currently supports:
- **scvi-tools models**: Various variational inference models (scVI, scANVI, etc.) mapping with scArches
- TODO: **scArches models**: VAE-based models (trVAE, scPoli, Expimap, etc.)
- TODO: **Foundation models**: Large-scale pre-trained models

## Configuration

```yaml
DATASETS:
  test:
    input:
      reference_mapping:
        query: test/input/pbmc68k.h5ad

    reference_mapping:
      scarches:
        layer: layers/counts  # or X
        var_key: feature_name  # optional
        model: test/input/model
        model_params:
          batch_key: sample_id
          labels_key: cell_type
          categorical_covariate_keys: [donor, condition]
          continuous_covariate_keys: [age]
        train_params:  # set to false for inference only (no query training)
          max_epochs: 10
          early_stopping: true
          check_val_every_n_epoch: 1
      neighbors:  # optional, passed to the kNN graph computation
        n_neighbors: 30
      umap_colors:  # optional, .obs columns to color the UMAP by
        - cell_type
```

### Input

The input AnnData object is the query dataset that should be mapped to the reference model.
The reference model should be defined as a Pytorch model directory under `scarches > model`.

### Configuration Parameters

- **`layer`** (default: `'X'`): Which slot of the query AnnData object to use as input counts, given as a path in the file
  - `'X'` uses the main expression matrix (`.X`)
  - `'layers/counts'` uses the counts layer from `.layers['counts']` (the `layers/` prefix is required)
  - `'raw/X'` uses `.raw.X`

- **`var_key`** (default: `None`): Column name in `.var` to use for gene matching between query and reference model
  - If `None`, uses the `.var` index (var_names)
  - Important for ensuring gene names are consistent between query and reference model

#### Model Parameters (`model_params`)
Parameters that align the query data structure with the reference model's expectations:

- **`batch_key`**: Column name in `.obs` containing batch/sample information (required)
- **`labels_key`**: Column name in `.obs` containing cell type labels (optional, can be `None` for unlabeled data)
- **`size_factor_key`**: Column name in `.obs` containing size factors (optional, only if the reference model uses one)
- **`categorical_covariate_keys`**: List of categorical covariate column names in `.obs` (optional, e.g., `["donor", "condition"]`)
- **`continuous_covariate_keys`**: List of continuous covariate column names in `.obs` (optional, e.g., `["age", "BMI"]`)
- Any remaining keys are passed to `load_query_data` of the model class.

#### Training Parameters (`train_params`)
Parameters that control the reference mapping training process (all optional). They are passed to the model's `train()` method on top of the defaults below.
Set `train_params: false` to skip training and only compute the latent representation of the query with the unchanged reference model.

- **`max_epochs`**: Maximum number of training epochs (default: 10)
- **`check_val_every_n_epoch`**: How often to run validation during training (default: 1, needed for loss curves)
- **`plan_kwargs`**: Training plan arguments (default: `{weight_decay: 0.0}`)
- Any other argument accepted by the model's `train()` method, e.g. **`early_stopping`**

> **Note:** The query data must have genes that overlap with the reference model's training data. Gene matching is performed using the `var_key` parameter or `.var` index. Only overlapping genes are used for mapping and missing genes from the query are padded with zeros.

## Output

The reference mapping workflow produces the following outputs:

* `<out_dir>/reference_mapping/dataset~<dataset>/file_id~<file_id>.zarr`: Mapped AnnData object containing:
  - **Latent embedding** (`obsm['X_emb']`): Low-dimensional representations in the reference space
  - **kNN graph** (`obsp['connectivities']`, `obsp['distances']`, `uns['neighbors']`) computed on `X_emb`
  - **UMAP** (`obsm['X_umap']`)

* `<out_dir>/reference_mapping/model/dataset~<dataset>/file_id~<file_id>/`: Updated model directory with query data integrated, ready for further reference mapping tasks
* `<images>/reference_mapping/model/dataset~<dataset>/file_id~<file_id>/`: training loss curves (only when training)
* `<images>/reference_mapping/umap/dataset~<dataset>/file_id~<file_id>/`: UMAP plots colored by `umap_colors`
