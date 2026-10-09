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
        layer: X  # or layer name like 'counts'
        model: test/input/model
        model_params:
          batch_key: sample_id
          labels_key: cell_type
          categorical_covariate: [donor, condition]
          continuous_covariate: [age]
        train_kwargs:
          max_epochs: 10
          early_stopping: true
          check_val_every_n_epoch: 1
```

### Input

The input AnnData object is the query dataset that should be mapped to the reference model.
The reference model should be defined as a Pytorch model directory under `scarches > model`.

### Configuration Parameters

- **`layer`** (default: `'X'`): Which data layer to use from the query AnnData object
  - `'X'` uses the main expression matrix (`.X`)
  - `'layers/counts'` uses the counts layer from `.layers['counts']`
  - Any other string uses the corresponding layer from `.layers[layer_name]`

- **`var_key`** (default: `None`): Column name in `.var` to use for gene matching between query and reference model
  - If `None`, uses the `.var` index (var_names)
  - Important for ensuring gene names are consistent between query and reference model

#### Model Parameters (`model_params`)
Parameters that align the query data structure with the reference model's expectations:

- **`batch_key`**: Column name in `.obs` containing batch/sample information (required)
- **`labels_key`**: Column name in `.obs` containing cell type labels (optional, can be `None` for unlabeled data)
- **`categorical_covariate`**: List of categorical covariate column names in `.obs` (optional, e.g., `["donor", "condition"]`)
- **`continuous_covariate`**: List of continuous covariate column names in `.obs` (optional, e.g., `["age", "BMI"]`)

#### Training Parameters (`train_params`)
Parameters that control the reference mapping training process (all optional):

- **`max_epochs`**: Maximum number of training epochs (default: 100)
- **`early_stopping`**: Whether to use early stopping to prevent overfitting (default: true)
- **`check_val_every_n_epoch`**: How often to run validation during training (default: 5)

> **Note:** The query data must have genes that overlap with the reference model's training data. Gene matching is performed using the `var_key` parameter or `.var` index. Only overlapping genes are used for mapping and missing genes from the query are padded with zeros.

## Output

The reference mapping workflow produces the following outputs:

* `<out_dir>/reference_mapping/dataset~<dataset>/file_id~<file_id>.zarr`: Mapped AnnData object containing:
  - **Latent embedding** (`obsm['X_emb']`): Low-dimensional representations in the reference space

* `<out_dir>/reference_mapping/model/dataset~<dataset>/file_id~<file_id>/`: Updated model directory with query data integrated, ready for further reference mapping tasks
