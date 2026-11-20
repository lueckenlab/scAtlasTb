# Cell Type Prediction

This module enables automated cell type annotation of single-cell RNA-seq datasets using pre-trained models. This allows for:

- Consistent cell type annotation across studies
- Probabilistic cell type assignments with confidence scores
- Majority voting and over-clustering analysis for robust predictions

The module currently supports:

- **CellTypist**: Automated cell type annotation using pre-trained models

## Environments

The following environments are needed for cell type prediction:

- [`celltypist`](https://github.com/HCA-integration/scAtlasTb/blob/main/envs/celltypist.yaml)

## Configuration

```yaml
DATASETS:
  test:
    input:
      celltype_prediction:
        preprocessed: test/input/preprocessing/dataset~all/file_id~pbmc/preprocessed.zarr
    celltype_prediction:
      reference_label: bulk_labels
      counts: layers/counts
      is_normalized: false
      celltypist:
        params:
          majority_voting: true
          over_clustering: bulk_labels
        models:
          - Healthy_COVID19_PBMC
          - Immune_All_Low
```

### Input

The input AnnData object should contain the single-cell RNA-seq data to be annotated.
Reference cell type labels can optionally be provided for evaluation and visualization purposes.

### Configuration Parameters

- **`counts`** (default: `'X'`): Which data layer to use from the AnnData object

- **`is_normalized`** (default: `false` for CellTypist): Boolean flag indicating whether the input data is already normalized
  - CellTypist expects log-normalized data; if `false`, the data are normalized to 10,000 counts per cell and log1p-transformed before prediction

- **`reference_label`** (optional): Column name in `.obs` containing reference cell type labels for comparison and visualization

#### CellTypist Parameters (`celltypist`)
Configuration for CellTypist cell type prediction:

- **`models`**: List of pre-trained model names to use for prediction (required)
  - Available models include `Healthy_COVID19_PBMC`, `COVID19_HumanChallenge_Blood`, `Immune_All_Low`, etc.
  - Multiple models can be applied to the same dataset

- **`params`**: Model parameters (optional)
  - **`majority_voting`**: Enable majority voting across over-clustering results (default: `false`)
  - **`over_clustering`**: Column name in `.obs` for over-clustering analysis (optional)

#### Sex prediction (`predict_sex`, optional)
If set, donor sex is predicted from X- and Y-linked gene expression:

- **`donor_key`** (required): `.obs` column with donor IDs
- **`donors`**: optional list of donors to restrict prediction to
- **`predict_column`** (default: `sex`): name of the output column
- **`reference_key`** (default: value of `predict_column`): `.obs` column with known sex for accuracy reporting
- **`x_genes`**, **`y_genes`**, **`x_threshold`** (default `0`), **`y_threshold`** (default `4`), **`imbalance_frac`** (default `0.1`): parameters of the X/Y expression rule
- **`y_nonpar_genes`**, **`y_par_genes`**, **`chrY_threshold`** (default `0.5`): parameters of the chrY non-PAR/PAR ratio rule; gene lists may be given as gene names, text files or URLs

> **Note:** CellTypist models are trained on specific tissue types and cell populations. Choose models appropriate for your data type (e.g., PBMC, immune cells, etc.).

## Output

### CellTypist
The cell type prediction workflow produces the following outputs:

* `<out_dir>/celltype_prediction/dataset~<dataset>/file_id~<file_id>.zarr`: Annotated AnnData object containing:
  - **Direct predictions** (`obs['celltypist_<model>:predicted_labels']`): Primary cell type predictions
  - **Majority voting results** (`obs['celltypist_<model>:majority_voting']`): Consensus predictions (if enabled)
  - **Over-clustering results** (`obs['celltypist_<model>:over_clustering']`): Fine-grained clustering results (if specified)
  - **Confidence scores** (`obs['celltypist_<model>:conf_score']`): Prediction confidence values

* `<images>/celltype_prediction/dataset~<dataset>/file_id~<file_id>/celltypist--<model>/`: Dot plots comparing predictions with reference labels (if `reference_label` is provided)
* If `predict_sex` is configured, the annotated AnnData object additionally contains `obs['<predict_column>']`, `obs['<predict_column>_chrY']` and `uns['predict_sex']`
