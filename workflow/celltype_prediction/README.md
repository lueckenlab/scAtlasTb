# Cell Type Prediction

This module adds cell type predictions from pre-trained
[CellTypist](https://www.celltypist.org/) models to a single-cell dataset. Use it
to annotate cells or compare predictions from several models. Predictions
include confidence scores; if you have existing labels, you can also generate
plots to compare them with the predictions.

## Quick start

Configure the input dataset and at least one CellTypist model:

```yaml
DATASETS:
  my_dataset:
    input:
      celltype_prediction:
        preprocessed: path/to/preprocessed.zarr
    celltype_prediction:
      celltypist:
        models:
          - Immune_All_Low
```

`preprocessed` is an example input name: point it to the AnnData/Zarr data you
want to annotate (or to the output of an upstream module). Choose a model
appropriate for your tissue and cell populations; browse available models at
[celltypist.org](https://www.celltypist.org/models). The `models` list is
required, and each listed model is run separately.

## Configuration

| Setting | Required? | Description |
|---|---|---|
| `input.celltype_prediction` | **Yes** | Input AnnData/Zarr object to annotate. |
| `celltypist.models` | **Yes** | One or more CellTypist model names. Models are downloaded automatically; see [Model download and environment](#model-download-and-environment). |
| `counts` | No | AnnData layer to use as expression data, such as `layers/counts`. Defaults to `X`. |
| `is_normalized` | No | Whether the selected data is already log1p-normalized to 10,000 counts per cell. Defaults to `false`. |
| `reference_label` | No | Name of an existing `.obs` column to compare against predictions. The column must exist if specified. |
| `celltypist.params` | No | Options passed through to `celltypist.annotate`. |
| `predict_sex` | No | Optionally run sex prediction; its results are included in the collected output. |

For example, this configuration uses raw counts, compares predictions with
existing labels, and enables majority voting:

```yaml
DATASETS:
  my_dataset:
    input:
      celltype_prediction:
        preprocessed: path/to/preprocessed.zarr
    celltype_prediction:
      counts: layers/counts
      is_normalized: false
      reference_label: bulk_labels
      celltypist:
        models:
          - Healthy_COVID19_PBMC
          - Immune_All_Low
        params:
          majority_voting: true
          over_clustering: leiden
```

### Expression data and normalization

CellTypist expects log1p-normalized expression, with each cell normalized to
10,000 total counts. By default, `is_normalized: false`: the module reads raw
counts from `X` (or the layer selected by `counts`), normalizes each cell to
10,000 counts, and applies `log1p` before prediction. Set
`is_normalized: true` only when the selected data is already normalized this
way; the module then uses it as-is. The flag does not detect or correct an
incorrectly described input.

If `var['feature_name']` is present, the module uses it for gene names rather
than `var_names` (often Ensembl IDs), to match the symbols used by CellTypist
models.

### CellTypist options

`celltypist.params` is passed directly to
[`celltypist.annotate`](https://celltypist.readthedocs.io/en/latest/celltypist.annotate.html).
In addition to the common options below, other supported annotate options can
be provided in the same mapping.

- **`majority_voting`**: Set to `true` to add labels smoothed by the majority
  prediction within clusters. This can make predictions more consistent within
  a cluster.
- **`over_clustering`**: Optionally name an `.obs` column containing cluster
  assignments to use for majority voting. If majority voting is enabled and
  this option is omitted, CellTypist can infer clusters internally from the
  input neighbor graph; the module loads `obsm`, `obsp`, and `uns` for this
  case.

`reference_label` is independent of voting: when set, it must name a column in
the input `.obs`. The module writes a `predicted_labels.png` comparison plot
for every model, and also a `majority_voting.png` plot when majority voting is
enabled.

## Model download and environment

The model-download rule retrieves each requested model using CellTypist and
keeps it in the module cache at
`<out_dir>/models/celltypist/<model>.pkl`. Downloading a model requires
internet access; the cached model is used by subsequent workflow runs.
Prediction runs in the module's `celltypist` Conda environment (Python 3.11;
see [`envs/celltypist.yaml`](../../envs/celltypist.yaml)).

## Outputs

For each model, the module writes an intermediate linked Zarr object at
`<out_dir>/celltypist/<wildcard_pattern>/<model>.zarr`. It contains the
CellTypist prediction columns in `.obs`, prefixed with `celltypist_<model>:`:
`predicted_labels` and `conf_score`, plus `majority_voting` when enabled.
Comparison plots are written under
`<image_dir>/<wildcard_pattern>/celltypist--<model>/`.

The default `all` target collects the intermediate results into
`<out_dir>/<wildcard_pattern>.zarr`, linked to the input object. This final
object retains the input data and adds the `.obs` columns from all configured
CellTypist models and, when enabled, `predict_sex`.
