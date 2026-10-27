import warnings
import numpy as np
from scipy import sparse

from utils.misc import ensure_sparse
from utils.annotate import add_wildcards


PCA_PARAMS = [
    'n_comps',
    'layer',
    'zero_center',
    'svd_solver',
    'return_info',
    'mask_var',
    'use_highly_variable',
    'dtype',
    'chunked',
    'chunk_size',
]

SCVI_MODEL_PARAMS = [
    'n_hidden',
    'n_latent',
    'n_layers',
    'dropout_rate',
    'dispersion',
    'gene_likelihood',
    'latent_distribution',
    # VAE parameters
    'n_batch',
    'n_labels',
    'n_hidden',
    'n_latent',
    'n_layers',
    'n_continuous_cov',
    'n_cats_per_cov',
    'dropout_rate',
    'dispersion',
    'log_variational',
    'gene_likelihood',
    'latent_distribution',
    'encode_covariates',
    'deeply_inject_covariates',
    'use_batch_norm',
    'use_layer_norm',
    'use_size_factor_key',
    'use_observed_lib_size',
    'library_log_means',
    'library_log_vars',
    'var_activation',
    'extra_encoder_kwargs',
    'extra_decoder_kwargs',
    'categorical_covariate_keys',
    'continuous_covariate_keys',
    'batch_representation',
    'batch_representation_kwargs',
    # vitkl parameters
    'library_n_hidden',
    'library_log_vars_weight',
    'scale_activation',
    'use_additive_background',
    'use_batch_in_decoder',
    'regularise_dispersion',
]

SCANVI_MODEL_PARAMS = SCVI_MODEL_PARAMS + ['linear_classifier']

SCVIX_MODEL_PARAMS = SCVI_MODEL_PARAMS + [
    # scvix assay
    'assay_key', # Mandatory 
    # VampPrior
    'prior',
    'pseudoinputs_data_indices',
    'n_prior_components',
    'adversarial_group_key',
]

SYSVI_MODEL_PARAMS = SCVI_MODEL_PARAMS + [
    'system_key',
    'prior',
    'n_prior_components',
    'pseudoinputs_data_indices',
]

DRVI_MODEL_PARAMS = [
    'categorical_covariate_keys',
    'continuous_covariate_keys',
    'n_latent',
    'encoder_dims',
    'decoder_dims',
    'prior',
    'prior_init_obs',
    'categorical_covariates',
    # DRVIModule VAE parameters
    'n_input',
    'n_split_latent',
    'split_aggregation',
    'split_method',
    'decoder_reuse_weights',
    'n_continuous_cov',
    'n_cats_per_cov',
    'encode_covariates',
    'deeply_inject_covariates',
    'categorical_covariate_dims',
    'covariate_modeling_strategy',
    'use_batch_norm',
    'affine_batch_norm',
    'use_layer_norm',
    'input_dropout_rate',
    'encoder_dropout_rate',
    'decoder_dropout_rate',
    'gene_likelihood',
    'prior',
    # 'prior_init_dataloader',
    'var_activation',
    # 'encoder_layer_factory',
    # 'decoder_layer_factory',
    'extra_encoder_kwargs',
    'extra_decoder_kwargs',
]


def add_metadata(adata, wildcards, params, **kwargs):
    """
    Add integration metatdata to integratd output
    :param adata:
    :param wildcards:
    :param params:
    :return:
    """
    adata.uns['integration'] = {
        'method': wildcards.method,
        'label_key': wildcards.label,
        'batch_key': wildcards.batch,
        'output_type': params['output_type'],
        'hyperparams': params['hyperparams'],
        **kwargs
    }
    add_wildcards(adata, wildcards, 'integration')
    # hyperparams = params['hyperparams']
    # if hyperparams is None:
    #     hyperparams = {}
    # add_wildcards(adata, dict(wildcards) | hyperparams, 'int')


def remove_slots(adata, output_type, keep_X=False):
    """
    Remove slots that are redundant to integration output
    """
    if isinstance(output_type, str):
        output_type = [output_type]
    
    del adata.layers
    if 'X_pca' in adata.obsm:
       del adata.obsm['X_pca']
    
    if 'full' in output_type:
        ensure_sparse(adata)
    elif 'embed' in output_type:
        if keep_X:
            assert adata.X is not None
            ensure_sparse(adata)
        else:
            del adata.X
    elif 'knn' in output_type:
        # del adata.X
        pass
    else:
        raise ValueError(f'Invalid output type {output_type}')
    return adata


def set_model_history_dtypes(model_history, dtype='float32'):
    """
    Quickfix to change encoding of the model history for saving in zarr file
    
    :param model_history: model.history from pytorch model (dictionary of pandas DataFrames)
    :param dtype: dtype to convert to
    """
    return {
        key: value.astype(dtype)
        for key, value in model_history.items()
    }


def standardize_training_logs(model):
    """
    Return standardized logs in format:

    {
        metric_name: {
            "train": [...],
            "validation": [...],
            "other": [...]
        }
    }
    """
    source = getattr(model, "history", None) or getattr(getattr(model, "trainer", None), "logs", None)

    if not source:
        warnings.warn("Model contains neither history nor logs.")
        return {}

    standardized = {}

    def parse_key(key):
        if key.endswith("_train"):
            return key[:-6], "train"
        if key.endswith("_validation"):
            return key[:-11], "validation"
        if key.startswith("train_"):
            return key[6:], "train"
        if key.startswith("val_"):
            return key[4:], "validation"
        if key in ("loss", "epoch_loss"):
            return "loss", "train"
        if key in ("val_loss", "validation_loss"):
            return "loss", "validation"
        return key, "train"

    for key, values in source.items():
        metric, split = parse_key(key)

        values = values.values if hasattr(values, "values") else values
        existing = standardized.get(metric, {}).get(split, None)
        if existing is not None:
            if existing is values or np.array_equal(existing, values):
                continue
            warnings.warn(f"Duplicate metric '{metric}' for split '{split}'. Keeping first.")
            continue

        standardized.setdefault(metric, {})[split] = values

    return standardized


def plot_model_history(
    model,
    output_dir,
    model_name="Model",
    prefix='',
    figsize_per_subplot=(4, 3),
    lw=1.5,
):
    """
    Plot standardized model training metrics grouped by metric type.
    Each metric gets its own subplot and its own ylim (automatic scaling).
    One figure per group.
    """
    import matplotlib.pyplot as plt
    import math
    from collections import defaultdict

    def classify_metric(metric_name):
        name = metric_name.lower()
        if any(x in name for x in ["accuracy", "f1"]):
            return "accuracy"
        if "calibration" in name:
            return "calibration"
        if any(x in name for x in ["loss", "elbo", "kl", "mse"]):
            return "loss"
        return "other"

    def best_grid(n_metrics, max_cols=6):
        """Find most square grid layout for n_metrics plots."""
        best_rows, best_cols = 1, 1
        best_score = float('inf')
        
        for cols in range(1, min(max_cols, n_metrics) + 1):
            rows = math.ceil(n_metrics / cols)
            score = abs(rows - cols) * 10 + (rows * cols - n_metrics)
            if score < best_score:
                best_score = score
                best_rows, best_cols = rows, cols
        
        return best_rows, best_cols

    def fig_size(nrows, ncols, min_figsize=(5, 5)):
        width = max(min_figsize[0], ncols * figsize_per_subplot[0])
        height = max(min_figsize[1], nrows * figsize_per_subplot[1])
        return (width, height)

    logs = standardize_training_logs(model)
    metrics = sorted(m for m in logs if logs[m])

    if not metrics:
        print("No metrics found to plot.", flush=True)
        return

    grouped = defaultdict(list)
    for metric in metrics:
        grouped[classify_metric(metric)].append(metric)
    print(grouped, flush=True)

    for group_name, group_metrics in grouped.items():
        n_metrics = len(group_metrics)
        nrows, ncols = best_grid(n_metrics)
        fig, axes = plt.subplots(nrows, ncols, figsize=fig_size(nrows, ncols), squeeze=False)
        axes_flat = axes.ravel()

        for ax, metric in zip(axes_flat, group_metrics):
            for split in logs[metric]:
                ax.plot(logs[metric][split], label=split, linewidth=lw)

            title = metric.replace("_", " ").title()
            ax.set_title(title)
            ax.set_xlabel("Epoch")
            ax.set_ylabel(title)
            ax.autoscale(enable=True, axis="y")
            ax.grid(False)
            ax.legend()

        for ax in axes_flat[n_metrics:]:
            fig.delaxes(ax)

        fig.suptitle(f"{model_name} — {group_name.capitalize()} Metrics", fontsize=16)
        fig.tight_layout(rect=[0, 0, 1, 0.96])

        output_path = f"{output_dir}/{prefix}training_metrics_{group_name}.png"
        plt.savefig(output_path)
        plt.close()
        print(f"Saved: {output_path}")


def check_output(adata, output_type):
    """
    Process data based on output type.
    If more than one output type is given, use the most processed output type: knn > embed > full
    :param adata: integrated anndata object
    :param output_type: string or list of output type
    :return: integrated anndata object with unintegrated anndata in .raw
    """
    if isinstance(output_type, str):
        output_type = [output_type]

    if 'full' in output_type:
        assert isinstance(adata.X, (np.ndarray, sparse.csr_matrix, sparse.csc_matrix))
    elif 'embed' in output_type:
        assert 'X_emb' in adata.obsm
    elif 'knn' in output_type:
        assert 'neighbors' in adata.uns
        assert 'connectivities' in adata.obsp
        assert 'distances' in adata.obsp
    else:
        raise ValueError(f'Invalid output type {output_type}')


def get_hyperparams(
    hyperparams: dict,
    model_params: list = None,
    train_params: list = None
):
    """
    Get hyperparameters and training parameters from hyperparameter dictionary
    :param hyperparams: dictionary of hyperparameters
    :param train_params: list of training parameters
    :return: hyperparams, train_params
    """
    if hyperparams is None:
        return {}, {}
    if model_params is None:
        model_params = []
        if train_params is not None:
            model_params = [x for x in hyperparams if x not in train_params]
    model_params = {k: v for k, v in hyperparams.items() if k in model_params}
    train_params = {k: v for k, v in hyperparams.items() if k not in model_params}
    return model_params, train_params


def clean_categorical_column(adata, column):
    assert column in adata.obs, f'Colmn {column} is missing'
    adata.obs[column] = adata.obs[column].astype(str).astype('category')
