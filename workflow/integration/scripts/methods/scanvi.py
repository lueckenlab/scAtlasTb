import warnings
warnings.filterwarnings("ignore", message="Can't initialize NVML")

import torch
import scvi
from pprint import pformat
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)

import dask
dask.config.set(scheduler='threads')

from integration_utils import add_metadata, get_hyperparams, remove_slots, set_model_history_dtypes, \
    SCVI_MODEL_PARAMS, plot_model_history, clean_categorical_column
from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg

input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_model = Path(snakemake.output.model)
output_plot_dir = Path(snakemake.output.plots)
output_model.mkdir(parents=True, exist_ok=True)
output_plot_dir.mkdir(parents=True, exist_ok=True)

wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch
label_key = wildcards.label

torch.set_float32_matmul_precision('medium')
scvi.settings.seed = params.get('seed', 0)
scvi.settings.progress_bar_style = 'tqdm'
scvi.settings.num_threads = snakemake.threads
dask.config.set({'num_workers': snakemake.threads})

model_params, train_params = get_hyperparams(
    hyperparams=params.get('hyperparams', {}),
    model_params=SCVI_MODEL_PARAMS,
)
categorical_covariate_keys = model_params.pop('categorical_covariate_keys', [])
if isinstance(categorical_covariate_keys, str):
    categorical_covariate_keys = [categorical_covariate_keys]
continuous_covariate_keys = model_params.pop('continuous_covariate_keys', [])
if isinstance(continuous_covariate_keys, str):
    continuous_covariate_keys = [continuous_covariate_keys]

train_scvi = {k: v for k, v in train_params.items() if not k.startswith('scanvi_')}
train_scanvi = {k.replace('scanvi_', ''): v for k, v in train_params.items() if k.startswith('scanvi_')}

train_scvi.setdefault('check_val_every_n_epoch', 1)
train_scanvi.setdefault('check_val_every_n_epoch', 1)

logging.info(
    f'model parameters:\n{pformat(model_params)}'
    f'\nscvi train parameters:\n{pformat(train_scvi)}'
    f'\nscanvi train parameters:\n{pformat(train_scanvi)}'
)

# check GPU
print(f'GPU available: {torch.cuda.is_available()}', flush=True)

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='layers/counts',
    var='var',
    obs='obs',
    uns='uns',
    dask=True,
    backed=True,
)

clean_categorical_column(adata, batch_key)
if label_key is not None and label_key != 'None':
    clean_categorical_column(adata, label_key)

# subset features
adata, _ = subset_hvg(adata, var_column='integration_features')

if isinstance(categorical_covariate_keys, list):
    for cov in categorical_covariate_keys:
        adata.obs[cov] = adata.obs[cov].astype('str')

scvi.model.SCVI.setup_anndata(
    adata,
    batch_key=batch_key,
    categorical_covariate_keys=categorical_covariate_keys, # Does not work with scArches
    continuous_covariate_keys=continuous_covariate_keys,
)

logging.info(f'Set up scVI with parameters:\n{pformat(model_params)}')
model = scvi.model.SCVI(
    adata,
    **model_params
)

logging.info(f'Train scVI with parameters:\n{pformat(train_scvi)}')
model.train(**train_scvi)

plot_model_history(
    model=model,
    output_dir=output_plot_dir,
    model_name="scVI",
    prefix='scvi_',
)

logging.info(f'Set up scANVI on top of scVI with parameters:\n{pformat(model_params)}')
model = scvi.model.SCANVI.from_scvi_model(
    model,
    labels_key=label_key,
    unlabeled_category='nan',
)

logging.info(f'Train scANVI with parameters:\n{pformat(train_scanvi)}')
# batch_size = train_scanvi.get("batch_size", 128)
# steps_per_epoch = ceil(adata.n_obs / batch_size)
checkpoint_path = output_model / f"scanvi_best.ckpt"
callbacks = [
    scvi.train.SaveCheckpoint(
        dirpath=output_model,
        filename=checkpoint_path.name,
        monitor="validation_accuracy",
        load_best_on_end=True,
        check_nan_gradients=True,
    )
]

try:
    model.train(
        callbacks=callbacks,
        **train_scanvi,
    )
except Exception as train_err:
    logging.exception("Training failed. Attempting checkpoint recovery...")
    
    if checkpoint_path.exists():
        try:
            model = scvi.model.SCANVI.load_from_checkpoint(
                checkpoint_path, 
                adata=adata,
            )
            logging.info(f"Successfully loaded SCANVI checkpoint from {checkpoint_path} after training failure.")
        except Exception:
            logging.exception("Checkpoint loading failed.")
            raise RuntimeError(
                "Training failed and checkpoint recovery also failed."
            ) from train_err
    else:
        logging.error(f"No checkpoint found at {checkpoint_path}")
        raise train_err

plot_model_history(
    model=model,
    output_dir=output_plot_dir,
    model_name="scANVI",
    prefix='scanvi_',
)

logging.info('Save model...')
model.save(output_model, overwrite=True)

# prepare output adata
adata.obsm["X_emb"] = model.get_latent_representation()
adata = remove_slots(adata=adata, output_type=params['output_type'], keep_X=True)
add_metadata(
    adata,
    wildcards,
    params,
    model_history=set_model_history_dtypes(model.history)
)

logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns'],
)
