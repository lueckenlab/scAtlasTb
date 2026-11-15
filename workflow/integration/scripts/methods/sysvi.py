import warnings
warnings.filterwarnings("ignore", message="Can't initialize NVML")

import torch
import scvi
from scvi.external import SysVI
from pathlib import Path
from pprint import pformat
from math import ceil
from pytorch_lightning.callbacks import TQDMProgressBar
import logging
logging.basicConfig(level=logging.INFO)

from integration_utils import add_metadata, get_hyperparams, remove_slots, set_model_history_dtypes, \
    SYSVI_MODEL_PARAMS, plot_model_history, clean_categorical_column
from utils.io import read_anndata, write_zarr_linked, to_memory
from utils.accessors import subset_hvg

input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_model = snakemake.output.model
output_plot_dir = snakemake.output.plots
Path(output_plot_dir).mkdir(parents=True, exist_ok=True)

wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch

torch.set_float32_matmul_precision('medium')
scvi.settings.seed = params.get('seed', 0)
scvi.settings.progress_bar_style = 'tqdm'
scvi.settings.num_threads = snakemake.threads

model_params, train_params = get_hyperparams(
    hyperparams=params.get('hyperparams', {}),
    model_params=SYSVI_MODEL_PARAMS,
)

logging.info(
    'User-provided parameters:\n'
    f'model parameters:\n{pformat(model_params)}\n'
    f'training parameters:\n{pformat(train_params)}'
)


# set correct early_stopping parameters
if train_params.pop('early_stopping', False):
    train_params |= dict(
        early_stopping=True,
        log_every_n_steps=train_params.get('log_every_n_steps', 1),
        check_val_every_n_epoch=train_params.get('check_val_every_n_epoch', 1),
        val_check_interval=train_params.get('val_check_interval', 1.0),
    )

# parse categorical and continuous covariate keys
categorical_covariate_keys = model_params.pop('categorical_covariate_keys', [])
if isinstance(categorical_covariate_keys, str):
    categorical_covariate_keys = [categorical_covariate_keys]
continuous_covariate_keys = model_params.pop('continuous_covariate_keys', [])
if isinstance(continuous_covariate_keys, str):
    continuous_covariate_keys = [continuous_covariate_keys]

# set system key
try:
    system_key = model_params.pop('system_key')
except KeyError:
    raise KeyError(
        'system_key is not configured, but is mandatory.'
        'Please make sure to set it in your config file.'
    )

# ensure that batch_key is included in categorical_covariate_keys
if batch_key not in categorical_covariate_keys:
    categorical_covariate_keys.append(batch_key)

# set defaults for training parameters
train_params.setdefault('check_val_every_n_epoch', 1) # needed to be able to plot loss curves

logging.info(
    'Parsed parameters:\n'
    f'model parameters:\n{pformat(model_params)}\n'
    f'training parameters:\n{pformat(train_params)}'
)

# check GPU
print(f'GPU available: {torch.cuda.is_available()}', flush=True)

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='layers/normcounts',
    var='var',
    obs='obs',
    uns='uns',
    dask=True,
    backed=True,
)

clean_categorical_column(adata, batch_key)
clean_categorical_column(adata, system_key)

# subset features
adata, _ = subset_hvg(adata, var_column='integration_features')

for cov in categorical_covariate_keys:
    adata.obs[cov] = adata.obs[cov].astype('str')

SysVI.setup_anndata(
    adata,
    batch_key=system_key,
    categorical_covariate_keys=categorical_covariate_keys, # Does not work with scArches
    continuous_covariate_keys=continuous_covariate_keys,
)

logging.info(f'Set up SysVI with parameters:\n{pformat(model_params)}')
model = SysVI(adata, **model_params)

logging.info(f'Train SysVI with parameters:\n{pformat(train_params)}')
batch_size = train_params.get("batch_size", 128)
steps_per_epoch = ceil(adata.n_obs / batch_size)
model.train(
    **train_params,
    simple_progress_bar=False,
    callbacks=[TQDMProgressBar(refresh_rate=steps_per_epoch)],
)

logging.info('Save model...')
model.save(output_model, overwrite=True)

# prepare output adata
adata.obsm["X_emb"] = model.get_latent_representation(adata=adata)
adata = remove_slots(adata=adata, output_type=params['output_type'], keep_X=True)
add_metadata(
    adata,
    wildcards,
    params,
    # history is not saved with standard model saving
    model_history=set_model_history_dtypes(model.history)
)

plot_model_history(
    model=model,
    output_dir=output_plot_dir,
    model_name="SysVI",
)

logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns'],
)
