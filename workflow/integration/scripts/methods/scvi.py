import warnings
warnings.filterwarnings("ignore", message="Can't initialize NVML")
warnings.filterwarnings(
    "ignore",
    message="The value argument must be within the support of the distribution",
    category=UserWarning,
)

import torch
import scvi
from pathlib import Path
from pprint import pformat
from math import ceil
from lightning.pytorch.callbacks import TQDMProgressBar
import logging
logging.basicConfig(level=logging.INFO)

import dask
dask.config.set(scheduler='threads')

from integration_utils import add_metadata, get_hyperparams, remove_slots, set_model_history_dtypes, \
    SCVI_MODEL_PARAMS, plot_model_history, clean_categorical_column
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
# scvi.settings.progress_bar_style = 'tqdm'
scvi.settings.num_threads = snakemake.threads
dask.config.set({'num_workers': snakemake.threads})

model_params, train_params = get_hyperparams(
    hyperparams=params.get('hyperparams', {}),
    model_params=SCVI_MODEL_PARAMS,
)

# parse categorical and continuous covariate keys
categorical_covariate_keys = model_params.pop('categorical_covariate_keys', [])
if isinstance(categorical_covariate_keys, str):
    categorical_covariate_keys = [categorical_covariate_keys]
continuous_covariate_keys = model_params.pop('continuous_covariate_keys', [])
if isinstance(continuous_covariate_keys, str):
    continuous_covariate_keys = [continuous_covariate_keys]

# set defaults for training parameters
train_params.setdefault('check_val_every_n_epoch', 1) # needed to be able to plot loss curves

logging.info(
    f'model parameters:\n{pformat(model_params)}\n'
    f'training parameters:\n{pformat(train_params)}'
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

# subset features
adata, subsetted = subset_hvg(adata, var_column='integration_features')

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
model = scvi.model.SCVI(adata, **model_params)

logging.info(f'Train scVI with parameters:\n{pformat(train_params)}')
batch_size = train_params.get("batch_size", 128)
steps_per_epoch = ceil(adata.n_obs / batch_size)
model.train(
    **train_params,
    # simple_progress_bar=True,
    # callbacks=[TQDMProgressBar(refresh_rate=steps_per_epoch)],
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
    # history is not saved with standard model saving
    model_history=set_model_history_dtypes(model.history)
)

plot_model_history(
    model=model,
    output_dir=output_plot_dir,
    model_name="scVI",
)

logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns'],
)
