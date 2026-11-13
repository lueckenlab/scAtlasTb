import logging
logging.basicConfig(level=logging.INFO)

import torch
from scarches.models.scpoli import scPoli
from pprint import pformat
from pathlib import Path

import dask
dask.config.set(scheduler='threads')

from integration_utils import add_metadata, get_hyperparams, remove_slots, \
    set_model_history_dtypes, plot_model_history, clean_categorical_column
from scArches_utils import SCPOLI_MODEL_PARAMS, SCPOLI_EARLY_STOPPING
from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg


input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_model = snakemake.output.model
output_plot_dir = snakemake.output.plots
Path(output_plot_dir).mkdir(parents=True, exist_ok=True)

wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch
label_key = wildcards.label
if label_key == 'None':
    label_key = None

torch.manual_seed(params.get('seed', 0))
torch.set_num_threads(snakemake.threads)
dask.config.set({'num_workers': snakemake.threads})

model_params, train_params = get_hyperparams(
    hyperparams=params.get('hyperparams', {}),
    model_params=SCPOLI_MODEL_PARAMS,
)

# set default model parameters
model_params = dict(
    condition_keys=[batch_key],
    cell_type_keys=[label_key],
    unknown_ct_names=['NA'],
) | model_params

# set defaults for training parameters
train_params.setdefault('check_val_every_n_epoch', 1) # needed to be able to plot loss curve
train_params.setdefault('early_stopping_kwargs', SCPOLI_EARLY_STOPPING)
train_params.setdefault('pretrain_epochs', int(0.8 * train_params.get('n_epochs', 0)))

logging.info(
    f'model parameters:\n{pformat(model_params)}\n'
    f'training parameters:\n{pformat(train_params)}',
)


# check GPU
logging.info(f'GPU available: {torch.cuda.is_available()}')

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='layers/counts',
    obs='obs',
    var='var',
    uns='uns',
    dask=True,
    backed=True,
)

clean_categorical_column(adata, batch_key)

cell_type_keys = model_params.get('cell_type_keys', [])
if isinstance(cell_type_keys, str):
    cell_type_keys = [cell_type_keys]
elif cell_type_keys is None:
    cell_type_keys = []
for key in cell_type_keys:
    if key is None:
        continue
    clean_categorical_column(adata, key)

# subset features
adata, _ = subset_hvg(adata, var_column='integration_features')

# prepare data for model
adata.X = adata.X.astype('float32')

logging.info(f'Set up scPoli with parameters:\n{pformat(model_params)}')
model = scPoli(adata=adata, **model_params)

logging.info(f'Train scPoli with parameters:\n{pformat(train_params)}')
model.train(**train_params)

logging.info('Plot model history...')
plot_model_history(
    model=model,
    output_dir=output_plot_dir,
    model_name="scPoli",
)

logging.info('Save model...')
model.save(output_model, overwrite=True)

# prepare output adata
adata.obsm["X_emb"] = model.get_latent(adata, mean=True)
adata.uns[f"scpoli_{batch_key}_embeddings"] = model.get_conditional_embeddings().X
adata = remove_slots(adata=adata, output_type=params['output_type'], keep_X=True)
add_metadata(
    adata,
    wildcards,
    params,
    model_history=dict(model.trainer.logs)
)

logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns']
)
