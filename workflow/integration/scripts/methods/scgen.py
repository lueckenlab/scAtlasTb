import logging
logging.basicConfig(level=logging.INFO)
from pprint import pformat
from pathlib import Path
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity

from integration_utils import add_metadata, get_hyperparams, remove_slots, \
    set_model_history_dtypes, plot_model_history, clean_categorical_column
from scArches_utils import SCGEN_MODEL_PARAMS
from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg

input_file = snakemake.input[0]
output_file = snakemake.output[0]
output_model = snakemake.output.model
output_plot_dir = snakemake.output.plots
Path(output_plot_dir).mkdir(parents=True, exist_ok=True)

wildcards = snakemake.wildcards
batch_key = wildcards.batch
label_key = wildcards.label
params = snakemake.params

torch.manual_seed(params.get('seed', 0))
torch.set_num_threads(snakemake.threads)

model_params, train_params = get_hyperparams(
    hyperparams=params.get('hyperparams', {}),
    model_params=SCGEN_MODEL_PARAMS,
)

# set defaults for training parameters
train_params.setdefault('check_val_every_n_epoch', 1) # needed to be able to plot loss curve

logging.info(
    f'model parameters:\n{pformat(model_params)}\n'
    f'training parameters:\n{pformat(train_params)}'
)

# check GPU
logging.info(f'GPU available: {torch.cuda.is_available()}')

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X='layers/normcounts',
    obs='obs',
    var='var',
    uns='uns',
    dask=True,
    backed=True,
)

clean_categorical_column(adata, batch_key)
if label_key is not None and label_key != 'None':
    clean_categorical_column(adata, label_key)

# subset features
all_var_names = adata.var_names.copy()
adata, _ = subset_hvg(adata, var_column='integration_features')

# prepare data for model
adata.X = adata.X.astype('float32')
adata = remove_sparsity(adata)

logging.info(f'Set up scGEN with parameters:\n{pformat(model_params)}')
model = sca.models.scgen(
    adata=adata,
    **model_params,
)

logging.info(f'Train scGEN with parameters:\n{pformat(train_params)}')
model.train(**train_params)

logging.info('Plot model history...')
plot_model_history(
    model=model,
    output_dir=output_plot_dir,
    model_name="scGEN",
)

logging.info('Save model...')
model.save(output_model, overwrite=True)

corrected_adata = model.batch_removal(
    adata,
    batch_key=batch_key,
    cell_label_key=label_key,
    return_latent=True,
)

# prepare output adata
adata.X = corrected_adata.X
adata.obsm["X_emb"] = corrected_adata.obsm["latent_corrected"]
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
    files_to_keep=['X', 'obsm', 'var', 'varm', 'varp', 'uns'],
    # subset linked feature slots (e.g. layers) to the integration features
    subset_mask=(None, all_var_names.isin(adata.var_names)),
)
