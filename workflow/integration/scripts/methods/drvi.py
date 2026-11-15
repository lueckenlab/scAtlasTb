import torch
import drvi
import scvi
from pathlib import Path
from pprint import pformat
from matplotlib import pyplot as plt
import anndata as ad
import logging
logging.basicConfig(level=logging.INFO)
import dask

from integration_utils import add_metadata, get_hyperparams, remove_slots, set_model_history_dtypes, \
    DRVI_MODEL_PARAMS, plot_model_history, clean_categorical_column
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
# drvi has no settings object; it trains through scvi-tools, whose global seed
# seeds torch, numpy and python
scvi.settings.seed = params.get('seed', 0)
# drvi.settings.progress_bar_style = 'tqdm'
# drvi.settings.num_threads = snakemake.threads
dask.config.set(scheduler='threads', num_workers=snakemake.threads)

model_params, train_params = get_hyperparams(
    hyperparams=params.get('hyperparams', {}),
    model_params=DRVI_MODEL_PARAMS,
)

# parse categorical and continuous covariate keys
categorical_covariate_keys = model_params.pop('categorical_covariate_keys', [])
if isinstance(categorical_covariate_keys, str):
    categorical_covariate_keys = [categorical_covariate_keys]
continuous_covariate_keys = model_params.pop('continuous_covariate_keys', [])
if isinstance(continuous_covariate_keys, str):
    continuous_covariate_keys = [continuous_covariate_keys]

# set defaults for training parameters
default_plan_kwargs = {
    "n_epochs_kl_warmup": train_params.get('max_epochs', 400),
}
train_params.setdefault('plan_kwargs', default_plan_kwargs)
train_params.setdefault('check_val_every_n_epoch', 1)

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
else:
    categorical_covariate_keys = []


drvi.model.DRVI.setup_anndata(
    adata,
    batch_key=batch_key,
    categorical_covariate_keys=categorical_covariate_keys,
    continuous_covariate_keys=continuous_covariate_keys,
    is_count_data=True,
)

logging.info(f'Set up DRVI with parameters:\n{pformat(model_params)}')
model = drvi.model.DRVI(adata, **model_params)

logging.info(f'Train DRVI with parameters:\n{pformat(train_params)}')
model.train(**train_params)

logging.info('Plot model history...')
plot_model_history(
    model=model,
    output_dir=output_plot_dir,
    model_name="DRVI",
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

# plot DRVI latent space
embed = ad.AnnData(X=adata.obsm["X_emb"], obs=adata.obs)
drvi.utils.tl.set_latent_dimension_stats(model, embed)

drvi.utils.pl.plot_latent_dimension_stats(embed, ncols=2)
fig = plt.gcf()
fig.savefig(Path(output_plot_dir) / 'latent_dimension_stats.png')
plt.close(fig)

drvi.utils.pl.plot_latent_dimension_stats(embed, ncols=2, remove_vanished=True)
fig = plt.gcf()
fig.savefig(Path(output_plot_dir) / 'latent_dimension_stats_no_vanished.png')
plt.close(fig)

logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns'],
)
