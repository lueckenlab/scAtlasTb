from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc
import logging
logging.basicConfig(level=logging.INFO)

from integration_utils import add_metadata, remove_slots, clean_categorical_column
from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg

input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

hyperparams = params.get('hyperparams', {})
hyperparams = {} if hyperparams is None else hyperparams

if 'covariates' in hyperparams:
    covariates = hyperparams['covariates']
    if isinstance(covariates, str):
        hyperparams['covariates'] = [covariates]

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

clean_categorical_column(adata, wildcards.batch)

# subset features
all_var_names = adata.var_names.copy()
adata, _ = subset_hvg(adata, var_column='integration_features')

# run method
logging.info(f'Run Combat with parameters {pformat(hyperparams)}...')
sc.pp.combat(
    adata,
    key=wildcards.batch,
    **hyperparams
)

# prepare output adata
adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['X', 'var', 'uns'],
    # subset linked feature slots (e.g. layers) to the integration features
    subset_mask=(None, all_var_names.isin(adata.var_names)),
)