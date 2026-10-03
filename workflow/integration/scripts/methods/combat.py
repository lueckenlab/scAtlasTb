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
    # keep obsm so that the input's X_pca, which remove_slots dropped, is not linked back in
    files_to_keep=['X', 'obsm', 'var', 'uns'],
)