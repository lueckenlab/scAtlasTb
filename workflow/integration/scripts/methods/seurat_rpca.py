"""
Seurat v5 RPCA integration: IntegrateLayers(method = RPCAIntegration) on one layer per batch
(https://satijalab.org/seurat/articles/seurat5_integration).

The method runs in R (seurat_rpca.R, called through seurat_utils). It starts from the raw counts
(layers/counts): each cell is log-normalized on all genes (Seurat's LogNormalize) before the
counts are subset to the integration features (var_mask), which are then scaled, reduced by PCA
and integrated. Output: embedding in obsm['X_emb'].

Hyperparameters:
    npcs: number of principal components computed per run (default: 50)
    dims: dimensions used to find anchors and returned as the embedding (default: 30)
    reference: batch, or list of batches, to integrate the other batches onto (default: none,
        i.e. all pairs of batches); pass a list inside a list, since lists are expanded into
        separate runs, e.g. reference: [[batch1, batch2]]
    scale_factor: LogNormalize scale factor (default: 10000)
    any other: passed to IntegrateLayers() and RPCAIntegration() under its R name, with
        underscores replaced by dots (e.g. k_weight: 100 -> k.weight = 100)
"""
import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
from pprint import pformat

from integration_utils import add_metadata, remove_slots, clean_categorical_column
from seurat_utils import read_counts, run_r_method, r_arguments
from utils.io import write_zarr_linked


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
batch_key = wildcards.batch

hyperparams = params.get('hyperparams', {})
hyperparams = {} if hyperparams is None else dict(hyperparams)

reference = hyperparams.pop('reference', None)
if reference is not None and not isinstance(reference, list):
    reference = [reference]
r_params = {
    'seed': params.get('seed', 0),
    'npcs': hyperparams.pop('npcs', 50),
    'dims': hyperparams.pop('dims', 30),
    'scale_factor': hyperparams.pop('scale_factor', 1e4),
    'reference': None if reference is None else [str(batch) for batch in reference],
    'integrate_args': r_arguments(hyperparams),
}

adata, library_size = read_counts(input_file, var_column='integration_features')
clean_categorical_column(adata, batch_key)

# run method
logging.info(f'Run Seurat RPCA with parameters {pformat(r_params)}...')
adata.obsm['X_emb'] = run_r_method(
    script=Path(__file__).with_suffix('.R'),
    adata=adata,
    library_size=library_size,
    obs={'batch': adata.obs[batch_key]},
    params=r_params,
    threads=snakemake.threads,
    tmpdir=snakemake.resources.get('tmpdir'),
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
    files_to_keep=['obsm', 'uns'],
)
