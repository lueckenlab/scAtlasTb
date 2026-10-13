import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import anndata as ad
import torch
from scimilarity.utils import lognorm_counts, align_dataset
from scimilarity import CellQuery

from utils.io import read_anndata, write_zarr_linked
from utils.misc import dask_compute


input_file = snakemake.input.zarr
model_path = snakemake.input.model
output_file = snakemake.output.zarr

layer = snakemake.params.get('layer', 'X')
var_key = snakemake.params.get('var_key')
model_params = snakemake.params.get('model_params', {})
gene_overlap_threshold = model_params.pop('gene_overlap_threshold', 5000)

logging.info(f'Read {input_file}...')
adata = read_anndata(
    input_file,
    X=layer,
    var='var',
    dask=True,
    backed=True
)
adata.var_names = adata.var_names if var_key is None else adata.var[var_key].astype(str).values

logging.info(f'Loading model: {model_path}...')
cq = CellQuery(model_path)

logging.info('Align data')
model_genes = cq.gene_order

# Handle duplicate var_names
if adata.var_names.duplicated().any():
    logging.warning(f'Found {adata.var_names.duplicated().sum()} duplicate var_names. Making unique...')
    adata.var_names_make_unique()

adata = adata[:, adata.var_names.isin(model_genes)].copy()
assert adata.n_vars > 0, 'No overlapping genes.'
logging.info(f'Found {adata.n_vars} overlapping genes.')

adata.layers['counts'] = adata.X.copy()
del adata.X
adata = dask_compute(adata, layers='counts')

adata = align_dataset(
    adata,
    target_gene_order=model_genes,
    gene_overlap_threshold=gene_overlap_threshold,
)

logging.info('Computing SCimilarity embeddings...')
adata = lognorm_counts(adata)
X_emb = cq.get_embeddings(adata.X)

if input_file.endswith('.h5ad'):
    scimilarity_shape = adata.shape[0]
    adata = read_anndata(input_file, dask=True, backed=True)
    assert adata.shape[0] == scimilarity_shape, \
        'Shape mismatch after reloading original data' \
        f'(adata shape: {adata.shape[0]}, expected: {scimilarity_shape})'

adata.obsm['X_emb'] = X_emb
logging.info(adata.__str__())

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata=adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['obsm']
)