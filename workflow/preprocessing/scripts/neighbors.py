"""
Build kNN graph on embedding
"""
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked
from utils.processing import assert_neighbors, sc, USE_GPU
from utils.misc import ensure_dense, ensure_sparse


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params.get('args', {})
extra_uns = snakemake.params.get('extra_uns', {})
files_to_overwrite = ['obsp', 'uns']

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, obs='obs', obsm='obsm', uns='uns', var='var')

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.write_zarr(output_file)
    exit(0)

# set defaults
if params == False:
    # use existing kNN graph
    adata.obsp = read_anndata(input_file, obs='obs', obsp='obsp').obsp
    assert_neighbors(adata, check_params=False)
    adata.uns['neighbors']['params'] = adata.uns['neighbors'].get('params', {})
    # keep the representation the graph was built on (e.g. X_emb, X_pca), so that
    # downstream steps like UMAP don't fall back to the full feature matrix
    use_rep = adata.uns['neighbors']['params'].get('use_rep')
    if use_rep is None:
        use_rep = 'X_pca' if 'X_pca' in adata.obsm else 'X'
    adata.uns['neighbors']['params'] |= dict(
        connectivities_key='connectivities',
        distances_key='distances',
        use_rep=use_rep,
    )
    assert_neighbors(adata)
else:
    # set representation for neighbors
    if 'use_rep' not in params: # determine use_rep if missing
        if 'X_pca' in adata.obsm:
            params['use_rep'] = 'X_pca'
            params['n_pcs'] = adata.obsm['X_pca'].shape[1]
        else:
            params['use_rep'] = 'X'
    if params['use_rep'] == 'X':
        adata.X = read_anndata(input_file, X='X').X
        ensure_dense(adata)
    elif params['use_rep'] == 'X_pca' and 'X_pca' not in adata.obsm:
        # compute PCA if PCA is missing
        adata.X = read_anndata(input_file, X='X').X
        ensure_sparse(adata)
        pca_params = adata.uns.get('preprocessing', {}).get('pca', {})
        logging.info(f'Compute PCA with parameters: {pca_params}...')
        files_to_overwrite.extend(['obsm', 'varm'])
        sc.pp.pca(adata, **pca_params)
    
    # set n_neighbors
    params['n_neighbors'] = min(params.get('n_neighbors', 15), adata.n_obs)
    
    # compute kNN graph
    logging.info(f'parameters: {params}')
    try:
        sc.pp.neighbors(adata, **params)
    except Exception as e:
        if not USE_GPU:
            raise e
        logging.error(f"Recompute with scanpy because of error: {e}")
        import scanpy
        scanpy.pp.neighbors(adata, **params)
    finally:
        assert_neighbors(adata)

# update .uns
adata.uns |= extra_uns

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=files_to_overwrite,
)
