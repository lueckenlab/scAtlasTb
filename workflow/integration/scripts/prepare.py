from anndata import AnnData
from dask import array as da
import numpy as np
from pathlib import Path
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, write_zarr_linked, get_store
from utils.accessors import subset_hvg
from utils.misc import dask_compute
from utils.processing import assert_neighbors, sc, _filter_genes


input_file = snakemake.input[0]
output_file = snakemake.output[0]
params = snakemake.params
norm_count_layer = params.get('norm_counts', 'X')
raw_count_layer = params.get('raw_counts', 'X')
var_mask = snakemake.wildcards.var_mask
save_subset = params.get('save_subset', False)
recompute_pca = params.get('recompute_pca', True)


def read_and_subset(
    input_file: [str, Path],
    in_layer: str,
    out_layer: str,
    var_column: str,
    files_to_keep: list,
    slot_map: dict,
    new_var_column: str = None,
    filter_zero_genes: bool = True,
    **kwargs,
):
    assert in_layer is not None, f'Please specify a layer key in the config under {in_layer}'
    
    logging.info(f'Read {input_file}...')
    adata = read_anndata(
        input_file,
        X=in_layer,
        obs='obs',
        var='raw/var' if 'raw/' in in_layer else 'var',
        uns='uns',
        backed=True,
        dask=True,
        **kwargs
    )
    
    if new_var_column is None:
        slot_map |= {f'layers/{out_layer}': in_layer}
        return adata, files_to_keep, slot_map
    
    logging.info('Determine var_mask...')
    subsetted = False
    if var_column is None or var_column == 'None':
        adata.var[new_var_column] = True
    else:
        adata.var[new_var_column] = adata.var[var_column]
    
    if filter_zero_genes:
        logging.info('Filter all zero genes...')
        # filter out which genes are all 0
        filtered_genes = _filter_genes(
            adata[:, adata.var[new_var_column]].copy(),
            min_cells=1,
            return_varnames=True,
        )
        adata.var[new_var_column] &= adata.var_names.isin(filtered_genes)
    
    if save_subset:
        logging.info('Subset HVG...')
        adata, subsetted = subset_hvg(
            adata,
            var_column=new_var_column,
            to_memory=False,
            compute_dask=False,
        )
        dask_compute(adata, layers=[in_layer], verbose=True)
    
    # determine output
    if subsetted and save_subset:
        files_to_keep.append(out_layer)
    else:
        slot_map |= {f'layers/{out_layer}': in_layer}
    
    return adata, files_to_keep, slot_map


files_to_keep = ['obs', 'var']
slot_map = {}

# check if var_mask exists in data
store = get_store(input_file)
if var_mask not in (None, 'None'):
    assert var_mask in store['var'], f'var_mask="{var_mask}" not found in {input_file}.'

adata_norm, files_to_link, slot_map = read_and_subset(
    input_file=input_file,
    in_layer=norm_count_layer,
    out_layer='normcounts',
    var_column=var_mask,
    new_var_column='integration_features',
    files_to_keep=files_to_keep,
    slot_map=slot_map,
    obsm='obsm',
    obsp='obsp',
    varm='varm',
    varp='varp',
)
adata_raw, files_to_link, slot_map = read_and_subset(
    input_file=input_file,
    in_layer=raw_count_layer,
    out_layer='counts',
    var_column=var_mask,
    files_to_keep=files_to_keep,
    slot_map=slot_map,
    filter_zero_genes=False,
)

assert adata_norm.var_names.equals(adata_raw.var_names), f'\nnorm:\n{adata_norm.var}\nraw:\n{adata_raw.var}'
assert adata_norm.n_obs == adata_raw.n_obs, f'\nnorm:\n{adata_norm.obs}\nraw:\n{adata_raw.obs}'

if input_file.endswith('.h5ad'):
    adata = read_anndata(
        input_file,
        obs='obs',
        obsm='obsm',
        obsp='obsp',
        uns='uns',
    )
    adata = AnnData(
        X=adata_norm.X,
        obs=adata.obs,
        var=adata_norm.var,
        obsm=adata.obsm,
        obsp=adata.obsp,
        uns=adata.uns,
        layers={
            'normcounts': adata_norm.X,
            'counts': adata_raw.X,
        },
    )
    files_to_keep.append('X')
elif input_file.endswith(('.zarr', '.zarr/')):
    if save_subset:
        files_to_keep.extend(['varm', 'varp'])
    adata = AnnData(
        obs=adata_norm.obs,
        obsm=adata_norm.obsm,
        obsp=adata_norm.obsp,
        var=adata_norm.var,
        layers={
            'normcounts': adata_norm.X,
            'counts': adata_raw.X,
        },
        uns=adata_norm.uns,
    )
else:
    raise ValueError(f'Invalid input file {input_file}')

# # preprocess if missing
# if recompute_pca or 'X_pca' not in adata.obsm:
#     dask_compute(adata, layers=['normcounts'], verbose=True)
#     logging.info('Compute PCA...')
#     import scanpy
#     scanpy.pp.pca(adata, layer='normcounts', mask_var='highly_variable')
#     files_to_keep.extend(['obsm', 'uns'])

# try:
#     assert not recompute_pca, 'PCA was recomputed'
#     assert_neighbors(adata)
#     logging.info(adata.uns['neighbors'].keys())
# except AssertionError as e:
#     logging.info(f'Compute neighbors due to Assertion Error: {e}...')
#     sc.pp.neighbors(adata, use_rep='X_pca')
#     files_to_keep.extend(['obsp', 'uns'])


logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
slot_map |= {'X': norm_count_layer}
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=files_to_keep,
    slot_map=slot_map,
)