import glob
import warnings
from pprint import pprint
import numpy as np
import anndata as ad
from anndata.experimental import read_elem
import zarr
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata


# check first prepare output linking
prep_file = glob.glob('test/out/integration/prepare/dataset~preprocessed/file_id~*/var_mask~*.zarr')[0]
print(prep_file)
ad.read_zarr(prep_file)
read_anndata(prep_file)

# direct integration outputs
outputs = glob.glob('test/out/integration/dataset~*/file_id~*/batch~*/var_mask~*/method~*/adata.zarr')
outputs = glob.glob('test/out/integration/postprocess/dataset~*/file_id~*/batch~*/var_mask~*/method~*--output_type~*.zarr')
if len(outputs) == 0:
    warnings.warn('No integration outputs found')

logging.info('Checking file outputs...')

for file in outputs:
    # logging.info(f'Checking {file}...')
    z = zarr.open(file)
    uns = read_elem(z["uns"])

    try:
        # Check Metadata
        assert 'wildcards' in uns
        assert 'integration_batch' in uns['wildcards']
        assert 'integration_hyperparams' in uns['wildcards']
        assert 'integration_label' in uns['wildcards']
        assert 'integration_method' in uns['wildcards']
        # assert 'dataset' in uns
        # assert 'methods' in uns
        assert 'integration' in uns
        for key in ['method', 'label_key', 'batch_key', 'output_type']:
            assert key in uns['integration']

        # Unintegrated data
        # adata_raw = raw.to_adata()
        # assert 'X_pca' in adata_raw.obsm
        # assert 'connectivities' in adata_raw.obsp
        # assert 'distances' in adata_raw.obsp

        # Output type specific outputs
        output_types = uns['integration']['output_type']
        output_types = [output_types] if isinstance(output_types, str) else output_types

        if uns['integration']['method'] == 'scanorama':
            assert len(output_types) == 2

        for ot in output_types:
            if ot not in ['knn', 'embed', 'full']:
                raise ValueError(f'Invalid output type {ot}')

        if 'knn' in output_types:
            assert 'obsp' in z, f'No obsp for output type knn\nfile: {file}'
            obsp = read_elem(z["obsp"])
            assert 'connectivities' in obsp
            assert 'distances' in obsp
        if 'embed' in output_types:
            assert 'obsm' in z, f'No obsm for output type embed\nfile: {file}'
            obsm = read_elem(z["obsm"])
            assert 'X_emb' in obsm
        if 'full' in output_types:
            assert 'X' in z, f'No X for output type full\nfile: {file}'
            X = read_elem(z['X'])
            assert X is not None
            # linked layers must match the (possibly feature-subsetted) var
            adata = read_anndata(file, X='X', var='var', layers='layers', dask=True, backed=True)
            assert X.shape[1] == adata.n_vars, \
                f'X has {X.shape[1]} features, but var has {adata.n_vars}\nfile: {file}'
            for key, layer in adata.layers.items():
                assert layer.shape == adata.shape, \
                    f'layer "{key}" has shape {layer.shape}, expected {adata.shape}\nfile: {file}'
            # assert 'corrected_counts' in layers
            # check that counts are different from input
            # if uns['integration']['method'] != 'unintegrated':
                # assert adata.X.data.shape != adata_raw.X.data.shape or np.any(X.data != adata_raw.X.data)

    except AssertionError as e:
        print('AssertionError for:', file)
        pprint(uns.keys())
        raise e
