from pathlib import Path
from tqdm import tqdm
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

import pandas as pd
from matplotlib import pyplot as plt
import anndata as ad
import scanpy as sc
import numpy as np
from dask import config as da_config

from load_data_utils import SCHEMAS, get_union
from utils.io import read_anndata, write_zarr
from utils.misc import ensure_sparse, dask_compute

in_file = snakemake.input.h5ad
schema_file = snakemake.input.get('schema')
annotation_file = snakemake.input.get('annotation_file')
out_file = snakemake.output.zarr
# out_plot = snakemake.output.plot

dask = snakemake.params.get('dask', True)
meta = snakemake.params.get('meta', {})
logging.info(f'meta:\n{pformat(meta)}')


# h5ad
logging.info(f'\033[0;36mread\033[0m {in_file}...')
try:
    adata = read_anndata(
        in_file,
        dask=True,
        backed=True,
        stride=500_000,
    )
    adata = ensure_sparse(adata)
except Exception as e:
    print(e)
    adata = sc.read_loom(in_file, sparse=True)
logging.info(adata)

# remove all unneeded uns entries
for key in list(adata.uns.keys()):
    if key in ['schema_version', 'batch_condition']:
        continue
    del adata.uns[key]

# ensure raw counts are kept in .X
adata.X = adata.raw.X if isinstance(adata.raw, ad._core.raw.Raw) else adata.X
del adata.raw

# # plot count distribution -> save to file
# x = to_memory(adata.X)
# plt.hist(x.data, bins=60)
# plt.savefig(out_plot)
# del x

# Adding general dataset info to uns and obs
adata.uns['meta'] = meta
for meta_i in ["organ", "study", "dataset"]:
    adata.obs[meta_i] = meta[meta_i]
    adata.uns[meta_i] = meta[meta_i]

# add annotation if available
if annotation_file is not None:
    logging.info(f'Add annotations from {annotation_file}...')
    barcode_column = meta['barcode_column']
    annotation = pd.read_csv(annotation_file, low_memory=False)

    # remove duplicates
    annotation = annotation.drop_duplicates(subset=barcode_column)
    print(annotation.head(), flush=True)

    # set index
    annotation.index = annotation[barcode_column].astype(str)
    adata.obs.index = adata.obs.index.astype(str)

    # remove column if existing to avoid conflict
    #if author_annotation in adata.obs.columns:
    #    logging.info(f'column {author_annotation} already exists, removing it')
    #    del adata.obs[author_annotation]

    # merge annotations
    for col in tqdm(annotation.columns, mininterval=5):
        adata.obs[col] = annotation[col]
        if adata.obs[col].isna().all():
            logging.warning(f'column {col} is empty or doesn\'t map, skipping')
            del adata.obs[col]


# assign sample and donor variables
adata.obs['donor'] = adata.obs[meta.pop('donor_column')]

tech_id_columns = [s.strip() for s in meta.pop('tech_id').split('+')]
logging.info(f'tech_id_columns: {tech_id_columns}')
adata.obs['tech_id'] = adata.obs[tech_id_columns].astype(str).apply(lambda x: '-'.join(x), axis=1)
adata.obs['sample'] = adata.obs['tech_id'] # keep for backwards compatibility

# CELLxGENE specific
if 'batch_condition' in adata.uns.keys():
    batch_columns = adata.uns['batch_condition']
    adata.obs['batch_condition'] = adata.obs[batch_columns].astype(str).apply(lambda x: '-'.join(x), axis=1)
else:
    adata.obs['batch_condition'] = meta['study']

# Checking schema version
if 'schema_version' not in adata.uns.keys():
    adata.uns['schema_version'] = float('nan')
if adata.uns['schema_version'] == '2.0.0':
    adata.obs['self_reported_ethnicity'] = adata.obs['ethnicity']
    adata.obs['self_reported_ethnicity_ontology_term_id'] = adata.obs['ethnicity_ontology_term_id']
    adata.obs['donor_id'] = adata.obs['donor']

# add author annotations column
author_annotation = meta.pop('author_annotation')
adata.obs['author_annotation'] = adata.obs[author_annotation]
logging.info(f'author_annotation: {author_annotation}')
print(pformat(list(adata.obs['author_annotation'].unique())), flush=True)
# use author annotations if no cell ontology available
if 'cell_type' not in adata.obs.columns:
    adata.obs['cell_type'] = 'nan'
if adata.obs['cell_type'].nunique() == 1:
    adata.obs['cell_type'] = adata.obs['author_annotation']

# Assigning other keys in meta to obs
for key, value in meta.items():
    if isinstance(value, list):
        continue
    adata.obs[key] = value

# save barcodes in separate column
if 'barcode' not in adata.obs.columns:
    adata.obs['barcode'] = adata.obs_names
adata.obs_names = adata.obs[['barcode', 'tech_id']].astype(str).agg('-'.join, axis=1)

# schemas translation
if schema_file is not None:
    schemas_df = pd.read_table(schema_file).dropna()
    logging.info(schemas_df)
    from_schema = meta['schema']
    assert from_schema in schemas_df.columns
    to_schema = 'cellxgene'
    SCHEMAS["NAMES"] = dict(zip(schemas_df[from_schema], schemas_df[to_schema]))
    adata.obs.rename(SCHEMAS["NAMES"], inplace=True)

# making sure all columns are in the object
all_columns = get_union(SCHEMAS["CELLxGENE_OBS"], SCHEMAS["TIER1"], SCHEMAS["EXTRA_COLUMNS"])

keep_covariates = meta.get('keep_covariates')
if isinstance(keep_covariates, str):
    keep_covariates = meta['keep_covariates'].split(',')
    all_columns = get_union(all_columns, keep_covariates)

for column in all_columns:
    if column not in adata.obs.columns:
        adata.obs[column] = np.nan
# keep only relevant columns
adata.obs = adata.obs[all_columns].copy()

# make sure all vars are present
if "feature_name" not in adata.var.columns:
    adata.var["feature_name"] = adata.var_names.tolist()
for column in SCHEMAS["CELLxGENE_VARS"]:
    if column not in adata.var.columns:
        adata.var[column] = np.nan
adata.var = adata.var[SCHEMAS["CELLxGENE_VARS"]]

if 'feature_id' in adata.var.columns:
    adata.var_names = adata.var['feature_id']
    del adata.var['feature_id']
adata.var.index.set_names('feature_id', inplace=True)

logging.info(f'\033[0;36mwrite\033[0m {out_file}...')
with da_config.set(num_workers=snakemake.threads):
    write_zarr(adata, out_file, compute=not dask)
