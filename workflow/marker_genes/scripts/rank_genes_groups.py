import numpy as np
import pandas as pd
from pprint import pprint
import warnings
warnings.filterwarnings('ignore')
import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc
rank_gene_func = sc.tl.rank_genes_groups

from utils.io import read_anndata, write_zarr_linked
from utils.accessors import subset_hvg
from utils.processing import filter_genes, get_pseudobulks #, sc, USE_GPU

# if USE_GPU:
#     rank_gene_func = sc.tl.rank_genes_groups_logreg
# else:
#     rank_gene_func = sc.tl.rank_genes_groups

def write_empty(adata, output_file, output_tsv):
    write_zarr_linked(
        adata,
        in_dir=input_file,
        out_dir=output_file,
        slot_map={'X': layer}
    )

    markers_df = pd.DataFrame(
        columns=[
            'gene',
            'z-score',
            'logfoldchange',
            'pval',
            'pval_adj',
            'cluster',
            '-log10 pvalue',
        ]
    )
    markers_df.to_csv(output_tsv, sep='\t')


input_file = snakemake.input[0]
output_file = snakemake.output.zarr
output_tsv = snakemake.output.tsv
group_key = snakemake.wildcards.group
sample_key = snakemake.params.sample
layer = snakemake.params.layer
args = snakemake.params.get('args', {})
pseudobulk = not sample_key is None
min_groups_size = snakemake.params.get('min_groups_size', 10)
min_cells = snakemake.params.get('min_cells', 100)
key = f'marker_genes_group={group_key}'

# determine how to read the data
obs_cols = [group_key, sample_key] if pseudobulk else [group_key]
obs = read_anndata(input_file, obs='obs').obs[obs_cols]
n_obs = obs.shape[0]
n_groups = obs.drop_duplicates().shape[0]
dask = pseudobulk or (n_obs > 1e6 or n_groups > 1e5)

logging.info(f'Reading {input_file}...')
adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var',
    uns='uns',
    dask=dask,
    backed=dask,
    stride=int(n_obs / 5),
)
logging.info(adata.__str__())

if adata.n_obs <= min_cells:
    logging.info('Too few cells in object, skipping marker genes')
    write_empty(adata, output_file, output_tsv)
    exit(0)

if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name'].astype(str)

duplicated = adata.var_names.duplicated(keep='first')
if duplicated.any():
    adata = adata[:, ~duplicated].copy()


# logging.info('Subset highly variable genes...')
# adata, _ = subset_hvg(adata, var_column='highly_variable')
# adata = filter_genes(adata, batch_key=sample_key, min_cells=1, min_counts=1)

# pseudobulk if sample key is set
if pseudobulk:
    adata.obs['pseudo_group'] = adata.obs[group_key].astype(str) + '_' + adata.obs[sample_key].astype(str)
    logging.info(f'Creating pseudobulks for sample={sample_key} and group={group_key}...')
    adata = get_pseudobulks(adata, group_key='pseudo_group', agg='sum')
    logging.info(f'Pseudobulk shape: {adata.shape}')

    # compute marker genes for all genes, because sample size is way more manageable
    args['n_genes'] = adata.n_vars

# filter groups with more than 1 observation
groups_counts = adata.obs[group_key].value_counts()
groups = groups_counts[groups_counts >= min_groups_size].index.tolist()
logging.info(groups_counts[groups])

if not groups:
    logging.info('No groups large enough for testing, writing empty output')
    write_empty(adata, output_file, output_tsv)
    exit(0)

logging.info(f'Running marker genes analysis for {group_key} and args={args}...')
rank_gene_func(
    adata,
    groupby=group_key,
    groups=groups,
    pts=True,
    use_raw=False,
    key_added=key,
    **args
)

logging.info(f'Writing {output_file}...')
if input_file.endswith('.h5ad'):
    logging.info('Retrieving original data from h5ad...')
    new_uns = adata.uns[key]
    adata = read_anndata(input_file, dask=True, backed=True)
    adata.uns[key] = new_uns
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    files_to_keep=['uns'],
    slot_map={'X': layer}
)


logging.info('Create DEG dataframe...')
result = adata.uns[key]
markers_df = []
for cluster in result['names'].dtype.names:
    # parse percent results
    pct_dict = {}
    for orig, new in dict(pts='pct_within', pts_rest='pct_outside').items():
        s = result[orig][cluster]
        # keep unique duplicates where values is max
        s = s.loc[s.groupby(level=0).idxmax()]
        pct_dict[new] = s
    pct_dict_df = pd.DataFrame(pct_dict).reset_index(names='gene')
    
    df = pd.DataFrame({
        'gene': result['names'][cluster],
        'z-score': result['scores'][cluster],
        'logfoldchange': result['logfoldchanges'][cluster],
        'pval': result['pvals'][cluster],
        'pval_adj': result['pvals_adj'][cluster], 
        'cluster': cluster
    })
    df = df.merge(pct_dict_df, on='gene', how='left')
    markers_df.append(df)
markers_df = pd.concat(markers_df).set_index('cluster')
markers_df['-log10 pvalue'] = -np.log10(markers_df['pval'])

markers_df.to_csv(output_tsv, sep='\t')