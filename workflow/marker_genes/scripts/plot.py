import numpy as np
import pandas as pd
import scanpy as sc
from pathlib import Path
import matplotlib
matplotlib.use('Agg')  # only save plots to file, never open windows
from matplotlib import pyplot as plt
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from utils.misc import dask_compute, ensure_dense, gene_labels

input_file = snakemake.input.zarr
input_markers = snakemake.input.markers
output_rankplot = snakemake.output.rankplot
output_dotplot = snakemake.output.dotplot
output_matrixplot = snakemake.output.matrixplot

Path(output_dotplot).mkdir(parents=True, exist_ok=True)
Path(output_matrixplot).mkdir(parents=True, exist_ok=True)

wildcards = snakemake.wildcards
group_key = wildcards.group
file_id = wildcards.file_id
title = f'Marker genes for file_id={file_id} group={group_key}'
marker_gene_key = f'marker_genes_group={group_key}'

args = snakemake.params.get('args', {})
group_order = args.pop('group_order', None) # TODO: unify group_order handling across all plotting rules
layer = args.pop('layer', snakemake.params.get('layer', 'X'))
args['key'] = marker_gene_key
n_genes = args.get('n_genes', 10)
n_groups_per_split = args.pop('n_groups_per_split', int(100 / n_genes))

markers = read_anndata(input_markers, uns='uns').uns
if marker_gene_key not in markers:
    logging.info(f'No marker genes found for group {group_key}, writing empty output')
    plt.savefig(output_rankplot)
    exit(0)

logging.info(f'Reading {input_file}...')
adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var',
    dask=True,
    backed=True,
)
adata.uns[marker_gene_key] = markers[marker_gene_key]

# get rank_genes_groups results
rank_genes_df = sc.get.rank_genes_groups_df(
    adata,
    key=marker_gene_key,
    group=None,
    log2fc_min=args.get('min_logfoldchange'),
)

if 'group' not in rank_genes_df.columns or rank_genes_df['group'].nunique() < 2:
    logging.info(f'Fewer than 2 groups found for group {group_key}, writing empty output')
    plt.savefig(output_rankplot)
    exit(0)

rank_genes_df = rank_genes_df.drop_duplicates(subset=['group', 'names'])

# partition groups to make plots visualisable
groups = rank_genes_df['group'].unique().tolist()
n_splits = max(1, int(len(groups) / n_groups_per_split))
splits = np.array_split(groups, n_splits)

# subset to genes for plotting
genes_to_plot = rank_genes_df.groupby('group').head(n_genes)['names'].unique().tolist()
labels = gene_labels(adata.var)
logging.info(f'Load {len(genes_to_plot)} genes to memory...')
adata = adata[
    adata.obs[group_key].isin(groups),
    adata.var_names.isin(genes_to_plot)
]
adata = dask_compute(adata.copy(), layers='X')

# show gene symbols instead of IDs, the gene ID is appended to duplicated symbols
adata.var_names = labels[adata.var_names].values
result = adata.uns[marker_gene_key]
result['names'] = pd.DataFrame(result['names']).apply(lambda names: names.map(labels)).to_records(index=False)
for pts_key in ['pts', 'pts_rest']:
    if pts_key in result:
        result[pts_key].index = result[pts_key].index.map(labels)
if adata.n_vars > 50:
    sc.pp.pca(adata)
else:
    ensure_dense(adata, layers='X')
logging.info(adata.__str__())

logging.info('Rankplot...')
try:
    sc.pl.rank_genes_groups(adata, **(args | {'show': False}))
    plt.suptitle(title)
    plt.savefig(output_rankplot, dpi=100, bbox_inches='tight')
except Exception as e:
    logging.error(e)
    logging.info('Error in rankplot, save empty plot...')
    plt.savefig(output_rankplot)

for i, partition_groups in enumerate(splits):
    logging.info(f'Dotplot for partition {i}...')
    dp = sc.pl.rank_genes_groups_dotplot(
        adata,
        **args,
        groups=partition_groups,
        return_fig=True,
    )
    dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5)
    dp.savefig(
        f'{output_dotplot}/split={i}_expression.png',
        dpi=100,
        bbox_inches='tight',
    )

    dp = sc.pl.rank_genes_groups_dotplot(
        adata,
        **args,
        groups=partition_groups,
        values_to_plot='logfoldchanges',
        cmap='bwr',
        vmin=-4,
        vmax=4,
        colorbar_title='log fold change',
        return_fig=True,
    )
    dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5)
    dp.savefig(
        f'{output_dotplot}/split={i}_logfoldchanges.png',
        dpi=100,
        bbox_inches='tight',
    )

    logging.info(f'Matrixplot for partition {i}...')
    mp = sc.pl.rank_genes_groups_matrixplot(
        adata,
        **args,
        groups=partition_groups,
        return_fig=True,
    )
    mp.add_totals().style(edge_color='black')
    mp.savefig(
        f'{output_matrixplot}/split={i}_expression.png',
        dpi=100,
        bbox_inches='tight'
    )

    mp = sc.pl.rank_genes_groups_matrixplot(
        adata,
        **args,
        groups=partition_groups,
        values_to_plot='logfoldchanges',
        cmap='bwr',
        vmin=-4,
        vmax=4,
        colorbar_title='log fold change',
        return_fig=True,
    )
    mp.add_totals().style(edge_color='black')
    mp.savefig(
        f'{output_matrixplot}/split={i}_logfoldchanges.png',
        dpi=100,
        bbox_inches='tight'
    )
