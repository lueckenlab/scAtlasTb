import scanpy as sc
import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from utils.misc import dask_compute, ensure_dense


sc.set_figure_params(frameon=False)
plt.rcParams['figure.figsize'] = 30, 25
input_file = snakemake.input[0]
output_png = snakemake.output.plots
Path(output_png).mkdir(parents=True, exist_ok=True)

wildcards = snakemake.wildcards
group_key = wildcards.group

markers = snakemake.params.markers
args = snakemake.params.get('args', {})
layer = args.pop('layer', snakemake.params.get('layer', 'X'))
args.pop('n_groups_per_split', None)
group_order = args.pop('group_order', None)


def parse_group_order(group_order):
    if isinstance(group_order, dict):
        group_order = group_order.get(group_key)

    if group_order is None:
        return None

    if isinstance(group_order, str):
        # If the provided string is a path to a text file, read the file
        # and treat each non-empty line as one group in order.
        p = Path(group_order)
        if p.exists() and p.is_file():
            try:
                order = [line.strip() for line in p.read_text().splitlines() if line.strip()]
            except Exception as exc:
                raise ValueError(f'Could not read group_order file {group_order}: {exc}')
        else:
            order = [group.strip() for group in group_order.split(',') if group.strip()]
    elif isinstance(group_order, (list, tuple)):
        order = [str(group) for group in group_order if group is not None and str(group) != '']
    else:
        raise ValueError(
            'Invalid group_order format. Expected a YAML list or comma-separated string '
            f'but got {type(group_order).__name__}: {group_order}'
        )

    if not order:
        logging.warning('group_order is empty. Continue without custom ordering.')
        return None
    return order


def unique_preserve_order(values):
    seen = set()
    ordered = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        ordered.append(value)
    return ordered


def filter_present_genes(genes, var_names):
    available = set(var_names)
    return [gene for gene in genes if gene in available]


def parse_marker_sets(markers):
    genes_to_plot = []
    gene_sets = {'markers': {}}

    if isinstance(markers, list):
        gene_sets['markers']['generic'] = markers
        genes_to_plot.extend(markers)
    elif isinstance(markers, dict):
        for key, value in markers.items():
            if isinstance(value, list):
                gene_sets['markers'][key] = value
                genes_to_plot.extend(value)
            elif isinstance(value, dict):
                gene_sets[key] = value
                genes_to_plot.extend(
                    gene
                    for genes in value.values()
                    for gene in genes
                )
            else:
                raise ValueError(f'Invalid markers format\n{key}: {value}')
    else:
        raise ValueError(f'Invalid markers format\n{markers}')

    if not gene_sets['markers']:
        del gene_sets['markers']

    return gene_sets, unique_preserve_order(genes_to_plot)


def apply_group_order(adata, group_key, group_order):
    group_order = parse_group_order(group_order)
    if group_order is None:
        return

    observed = adata.obs[group_key].dropna().astype(str).unique().tolist()
    ordered = [group for group in group_order if group in observed]
    if not ordered:
        logging.warning(
            f'No overlap between configured group_order and observed groups in {group_key}. '
            'Continue with default ordering.'
        )
        return

    remaining = [group for group in observed if group not in ordered]
    adata.obs[group_key] = pd.Categorical(
        adata.obs[group_key].astype(str),
        categories=ordered + remaining,
    )
    logging.info(f'Applied group_order: {len(ordered)} configured + {len(remaining)} remaining groups.')

logging.info(f'Reading {input_file}...')
adata = read_anndata(
    input_file,
    X=layer,
    obs='obs',
    var='var',
    uns='uns',
    dask=True,
    backed=True,
)

if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name']

if group_key not in adata.obs.columns:
    raise KeyError(
        f'group_key "{group_key}" is not present in adata.obs. '
        f'Available columns: {pformat(adata.obs.columns.tolist())}'
    )

gene_sets, genes_to_plot = parse_marker_sets(markers)

# subset to genes for plotting
logging.info(f'{len(genes_to_plot)} genes requested...')
gene_mask = adata.var_names.isin(genes_to_plot)
logging.info(f'{gene_mask.sum()} genes found overlapping...')
adata = adata[:, gene_mask].copy()

logging.info('Load genes to memory...')
adata = dask_compute(adata)
adata = ensure_dense(adata)

apply_group_order(adata, group_key, group_order)

logging.info(adata.__str__())

if min(adata.shape) == 0:
    logging.info('No data, exit...')
    print(adata, flush=True)
    exit()

for set_name, gene_set in gene_sets.items():
    logging.info(f'Processing gene set "{set_name}"...')
    # match marker genes and var_names
    print(
        pformat({k: len(v) for k, v in gene_set.items()}),
        flush=True
    )
    gene_set = {
        k: filter_present_genes(v, adata.var_names)
        for k, v in gene_set.items()
    }
    
    # remove groups without matching genes
    gene_set = {k: v for k, v in gene_set.items() if len(v) > 0}
    print(
        f'matching genes:\n',
        pformat({k: len(v) for k, v in gene_set.items()}),
        flush=True
    )
    
    if not gene_set:
        logging.info(f'No genes in marker gene set {set_name}, skipping...\n{gene_set}')
        continue

    logging.info('Dotplot...')
    dp = sc.pl.dotplot(
        adata,
        gene_set,
        groupby=group_key,
        use_raw=False,
        standard_scale='var',
        # title=title,
        show=False,
        swap_axes=False,
        return_fig=True,
    )
    dp.add_totals().style(dot_edge_color='black', dot_edge_lw=0.5)
    dp.savefig(f'{output_png}/{set_name}.png', dpi=100, bbox_inches='tight')    
    
    logging.info('Matrixplot...')
    mp = sc.pl.matrixplot(
        adata,
        gene_set,
        groupby=group_key,
        use_raw=False,
        standard_scale='var',
        show=False,
        swap_axes=False,
        return_fig=True,
    )
    mp.savefig(f'{output_png}/{set_name}_matrix.png', dpi=100, bbox_inches='tight')