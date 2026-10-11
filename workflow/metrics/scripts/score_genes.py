from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
import numpy as np
import anndata as ad
from tqdm import tqdm

from metrics.bootstrap import get_bootstrap_adata
from utils.io import read_anndata, write_zarr_linked
from utils.processing import sc
from utils.misc import dask_compute
from utils.accessors import parse_gene_names


def bins_by_quantiles(matrix: [np.array, np.matrix], n_quantiles: int):
    assert n_quantiles > 1
    
    def bin_values(x):
        bins = np.linspace(0, 1, n_quantiles + 1)
        quantiles = np.quantile(x, bins)
        binned_idx = np.digitize(x, quantiles) - 1
        return bins[binned_idx]
    
    return np.apply_along_axis(bin_values, axis=0, arr=matrix)


input_file = snakemake.input[0]
output_file = snakemake.output[0]

params = snakemake.params
unintegrated_layer = params.get('unintegrated_layer', 'X')
raw_counts_layer = params.get('raw_counts_layer', unintegrated_layer)
gene_sets = params['gene_sets']
n_random_permutations = params.get('n_permutations', 100)
n_random_genes = params.get('n_random_genes', 50)
ctrl_size = params.get('ctrl_size', 50)
n_quantiles = params.get('n_quantiles', 2)
MAX_OBS = int(params.get('MAX_OBS', 4e6))

files_to_keep = ['obs', 'obsm']   # TODO: only overwrite specific obsm slots


logger.info(f'Read {input_file} ...')
adata = read_anndata(
    input_file,
    X=unintegrated_layer,
    obs='obs',
    obsm='obsm',
    var='var',
    dask=True,
    backed=True,
)

# adata.layers['raw_counts'] = read_anndata(
#     input_file,
#     X=raw_counts_layer,
#     dask=True,
#     backed=True,
# ).X
 
if 'feature_name' in adata.var.columns:
    adata.var['feature_id'] = adata.var_names
    adata.var_names = adata.var['feature_name'].astype(str).values

all_obs_names = adata.obs_names.copy()
# feature_name can be non-unique (CxG feature collisions); use unique feature_id for the
# var subset_mask so duplicate columns dropped during dedup aren't double-counted.
all_var_ids = adata.var['feature_id'].copy() if 'feature_id' in adata.var.columns else adata.var_names.copy()

logging.debug('Parse gene sets...')
logging.debug(f'Initial gene sets: {gene_sets}')
logging.debug(f'adata.var_names: {adata.var_names}')

# filter all gene sets to genes in adata
for set_name, gene_list in gene_sets.items():
    gene_sets[set_name] = parse_gene_names(adata, gene_list)

# get all genes from gene sets
genes = list(set().union(*gene_sets.values()))

if len(genes) == 0:
    logging.warning('No genes of interest found in dataset')
    write_zarr_linked(
        adata,
        in_dir=input_file,
        out_dir=output_file,
    )
    exit(0)

# deal with duplicate gene names
if adata.var_names.duplicated().sum() > 0:
    duplicated_gene_mask = adata.var_names.duplicated()
    duplicated_gene_names = adata.var_names[duplicated_gene_mask]
    logger.warning(f"Found and removed {duplicated_gene_mask.sum()} duplicate gene names: {list(duplicated_gene_names)}")
    adata = adata[:, ~duplicated_gene_mask]

# subset adata if dataset too large TODO: move to prepare script?
if adata.n_obs > MAX_OBS:
    adata = get_bootstrap_adata(adata, size=MAX_OBS)

# Extract control genes and precompute random sets
control_genes = adata.var_names[~adata.var_names.isin(genes)].tolist()
effective_n_random_genes = min(n_random_genes, len(control_genes))
if effective_n_random_genes < n_random_genes:
    logger.warning(
        f'Only {len(control_genes)} control genes available, '
        f'reducing n_random_genes from {n_random_genes} to {effective_n_random_genes}'
    )
random_gene_sets = [
    np.random.choice(control_genes, size=effective_n_random_genes, replace=False)
    for _ in range(n_random_permutations)
]
all_random_genes = set().union(*random_gene_sets) if random_gene_sets else set()

# Include additional control genes for validation
n_control = min(ctrl_size, sum(1 for g in control_genes if g not in all_random_genes))
control_reserve = [g for g in control_genes if g not in all_random_genes][:n_control]

# Subset adata and prepare for scoring
genes_for_scoring = set(genes) | all_random_genes | set(control_reserve)
adata = adata[:, adata.var_names.isin(genes_for_scoring)].copy()
logging.info(f'Subset to {adata.n_vars} genes ({len(genes)} from user + {len(all_random_genes)} random + {n_control} control)')
adata = dask_compute(adata, layers='X')

# Score random gene sets
logging.info('Computing random gene scores...')
random_scores = []
for random_genes in random_gene_sets:
    sc.tl.score_genes(
        adata,
        gene_list=random_genes,
        ctrl_size=ctrl_size
    )
    random_scores.append(adata.obs['score'].values)
if 'score' in adata.obs.columns:
    del adata.obs['score']
adata.obsm['random_gene_scores'] = np.stack(random_scores, axis=1) if random_scores else np.empty((adata.n_obs, 0))
logging.info(f'Stored {len(random_scores)} random gene score sets in obsm')

# Score gene sets of interest
for set_name, gene_list in tqdm(gene_sets.items(), desc='Gene scores'):
    if len(gene_list) == 0:
        logging.warning(f'Gene set "{set_name}" is empty, skipping')
        continue
    sc.tl.score_genes(
        adata,
        gene_list=gene_list,
        ctrl_size=ctrl_size,
        score_name=f'gene_score:{set_name}'
    )

# Final subset to genes of interest only
adata = adata[:, genes].copy()

logging.info(f'Write to {output_file}...')
logging.info(adata.__str__())
obs_mask = all_obs_names.isin(adata.obs_names)
final_var_ids = adata.var['feature_id'] if 'feature_id' in adata.var.columns else adata.var_names
var_mask = np.asarray(all_var_ids.isin(final_var_ids))
write_zarr_linked(
    adata,
    in_dir=input_file,
    out_dir=output_file,
    subset_mask=(obs_mask, var_mask),
    files_to_keep=files_to_keep,
)