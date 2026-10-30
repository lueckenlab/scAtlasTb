import pandas as pd
from pprint import pformat
from utils.io import read_anndata, write_zarr_linked

input_file = snakemake.input[0]
output_file = snakemake.output[0]
gene_scores_available = 'gene_scores' in snakemake.input.keys()

files_to_keep = ['uns', 'obs'] if gene_scores_available else ['uns']

adata = read_anndata(input_file, **dict.fromkeys(files_to_keep), verbose=False)

metrics_df = pd.read_table(snakemake.input.metrics, sep='\t')
metrics_df = metrics_df.apply(lambda col: col.astype(str) if col.dtype in ('object', 'category') else col)
adata.uns['metrics'] = metrics_df

if gene_scores_available:
    gs_obs = read_anndata(str(snakemake.input.gene_scores), obs='obs', verbose=False).obs
    gene_score_cols = [c for c in gs_obs.columns if c.startswith('gene_score:')]
    print(f'Add {len(gene_score_cols)} gene score columns to adata.obs:\n{pformat(gene_score_cols)}')
    adata.obs[gene_score_cols] = gs_obs[gene_score_cols]

write_zarr_linked(adata, input_file, output_file, files_to_keep=files_to_keep)