import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import pandas.api.types as ptypes
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg') # non-interactive backend for parallel plotting
import seaborn as sns
from tqdm import tqdm
import traceback
import concurrent.futures
from joblib import Parallel, delayed

from utils.io import read_anndata
from utils.plots import plot_stacked_bar, plot_violin


def nmi(df, x, y):
    from sklearn.metrics.cluster import normalized_mutual_info_score

    df = df[[x, y]]
    df = df[df[x].notna() & df[y].notna()]
    return normalized_mutual_info_score(df[x], df[y])


input_file = snakemake.input.zarr
output_dir = Path(snakemake.output.plots)
output_dir.mkdir(exist_ok=True)
cluster_keys = snakemake.params.get('cluster_keys')
covariates = snakemake.params.get('covariates')
max_categories = snakemake.params.get('max_categories', 100)
threads = snakemake.threads

logging.info(f'Read {input_file}...')
obs = read_anndata(input_file, obs='obs').obs

# filter covariates
logging.info(f'Categories: {covariates}')
covariates = set([c for c in covariates if c in obs.columns])
for col in list(covariates):
    if obs[col].dtype == 'object':
        obs[col] = obs[col].astype('category')
    if obs[col].dtype == 'category' and obs[col].nunique() > max_categories:
        logging.info(f'Skipping {col}, too many values')
        covariates.remove(col)
logging.info(f'Categories after filtering: {covariates}')

# Make obs read-only to prevent accidental modifications
obs.flags.writeable = False


def call_plot(covariate, cluster_key, output_dir):
    try:
        if obs[covariate].dtype == 'category':
            score = nmi(obs, cluster_key, covariate)
            fig = plot_stacked_bar(
                obs,
                category_key=cluster_key,
                covariate_key=covariate,
                count_type='proportion',
                title_suffix=f' (NMI: {score:.3f})',
                return_fig=True,
            )
        elif ptypes.is_numeric_dtype(obs[covariate]):
            fig = plot_violin(
                obs,
                category_key=cluster_key,
                covariate_key=covariate,
                return_fig=True,
            )
        else:
            raise ValueError(f'Invalid covariate type: {obs[covariate].dtype}, must be category or numeric')
        
        fig.savefig(output_dir / f'{cluster_key}--{covariate}.png', bbox_inches='tight')
        plt.close(fig)
        
    except Exception as e:
        logging.error(f'Error plotting {cluster_key}--{covariate}: {e}')
        traceback.print_exc()
    finally:
        # Ensure all figures are closed
        plt.close('all')


tasks = [(covariate, cluster_key) for covariate in covariates for cluster_key in cluster_keys]

results = Parallel(n_jobs=threads, backend='threading', return_as='generator')(
    delayed(call_plot)(covariate, cluster_key, output_dir)
    for covariate, cluster_key in tasks
)

for _ in tqdm(results, total=len(tasks), desc=f'Plotting with {threads} threads'):
    pass
