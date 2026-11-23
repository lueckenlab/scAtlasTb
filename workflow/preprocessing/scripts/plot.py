import logging
logging.basicConfig(level=logging.INFO)
from pathlib import Path
import traceback
import warnings
warnings.filterwarnings("ignore")

from matplotlib import pyplot as plt
import matplotlib as mpl
import scanpy as sc
import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype, is_string_dtype, is_categorical_dtype
from pprint import pformat
from tqdm import tqdm
from joblib import Parallel, delayed

from utils.io import read_anndata, get_file_reader
from utils.misc import ensure_dense, remove_outliers, dask_compute
from utils.accessors import parse_gene_names


input_file = snakemake.input[0]
output_dir = Path(snakemake.output.plots)
output_dir.mkdir(exist_ok=True)

params = dict(snakemake.params.items())
basis = params.pop('basis', 'X_umap')
gene_chunk_size = params.pop('gene_chunk_size', 12)
plot_centroids = params.pop('plot_centroids', [])

wildcards_string = '\n'.join([f'{k}: {v}' for k, v in snakemake.wildcards.items()])
logging.info(f'Wildcard string: {wildcards_string}')

logging.info(f'Read file {input_file}...')
kwargs = dict(
    obs='obs',
    obsm='obsm',
    var='var',
)

# check if .X exists
read_func, _ = get_file_reader(input_file)
if 'X' in read_func(input_file, 'r'):
    kwargs |= dict(X='X', dask=True, backed=True)

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file, **kwargs)
assert basis in adata.obsm.keys(), f'"{basis}" not in adata.obsm'
n_cells = adata.n_obs # save original number of cells
obs_columns = adata.obs.columns.tolist()
ensure_dense(adata, basis)

if adata.n_obs == 0:
    logging.info('No cells, skip...')
    exit()

if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name'].astype(str)

# set minimum cells per category
min_cells_per_category = params.pop('min_cells_per_category', 1e-4)
if min_cells_per_category < 1:
    min_cells_per_category *= n_cells
logging.info(f'Remove categories with fewer than {min_cells_per_category:.1f} cells')

# parse colors
colors = params.pop('color', None)
logging.info(f'Configured colors:\n{pformat(colors)}')
colors = colors if isinstance(colors, list) else [colors]
colors += plot_centroids
colors = list(dict.fromkeys(colors))

# get gene colors
gene_colors = [col for col in colors if col not in obs_columns]
gene_colors = parse_gene_names(adata, gene_colors)
gene_colors.sort()

# filter colors that aren't in object
colors = [
    color for color in colors
    if color in obs_columns
    and adata.obs[color].nunique() > 1
]
logging.info(f'Colors from obs after filtering:\n{pformat(colors)}')

for color in colors:
    column = adata.obs[color]

    if is_categorical_dtype(column) or is_string_dtype(column):
        # parse data types
        column = column.astype(object) \
            .replace(['NaN', 'None', '', 'nan', 'unknown'], float('nan')) \
            .astype('category')

        # remove categories with too few cells
        value_counts = column.value_counts()
        categories_to_remove = value_counts[value_counts <= min_cells_per_category].index
        column = column.cat.remove_categories(categories_to_remove)
        
        # update column
        adata.obs[color] = column
        # adata.obs[color] = column.codes if len(column.categories) > 102 else column

if len(colors) == 0:
    logging.info('No valid colors, skip...')
    colors = [None]


logging.info('Remove outliers...')
outlier_factor = params.pop('outlier_factor', 0)
adata = remove_outliers(adata, 'max', factor=outlier_factor, rep=basis)
adata = remove_outliers(adata, 'min', factor=outlier_factor, rep=basis)

# match gene patterns
if gene_colors:
    logging.info(f'Subset to {len(gene_colors)} requested genes...')
    adata = adata[:, adata.var_names.isin(gene_colors)].copy()
    logging.info(adata.__str__())
    adata = dask_compute(adata, layers='X')
    print(adata.var, flush=True)
else:
    del adata.X
    del adata.var

logging.info('Shuffle cells...')
if n_cells > 1e6:
    adata = adata[adata.obs.sample(frac=0.7).index]

if adata.is_view:
    logging.info('Convert view to copy...')
    adata = adata.copy()

# set minimum point size
min_size = 0.5 if adata.n_obs > 500_000 else 1
default_size = max(min_size, 200_000 / adata.n_obs)
size = params.get('size', default_size)
if size is None:
    size = default_size
params['size'] = min(200, max(default_size, size))


def plot_centroids_on_embedding(ax, adata, color, basis, legend, category_numbers, legend_fontsize=10):
    """
    Plot category numbers at centroid positions on embedding.
    
    Args:
        ax: Matplotlib axes object
        adata: AnnData object
        color: Column name in adata.obs for grouping
        basis: Key in adata.obsm for coordinates
        legend: Matplotlib legend object
        category_numbers: Dict mapping categories to their numbers
        legend_fontsize: Font size for category labels
    """
    categories = [cat for cat in adata.obs[color].cat.categories if cat in adata.obs[color].unique()]
    
    # Compute centroids for each category
    coords = adata.obsm[basis][:, :2]
    centroids = pd.DataFrame(coords, index=adata.obs[color]) \
        .groupby(level=0, dropna=False) \
        .median() \
        .reindex(categories)
    
    # Extract colors from legend handles
    color_map = {
        text.get_text(): handle.get_facecolor()[0]
        for handle, text in zip(legend.legend_handles, legend.get_texts())
        if hasattr(handle, 'get_facecolor')
    }
    
    # Plot category numbers at centroids with matching colors
    for cat, row in centroids.iterrows():
        bg_color = color_map.get(cat, 'white')
        
        # Convert color to RGB for luminance calculation
        try:
            # Use matplotlib to convert any color format to RGBA
            r, g, b = mpl.colors.to_rgba(bg_color)[:3]
            luminance = 0.299 * r + 0.587 * g + 0.114 * b
            text_color = 'white' if luminance < 0.4 else 'black'
        except (ValueError, TypeError):
            text_color = 'white'
        
        label_value = category_numbers[cat]
        ax.text(
            row.iloc[0], row.iloc[1],
            s=str(label_value),
            fontsize=legend_fontsize,
            fontweight="bold",
            ha="center",
            va="center",
            color=text_color,
            bbox=dict(
                boxstyle="circle",
                facecolor='none',
                alpha=0.4,
                edgecolor='none',
            )
        ).set_path_effects([
            mpl.patheffects.Stroke(
                linewidth=1.5,
                foreground=color_map.get(cat, 'white')
            ),
            mpl.patheffects.Normal()
        ])

    # Add category numbers and group sizes to legend labels.
    category_numbers_str = {str(k): v for k, v in category_numbers.items()}

    def formatter(label):
        mapped = category_numbers_str.get(label)
        return f"{mapped}: {label}" if isinstance(mapped, int) else label

    add_group_sizes_to_legend(legend=legend, adata=adata, color=color, label_formatter=formatter)


def add_group_sizes_to_legend(legend, adata, color, label_formatter=None):
    """Append group sizes to legend labels for categorical columns."""
    category_counts_str = {
        str(category): int(count)
        for category, count in adata.obs[color].value_counts(dropna=False).items()
    }
    for text in legend.get_texts():
        label = text.get_text()
        count = category_counts_str.get(label)
        if count is not None:
            prefix = label if label_formatter is None else label_formatter(label)
            text.set_text(f"{prefix} (n={count})")


def plot_color(
    adata,
    color,
    basis,
    n_cells=None,
    plot_centroids=False,
    verbose=True,
    file_name=None,
    title='',
    output_dir=Path(),
    **kwargs
):

    if n_cells is None:
        n_cells = adata.n_obs
    
    palette = None
    if file_name is None:
        file_name = str(color)
    colors = color if isinstance(color, list) else [color]

    fig_params = dict(
        frameon=False,
        vector_friendly=True,
        fontsize=9,
        figsize=(6,6),
        dpi=200,
        dpi_save=200,
        format='png',
    )
    
    if len(colors) > 4:
        fig_params = dict(frameon=False)
    sc.set_figure_params(**fig_params)
    mpl.rcParams['figure.constrained_layout.use'] = True

    # select palette according to obs column
    for color in colors:
        if color in adata.obs.columns:
            color_vec = adata.obs[color]
            if is_categorical_dtype(color_vec):
                if color_vec.nunique() > 102:
                    palette = 'turbo'
                elif color_vec.nunique() > 20:
                    palette = sc.pl.palettes.godsnot_102
            elif is_numeric_dtype(color_vec):
                if color_vec.min() < 0:
                    palette = 'coolwarm'
                else:
                    palette = 'plasma'

    # centroid plotting setup (sctk-like numbering)
    color = colors[0] if colors else None
    plot_centroids = (
        plot_centroids
        and len(colors) == 1
        and color in adata.obs.columns
        and is_categorical_dtype(adata.obs[color])
        and adata.obs[color].nunique() <= 102
    )

    try:
        fig = sc.pl.embedding(
            adata,
            basis=basis,
            color=colors,
            show=False,
            return_fig=True,
            palette=palette,
            **kwargs,
        )
        suptitle_text = f'{title}\nn={n_cells}'
        suptitle_fontsize = 12
        fig.suptitle(suptitle_text, fontsize=suptitle_fontsize)

        # fix legend
        ax = fig.get_axes()[0]
        legend = ax.get_legend()

        if palette == 'turbo':
            if legend:
                legend.remove()
                legend = None
        
        elif legend and plot_centroids:
            categories = [cat for cat in adata.obs[color].cat.categories if cat in adata.obs[color].unique()]
            category_numbers = {
                cat: idx + 1 if len(str(cat)) > 3 else cat
                for idx, cat in enumerate(categories)
            }
            fontsize = kwargs.get('legend_fontsize', 10)
            plot_centroids_on_embedding(
                ax=ax,
                adata=adata,
                color=color,
                basis=basis,
                legend=legend,
                category_numbers=category_numbers,
                legend_fontsize=fontsize,
            )

        elif legend and len(colors) == 1 and is_categorical_dtype(adata.obs[color]):
            add_group_sizes_to_legend(
                legend=legend,
                adata=adata,
                color=color,
            )


        if legend:
            # adjust legend
            legend_bbox = legend.get_window_extent()
            fig_width, fig_height = fig.get_size_inches()
            fig_width = fig_width + (legend_bbox.width / fig.dpi)
            fig.set_size_inches((fig_width, fig_height))

        if verbose:
            logging.info(f'Plotting color "{file_name}" successful.')
    
    except Exception as e:
        traceback.print_exc()
        logging.error(f'Failed to plot "{file_name}": {e}')
        plt.plot([])
    
    out_path = output_dir / f'{file_name}.png'
    
    try:        
        plt.savefig(out_path, bbox_inches='tight')
    except Exception as e:
        logging.error(f'Failed to save plot "{file_name}" to {out_path}: {e}')
        traceback.print_exc()
    finally:
        plt.close('all')  # Free memory from figure


logging.info('Parameters:\n' + pformat(params))
# Run plotting in parallel
list(tqdm(
    Parallel(return_as='generator', backend='threading')(delayed(plot_color)(
        adata=adata,
        color=color,
        basis=basis,
        n_cells=n_cells,
        plot_centroids=color in plot_centroids,
        **params,
        title=wildcards_string,
        file_name=color,
        output_dir=output_dir,
    ) for color in colors),
    desc="Plotting colors",
    total=len(colors),
    miniters=1,
))

gene_colors = {
    f'genes_group={idx}': gene_colors[i:i + gene_chunk_size]
    for idx, i in enumerate(range(0, len(gene_colors), gene_chunk_size))
}
if gene_colors:
    params['ncols'] = params.get('ncols', 4)
    list(tqdm(
        Parallel(return_as='generator', backend='threading')(delayed(plot_color)(
            adata=adata,
            color=color,
            basis=basis,
            n_cells=n_cells,
            verbose=False,
            **params,
            title=wildcards_string,
            file_name=title,
            output_dir=output_dir,
        ) for title, color in gene_colors.items()),
        desc="Plotting genes",
        total=len(gene_colors),
        miniters=1,
    ))