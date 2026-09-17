from pathlib import Path
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
import pprint
pp = pprint.PrettyPrinter(indent=2)
import logging
from tqdm import tqdm
import traceback
import gc
from joblib import Parallel, delayed
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend to reduce memory overhead

plt.rcParams['svg.fonttype'] = 'none'
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata, parse_args_, parse_set_nested, is_interactive
from qc_utils import parse_parameters, QC_FLAGS, log_auto

def plot_bar_all(adata, dataset, output_plots, dpi: int = 150):
    plt.figure(figsize=(4, 5))
    plt.grid(False)
    sns.countplot(
        x='qc_status',
        data=adata.obs,
        hue='qc_status',
    )
    ax = plt.gca()
    for pos in ['right', 'top']: 
        ax.spines[pos].set_visible(False)
    for container in ax.containers:
        ax.bar_label(container)
    plt.xlabel('Cell QC Status')
    plt.ylabel('Count')
    plt.title(f'Counts of cells QC\'d\n{dataset}')
    plt.tight_layout()
    plt.savefig(output_plots / 'bar_cells_passed_all.png', bbox_inches='tight', dpi=dpi)
    plt.close()


def plot_composition(
    df,
    group_key,
    plot_dir,
    dataset: str = "dataset",
    file_id: str = "file_id",
    dpi: int = 150,
):
    df[group_key] = pd.Categorical(df[group_key])
    n_hues = df[group_key].nunique()

    # Get counts by group and QC status
    counts = (
        df.groupby([group_key, "qc_status"], observed=False)
        .size()
        .unstack(fill_value=0)
    )

    # Sort by fraction failed (descending) for consistent ordering
    order = (
        pd.DataFrame({
            "failed": counts.get("failed", 0),
            "ambiguous": counts.get("ambiguous", 0),
        })
        .div(counts.sum(axis=1), axis=0)
        .sort_values(["failed", "ambiguous"], ascending=False)
        .index
    )
    counts_ordered = counts.loc[order, :]

    # Calculate proportions and totals
    proportions = counts_ordered.div(counts_ordered.sum(axis=1), axis=0) * 100
    total_counts = counts_ordered.sum(axis=1)

    # Create JointGrid-like figure with space for legend
    from matplotlib.gridspec import GridSpec
    # Calculate dynamic height based on number of categories with better scaling
    # Use 0.3 inches per category with a minimum of 6 inches
    fig_height = max(6, n_hues * 0.3)
    fig_width = 11
    # make the width the height minus 10% when > 10 categories
    # this ensures a square-ish plot for many categories
    if n_hues > 10:
        fig_width = fig_height - (fig_height * 0.1)
    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = GridSpec(
        nrows=1,
        ncols=3,
        width_ratios=[0.75, 0.2, 0.05],
        wspace=0.1,
        hspace=0,
        # Reserve vertical room for the figure-level title so it never overlaps bars.
        top=0.86,
        bottom=0.02
    )

    ax_main = fig.add_subplot(gs[0, 0])
    ax_margin = fig.add_subplot(gs[0, 1], sharey=ax_main)

    # Main plot: Relative abundance stacked barplot
    bottom = pd.Series(0, index=proportions.index)
    bars_per_group = {idx: [] for idx in proportions.index}

    for status in proportions.columns:
        widths = proportions[status]
        bars = ax_main.barh(
            proportions.index,
            widths,
            left=bottom,
            edgecolor='white',
            linewidth=0.5,
            label=status,
            alpha=0.9,
        )
        bottom += widths
        for idx, bar in zip(proportions.index, bars):
            bars_per_group[idx].append((bar, status))

    ax_main.set_xlabel('% Cells', fontsize="large")
    ax_main.set_ylabel(group_key, fontsize="large")
    ax_main.set_xlim([0, 100])
    ax_main.margins(y=0)
    
    # Add detailed labels for main plot
    # Ensure the figure is drawn before querying the renderer/text extents
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for idx, bar_status_pairs in bars_per_group.items():
        y_pos = (
            bar_status_pairs[0][0].get_y()
            + bar_status_pairs[0][0].get_height() / 2
        )

        for bar, status in bar_status_pairs:
            count_value = counts_ordered.loc[idx, status]
            label_text = str(int(count_value))

            temp_text = ax_main.text(
                0, 0, label_text, fontsize="medium", fontweight="bold"
            )
            text_width_px = temp_text.get_window_extent(renderer=renderer).width
            temp_text.remove()

            left_px = ax_main.transData.transform((bar.get_x(), 0))[0]
            right_px = ax_main.transData.transform(
                (bar.get_x() + bar.get_width(), 0)
            )[0]
            bar_width_px = right_px - left_px

            # Place label inside bar if it fits, otherwise outside
            if (bar_width_px >= text_width_px + 3) or (proportions.shape[0] < 10):
                x_pos = bar.get_x() + bar.get_width() / 2
                color = 'white'
                ha = 'center'
            else:
                x_pos = (
                    bar.get_x()
                    + bar.get_width()
                    + proportions.sum(1).max() / 100
                )
                color = bar.get_facecolor()
                ha = 'left'

            ax_main.text(
                x_pos,
                y_pos,
                label_text,
                va='center',
                ha=ha,
                fontsize="medium",
                fontweight='bold',
                color=color,
                rotation=270 if proportions.shape[0] < 10 else 0,
            )

    # Margin plot: stacked barplot of total counts on ax_margin
   # Margin plot: Total cell counts as histogram with annotations
    ax_margin.barh(
        total_counts.index,
        total_counts.values,
        edgecolor='white',
        linewidth=0.5,
        alpha=0.9,
        color='#555555',  # Dark gray
    )

    ax_margin.tick_params(axis='y', labelleft=False, left=False)
    ax_margin.tick_params(
        axis='x', labelbottom=True, labelsize="medium", rotation=270
    )
    ax_margin.set_xlabel('# Cells', fontsize="large", labelpad=10)
    ax_margin.locator_params(axis='x', nbins=5)
    ax_margin.ticklabel_format(axis='x', style='plain')

    # Format axes
    for pos in ['right', 'top']:
        ax_main.spines[pos].set_visible(False)
        ax_margin.spines[pos].set_visible(False)

    # Add legend at figure level in the reserved space
    handles, labels = ax_main.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        title="QC Status",
        loc="center left",
        bbox_to_anchor=(0.88, 0.5),
        frameon=False,
        fontsize="large",
    )

    fig.suptitle(
        f"Cells that passed QC\n{dataset=}\n{file_id=}",
        fontsize="large" if proportions.shape[0] < 10 else "xx-large",
        y=0.98,
    )
    fig.savefig(plot_dir / f'bar_by={group_key}.png', bbox_inches='tight', dpi=dpi)
    plt.close()
    del fig, ax_main, ax_margin


def safe_call_plot(*args, **kwargs):
    try:
        plot_composition(*args, **kwargs)
        return None
    except Exception as e:
        logging.error(f"Error in plot job: {e}")
        traceback.print_exc()
        return e


def plot_violin(
        df: pd.DataFrame,
        hue: str | list[str],
        metrics: str | list[str],
        dataset: str = "dataset",
        facet: str | None = None,
        output_plots: Path = Path("."),
        suffix: str = "",
        dpi: int = 150,
):
    if isinstance(metrics, str):
        metrics = [metrics]
    n_panels = len(metrics)
    figsize_metrics = 4 * n_panels
    if facet is None:
        fig, axes = plt.subplots(1, n_panels, figsize=(figsize_metrics, 6))
    else:
        fig_width = max(figsize_metrics * 2, df[facet].nunique() * 1)
        fig, axes = plt.subplots(n_panels, 1, figsize=(fig_width, figsize_metrics))
    axes = np.atleast_1d(axes)
    plt.grid(False)
    for i, qc_metric in enumerate(metrics):
        sns.violinplot(
            data=df.reset_index(drop=True),
            x=hue if facet is None else facet,
            y=qc_metric,
            hue=hue,
            fill=False,
            split=facet is None,
            inner='quartile',
            legend=False,
            ax=axes[i]
        )
        axes[i].set_xlabel(hue)
        axes[i].set_ylabel(qc_metric)
        axes[i].set_title(f'{qc_metric} distribution')
        if facet is not None:
            axes[i].tick_params(axis='x', labelrotation=90)
        for pos in ['right', 'top']: 
            axes[i].spines[pos].set_visible(False) 
    plt.suptitle(f'Cells that passed QC\n{dataset}', fontsize="xx-large")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    filename = f'violin_{hue}_{facet}_{len(metrics)}metrics{suffix}.png'
    plt.savefig(output_plots / filename, bbox_inches='tight', dpi=dpi)


def plot_removed(
    snakemake,
    **kwargs,
):
    ## Variables ## ------------------------------------------------------------
    for k, v in kwargs.items(): parse_set_nested(snakemake, k, v)
    pp.pprint(snakemake)
    input_zarr = snakemake.input.zarr
    output_plots = Path(snakemake.output.plots)
    output_plots.mkdir(parents=True, exist_ok=True)
    threads = snakemake.threads
    dpi = snakemake.params.get('dpi', 150)
    # get parameters
    file_id = snakemake.wildcards.file_id
    scautoqc_metrics = snakemake.params.get('scautoqc_metrics', QC_FLAGS)

    ## Loading data ## ---------------------------------------------------------
    logging.info(f'Read {input_zarr}...')
    adata = read_anndata(input_zarr, obs='obs', uns='uns')
    # If no cells filtered out, save empty plots
    if adata.obs.shape[0] == 0:
        logging.info('Empty data, skipping plots...')
        return
    dataset, hues = parse_parameters(adata, snakemake.params, filter_hues=True)

    ## Pre-processing ## -------------------------------------------------------
    scautoqc_metrics = [m for m in scautoqc_metrics if m in adata.obs.columns]

    ## Main code ## ------------------------------------------------------------
    logging.info(f'{hues=}')
    logging.info('Plot removed cells...')
    plot_bar_all(adata, dataset, output_plots, dpi)

    # Filter out None (and non-categorical) hues before creating composition plots
    hues_for_composition = [
        h for h in hues
        if h is not None and pd.api.types.is_categorical_dtype(adata.obs[h])
    ]

    # run in parallel using joblib; results is a list of return values (None or Exception)
    # Use threading backend to avoid copying data across processes (memory-efficient)
    results = Parallel(n_jobs=threads, backend='threading', pre_dispatch='2*n_jobs')(
        delayed(safe_call_plot)(
            adata.obs.copy(),
            dataset=dataset,
            file_id=file_id,
            group_key=g,
            plot_dir=output_plots,
            dpi=dpi,
        ) for g in tqdm(hues_for_composition)
    )

    gc.collect()
    plt.close('all')
    plt.rcdefaults()  # Reset matplotlib state

    # log any errors
    errors = [(g, r) for g, r in zip(hues_for_composition, results) if r is not None]
    if errors:
        logging.error(f"{len(errors)} group(s) failed during plotting.")
        for g, e in errors:
            logging.error(f"Group {g} failed: {e}")

    logging.info('Plot violin plots per QC metric...')
    logging.info(f'{scautoqc_metrics=}')
    df = adata.obs.copy()
    df = df.reset_index(drop=True) # when many cells and no deduplication
    plot_violin(
        df=df,
        hue="qc_status",
        metrics=scautoqc_metrics,
        dataset=dataset,
        output_plots=output_plots,
        dpi=dpi,
    )

    logging.info('Plot violin plots per QC metric for log scale...')
    qc_metric_log_list = []
    for i, qc_metric in enumerate(scautoqc_metrics):
        temp, metric_base = log_auto(df[qc_metric])
        if metric_base > 1:
            qc_metric_log = f"log{metric_base}_{qc_metric}"
            df[qc_metric_log] = temp
            qc_metric_log_list.append(qc_metric_log)

    if len(qc_metric_log_list) > 0:
        plot_violin(
            df=df,
            hue="qc_status",
            metrics=qc_metric_log_list,
            dataset=dataset,
            output_plots=output_plots,
            suffix="_log",
            dpi=dpi,
        )

if __name__ == "__main__" and not is_interactive():
    if "snakemake" in globals():
        args = snakemake
    else:
        args = parse_args_()
    logging.info("Starting plot_removed.")
    plot_removed(args)
    logging.info("Run completed successfully.")
