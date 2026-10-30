from pathlib import Path
import re
import warnings
import logging
from tqdm import tqdm

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from qc_utils import QC_FLAGS

plt.rcParams['svg.fonttype'] = 'none'
warnings.simplefilter("ignore", UserWarning)
logging.basicConfig(level=logging.INFO)
sns.set_theme(style='white', context='paper')


def infer_file_id(path, fallback):
    match = re.search(r'file_id~([^/]+)', str(path))
    if match:
        return match.group(1)
    return fallback


def parse_file_id_parts(file_id):
    file_id = str(file_id)

    # Keep only the part after the last ':' before the first '--'.
    prefix = file_id.split('--', 1)[0]
    file_name = prefix.rsplit(':', 1)[-1]

    if file_name == '':
        file_name = file_id

    split_data_value = ''
    for token in file_id.split('--')[1:]:
        if token.startswith('split_data_value='):
            split_data_value = token.split('=', 1)[1] or ''
            break

    return file_name, split_data_value


def safe_name(name):
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', str(name)).strip('_')


def _extract_bounds(df, study, lineage, lineage_key, thresh_type, metrics):
    """Helper to cleanly extract (min, max) threshold tuples for a framework."""
    match = df[(df["study"] == study) & (df[lineage_key] == lineage) & (df["threshold_type"] == thresh_type)]
    if match.empty:
        return {m: (None, None) for m in metrics}

    row = match.iloc[0]
    return {m: (row.get(f"{m}_min"), row.get(f"{m}_max")) for m in metrics}


def _infer_right_tail_clip_max(values, clip_quantile=0.999, tail_ratio=10.0):
    """Return clipped max for long right tails, else full max."""
    values = pd.to_numeric(values, errors='coerce')
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return None, False

    x_max = float(values.max())
    if len(values) < 20:
        return x_max, False

    q50 = float(values.quantile(0.50))
    q75 = float(values.quantile(0.75))
    q99 = float(values.quantile(0.99))
    q_clip = float(values.quantile(clip_quantile))

    body_scale = q75 - q50
    if body_scale <= 0:
        body_scale = max(float(values.quantile(0.95) - q50), 1e-12)
    tail_span = x_max - q99

    use_tail_clip = (tail_span > tail_ratio * body_scale) and (q_clip < x_max)
    return (q_clip if use_tail_clip else x_max), use_tail_clip


def ridge_plot(
    df,
    row_group,
    col_group,
    metric,
    out_file=None,
    dpi=100,
    thresholds_df=None,
    threshold_color='black',
    hspace=-0.7,
    linewidth=1.5,
    long_tail_clip_quantile=0.999,
    long_tail_ratio=3.0,
):
    THRESHOLD_TYPES = [
        ('sctk_autoqc', '--', 'sctk scAutoQC'),
        ('updated',     '-',  'Applied threshold'),
        ('alternative', ':',  'Alternative threshold'),
    ]

    # ── Grid setup ───────────────────────────────────────────────────────────
    sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

    grid = sns.FacetGrid(
        df,
        row=row_group,
        col=col_group,
        hue=row_group,
        aspect=4,
        height=0.5,
        palette=sns.cubehelix_palette(df[row_group].nunique(), rot=-.25, light=.7),
        margin_titles=False,
    )

    # ── KDE layers ───────────────────────────────────────────────────────────
    x_clip_max, use_tail_clip = _infer_right_tail_clip_max(
        df[metric],
        clip_quantile=long_tail_clip_quantile,
        tail_ratio=long_tail_ratio,
    )

    # Only clip when a long right tail is detected.
    kde_kwargs = dict(bw_adjust=1)
    if use_tail_clip and x_clip_max is not None:
        kde_kwargs['cut'] = 0
        # Keep lower side open; only constrain the right side.
        kde_kwargs['clip'] = (-np.inf, x_clip_max)

    grid.map(sns.kdeplot, metric, fill=True, alpha=1, **kde_kwargs)
    grid.map(sns.kdeplot, metric, color='white', linewidth=linewidth, **kde_kwargs)
    grid.refline(y=0, linewidth=linewidth, linestyle="-", color=None, clip_on=False)

    if use_tail_clip and x_clip_max is not None:
        xlim_min, _ = grid.axes.flat[0].get_xlim()
        if x_clip_max > xlim_min:
            grid.set(xlim=(xlim_min, x_clip_max))

    # ── Layout ───────────────────────────────────────────────────────────────
    grid.set(yticks=[], ylabel='')
    grid.set_xlabels(metric)
    grid.despine(bottom=True, left=True)
    grid.figure.subplots_adjust(left=0.18, right=0.98, top=0.98, hspace=hspace, wspace=0.12)

    # ── Column titles (workaround for seaborn bug) ────────────────────────────
    for ax in grid.axes.flat:
        ax.set_title("")
    if df[col_group].nunique() > 1:
        for ax, title in zip(grid.axes[0], grid.col_names):
            ax.set_title(title.split(" | ")[-1].split(" = ")[-1])

    # ── Row labels ───────────────────────────────────────────────────────────
    for row_idx, row_name in enumerate(grid.row_names):
        grid.axes[row_idx, 0].text(
            -0.08, 0.04, str(row_name),
            transform=grid.axes[row_idx, 0].transAxes,
            ha='right', va='bottom', fontsize=12
        )

    if thresholds_df is None or thresholds_df.empty:
        return grid

    # ── Threshold lines ───────────────────────────────────────────────────────
    study_first = (row_group == 'file_name' and col_group == 'split_data_value')
    metric_facets = (col_group == 'metric_name' and 'metric_name' in df.columns)
    fixed_study = None
    fixed_lineage = None
    if metric_facets:
        if row_group == 'file_name' and 'split_data_value' in df.columns:
            lineage_values = [str(v) for v in df['split_data_value'].dropna().unique()]
            if len(lineage_values) == 1:
                fixed_lineage = lineage_values[0]
        elif row_group == 'split_data_value' and 'file_name' in df.columns:
            study_values = [str(v) for v in df['file_name'].dropna().unique()]
            if len(study_values) == 1:
                fixed_study = study_values[0]

    def _metric_parts(metric_name):
        metric_name = str(metric_name)
        if metric_name.startswith('log1p_'):
            return metric_name.split('_', 1)[1], True
        return metric_name, False

    seen_types = set()

    # draw threshold lines in a separate overlay axis to ensure they are visible above the KDE layers
    overlay = grid.figure.add_axes([0, 0, 1, 1], facecolor='none', zorder=10)
    overlay.set_xlim(0, 1)
    overlay.set_ylim(0, 1)
    overlay.axis('off')

    to_display = grid.figure.transFigure.inverted().transform

    def to_fig(ax, x_data, y_axes):
        x_fig, _ = to_display(ax.transData.transform((x_data, 0)))
        _, y_fig = to_display(ax.transAxes.transform((0, y_axes)))
        return x_fig, y_fig

    def _same_value(x, y):
        if pd.isna(x) and pd.isna(y):
            return True
        if pd.isna(x) or pd.isna(y):
            return False
        return np.isclose(float(x), float(y), rtol=0, atol=1e-12)

    def draw_threshold(
        ax,
        study_name,
        lineage_name,
        base_metric,
        log_transform,
        thresh_type,
        linestyle,
        skip_bounds=None,
    ):
        bounds = _extract_bounds(thresholds_df, study_name, lineage_name, 'lineage', thresh_type, [base_metric])
        vmin, vmax = bounds.get(base_metric, (None, None))
        drawn = False
        for bound_name, v in (('min', vmin), ('max', vmax)):
            if skip_bounds is not None and bound_name in skip_bounds:
                continue
            if pd.notna(v):
                v = np.log1p(float(v)) if log_transform else float(v)
                xlim = ax.get_xlim()
                # Only draw threshold if it falls within the current plot x-limits
                if v < xlim[0] or v > xlim[1]:
                    continue
                x, y_bot = to_fig(ax, v, 0)
                _, y_top = to_fig(ax, v, 1 - abs(hspace))
                overlay.vlines(x, y_bot, y_top, colors=threshold_color, linestyles=linestyle, linewidth=1)
                drawn = True
        return drawn

    for i, r_name in enumerate(grid.row_names):
        for j, c_name in enumerate(grid.col_names):
            if metric_facets:
                base_metric, log_transform = _metric_parts(c_name)
                if row_group == 'file_name' and fixed_lineage is not None:
                    study, lineage = str(r_name), fixed_lineage
                elif row_group == 'split_data_value' and fixed_study is not None:
                    study, lineage = fixed_study, str(r_name)
                else:
                    continue
            else:
                base_metric, log_transform = _metric_parts(metric)
                study, lineage = (str(r_name), str(c_name)) if study_first else (str(c_name), str(r_name))

            for thresh_type, linestyle, _ in THRESHOLD_TYPES:
                skip_bounds = None
                if thresh_type == 'updated':
                    updated_bounds = _extract_bounds(thresholds_df, study, lineage, 'lineage', 'updated', [base_metric])
                    autoqc_bounds = _extract_bounds(thresholds_df, study, lineage, 'lineage', 'sctk_autoqc', [base_metric])
                    updated_min, updated_max = updated_bounds.get(base_metric, (None, None))
                    autoqc_min, autoqc_max = autoqc_bounds.get(base_metric, (None, None))
                    skip_bounds = {
                        bound_name
                        for bound_name, updated_value, autoqc_value in (
                            ('min', updated_min, autoqc_min),
                            ('max', updated_max, autoqc_max),
                        )
                        if _same_value(updated_value, autoqc_value)
                    }
                    if len(skip_bounds) == 2:
                        continue
                if draw_threshold(
                    grid.axes[i, j],
                    study,
                    lineage,
                    base_metric,
                    log_transform,
                    thresh_type,
                    linestyle,
                    skip_bounds=skip_bounds,
                ):
                    seen_types.add(thresh_type)

    # ── Threshold legend ──────────────────────────────────────────────────────
    legend_handles = [
        Line2D([], [], color=threshold_color, linestyle=ls, lw=linewidth, label=label)
        for thresh_type, ls, label in THRESHOLD_TYPES
        if thresh_type in seen_types
    ]
    if legend_handles:
        grid.figure.legend(
            handles=legend_handles,
            labels=[h.get_label() for h in legend_handles],
            loc='center left',
            bbox_to_anchor=(1.02, 0.5),
            bbox_transform=grid.figure.transFigure,
            frameon=False,
            fontsize=10,
        )
    
    if out_file:
        grid.savefig(out_file, dpi=dpi, bbox_inches='tight')
        plt.close(grid.figure)


def plot_heatmap(data_df, title, filename, out_dir, annot=False, dpi=200):
    import textwrap
    import matplotlib.cm as cm

    plot_df = data_df.copy()
    if plot_df.shape[0] < plot_df.shape[1]:
        plot_df = plot_df.T
    
    n_rows, n_cols = plot_df.shape
    base_cell_size = 0.4
    label_scale = 0.08
    min_dim = 2.0
    col_label_pad = max(str(v) for v in plot_df.columns)
    row_label_pad = max((str(v) for v in plot_df.index), key=len)
    row_label_width = len(row_label_pad) * label_scale
    col_label_height = len(col_label_pad) * label_scale
    cell_size = max(
        base_cell_size,
        max(0, min_dim - row_label_width) / max(n_cols, 1),
        max(0, min_dim - col_label_height) / max(n_rows, 1),
    )
    fig, (ax, cb_ax) = plt.subplots(
        1, 2,
        figsize=(
            row_label_width + n_cols * cell_size,
            col_label_height + n_rows * cell_size
        ),
        gridspec_kw={'width_ratios': [1, 0.05]},
    )
    sns.heatmap(
        plot_df,
        cmap='Blues',
        annot=annot,
        fmt='.2f' if annot else '',
        ax=ax,
        square=True,
        linewidths=0.5,
        linecolor=cm.Blues(0.3),
        cbar_ax=cb_ax,
    )
    
    # adjust ticks
    ax.xaxis.tick_top()
    ax.xaxis.set_label_position('top')
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.tick_params(axis='x', rotation=90, length=0)
    ax.tick_params(axis='y', rotation=0)

    # adjust legend
    colorbar = ax.collections[0].colorbar
    colorbar.ax.yaxis.set_label_position('right')
    cb_pos = colorbar.ax.get_position()
    colorbar.ax.set_position([
        cb_pos.x0,
        cb_pos.y0 + cb_pos.height * 0.25,
        cb_pos.width,
        cb_pos.height * 0.5,
    ])
    
    # adjust colorbar label to wrap based on its width
    fig.canvas.draw()
    cbar_height = colorbar.ax.get_window_extent().height
    chars_per_line = max(1, int(cbar_height / fig.dpi / 0.09))
    wrapped = "\n".join(textwrap.wrap(title, chars_per_line))
    n_lines = wrapped.count('\n') + 1
    colorbar.set_label("\n".join(
        textwrap.wrap(title, chars_per_line)),
        rotation=270,
        labelpad=n_lines * 10,
        fontsize=10,
    )
    
    out_file = out_dir / filename
    fig.savefig(out_file, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    return out_file.name


# removed wrapper: heatmap creation inlined where needed below


input_files = list(snakemake.input.qc_metrics)
threshold_files = list(snakemake.input.tsv)
output_dir = Path(snakemake.output.plots)
output_dir.mkdir(parents=True, exist_ok=True)
dpi = snakemake.params.get('dpi', 150)
qc_metric_columns = snakemake.params.get('scautoqc_metrics', QC_FLAGS)
qc_metric_columns += [
    f'log1p_{col}' for col in qc_metric_columns
    if col in ['n_counts', 'n_genes']
]
threshold_color = snakemake.params.get('threshold_color', '#E66101')
# Hardcoded linestyles for thresholds
threshold_linestyle = '-'
threshold_linestyle2 = '--'
output_ridge = output_dir / 'ridge_plots'
output_ridge.mkdir(parents=True, exist_ok=True)
written_files = []

## Read data

frames = []
if len(input_files) == 0:
    logging.info('No qc_metrics input files found. Skipping ridge plots.')

for idx, path in tqdm(
    enumerate(input_files),
    total=len(input_files),
    desc='Reading qc_metrics files',
):
    try:
        df = pd.read_parquet(path)
    except Exception as exc:
        logging.warning(f'Could not read qc metrics parquet: {path} ({exc})')
        continue

    if df.shape[0] == 0:
        continue
    file_id = infer_file_id(path, fallback=f'file_{idx + 1}')
    file_name, split_data_value = parse_file_id_parts(file_id)
    df = df.copy()
    df['file_id'] = file_id
    df['file_name'] = file_name
    df['split_data_value'] = split_data_value
    frames.append(df)

qc_df = pd.concat(frames, ignore_index=True) if len(frames) > 0 else None
if qc_df is None:
    logging.info('All qc_metrics files are empty or unreadable. Skipping ridge plots.')

threshold_frames = []
for path in threshold_files:
    try:
        thresholds = pd.read_table(path)
    except Exception as exc:
        logging.warning(f'Could not read threshold TSV: {path} ({exc})')
        continue

    if thresholds.shape[0] == 0:
        continue
    thresholds = thresholds.copy()
    thresholds['file_id'] = thresholds['file_id'].astype(str)
    parsed = thresholds['file_id'].apply(parse_file_id_parts)
    thresholds['study'] = parsed.apply(lambda x: x[0])
    thresholds['lineage'] = parsed.apply(lambda x: x[1])
    threshold_frames.append(thresholds)

# If threshold files were provided, build thresholds DataFrame now
if len(threshold_frames) > 0:
    thresholds_df = pd.concat(threshold_frames, ignore_index=True)
else:
    thresholds_df = None

## Removed-cell heatmaps

logging.info('Creating removed-cell heatmaps...')
if thresholds_df is not None:
    thresholds = thresholds_df.query('threshold_type == "updated"')
    if thresholds.shape[0] == 0:
        logging.info('No updated thresholds found. Skipping removed-cell heatmaps.')
    else:
        thresholds = thresholds.drop_duplicates(subset=['study', 'lineage'])
        removed_df = thresholds.pivot(index=['study'], columns='lineage', values='n_removed')
        removed_frac_df = thresholds.pivot(index=['study'], columns='lineage', values='removed_frac')

        plot_heatmap(
            removed_df,
            'Cells removed',
            'n_removed_heatmap.svg',
            output_dir,
            annot=False,
            dpi=dpi
        )
        
        plot_heatmap(
            removed_frac_df,
            'Fraction of cells removed',
            'removed_frac_heatmap.svg',
            output_dir,
            annot=True,
            dpi=dpi
        )
        logging.info(f'Created removed-cell heatmaps in {output_dir}')
else:
    logging.info('No thresholds TSV rows found, skipping removed-cell heatmaps.')


## Ridge plots for numeric QC metrics

if qc_df is not None:
    metric_columns = [
        col for col in qc_metric_columns
        if col in qc_df.columns and pd.api.types.is_numeric_dtype(qc_df[col])
    ]

    if len(metric_columns) == 0:
        logging.info('No numeric QC metric columns found. Skipping ridge plots.')
    else:
        logging.info(f'Found {len(metric_columns)} numeric QC metric columns for plotting: {metric_columns}')
        n_file_names = qc_df['file_name'].nunique()
        n_split_values = qc_df['split_data_value'].nunique()
        if n_file_names < n_split_values:
            row_group = 'split_data_value'
            col_group = 'file_name'
        else:
            row_group = 'file_name'
            col_group = 'split_data_value'

        # If there would only be a single facet column, generate one combined plot
        # with one column per metric instead of one file per metric.
        if qc_df[col_group].nunique() == 1 and len(metric_columns) > 1:
            combined_frames = []
            for metric in metric_columns:
                tmp = qc_df[[row_group, 'file_name', 'split_data_value', metric]].copy()
                tmp[metric] = pd.to_numeric(tmp[metric], errors='coerce')
                tmp = tmp[np.isfinite(tmp[metric])]
                if tmp.shape[0] < 2 or tmp[metric].nunique() < 2:
                    continue
                tmp = tmp.rename(columns={metric: 'metric_value'})
                tmp['metric_name'] = metric
                combined_frames.append(tmp)

            if len(combined_frames) == 0:
                logging.info('No plottable metric values found for combined ridge plot.')
            else:
                combined_df = pd.concat(combined_frames, ignore_index=True)
                out_file = output_ridge / 'all_metrics_ridge.svg'
                ridge_plot(
                    combined_df,
                    row_group=row_group,
                    col_group='metric_name',
                    metric='metric_value',
                    out_file=out_file,
                    dpi=dpi,
                    thresholds_df=thresholds_df,
                    threshold_color=threshold_color,
                )
                written_files.append(out_file.name)
                logging.info('Created combined ridge plot across metrics with threshold overlays.')
        else:
            for metric in tqdm(metric_columns, desc='Creating ridge plots'):
                metric_df = qc_df[['file_id', 'file_name', 'split_data_value', metric]].copy()
                metric_df[metric] = pd.to_numeric(metric_df[metric], errors='coerce')
                metric_df = metric_df[np.isfinite(metric_df[metric])]

                if metric_df.shape[0] < 2 or metric_df['file_name'].nunique() == 0:
                    continue
                if metric_df[metric].nunique() < 2:
                    continue

                out_file = output_ridge / f'{safe_name(metric)}.svg'
                ridge_plot(
                    metric_df,
                    row_group=row_group,
                    col_group=col_group,
                    metric=metric,
                    out_file=out_file,
                    dpi=dpi,
                    thresholds_df=thresholds_df,
                    threshold_color=threshold_color,
                )
                written_files.append(out_file.name)

if len(written_files) > 0:
    logging.info(f'Created {len(written_files)} ridge plots in {output_ridge}')
elif qc_df is not None:
    logging.info('No ridge plots were created after metric filtering.')
