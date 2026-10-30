import logging
logging.basicConfig(level=logging.INFO)
import textwrap
import pandas as pd

from dataclasses import dataclass
from typing import Literal
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
from matplotlib.gridspec import GridSpec
import seaborn as sns

plt.rcParams["svg.fonttype"] = "none"


# ============================================================================
# Helper
# ============================================================================

def round_values(x, prefix='', n_digits=3):
    """Format numeric values for plot annotations with fixed or scientific notation."""
    if x < 10 ** (-n_digits):
        return f'{prefix}{x:.2e}'
    elif pd.notna(x):
        return f'{prefix}{x:.{n_digits}f}'
    return ''


# ============================================================================
# PCR Plotter
# ============================================================================

PlotKind = Literal['barplot', 'violin']


@dataclass
class PCRPlotConfig:
    """Configuration for PCR score plots."""
    # Shared sizing for all plot kinds: (height_base, height_per_covariate, width_scale)
    plot_sizing: tuple[float, float, float] = 4.0, 0.3, 1.0
    min_width: float = 4.0
    panel_fraction: float = 0.14
    min_gap_inches: float = 0.7
    label_char_width_inches: float = 0.08
    label_gap_padding_inches: float = 0.2
    main_xlabel: str = 'Principal component regression score'
    left_margin: float = 0.04
    right_margin: float = 0.02
    bottom_margin: float = 0.11
    top_margin: float = 0.93
    dpi: int = 150
    # Font sizes
    title_fontsize: int = 16
    legend_title_fontsize: int = 11
    legend_fontsize: int = 11
    tick_fontsize: int = 11
    label_fontsize: int = 12
    bar_label_fontsize: int = 11
    perm_text_fontsize: int = 11


class PCRPlotter:
    """Plots PCR scores as a barplot and violin plot with a left-side n_covariates panel."""

    def __init__(self, df: pd.DataFrame, config: PCRPlotConfig | None = None):
        """Prepare sorted plotting data and precomputed label strings."""
        self.config = config or PCRPlotConfig()

        df = df.sort_values(['pcr', 'n_covariates', 'covariate'], ascending=[False, True, False])
        self._df = df
        self._covariate_order = list(df['covariate'].drop_duplicates())

        stats = df.groupby('covariate', sort=False).first().reindex(self._covariate_order).copy()
        stats['pcr_string'] = stats['pcr'].apply(round_values, prefix='pcr=')
        stats['signif'] = stats['p-val'].apply(
            lambda x: '**' if x <= 0.01 else '*' if x <= 0.05 else ''
        )
        stats['z_score'] = stats['z_score'].apply(round_values, prefix='z=', n_digits=2)
        stats['p-val'] = stats['p-val'].apply(round_values, prefix='p-val=', n_digits=3)
        self._covariate_stats = stats

        self._covariate_bar_labels = stats[
            ['z_score', 'signif']
        ].astype(str).agg(lambda x: ', '.join([s for s in x if s]), axis=1)

        self._n_covariates_values = stats.loc[self._covariate_order, 'n_covariates'].values
        self._covariate_is_categorical = self._resolve_categorical_flags(stats)
        self._num_covariates = len(self._covariate_order)

    def _resolve_categorical_flags(self, stats: pd.DataFrame) -> pd.Series:
        """Resolve per-covariate categorical status from optional metadata columns."""
        if 'covariate_type' in stats.columns:
            return stats['covariate_type'].astype(str).str.strip().str.lower().eq('categorical')
        # Backward-compatible default: treat all covariates as categorical.
        return pd.Series(True, index=stats.index, dtype=bool)

    def _compute_wrapped_labels(self, panel_width_inches: float) -> list[str]:
        """Wrap covariate labels to fit within the label panel width and row height.

        Derives the maximum character wrap width from available horizontal space,
        then tightens it until all labels fit within the vertical row height.
        """
        cfg = self.config
        _, height_per_cov, _ = cfg.plot_sizing

        # Horizontal limit: how many chars fit in the available panel width.
        available_inches = max(panel_width_inches - cfg.label_gap_padding_inches, 0.5)
        max_chars_by_width = max(6, int(available_inches / cfg.label_char_width_inches))

        # Vertical limit: how many lines fit in one covariate row.
        char_height_inches = cfg.label_char_width_inches * 1.6
        max_lines = max(1, int(height_per_cov / char_height_inches))

        # Tighten wrap until all labels respect the row height.
        for wrap_at in range(max_chars_by_width, 5, -1):
            wrapped = [textwrap.fill(str(c), wrap_at) for c in self._covariate_order]
            if all(len(w.split('\n')) <= max_lines for w in wrapped):
                return wrapped

        # Fallback: wrap at 6 chars regardless of overflow.
        return [textwrap.fill(str(c), 6) for c in self._covariate_order]

    def _estimate_label_panel_width_inches(self) -> float:
        """Estimate required label panel width from a preliminary label wrap.

        Uses the full min_gap_inches as an initial panel width estimate to seed
        the wrap, then measures the longest resulting line.
        """
        cfg = self.config
        wrapped = self._compute_wrapped_labels(cfg.min_gap_inches * 4)
        longest_len = max(
            (max(len(line) for line in label.split('\n')) for label in wrapped),
            default=0,
        )
        return max(
            cfg.min_gap_inches,
            longest_len * cfg.label_char_width_inches + cfg.label_gap_padding_inches,
        )

    def _compute_column_width_ratios(self, figure_width_inches: float) -> tuple[float, float, float]:
        """Compute left-label-main width ratios for GridSpec."""
        cfg = self.config
        left_panel_ratio = max(cfg.panel_fraction, 0.01)
        label_panel_ratio = max(self._estimate_label_panel_width_inches() / figure_width_inches, 0.01)
        main_ratio = max(1.0 - left_panel_ratio - label_panel_ratio, 0.2)
        return left_panel_ratio, label_panel_ratio, main_ratio

    def _style_main_axis(self, ax_main: Axes, title: str) -> None:
        """Apply shared styling for the main plotting axis."""
        cfg = self.config
        ax_main.set_title(title, fontsize=cfg.title_fontsize)
        ax_main.set_ylabel('', fontsize=cfg.label_fontsize)
        ax_main.set_xlabel(cfg.main_xlabel, fontsize=cfg.label_fontsize)
        ax_main.tick_params(axis='x', labelrotation=90, labelsize=cfg.tick_fontsize)
        ax_main.tick_params(axis='y', labelsize=cfg.tick_fontsize, labelleft=False)
        sns.despine(ax=ax_main, top=True, right=True, left=False, bottom=False)

    def _annotate_permuted_std(self, ax_main: Axes) -> None:
        """Annotate permuted bars with precomputed standard deviation labels."""
        if len(ax_main.containers) < 2:
            return
        perm_std_values = pd.to_numeric(self._covariate_stats['perm_std'], errors='coerce').fillna(0.0)
        perm_labels = perm_std_values.apply(round_values, prefix='std=')
        x_min, x_max = ax_main.get_xlim()
        x_pad = (x_max - x_min) * 0.02
        for rect, label, std in zip(ax_main.containers[1], perm_labels, perm_std_values):
            x_tip = rect.get_x() + rect.get_width()
            y_mid = rect.get_y() + rect.get_height() / 2
            ax_main.text(
                x_tip + max(float(std), 0.0) + x_pad,
                y_mid,
                label,
                va='center',
                ha='left',
                fontsize=self.config.perm_text_fontsize,
                clip_on=False,
            )

    def _add_right_summary_twin_axis(self, ax_main: Axes) -> None:
        """Attach right-side summary labels aligned to y-ticks."""
        cfg = self.config
        ax_right = ax_main.twinx()
        ax_right.set_ylim(ax_main.get_ylim())
        ax_right.set_yticks(ax_main.get_yticks())
        ax_right.set_yticklabels(
            [self._covariate_bar_labels.get(c, '') for c in self._covariate_order],
            fontsize=cfg.label_fontsize,
        )
        ax_right.tick_params(left=False, right=False, length=0, labelsize=cfg.tick_fontsize)
        ax_right.spines['top'].set_visible(False)
        ax_right.spines['right'].set_visible(False)

    def _draw_left_panel(self, ax: Axes, y_lim: tuple) -> None:
        """Draw the left marginal panel with inverted horizontal bar labels."""
        cfg = self.config
        values = pd.Series(self._n_covariates_values, dtype=float)
        is_categorical = self._covariate_is_categorical.reindex(self._covariate_order).fillna(True).to_numpy()
        y_positions = list(range(len(values)))

        categorical_values = [v for v, is_cat in zip(values.values, is_categorical) if is_cat]
        max_value = max(max(categorical_values, default=0.0), 1.0)

        if categorical_values:
            categorical_y = [y for y, is_cat in zip(y_positions, is_categorical) if is_cat]
            ax.barh(categorical_y, categorical_values, color='black', alpha=0.75, height=0.8)

        ax.set_xlim(max_value * 1.2, 0)
        ax.set_ylim(y_lim)
        ax.set_xlabel('# covariates', fontsize=cfg.label_fontsize)
        ax.set_yticks([])

        for y, v, is_cat in zip(y_positions, values.values, is_categorical):
            label = str(int(v)) if float(v).is_integer() else f'{v:.2f}'
            x_pos = v + max_value * 0.03 if is_cat else max_value * 0.06
            ax.text(x_pos, y, label, va='center', ha='right', color='black', clip_on=False, fontsize=cfg.label_fontsize)

        sns.despine(ax=ax, top=True, right=False, left=True, bottom=False)

    def _draw_covariate_label_panel(self, ax: Axes, y_lim: tuple, panel_width_inches: float) -> None:
        """Draw wrapped covariate labels in a dedicated middle panel."""
        wrapped = self._compute_wrapped_labels(panel_width_inches)

        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(y_lim)
        ax.set_xticks([])
        ax.set_yticks([])

        for y, label in enumerate(wrapped):
            ax.text(1.0, y, label, va='center', ha='right', fontsize=self.config.tick_fontsize)

        for spine in ax.spines.values():
            spine.set_visible(False)

    def _create_barplot(self, title: str, ax_left: Axes, ax_labels: Axes, ax_main: Axes, label_panel_width_inches: float) -> None:
        """Draw the stratified PCR barplot on the main axis."""
        sns.barplot(
            data=self._df,
            x='pcr', y='covariate', hue='permuted',
            order=self._covariate_order,
            errorbar='sd', dodge=True,
            err_kws={'linewidth': 1}, capsize=0.1,
            ax=ax_main,
        )
        cfg = self.config
        _, x_max = ax_main.get_xlim()
        ax_main.set_xlim(0, x_max * 1.25)
        logging.info(self._covariate_bar_labels)
        ax_main.bar_label(ax_main.containers[0], labels=self._covariate_bar_labels, padding=10, fontsize=cfg.bar_label_fontsize)
        legend = ax_main.get_legend()
        if legend is not None:
            try:
                legend.get_title().set_fontsize(cfg.legend_title_fontsize)
            except Exception:
                pass
            for text in legend.get_texts():
                text.set_fontsize(cfg.legend_fontsize)
        self._annotate_permuted_std(ax_main)

        self._draw_left_panel(ax_left, ax_main.get_ylim())
        self._draw_covariate_label_panel(ax_labels, ax_main.get_ylim(), label_panel_width_inches)
        self._style_main_axis(ax_main, title)

    def _create_violin(self, title: str, ax_left: Axes, ax_labels: Axes, ax_main: Axes, label_panel_width_inches: float) -> None:
        """Draw violin and strip overlays for permuted and observed PCR values."""
        sns.boxenplot(
            data=self._df[self._df['permuted']],
            x='pcr', y='covariate',
            order=self._covariate_order,
            fill=False, dodge=False,
            ax=ax_main,
        )
        sns.stripplot(
            data=self._df[~self._df['permuted']],
            x='pcr', y='covariate',
            order=self._covariate_order,
            color='red', size=8, marker='o', dodge=False,
            ax=ax_main,
        )

        cfg = self.config
        legend = ax_main.get_legend()
        if legend is not None:
            try:
                legend.get_title().set_fontsize(cfg.legend_title_fontsize + 2)
            except Exception:
                pass
            for text in legend.get_texts():
                text.set_fontsize(cfg.legend_fontsize + 2)

        self._add_right_summary_twin_axis(ax_main)
        self._draw_left_panel(ax_left, ax_main.get_ylim())
        self._draw_covariate_label_panel(ax_labels, ax_main.get_ylim(), label_panel_width_inches)
        self._style_main_axis(ax_main, title)

    def plot(self, title: str, output_path: str, kind: PlotKind) -> None:
        """Create and save a PCR plot of the requested kind."""
        cfg = self.config
        height_base, height_per_cov, width_scale = cfg.plot_sizing
        height = height_base + self._num_covariates * height_per_cov
        width = max(cfg.min_width, height * width_scale)

        fig = plt.figure(figsize=(width, height), dpi=cfg.dpi)
        fig.subplots_adjust(
            left=cfg.left_margin,
            right=1.0 - cfg.right_margin,
            bottom=cfg.bottom_margin,
            top=cfg.top_margin,
        )

        left_ratio, label_ratio, main_ratio = self._compute_column_width_ratios(width)
        label_panel_width_inches = width * label_ratio

        gs = GridSpec(
            nrows=1,
            ncols=3,
            figure=fig,
            width_ratios=[left_ratio, label_ratio, main_ratio],
            wspace=0.02,
        )

        ax_left = fig.add_subplot(gs[0, 0])
        ax_labels = fig.add_subplot(gs[0, 1], sharey=ax_left)
        ax_main = fig.add_subplot(gs[0, 2], sharey=ax_left)

        logging.info(f'{kind}...')
        if kind == 'barplot':
            self._create_barplot(title, ax_left, ax_labels, ax_main, label_panel_width_inches)
        else:
            self._create_violin(title, ax_left, ax_labels, ax_main, label_panel_width_inches)

        logging.info(f'Save to {output_path}...')
        fig.savefig(output_path, bbox_inches='tight', dpi=cfg.dpi)
        plt.close(fig)


# ============================================================================
# Snakemake script
# ============================================================================

input_file = snakemake.input.tsv
output_bar = snakemake.output.barplot
output_violin = snakemake.output.violinplot
dataset = snakemake.wildcards.dataset
file_id = snakemake.wildcards.file_id
n_permute = snakemake.params.n_permute

logging.info('Read TSV...')
df = pd.read_table(input_file)

if df.shape[0] == 0:
    logging.info('Empty TSV, skip plotting')
    for output in [output_bar, output_violin]:
        fig, ax = plt.subplots()
        ax.axis('off')
        fig.savefig(output, bbox_inches='tight')
        plt.close(fig)
    exit(0)

title = 'Linear variance contribution of covariates'
title += f'\ndataset={dataset}\nfile_id={file_id}\n{n_permute} permutations'

PCRPlotter(df).plot(title, output_bar, kind='barplot')
PCRPlotter(df).plot(title, output_violin, kind='violin')