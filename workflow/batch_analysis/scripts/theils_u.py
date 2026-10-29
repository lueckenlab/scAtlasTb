import logging
logging.basicConfig(level=logging.INFO)
from pprint import pformat
from dataclasses import dataclass, field
from typing import Optional, List, Callable, Any
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from copy import copy
from pandas.api.types import is_numeric_dtype
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

plt.rcParams["svg.fonttype"] = "none"

"""Compute and visualize Theil's U association matrices for AnnData obs covariates.

This script performs three high-level steps:
1. Validate and preprocess selected covariates from `obs`.
2. Compute pairwise Theil's U on sample-level aggregated covariate values.
3. Render either a clustered heatmap (with optional row-side bar panel) or a
    standard heatmap fallback.

The script is designed to run as a Snakemake rule script.
"""

# ============================================================================
# Theil's U Computation
# ============================================================================

class TheilsUAnalyzer:
    """Build a Theil's U matrix from observation-level metadata.

    The class encapsulates covariate filtering, sample-level aggregation, and
    pairwise Theil's U computation. The output matrix is square with both index
    and columns corresponding to the retained covariates.
    """
    
    def __init__(
        self,
        obs: pd.DataFrame,
        sample_key: str,
        covariates: List[str],
        max_unique_continuous: int = 10
    ):
        """
        Initialize Theil's U analyzer
        
        Args:
            obs: Observations dataframe
            sample_key: Column name(s) to use as sample identifier
            covariates: List of covariate column names to analyze
            max_unique_continuous: Skip numeric covariates with more unique values
        """
        self.obs = obs
        self.sample_key = sample_key
        self.covariates = covariates
        self.max_unique_continuous = max_unique_continuous
        
        self._valid_covariates: List[str] = []
        self._aggregated_df: Optional[pd.DataFrame] = None
        self._theils_u_matrix: Optional[pd.DataFrame] = None
    
    def prepare_sample_key(self) -> str:
        """Resolve the sample identifier column used for aggregation.

        Returns:
            Name of the sample key column to use downstream. If multiple sample
            keys are provided, a composite key is created. If missing/None, the
            observation index is used.
        """
        if not self.sample_key or self.sample_key == 'None':
            logging.info('Using index as sample key...')
            sample_key = 'index'
            self.obs[sample_key] = self.obs.index
            return sample_key
        
        sample_keys = [x.strip() for x in self.sample_key.split(',')]
        if len(sample_keys) > 1:
            sample_key = '--'.join(sample_keys)
            self.obs[sample_key] = (
                self.obs[sample_keys]
                .astype(str)
                .agg('-'.join, axis=1)
                .astype('category')
            )
            return sample_key
        
        return self.sample_key
    
    def filter_covariates(self) -> List[str]:
        """Filter candidate covariates to those suitable for Theil's U.

        Filtering criteria:
        - column must exist in `obs`
        - at least two non-null unique values
        - numeric columns with too many unique values are excluded

        Returns:
            List of retained covariate names.
        """
        logging.info(f'Filter covariates:\n{pformat(self.covariates)}')
        valid_covariates = []
        
        for covariate in self.covariates:
            if covariate not in self.obs.columns:
                logging.info(f'Skipping covariate not in obs: {covariate}')
                continue
            
            if self.obs[covariate].dropna().nunique() < 2:
                logging.info(f'Skipping covariate with fewer than 2 non-NA unique values: {covariate}')
                continue
            
            if (is_numeric_dtype(self.obs[covariate]) and 
                self.obs[covariate].nunique() > self.max_unique_continuous):
                logging.info(
                    f'Skipping numeric covariate with > {self.max_unique_continuous} '
                    f'unique values: {covariate}'
                )
                continue
            
            valid_covariates.append(covariate)
        
        self._valid_covariates = valid_covariates
        logging.info(f'Valid covariates:\n{pformat(valid_covariates)}')
        return valid_covariates
    
    def aggregate_by_sample(self, sample_key: str) -> pd.DataFrame:
        """Aggregate covariates per sample using mode.

        Args:
            sample_key: Name of the sample grouping column.

        Returns:
            Aggregated dataframe with one row per sample and one column per
            retained covariate.

        Raises:
            AssertionError: If no covariates remain after filtering.
        """
        logging.info(f'Aggregate covariates by {sample_key}...')
        
        # Remove sample_key from covariates if present
        covariates = [c for c in self._valid_covariates if c != sample_key]
        assert len(covariates) > 0, 'No valid covariates found in obs'
        
        self._aggregated_df = self.obs[[sample_key] + covariates].groupby(
            sample_key, observed=True
        ).agg(
            lambda x: x.mode().iloc[0] if not x.mode().empty else np.nan
        )
        
        # Drop covariates that become constant after aggregation
        constant_cols = [
            col for col in self._aggregated_df.columns
            if self._aggregated_df[col].dropna().nunique() < 2
        ]
        if constant_cols:
            logging.info(f'Dropping covariates constant after aggregation: {constant_cols}')
            self._aggregated_df = self._aggregated_df.drop(columns=constant_cols)
            self._valid_covariates = [c for c in self._valid_covariates if c not in constant_cols]
        
        return self._aggregated_df
    
    def compute_theils_u(self) -> pd.DataFrame:
        """Compute the pairwise Theil's U matrix for aggregated covariates.

        Returns:
            Square dataframe where each entry `(x, y)` is `U(x|y)`.
        """
        logging.info(f'Compute Theil\'s U matrix...')
        
        df = self._aggregated_df
        cols = df.columns
        n = len(cols)
        
        # Precompute entropies
        entropies = {col: self._entropy(df[col]) for col in cols}
        u = np.zeros((n, n), dtype=float)
        
        for i, x in enumerate(cols):
            h_x = entropies[x]
            for j, y in enumerate(cols):
                if i == j:
                    u[i, j] = 1.0
                else:
                    u[i, j] = (h_x - self._conditional_entropy(df[x], df[y])) / h_x
        
        self._theils_u_matrix = pd.DataFrame(u, index=cols, columns=cols)
        return self._theils_u_matrix
    
    @staticmethod
    def _entropy(x: pd.Series) -> float:
        """Compute Shannon entropy $H(X)$ for a categorical series."""
        p_x = x.value_counts(normalize=True)
        p_x = p_x[p_x > 0]  # filter zeros
        return -(p_x * np.log2(p_x)).sum()
    
    @staticmethod
    def _conditional_entropy(x: pd.Series, y: pd.Series) -> float:
        """Compute conditional entropy $H(X \mid Y)$ for categorical series."""
        p_xy = pd.crosstab(x, y, normalize=True)
        p_x_given_y = p_xy.div(p_xy.sum(axis=0), axis=1)
        probs = p_x_given_y.values
        log_p = np.zeros_like(probs, dtype=float)
        np.log2(probs, out=log_p, where=probs > 0)
        return -np.nansum(p_xy.values * log_p)
    
    def run(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Run the full analysis pipeline.

        Returns:
            Tuple of:
            - Theil's U matrix dataframe
            - sample-level aggregated covariate dataframe
        """
        sample_key = self.prepare_sample_key()
        self.filter_covariates()
        self.aggregate_by_sample(sample_key)
        theils_u_matrix = self.compute_theils_u()
        
        return theils_u_matrix, self._aggregated_df


# ============================================================================
# Generalized Clustermap Plotter
# ============================================================================

@dataclass
class ClustermapConfig:
    """Configuration container for heatmap/clustermap rendering.

    Attributes group naturally into:
    - sizing/typography (`min_fig_size`, `size_per_item`, font settings)
    - color/value mapping (`cmap`, `vmin`, `vmax`, `nan_color`)
    - panel/layout controls (margins, gaps, dendrogram and bar panel scaling)
    - seaborn kwargs passthrough (`*_kwargs` dictionaries)
    """
    min_fig_size: float = 6.0
    size_per_item: float = 0.6
    min_fontsize: int = 9
    max_fontsize: int = 12
    fontsize_multiplier: float = 12.0
    bar_width: float = 0.1
    bar_gap: float = 0.02
    dpi: int = 300
    cmap: str = "coolwarm"
    vmin: float = 0
    vmax: float = 1
    dendrogram_ratio: float = 0.1
    show_dendrogram: bool = True
    show_barplot: bool = True
    nan_color: str = "#d9d9d9"
    layout_left: float = 0.06
    layout_right: float = 0.98
    layout_bottom: float = 0.10
    layout_top: float = 0.95
    layout_gap: float = 0.012
    title_y: float = 0.975
    title_gap: float = 0.015
    cbar_width: float = 0.018
    cbar_gap: float = 0.02
    cbar_height_scale: float = 0.55
    outer_right_margin: float = 0.015
    bar_panel_width_scale: float = 1.5
    dendrogram_height_scale: float = 0.7
    dendrogram_min_height: float = 0.04
    dendrogram_max_height: float = 0.08
    title_fontsize: Optional[int] = None
    default_heatmap_kwargs: dict[str, Any] = field(
        default_factory=lambda: {
            "annot": True,
            "fmt": ".2f",
            "cbar_kws": {},
            "xticklabels": True,
            "yticklabels": True,
        }
    )
    default_clustermap_kwargs: dict[str, Any] = field(
        default_factory=lambda: {
            "col_cluster": True,
            "row_cluster": True,
            "cbar_pos": (0.02, 0.64, 0.02, 0.24),
        }
    )
    clustermap_kwargs: dict[str, Any] = field(default_factory=dict)
    heatmap_kwargs: dict[str, Any] = field(default_factory=dict)


class ClustermapPlotter:
    """Render Theil's U matrices as heatmaps with optional clustering panels.

    Supports:
    - clustered heatmap with top/left dendrograms
    - replacement of row dendrogram by a metadata bar panel
    - automatic layout adjustments to avoid overlap with right-side labels and
        colorbar
    """
    
    def __init__(
        self,
        data: pd.DataFrame,
        metadata_df: Optional[pd.DataFrame] = None,
        metadata_agg_func: Optional[Callable] = None,
        config: Optional[ClustermapConfig] = None
    ):
        """
        Initialize clustermap plotter
        
        Args:
            data: Data matrix to plot (rows and columns will be clustered)
            metadata_df: Optional dataframe for computing barplot metadata
            metadata_agg_func: Function to aggregate metadata (e.g., lambda col: col.nunique())
            config: Plot configuration
        """
        self.data = data
        self.metadata_df = metadata_df if metadata_df is not None else data
        self.metadata_agg_func = metadata_agg_func or (lambda col: col.nunique())
        self.config = config or ClustermapConfig()
        
        self._fig: Optional[Figure] = None
        self._heatmap_ax: Optional[Axes] = None
        self._clustermap = None  # Store clustermap object for barplot access
        self._reordered_labels: Optional[pd.Index] = None
        self._cbar_label: str = ""
        self._fig_size: float = 0
        self._fontsize: int = 0
        self._title_fontsize: int = 0

    def _get_plot_cmap(self):
        """Return a matplotlib colormap instance with explicit NaN color."""
        cmap = copy(plt.get_cmap(self.config.cmap))
        cmap.set_bad(self.config.nan_color)
        return cmap

    @staticmethod
    def _merge_plot_kwargs(defaults: dict[str, Any], overrides: dict[str, Any]) -> dict[str, Any]:
        """Merge plot kwargs, including nested dicts for selected keys.

        Args:
            defaults: Base keyword arguments.
            overrides: User-provided keyword arguments.

        Returns:
            Merged keyword arguments with nested merge for `annot_kws` and
            `cbar_kws`.
        """
        merged = {**defaults, **overrides}
        for nested_key in ("annot_kws", "cbar_kws"):
            if nested_key in defaults or nested_key in overrides:
                merged[nested_key] = {
                    **defaults.get(nested_key, {}),
                    **overrides.get(nested_key, {}),
                }
        return merged

    def _base_heatmap_kwargs(self) -> dict[str, Any]:
        """Build shared seaborn heatmap kwargs used in both rendering paths."""
        return {
            **self.config.default_heatmap_kwargs,
            "cmap": self._get_plot_cmap(),
            "vmin": self.config.vmin,
            "vmax": self.config.vmax,
            "annot_kws": {"size": self._fontsize},
        }

    def _style_colorbar_label(self, cbar_ax: Axes) -> None:
        """Apply consistent colorbar label position and typography.

        Args:
            cbar_ax: Colorbar axes to style.
        """
        cbar_ax.yaxis.set_label_position('left')
        cbar_ax.set_ylabel(
            self._cbar_label,
            rotation=90,
            ha='center',
            va='bottom',
            fontsize=self._fontsize,
            labelpad=8,
        )

    def _get_right_label_margin(self) -> float:
        """Estimate required figure margin for right-side y tick labels.

        Returns:
            Fractional figure width to reserve on the right.
        """
        self._fig.canvas.draw()
        renderer = self._fig.canvas.get_renderer()
        yticklabels = [lab for lab in self._heatmap_ax.get_yticklabels() if lab.get_visible()]

        fig_w_px = self._fig.get_size_inches()[0] * self._fig.dpi
        if not yticklabels or fig_w_px <= 0:
            return 0.08

        max_label_w = max(lab.get_window_extent(renderer=renderer).width for lab in yticklabels)
        return max_label_w / fig_w_px + 0.03

    def _get_bottom_label_margin(self) -> float:
        """Estimate required figure margin for rotated x tick labels.

        Returns:
            Fractional figure height to reserve below the heatmap.
        """
        self._fig.canvas.draw()
        renderer = self._fig.canvas.get_renderer()
        xticklabels = [lab for lab in self._heatmap_ax.get_xticklabels() if lab.get_visible()]

        fig_h_px = self._fig.get_size_inches()[1] * self._fig.dpi
        if not xticklabels or fig_h_px <= 0:
            return self.config.layout_bottom

        max_label_h = max(lab.get_window_extent(renderer=renderer).height for lab in xticklabels)
        return max(self.config.layout_bottom, max_label_h / fig_h_px + 0.02)

    def _get_row_fraction(self, row_ax: Axes) -> float:
        """Derive row-panel width share from seaborn layout and config scaling.

        Args:
            row_ax: Axes currently used as row dendrogram/bar panel.

        Returns:
            Width fraction for the row-side panel, clamped to sensible bounds.
        """
        original_row_pos = row_ax.get_position()
        original_heatmap_pos = self._heatmap_ax.get_position()
        total_original_w = original_row_pos.width + original_heatmap_pos.width
        if total_original_w <= 0:
            base_fraction = 0.18
        else:
            base_fraction = original_row_pos.width / total_original_w
        return min(0.42, max(0.12, base_fraction * self.config.bar_panel_width_scale))

    def _get_dendrogram_height(self, total_h: float) -> float:
        """Compute dendrogram height from total panel height with clipping."""
        scaled = total_h * self.config.dendrogram_ratio * self.config.dendrogram_height_scale
        return min(self.config.dendrogram_max_height, max(self.config.dendrogram_min_height, scaled))

    def plot(
        self,
        title: str,
        output_path: str,
        cbar_label: str = "Value",
        bar_xlabel: str = "# unique"
    ) -> None:
        """
        Create and save the complete plot.
        
        Args:
            title: Plot title
            output_path: Path to save figure
            cbar_label: Label for colorbar
            bar_xlabel: Label for barplot x-axis
        """
        self._calculate_plot_parameters()
        
        if self.config.show_dendrogram and self.data.shape[0] > 2:
            self._create_clustermap(cbar_label)
        else:
            self._create_simple_heatmap(cbar_label)
        
        if self.config.show_barplot and self._clustermap is not None:
            self._add_barplot(bar_xlabel)

        self._configure_heatmap_axis()

        if self._clustermap is not None:
            self._apply_clustermap_layout()
        
        # Avoid tight_layout for clustermap: it can distort seaborn's custom axis geometry.
        self._fig.suptitle(title, y=self.config.title_y, fontsize=self._title_fontsize)
        if self._clustermap is None:
            self._fig.tight_layout(rect=(0, 0, 1, 0.95))
        plt.savefig(output_path, dpi=self.config.dpi, pad_inches=0.2, bbox_inches='tight')
        plt.close(self._fig)
    
    def _calculate_plot_parameters(self) -> None:
        """Derive figure and text sizes from matrix dimensions and config."""
        n_items = max(self.data.shape)

        self._fig_size = max(self.config.min_fig_size, n_items * self.config.size_per_item)
        self._fontsize = max(
            self.config.min_fontsize,
            min(
                self.config.max_fontsize,
                int((self._fig_size / n_items) * self.config.fontsize_multiplier),
            ),
        )
        self._title_fontsize = self.config.title_fontsize if self.config.title_fontsize is not None else self.config.max_fontsize
    
    def _create_clustermap(self, cbar_label: str) -> None:
        """Create seaborn clustermap and cache core axes/ordering metadata."""
        self._cbar_label = cbar_label
        # Keep NaNs for plotting, but use a filled copy for linkage computation.
        cluster_data = self.data.fillna(0.0)

        # Theil's U is directional, so symmetrize before deriving distances for clustering.
        symmetric_similarity = (cluster_data.values + cluster_data.values.T) / 2.0
        np.fill_diagonal(symmetric_similarity, 1.0)

        # Convert symmetric similarity to distance for hierarchical clustering.
        distance_matrix = 1 - symmetric_similarity
        np.fill_diagonal(distance_matrix, 0.0)
        distances = squareform(distance_matrix, checks=True)
        col_linkage = linkage(distances, method='average')

        clustermap_defaults: dict[str, Any] = {
            "col_linkage": col_linkage,
            "row_linkage": col_linkage,  # Use same linkage for rows
            "figsize": (self._fig_size, self._fig_size),
            "dendrogram_ratio": self.config.dendrogram_ratio,
        }
        clustermap_defaults = {**self.config.default_clustermap_kwargs, **clustermap_defaults}
        clustermap_defaults = {**self._base_heatmap_kwargs(), **clustermap_defaults}
        clustermap_kwargs = self._merge_plot_kwargs(clustermap_defaults, self.config.clustermap_kwargs)
        self._clustermap = sns.clustermap(self.data, **clustermap_kwargs)
        
        self._fig = self._clustermap.fig
        self._heatmap_ax = self._clustermap.ax_heatmap
        # Get reordered labels from row dendrogram for consistent ordering
        self._reordered_labels = self.data.index[self._clustermap.dendrogram_row.reordered_ind]

    def _apply_clustermap_layout(self) -> None:
        """Reposition clustermap axes to avoid overlaps and preserve readability.

        This method places row panel, heatmap, top dendrogram, and colorbar
        explicitly using figure-relative coordinates.
        """
        row_ax = self._clustermap.ax_row_dendrogram
        col_ax = self._clustermap.ax_col_dendrogram
        cbar_ax = self._clustermap.cax

        right_label_margin = self._get_right_label_margin()
        bottom_label_margin = self._get_bottom_label_margin()

        left = self.config.layout_left
        right = self.config.layout_right
        bottom = bottom_label_margin
        top = min(self.config.layout_top, self.config.title_y - self.config.title_gap)
        gap = self.config.layout_gap
        cbar_w = self.config.cbar_width
        cbar_gap = self.config.cbar_gap
        right_outer_margin = self.config.outer_right_margin

        # Keep the colorbar inside the figure and reserve label space between it and the heatmap.
        cbar_x = right - right_outer_margin - cbar_w
        content_right = cbar_x - cbar_gap - right_label_margin
        content_right = max(content_right, left + 0.35)

        # Preserve relative width allocated by seaborn between row panel and heatmap.
        row_fraction = self._get_row_fraction(row_ax)

        total_w = content_right - left
        row_w = total_w * row_fraction
        available_heatmap_w = total_w - row_w - gap

        # Keep dendrogram compact and reserve only a narrow band for the title.
        total_h = top - bottom
        dendrogram_h = self._get_dendrogram_height(total_h)
        available_heatmap_h = total_h - dendrogram_h - gap

        # Enforce a square heatmap panel by using one side length for width and height.
        heatmap_side = min(available_heatmap_w, available_heatmap_h)

        used_w = row_w + gap + heatmap_side
        x_offset = max(0.0, (total_w - used_w) / 2)
        used_h = heatmap_side + gap + dendrogram_h
        y_offset = max(0.0, (total_h - used_h) / 2)

        row_x = left + x_offset
        heatmap_x = row_x + row_w + gap
        heatmap_y = bottom + y_offset

        row_ax.set_position([row_x, heatmap_y, row_w, heatmap_side])
        self._heatmap_ax.set_position([heatmap_x, heatmap_y, heatmap_side, heatmap_side])
        col_ax.set_position([heatmap_x, heatmap_y + heatmap_side + gap, heatmap_side, dendrogram_h])

        # Vertically center the colorbar next to the square heatmap while keeping it inside the figure.
        cbar_h = min(heatmap_side, max(0.18, heatmap_side * self.config.cbar_height_scale))
        cbar_y = heatmap_y + (heatmap_side - cbar_h) / 2
        cbar_ax.set_position([cbar_x, cbar_y, cbar_w, cbar_h])
        self._style_colorbar_label(cbar_ax)
    
    def _create_simple_heatmap(self, cbar_label: str) -> None:
        """Create standard seaborn heatmap when clustering is disabled/unneeded."""
        self._cbar_label = cbar_label
        self._fig, self._heatmap_ax = plt.subplots(figsize=(self._fig_size, self._fig_size))

        heatmap_defaults: dict[str, Any] = {
            "ax": self._heatmap_ax,
        }
        heatmap_defaults = {**self._base_heatmap_kwargs(), **heatmap_defaults}
        heatmap_kwargs = self._merge_plot_kwargs(heatmap_defaults, self.config.heatmap_kwargs)
        sns.heatmap(self.data, **heatmap_kwargs)
        
        self._reordered_labels = self.data.index
        # Place label to the left of the colorbar.
        cbar_ax = self._heatmap_ax.collections[0].colorbar.ax
        self._style_colorbar_label(cbar_ax)
    
    def _add_barplot(self, bar_xlabel: str) -> None:
        """Draw metadata bar panel in the row-dendrogram axis.

        Args:
            bar_xlabel: X-axis label for the metadata bar panel.
        """
        # Clear the row dendrogram axis and replace with barplot
        bar_ax = self._clustermap.ax_row_dendrogram
        bar_ax.clear()
        
        # Compute metadata values in reordered labels order
        metadata_values = [
            self.metadata_agg_func(self.metadata_df[col]) 
            for col in self._reordered_labels
        ]

        if len(metadata_values) == 0:
            return

        max_value = float(max(metadata_values))
        if max_value <= 0:
            max_value = 1.0
        y_positions = np.arange(len(metadata_values), dtype=float) + 0.5
        
        # Create horizontal barplot
        bar_ax.barh(y_positions, metadata_values, color='black', alpha=0.75, height=0.8)
        
        # Match heatmap alignment
        bar_ax.set_ylim(self._heatmap_ax.get_ylim())
        bar_ax.set_xlim(max_value * 1.2, 0)
        bar_ax.set_xlabel(bar_xlabel, fontsize=self._fontsize)
        bar_ax.set_yticks([])
        
        # Add value labels just beyond the bar tip (leftward, since x-axis is inverted).
        for y, v in zip(y_positions, metadata_values):
            label = v if v >= 1 else f"{v:.2f}"
            bar_ax.text(
                v + max_value * 0.03, y, str(label),
                va='center', ha='right',
                color='black',
                clip_on=False,
            )
        
        sns.despine(ax=bar_ax, left=True, bottom=False)
    
    def _configure_heatmap_axis(self) -> None:
        """Apply consistent tick label rotation and font sizing for heatmap axes."""
        self._heatmap_ax.set_ylabel('')
        self._heatmap_ax.tick_params(axis='x', labelrotation=90, labelsize=self._fontsize, pad=2)
        self._heatmap_ax.tick_params(axis='y', labelrotation=0, labelsize=self._fontsize)

        for label in self._heatmap_ax.get_xticklabels():
            label.set_ha('center')
            label.set_va('top')

        for label in self._heatmap_ax.get_yticklabels():
            label.set_rotation(0)
            label.set_fontsize(self._fontsize)


# ============================================================================
# Snakemake Script
# ============================================================================
from utils.io import read_anndata


# Extract snakemake parameters
input_file = snakemake.input[0]
output_tsv = snakemake.output.tsv
output_png = snakemake.output.plot

sample_key = snakemake.params.get('sample_key')
covariates = snakemake.params.get('covariates', [])
max_unique_continuous = snakemake.params.get('max_unique_continuous', 10)
wildcards = snakemake.wildcards

title = f"Theil's U Heatmap\ndataset: {wildcards['dataset']}\nfile_id: {wildcards['file_id']}"

# Read data
logging.info(f'Read {input_file}...')
obs = read_anndata(input_file, obs='obs', verbose=False).obs

# Run Theil's U analysis
analyzer = TheilsUAnalyzer(
    obs=obs,
    sample_key=sample_key,
    covariates=covariates,
    max_unique_continuous=max_unique_continuous
)
theils_u_matrix, aggregated_df = analyzer.run()

logging.info(f"Saving Theil's U matrix to {output_tsv} ...")
theils_u_matrix.to_csv(output_tsv, sep='\t')

# Plot
logging.info(f'Plot Theil\'s U heatmap to {output_png} ...')
plotter = ClustermapPlotter(
    data=theils_u_matrix,
    metadata_df=aggregated_df,
    metadata_agg_func=lambda col: col.nunique(),
    config=ClustermapConfig(
        show_dendrogram=True,
        show_barplot=True,
    ),
)
plotter.plot(
    title=title,
    output_path=output_png,
    cbar_label="Theil's U",
    bar_xlabel="# unique"
)
logging.info('Done.')
