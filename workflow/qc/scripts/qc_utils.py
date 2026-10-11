import ast
from types import SimpleNamespace

import numpy as np
import pandas as pd
import anndata as ad


QC_FLAGS = [
    'n_counts',
    'log1p_n_counts',
    'n_genes',
    'log1p_n_genes',
    'percent_mito',
    'n_counts_mito',
    'percent_ribo',
    'n_counts_ribo',
    'percent_hb',
    'n_counts_hb',
]


def parse_parameters(adata: ad.AnnData, params: dict, filter_hues: bool = False):
    dataset = params.get('dataset', 'None')
    hues = params.get('hue', [])

    split_datasets = dataset.split('--')
    if len(split_datasets) > 1:
        dataset = ' '.join([split_datasets[0], split_datasets[-1]])

    if isinstance(hues, str):
        hues = [hues]
    hues = [hue for hue in hues if hue in adata.obs.columns]
    if filter_hues:
        max_groups = min(params.get('max_groups', 100), 102)
        hues = {
            hue for hue in hues
            if (
                # Keep numeric columns regardless of unique value count (for continuous colormaps)
                pd.api.types.is_numeric_dtype(adata.obs[hue])
                # For categorical/non-numeric, apply max_groups filter
                or (1 < adata.obs[hue].nunique() < max_groups)
            )
        }
        hues = list(hues)
    if len(hues) == 0:
        hues = [None]

    return dataset, hues

def parse_autoqc(autoqc_thresholds: pd.DataFrame):
    if autoqc_thresholds is None or autoqc_thresholds.empty or not 'side' in autoqc_thresholds.columns:
        return autoqc_thresholds
    # side has precedence over low/high values
    for key in autoqc_thresholds.index:
        side = autoqc_thresholds.loc[key, 'side'].lower()
        if side == 'max_only':
            autoqc_thresholds.loc[key, 'low'] = None
        elif side == 'min_only':
            autoqc_thresholds.loc[key, 'high'] = None

    return autoqc_thresholds


def update_thresholds(thresholds: dict, autoqc_thresholds: pd.DataFrame):
    if autoqc_thresholds.empty:
        return thresholds
    
    for key in autoqc_thresholds.index:
        thresholds[f'{key}_min'] = autoqc_thresholds.loc[key, 'low']
        thresholds[f'{key}_max'] = autoqc_thresholds.loc[key, 'high']
    return thresholds

def fill_autoqc_thresholds_from_df(metrics, df):
    for metric in metrics.index:
        if metric not in df.columns:
            min_val = -np.inf
            max_val = np.inf
        else:
            min_val = df[metric].min()
            max_val = df[metric].max()
        if pd.isna(metrics.loc[metric, 'min']):
            metrics.loc[metric, 'min'] = min_val if np.isfinite(min_val) else -np.inf
        if pd.isna(metrics.loc[metric, 'max']):
            metrics.loc[metric, 'max'] = max_val if np.isfinite(max_val) else np.inf
    return metrics


def get_thresholds(
    threshold_keys: list = None,
    autoqc_thresholds: pd.DataFrame = None,
    user_thresholds: [str, dict] = None,
    transform: bool = True,
    init_nan: bool = False,
    df: pd.DataFrame = None,
):
    """
    :param threshold_keys: keys in mappings that define different QC paramters
    :param autoqc_thresholds: autoQC thresholds as provided by sctk under adata.uns['scautoqc_ranges']
    :param user_thresholds: user defined thresholds
    :return: thresholds: dict of key: thresholds tuple
    """
    import ast
    
    # set default thresholds
    if threshold_keys is None:
        threshold_keys = ['n_counts', 'n_genes', 'percent_mito']
    
    # initialise thresholds
    if init_nan:
        if df is not None:
            thresholds = {f'{key}_min': df[key].min() for key in threshold_keys}
            thresholds |= {f'{key}_max': df[key].max() for key in threshold_keys}
        else:
            thresholds = {f'{key}_min': np.nan for key in threshold_keys}
            thresholds |= {f'{key}_max': np.nan for key in threshold_keys}
    else:
        thresholds = {f'{key}_min': 0 for key in threshold_keys}
        thresholds |= {f'{key}_max': np.inf for key in threshold_keys}

    if autoqc_thresholds is not None:
        thresholds = update_thresholds(thresholds, autoqc_thresholds)

    # update to user thresholds
    if user_thresholds is None:
        user_thresholds = {}
    elif isinstance(user_thresholds, str):
        user_thresholds = ast.literal_eval(user_thresholds)
    elif isinstance(user_thresholds, dict):
        pass
    else:
        ValueError('thresholds must be a dict or string')
    thresholds |= user_thresholds
    
    if transform:
        # transform to shape expected by plot_qc_joint
        return {
            key: (thresholds[f'{key}_min'], thresholds[f'{key}_max'])
            for key in threshold_keys
        }
    return {
        key: value
        for key, value in thresholds.items()
        if any([key.startswith(x) for x in threshold_keys])
    }


def apply_thresholds(
    data,
    thresholds: dict,
    threshold_keys: list,
    column_name='passed_qc',
    inplace=True,
):
    """
    :param data: AnnData object or obs DataFrame
    :param thresholds: dict of key: thresholds tuple as returned by get_thresholds
    :param inplace: if True, write the result back to the provided object/Frame
    :return: boolean Series with pass/fail status
    """
    if isinstance(data, pd.DataFrame):
        obs = data if inplace else data.copy()
    elif hasattr(data, 'obs'):
        obs = data.obs if inplace else data.obs.copy()
    else:
        raise TypeError('data must be an AnnData-like object with .obs or a pandas DataFrame')

    obs[column_name] = True
    if obs.shape[0] == 0:
        return obs[column_name]
    for key in threshold_keys:
        lower, upper = thresholds[key]
        if pd.isna(lower) and pd.isna(upper):
            continue
        if pd.isna(lower):
            passed = obs[key] <= upper
        elif pd.isna(upper):
            passed = obs[key] >= lower
        else:
            passed = obs[key].between(lower, upper)
        obs[column_name] = obs[column_name] & passed.fillna(False)
    return obs[column_name]

def log_auto_base(values: list, count_shift: int = 1) -> int:
    """
    Decide whether to use log base 2 or log base 10 for a numeric list.

    Decision hierarchy
    ------------------
    1. If all values are exact powers of 2 use log2.
       Good for: read-depths, bit-widths, genomic bin sizes
    2. If all values are exact powers of 10 use log10.
       Good for: dilution series, pH, concentration, decades
    3. Span heuristic: compare the dynamic range expressed in each base.
       log2(max/min) < log10(max/min) * log2(10) is always true by definition,
       so check data spans < 10 doublings vs < 3 decimals places
       if log2 span <= 10 use log2 (fold-changes), otherwise fall back to log10

    Parameters
    ----------
    values : sequence of positive floats

    Returns
    -------
    2  if log base 2 is recommended
    10 if log base 10 is recommended

    Raises
    ------
    ValueError if the list is empty or contains non-positive values.
    """
    import math

    # Drop NA values and validate input before computing min/max
    data = [float(v + count_shift) for v in values if pd.notna(v)]
    if not data:
        raise ValueError("Input list is empty.")
    if any(v <= 0 for v in data):
        raise ValueError("All values must be strictly positive for log scaling.")

    # rule: dynamic-range
    vmin, vmax = min(data), max(data)

    if vmin == vmax or (
        (0 + count_shift) <= vmin <= (1 + count_shift)
        and (0 + count_shift) <= vmax < (1 + count_shift)
    ):
        # No spread so either base works; or all within 0 and 1
        # no point in log transforming
        return 1

    def _is_integer_close(x: float, tol: float = 1e-3) -> bool:
        return abs(x - round(x)) < tol

    def _all_powers_of(base: float) -> bool:
        return all(
            _is_integer_close(
                math.log(v) / math.log(base)
            )
            for v in data
        )

    # rule: exact powers of 2
    if _all_powers_of(2):
        return 2.0

    # rule: exact powers of 10
    if _all_powers_of(10):
        return 10.0

    log2_span  = math.log2(vmax / vmin)   # number of doublings
    log10_span = math.log10(vmax / vmin)  # decimal order

    # Prefer log2 when the spread fits comfortably within ~10 doublings
    # (~3 orders of magnitude); beyond that, log10 reads more naturally
    LOG2_SPAN_THRESHOLD = 10.0

    if log2_span <= LOG2_SPAN_THRESHOLD and log10_span <= 3:
        return 2.0
    else:
        return 10.0

def log1p_base(_x, base):
    return np.log1p(_x) / np.log(base)

def log_auto(_x):
    base = log_auto_base(_x)
    return log1p_base(_x, base), base

def plot_qc_joint(
    df: pd.DataFrame,
    x: str,
    y: str,
    log_x: bool | int = 1,
    log_y: bool | int = 1,
    hue: str = None,
    main_plot_function=None,
    marginal_hue=None,
    x_threshold=None,
    y_threshold=None,
    x_threshold2=None,
    y_threshold2=None,
    threshold_color='black',
    threshold_color2='black',
    threshold_linestyle='-',
    threshold_linestyle2='--',
    title='',
    return_df=False,
    marginal_kwargs: dict = None,
    fig=None,
    subplot_spec=None,
    sharex=None,
    sharey=None,
    x_max=None,
    y_max=None,
    lim_pad = 0.05,
    **kwargs,
):
    """
    Plot scatter plot with marginal histograms from df columns.

    :param df: observation dataframe
    :param x: df column for x axis
    :param y: df column for y axis
    :param log_x: Whether to log-transform x (No transformation = False | 0 | 1)
    :param log_y: Whether to log-transform y (No transformation = False | 0 | 1)
    :param hue: df column with annotations for color coding scatter plot points
    :param marginal_hue: df column with annotations for color coding marginal plot distributions
    :param x_threshold: tuple of (min, max) filter thresholds for x axis
    :param y_threshold: tuple of (min, max) filter thresholds for y axis
    :param fig: existing matplotlib Figure used when `subplot_spec` is provided
    :param subplot_spec: optional matplotlib SubplotSpec to draw into an existing grid layout.
        When provided, creates joint/marginal axes inside this slot instead of creating a new JointGrid figure.
    :param sharex: existing Axes to share x axis with (subplot_spec mode only)
    :param sharey: existing Axes to share y axis with (subplot_spec mode only)
    :param title: title text for plot
    :return: seaborn plot (and df dataframe with updated values, if `return_df=True`)
    """
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, FixedLocator, FixedFormatter
    import seaborn as sns

    df = df.copy()

    if main_plot_function is None:
        main_plot_function = sns.scatterplot

    def _normalize_threshold_pair(threshold):
        if threshold is None:
            return (0, np.inf)
        lower, upper = threshold
        lower = 0 if lower is None or pd.isna(lower) else lower
        upper = np.inf if upper is None or pd.isna(upper) else upper
        return (lower, upper)

    x_threshold  = _normalize_threshold_pair(x_threshold)
    y_threshold  = _normalize_threshold_pair(y_threshold)
    x_threshold2 = _normalize_threshold_pair(x_threshold2)
    y_threshold2 = _normalize_threshold_pair(y_threshold2)

    orig_x = df[x].copy() if log_x > 1 else None
    orig_y = df[y].copy() if log_y > 1 else None

    if log_x > 1:
        df[x] = log1p_base(orig_x, log_x)
        x_threshold, x_threshold2 = log1p_base(x_threshold, log_x), log1p_base(x_threshold2, log_x)

    if log_y > 1:
        df[y] = log1p_base(orig_y, log_y)
        y_threshold, y_threshold2 = log1p_base(y_threshold, log_y), log1p_base(y_threshold2, log_y)

    def thresholds_equal(a, b):
        a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
        return a.shape == b.shape and np.allclose(a, b, equal_nan=True)

    if thresholds_equal(x_threshold, x_threshold2): x_threshold2 = (0, np.inf)
    if thresholds_equal(y_threshold, y_threshold2): y_threshold2 = (0, np.inf)

    if marginal_kwargs is None:
        marginal_kwargs = dict(legend=False, palette=kwargs.get('palette'))

    if marginal_hue in df.columns and df[marginal_hue].nunique() > 100:
        marginal_hue = None
    use_marg_hue = marginal_hue is not None

    if not use_marg_hue:
        marginal_kwargs.pop('palette', None)

    hue_order = df[hue].value_counts(ascending=False).index if (hue in df.columns and hasattr(df[hue], 'cat')) else None

    def _draw_thresholds(ax_joint, ax_marg_x, ax_marg_y):
        for thresh, color, ls, axline, marg_ax in [
            (x_threshold,  threshold_color,  threshold_linestyle,  'axvline', ax_marg_x),
            (x_threshold2, threshold_color2, threshold_linestyle2, 'axvline', ax_marg_x),
            (y_threshold,  threshold_color,  threshold_linestyle,  'axhline', ax_marg_y),
            (y_threshold2, threshold_color2, threshold_linestyle2, 'axhline', ax_marg_y),
        ]:
            for t, t_def in zip(thresh, (0, np.inf)):
                if t != t_def:
                    kw = dict(color=color, linestyle=ls)
                    getattr(ax_joint, axline)(**{'x' if axline == 'axvline' else 'y': t}, **kw)
                    getattr(marg_ax,  axline)(**{'x' if axline == 'axvline' else 'y': t}, **kw)

    def _plot_on_axes(ax_joint, ax_marg_x, ax_marg_y, x_max=None, y_max=None):
        main_plot_function(data=df, x=x, y=y, hue=hue, hue_order=hue_order, ax=ax_joint, **kwargs)

        if hue is not None and kwargs.get('legend', True):
            handles, labels = ax_joint.get_legend_handles_labels()
            if handles:
                ax_joint.legend(
                    handles=handles, labels=labels,
                    markerscale=(60 / kwargs.get('s', 20)) ** 0.5,
                    fontsize=ax_joint.xaxis.label.get_size(),
                )

        hist_kwargs = dict(
            data=df, hue=marginal_hue,
            hue_order=hue_order if use_marg_hue else None,
            element='step' if use_marg_hue else 'bars',
            fill=False, bins=100, **marginal_kwargs,
        )
        sns.histplot(x=x, ax=ax_marg_x, **hist_kwargs)
        sns.histplot(y=y, ax=ax_marg_y, **hist_kwargs)
        _draw_thresholds(ax_joint, ax_marg_x, ax_marg_y)

        # hide from marginal axes
        for ax in (ax_marg_x, ax_marg_y):
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_xlabel('')
            ax.set_ylabel('')
        ax_marg_x.tick_params(axis='x', bottom=False, labelbottom=False)
        ax_marg_y.tick_params(axis='y', left=False, labelleft=False)

        def human_format(v, _):
            if v >= 1e6: return f"{v/1e6:.3g}M"
            if v >= 1e3: return f"{v/1e3:.3g}K"
            return f"{v:.3g}"

        def set_log_ticks(ax_main, ax_marg, axis, raw_vals, log_base):
            orig_max = np.nanmax(raw_vals) if len(raw_vals) else 1.0
            max_exp = int(np.floor(np.log10(max(orig_max, 1))))
            tick_vals = [0] + [10 ** e for e in range(0, max_exp + 1)]
            tick_pos = [log1p_base(v, log_base) for v in tick_vals]
            for ax in (ax_main, ax_marg):
                getattr(ax, f'{axis}axis').set_major_locator(FixedLocator(tick_pos))
                getattr(ax, f'{axis}axis').set_major_formatter(FixedFormatter([human_format(v, None) for v in tick_vals]))

        def _col_max(col, raw_vals, log_base):
            raw = np.nanmax(raw_vals) if raw_vals is not None and len(raw_vals) else np.nanmax(df[col])
            return (log1p_base(raw, log_base) if log_base > 1 else raw) or 1.0

        if x_max is None: x_max = _col_max(x, orig_x, log_x) * (1 + lim_pad)
        if y_max is None: y_max = _col_max(y, orig_y, log_y) * (1 + lim_pad)

        if log_x > 1: set_log_ticks(ax_joint, ax_marg_x, 'x', orig_x, log_x)
        else:
            for ax in (ax_joint, ax_marg_x): ax.xaxis.set_major_formatter(FuncFormatter(human_format))

        if log_y > 1: set_log_ticks(ax_joint, ax_marg_y, 'y', orig_y, log_y)
        else:
            for ax in (ax_joint, ax_marg_y): ax.yaxis.set_major_formatter(FuncFormatter(human_format))
        
        # format non-shared marginal axis ticks
        ax_marg_x.yaxis.set_major_formatter(FuncFormatter(human_format))
        ax_marg_y.xaxis.set_major_formatter(FuncFormatter(human_format))

        ax_joint.set_xlim(0, x_max)
        ax_joint.set_ylim(0, y_max)
        ax_joint.set_xlabel(f"log {x}" if log_x > 1 else x)
        ax_joint.set_ylabel(f"log {y}" if log_y > 1 else y)
        ax_joint.spines['top'].set_visible(False)
        ax_joint.spines['right'].set_visible(False)

        return SimpleNamespace(
            fig=ax_joint.figure, ax_joint=ax_joint,
            ax_marg_x=ax_marg_x, ax_marg_y=ax_marg_y,
            _figsize=ax_joint.figure.get_size_inches(),
        )

    if subplot_spec is not None:
        if fig is None:
            fig = getattr(subplot_spec, 'figure', None) or subplot_spec.get_gridspec().figure
        inner = subplot_spec.subgridspec(2, 2, height_ratios=[1, 4], width_ratios=[4, 1], hspace=0, wspace=0)
        ax_joint  = fig.add_subplot(inner[1, 0])
        ax_marg_x = fig.add_subplot(inner[0, 0], sharex=ax_joint)
        ax_marg_y = fig.add_subplot(inner[1, 1], sharey=ax_joint)
        fig.add_subplot(inner[0, 1]).set_axis_off()

        g = _plot_on_axes(ax_joint, ax_marg_x, ax_marg_y, x_max=x_max, y_max=y_max)
        if title:
            ax_marg_x.set_title(title, fontsize=12, pad=10)
        if return_df:
            return g, df
        return g

    g = sns.JointGrid(data=df, x=x, y=y)
    g = _plot_on_axes(g.ax_joint, g.ax_marg_x, g.ax_marg_y, x_max=x_max, y_max=y_max)
    g.fig.suptitle(title, fontsize=14)
    if return_df:
        return g, df
    return g


class DatashaderKDE:
    """Fast KDE for JointGrid with outlier sensitivity and scale awareness."""

    def __init__(
        self,
        sigma=5,
        plot_width=400,
        plot_height=400,
        outlier_boost=1.0,
        outlier_quantile=0.95,
    ):
        self.sigma = sigma
        self.plot_width = plot_width
        self.plot_height = plot_height
        self.outlier_boost = outlier_boost
        self.outlier_quantile = outlier_quantile

    def plot(
        self,
        x_data,
        y_data,
        ax=None,
        cmap='plasma',
        alpha=0.8,
        levels=5,
        thresh=0.05,
        padding=0.05,
        **kwargs,
    ):
        """Plot KDE density. Usage: g.plot_joint(kde.plot, ...)"""
        import matplotlib.pyplot as plt

        if ax is None:
            ax = plt.gca()

        # Process data pipeline
        x_clean, y_clean = self._clean_data(x_data, y_data)
        if len(x_clean) == 0:
            return ax

        if self.outlier_boost > 1.0:
            x_clean, y_clean = self._boost_outliers(x_clean, y_clean)

        x_range, y_range = self._compute_ranges(x_clean, y_clean, padding)
        agg = self._rasterize(x_clean, y_clean, x_range, y_range)

        sigma_data = self._compute_sigma(ax, x_range)
        smoothed = self._smooth(agg.values, sigma_data)

        # Render
        self._render(ax, smoothed, x_range, y_range, cmap, alpha, thresh, levels)
        return ax

    def _clean_data(self, x, y):
        """Remove NaN/Inf values."""
        x, y = np.asarray(x), np.asarray(y)
        mask = np.isfinite(x) & np.isfinite(y)
        return x[mask], y[mask]

    def _boost_outliers(self, x, y):
        """Amplify low-density outliers by replication."""
        from sklearn.neighbors import KDTree

        xy = np.column_stack([x, y])
        tree = KDTree(xy)
        densities = tree.query(xy, k=20, return_distance=False).shape[1] / 20

        outlier_mask = densities < np.quantile(densities, self.outlier_quantile)
        n_outliers = outlier_mask.sum()

        if n_outliers == 0:
            return x, y

        boost_factor = int(self.outlier_boost)
        outlier_points = np.repeat(xy[outlier_mask], boost_factor, axis=0)
        main_points = xy[~outlier_mask]
        xy_boosted = np.vstack([main_points, outlier_points])

        return xy_boosted[:, 0], xy_boosted[:, 1]

    def _compute_ranges(self, x, y, padding):
        """Compute padded data limits."""
        x_range = [
            x.min() - padding * (x.max() - x.min()),
            x.max() + padding * (x.max() - x.min()),
        ]
        y_range = [
            y.min() - padding * (y.max() - y.min()),
            y.max() + padding * (y.max() - y.min()),
        ]
        return x_range, y_range

    def _rasterize(self, x, y, x_range, y_range):
        """Rasterize points to density aggregate."""
        import datashader as ds

        df = pd.DataFrame({"x": x, "y": y})
        canvas = ds.Canvas(self.plot_width, self.plot_height, x_range, y_range)
        return canvas.points(df, "x", "y")

    def _compute_sigma(self, ax, x_range):
        """Scale-aware sigma calculation."""
        x_is_log = ax.get_xscale() == 'log'
        y_is_log = ax.get_yscale() == 'log'

        pixel_size = (x_range[1] - x_range[0]) / self.plot_width
        base_sigma = max(0.1, min(self.sigma * pixel_size, 50))

        log_factor = (
            0.3
            if (x_is_log and y_is_log)
            else (0.5 if (x_is_log or y_is_log) else 1.0)
        )
        outlier_factor = 0.8 if self.outlier_boost > 1.0 else 1.0

        return base_sigma * log_factor * outlier_factor

    def _smooth(self, data, sigma):
        """Apply gaussian smoothing with error handling."""
        from scipy.ndimage import gaussian_filter

        try:
            return gaussian_filter(data, sigma=sigma)
        except OverflowError:
            return gaussian_filter(data, sigma=2.0)

    def _render(
        self,
        ax,
        smoothed,
        x_range,
        y_range,
        cmap,
        alpha,
        thresh,
        levels,
        **kwargs,
    ):
        """Render smoothed density with masking and contours."""
        smoothed_norm = (
            smoothed / smoothed.max() if smoothed.max() > 0 else smoothed
        )
        masked = np.ma.masked_where(smoothed_norm < thresh, smoothed_norm)

        extent = [*x_range, *y_range]
        ax.imshow(
            masked,
            origin="lower",
            extent=extent,
            aspect="auto",
            cmap=cmap,
            alpha=alpha,
            interpolation="bilinear",
            **kwargs,
        )


def plot_density(
    df,
    x,
    y,
    thresholds,
    auto_thresholds=None,
    log_x=1,
    log_y=1,
    large_threshold=500_000,
    kde_plot_kwargs: dict=None,
    datashader_kwargs: dict=None,
    force_kde=False,
    threshold_color='black',
    threshold_color2='black',
    threshold_linestyle='-',
    threshold_linestyle2='--',
    title='',
    fig=None,
    subplot_spec=None,
    sharex=None,
    sharey=None,
):
    """
    Plots bivariate density using seaborn KDE for small datasets
    and Datashader for large datasets, integrated with plot_qc_joint.
    
    Parameters
    ----------
    df : pd.DataFrame
    x, y : str
        Columns to plot
    thresholds : dict
        Optional thresholds for joint QC plotting
    log_x, log_y : bool
        Whether to log-transform x or y
    large_threshold : int
        Dataset size above which Datashader is used
    """
    import seaborn as sns
    from matplotlib import pyplot as plt

    if kde_plot_kwargs is None:
        kde_plot_kwargs = {}
    kde_plot_kwargs = dict(fill=True, cmap='plasma', alpha=0.8, **kde_plot_kwargs)

    force_kde |= 100_000 < df.shape[0] < large_threshold
    if force_kde:
        # subset to max of 300k cells due to high computational cost
        df = df.sample(n=int(min(300_000, df.shape[0])), random_state=42)

    df_plot = df.copy()

    if len(df_plot) > large_threshold:
        if datashader_kwargs is None:
            datashader_kwargs = {}
        datashader_kwargs = dict(
            plot_width=800,
            plot_height=800,
            outlier_boost=2,
            outlier_quantile=0.95,
        ) | datashader_kwargs
        kde = DatashaderKDE(**datashader_kwargs)
        kde_plot_func = kde.plot
    else:
        kde_plot_func = sns.kdeplot
        if len(df_plot) > 5e4:
            # heuristic for faster fitting
            kde_plot_kwargs.update(bw_adjust=4, gridsize=75, thresh=0.2, levels=5)

    plot_qc_joint(
        df_plot,
        x=x,
        y=y,
        main_plot_function=kde_plot_func,
        log_x=log_x,
        log_y=log_y,
        x_threshold=thresholds.get(x),
        y_threshold=thresholds.get(y),
        x_threshold2=(auto_thresholds.get(x) if auto_thresholds else None),
        y_threshold2=(auto_thresholds.get(y) if auto_thresholds else None),
        threshold_color=threshold_color,
        threshold_color2=threshold_color2,
        threshold_linestyle=threshold_linestyle,
        threshold_linestyle2=threshold_linestyle2,
        title=title,
        fig=fig,
        subplot_spec=subplot_spec,
        sharex=sharex,
        sharey=sharey,
        **kde_plot_kwargs,
    )
