import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns


def human_format(num):  # https://stackoverflow.com/a/45846841
    num = float('{:.3g}'.format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return '{}{}'.format('{:f}'.format(num).rstrip('0').rstrip('.'), ['', 'K', 'M', 'B', 'T'][magnitude])


def plot_stacked_bar(
    df,
    category_key,
    covariate_key,
    count_type='count',
    palette='tab20',
    title_suffix='',
    fig_width=8,
    fig_height_min=5,
    fig_height_factor=0.3,
    n_legend_categories=50,
    return_fig=False,
):
    """
    Plot a stacked bar chart of the counts of a covariate by a category.
    
    :param df: DataFrame containing the data
    :param category_key: Column name of the category to stratify by, ends up in y-axis
    :param covariate_key: Column name of the covariate to show counts for, ends up in x-axis
    :param count_type: Type of count to show, either 'count' or 'proportion'
    :param palette: Name of the color palette to use, or a list of colors
    :param title_suffix: Suffix to add to the title
    :param fig_width: Width of the figure
    :param fig_height_min: Minimum height of the figure
    :param fig_height_factor: Factor to multiply the number of categories with to determine the height of the figure
    :param n_legend_categories: Maximum number of categories to show in the legend
    """
    colors = palette
    if palette is not None:
        colors = plt.get_cmap(palette).colors
    
    # determine number of cells per group and category
    df = df[[category_key, covariate_key]].copy()
    df[category_key] = df[category_key].astype(str)
    counts_df = df[[category_key, covariate_key]].value_counts(dropna=False).to_frame()
    counts_df['total_count'] = counts_df.groupby(category_key, observed=True).transform('sum')
    
    if count_type == 'proportion':
        counts_df['proportion'] = counts_df['count'] / counts_df['total_count']
    counts_df = counts_df.reset_index(drop=False)
    
    # determine sort order of categories
    total_counts = counts_df.drop_duplicates(category_key) \
        .set_index(category_key)['total_count'].sort_values(ascending=True)
    counts_df[category_key] = pd.Categorical(
        counts_df[category_key].astype(str),
        categories=total_counts.index,
        ordered=True
    )

    # create plot
    num_categories = counts_df[category_key].nunique()
    fig_height = max(fig_height_min, num_categories * fig_height_factor)
    fig = plt.figure(figsize=(fig_width, fig_height))
    
    ax = counts_df.pivot(
        index=category_key,
        columns=covariate_key,
        values=count_type
    ).fillna(0) \
    .sort_index() \
    .plot.barh(
        stacked=True,
        color=colors,
        ax=plt.gca()
    )
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    longest_label = 0
    for y_position, count in enumerate(total_counts):
        x_position=plt.xlim()[1]
        label=f'n={human_format(count)}'
        ax.text(
            x_position,
            y_position,
            label,
            ha='left',
            va='center',
            color='black'
        )
        longest_label = max(longest_label, len(label))
        
    if df[covariate_key].nunique() < n_legend_categories:
        # adjust figure size to avoid overlapping labels with legend
        label_length_factor = longest_label * 0.13
        fig.set_size_inches(fig_width + label_length_factor, fig_height)
        ax.legend(loc='center left', bbox_to_anchor=(1 + label_length_factor / fig_width, 0.5))
    else:
        ax.legend_ = None
    
    title = f'{count_type.capitalize()} of {covariate_key} by {category_key}{title_suffix}'
    ax.set_title(title)
    
    if return_fig:
        return fig


def plot_violin(
    df: pd.DataFrame,
    category_key: str,
    covariate_key: str,
    fig_width: float = 8,
    fig_height_min: float = 5,
    fig_height_factor: float = 0.3,
    return_fig: bool = False,
):
    """
    Plot a violin plot of the distribution of a covariate by a category.
    
    :param df: DataFrame containing the data
    :param category_key: Column name of the category to stratify by, ends up in y-axis
    :param covariate_key: Column name of the covariate to show counts for, ends up in x-axis
    :param fig_width: Width of the figure
    :param fig_height_min: Minimum height of the figure
    :param fig_height_factor: Factor to multiply the number of categories with to determine the height of the figure
    """
    # prepare data
    df = df[[category_key, covariate_key]].copy()
    df = df.dropna(subset=[covariate_key])
    df[category_key] = df[category_key].astype(str)
    
    # get stats
    total_counts = df[category_key].value_counts().sort_values(ascending=True)
    clusters = total_counts.index.values

    grouped_data = [
        df[df[category_key] == cluster][covariate_key]
        for cluster in clusters
    ]
    
    # build plot
    num_categories = len(grouped_data)
    fig_height = max(fig_height_min, num_categories * fig_height_factor)
    plt.figure(figsize=(fig_width, fig_height))

    plt.violinplot(
        grouped_data,
        showmeans=False,
        showmedians=True,
        vert=False,
    )
    plt.yticks(ticks=range(1, len(grouped_data) + 1), labels=clusters)
    plt.title(f'Distribution of {covariate_key} by {category_key}')
    plt.xlabel(covariate_key)
    plt.ylabel(category_key)
    
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # add labels for cluster sizes
    for i, count in enumerate(total_counts):
        plt.text(
            plt.xlim()[1],
            i + 1,
            f'n={human_format(count)}',
            va='center',
            ha='left',
        )
    
    if return_fig:
        return plt.gcf()
