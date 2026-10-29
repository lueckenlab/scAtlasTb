input_file <- snakemake@input$tsv
output_tsv <- snakemake@output$tsv
output_svg <- snakemake@output$svg
output_svg_overall <- snakemake@output$svg_overall

id_vars <- snakemake@params$id_vars
value_var <- snakemake@params$value_var
variable_var <- snakemake@params$variable_var
weight_batch <- snakemake@params$weight_batch
weight_batch <- ifelse(is.null(weight_batch), 0.4, weight_batch)
group_col <- snakemake@params$group_col
group_col_order <- snakemake@params$group_col_order
metric_order <- snakemake@params$metric_order

scale <- snakemake@params$scale
scale <- ifelse(is.null(scale), FALSE, scale)
scale_by_group <- snakemake@params$scale_by_group

n_top <- snakemake@params$n_top
n_top <- if (is.null(n_top) || length(n_top) == 0 || is.na(n_top)) Inf else as.integer(n_top)
dpi <- snakemake@params$dpi
dpi <- ifelse(is.null(dpi), 300, dpi)

parse_order <- function(x, name = 'metric_order') {
  if (is.null(x)) return(NULL)
  if (is.list(x)) x <- unlist(x)
  if (!is.character(x)) stop(paste0(name, " must be a character vector or list of metric ids"))
  # disallow comma-separated single strings to avoid implicit parsing
  if (length(x) == 1 && grepl(',', x)) {
    stop(paste0(name, " must be provided as a list/character vector, not a comma-separated string"))
  }
  trimws(x)
}

# Single-list ordering: user provides `metric_order` as a character vector/list.
metric_order <- parse_order(metric_order, 'metric_order')
group_col_order <- parse_order(group_col_order, 'group_col_order')

tryCatch({
 library(funkyheatmap)
  }, error = function(e) {
    install.packages('funkyheatmap', repos = snakemake@params$cran_url)
  }
)
suppressPackageStartupMessages({
  library(funkyheatmap)
  library(data.table)
  library(tidyverse)
  library(dynutils)
  library(RColorBrewer)
})
source(snakemake@input$r_utils)

# read files
dt <- fread(input_file)
dt[, (id_vars) := lapply(.SD, as.character), .SDcols = id_vars]
print(head(dt))

extra_columns <- readLines(snakemake@input$extra_columns)
print(extra_columns)
id_vars <- unique(c(id_vars, extra_columns))

# remove unintegrated output types without corresponding method
if ('unintegrated' %in% dt$method & uniqueN(dt$method) > 1) {
  ot_count <- dt[method != 'unintegrated', .(output_type, method)] %>%
    unique %>%
    .[, output_type] %>%
    table
  keep_ot <- ot_count[ot_count > 1]
  df <- dt[, output_type %in% keep_ot]
}

# define groups of interest to be plotted in funkyheatmap
integration_setup <- id_vars
bio_metrics <- unique(dt[metric_type == 'bio_conservation' & !is.na(score), metric]) 
batch_metrics <- unique(dt[metric_type == 'batch_correction' & !is.na(score), metric])
metrics <- c(batch_metrics, bio_metrics)

# If user provided a combined `metric_order`, split it into batch/bio parts
metric_order_batch <- NULL
metric_order_bio <- NULL
if (!is.null(metric_order)) {
  metric_order_batch <- metric_order[metric_order %in% batch_metrics]
  metric_order_bio <- metric_order[metric_order %in% bio_metrics]
}

# do not override detected metric lists here; consolidation with automatic
# ordering happens later (after order_metric_block is defined)
metrics <- c(batch_metrics, bio_metrics)

# subset data.table to columns of interest & transform
all_columns <- c(id_vars, variable_var, value_var, 'metric_type')
id_formula <- paste(paste(id_vars, collapse='+'), variable_var, sep = '~')
cat(paste(id_formula, '\n'))
metrics_tab <- dcast(
  subset(dt, select = all_columns),
  as.formula(id_formula),
  value.var = value_var,
  fun.aggregate = function(x) mean(x, na.rm = TRUE)
)
print(head(metrics_tab))

# remove columns that are all empty
metrics_tab <- metrics_tab[, .SD, .SDcols = colSums(is.na(metrics_tab)) < nrow(metrics_tab)]

# scores should be already scaled [0,1] - however, we aim to compute the scores based on the min-max scaled metrics
safe_scale_minmax <- function(x) {
  if (all(is.na(x)) || length(unique(x[!is.na(x)])) <= 1) {
    return(x)
  }
  scale_minmax(x)
}

scaled_metrics_tab <- metrics_tab[, metrics, with=FALSE]
if (nrow(scaled_metrics_tab) > 1) {
  if (scale_by_group && !is.null(group_col) && group_col %in% colnames(metrics_tab)) {
    grp <- as.character(metrics_tab[[group_col]])
    grp[is.na(grp)] <- "NA"
    scaled_metrics_tab[, .scale_group := grp]
    scaled_metrics_tab[, (metrics) := lapply(.SD, safe_scale_minmax), by = .scale_group, .SDcols = metrics]
    scaled_metrics_tab[, .scale_group := NULL]
  } else {
    # default scaling: global
    scaled_metrics_tab <- as.data.table(apply(as.matrix(scaled_metrics_tab), 2, safe_scale_minmax))
  }
}


if (scale) {
  # update original metrics_tab with scaled values for the metrics columns
  metrics_tab[, (metrics) := scaled_metrics_tab]
}

cat('Scaled metrics table: \n')
print(head(scaled_metrics_tab))

# calculate average score by metric type
summarize_prefixed_metrics <- function(scaled_tab, metric_names, dataset_values = NULL) {
  # Summary metrics are named as: <prefix><sep><dataset>
  summary_metric_sep <- "__"
  prefixed <- metric_names[grepl(":", metric_names, fixed = TRUE)]
  if (length(prefixed) == 0) {
    return(list(scaled_tab = scaled_tab, score_metrics = metric_names, summary_metrics = character(0)))
  }

  prefixes <- unique(sub(":.*", "", prefixed))
  summary_metrics <- character(0)

  dataset_values <- if (is.null(dataset_values)) NULL else as.character(dataset_values)
  if (!is.null(dataset_values)) {
    # Keep missing dataset labels explicit so rows are still grouped consistently.
    dataset_values[is.na(dataset_values)] <- "NA"
    dataset_levels <- unique(dataset_values)
  }

  for (prefix in prefixes) {
    cols <- prefixed[startsWith(prefixed, paste0(prefix, ":"))]
    summary_col <- rowMeans(scaled_tab[, cols, with = FALSE], na.rm = TRUE)

    if (is.null(dataset_values)) {
      # Without dataset labels, create one summary metric per prefix.
      scaled_tab[, (prefix) := summary_col]
      summary_metrics <- c(summary_metrics, prefix)
    } else {
      # With dataset labels, create one summary metric per prefix and dataset.
      for (ds in dataset_levels) {
        metric_name <- paste(prefix, ds, sep = summary_metric_sep)
        scaled_tab[, (metric_name) := NA_real_]
        scaled_tab[dataset_values == ds, (metric_name) := summary_col[dataset_values == ds]]
        summary_metrics <- c(summary_metrics, metric_name)
      }
    }
  }

  score_metrics <- c(metric_names[!metric_names %in% prefixed], summary_metrics)
  list(scaled_tab = scaled_tab, score_metrics = score_metrics, summary_metrics = summary_metrics)
}

dataset_values <- if ("dataset" %in% colnames(metrics_tab)) metrics_tab$dataset else NULL

batch_summary <- summarize_prefixed_metrics(scaled_metrics_tab, batch_metrics, dataset_values)
scaled_metrics_tab <- batch_summary$scaled_tab
batch_metrics_score <- batch_summary$score_metrics
batch_metrics_overall <- batch_summary$summary_metrics

bio_summary <- summarize_prefixed_metrics(scaled_metrics_tab, bio_metrics, dataset_values)
scaled_metrics_tab <- bio_summary$scaled_tab
bio_metrics_score <- bio_summary$score_metrics
bio_metrics_overall <- bio_summary$summary_metrics

summary_metrics_overall <- unique(c(batch_metrics_overall, bio_metrics_overall))
if (length(summary_metrics_overall) > 0) {
  metrics_tab[, (summary_metrics_overall) := scaled_metrics_tab[, summary_metrics_overall, with = FALSE]]
}

score_group_batch <- rowMeans(scaled_metrics_tab[, batch_metrics_score, with=FALSE], na.rm = TRUE)
if(length(score_group_batch) > 0){
  metrics_tab <- add_column(metrics_tab, "Batch Correction" = score_group_batch, .before = batch_metrics[1])
}

score_group_bio <- rowMeans(scaled_metrics_tab[, bio_metrics_score, with=FALSE], na.rm = TRUE)
if(length(score_group_bio) > 0){
  metrics_tab <- add_column(metrics_tab, "Bio Conservation" = score_group_bio, .before = bio_metrics[1])
}

# weighted overall score
score_all <- (weight_batch * score_group_batch + (1 - weight_batch) * score_group_bio)
if(length(score_all) > 0){
  metrics_tab <- add_column(metrics_tab, "Overall Score" = score_all, .before = "Batch Correction")
  metrics_tab <- metrics_tab[order(-`Overall Score`, na.last = TRUE)]
}

metric_columns <- c(
  "Overall Score",
  "Batch Correction",
  batch_metrics,
  batch_metrics_overall,
  "Bio Conservation",
  bio_metrics,
  bio_metrics_overall
)
metric_columns <- metric_columns[metric_columns %in% colnames(metrics_tab)]
non_metric_columns <- setdiff(colnames(metrics_tab), metric_columns)
metrics_tab <- metrics_tab[, c(non_metric_columns, metric_columns), with = FALSE]

# write output summary.csv file before plotting ? Michaela to check the output dir?
cat('Writing output file: ', output_tsv, '\n')
fwrite(metrics_tab, output_tsv, sep='\t')

# Check if we have enough data to create a meaningful heatmap
if (nrow(metrics_tab) < 2) {
  cat("Warning: Only", nrow(metrics_tab), "row(s) available. Skipping funkyheatmap generation.\n")
  cat("Creating empty SVG placeholder...\n")
  svg(output_svg, width = 8, height = 4)
  plot.new()
  text(
    0.5, 0.5,
    paste(
      "Insufficient data for heatmap visualization\n",
      "Only", nrow(metrics_tab), "method(s) available\n",
      "See", basename(output_tsv), "for results"
    ),
    cex = 1.2, col = "gray40"
  )
  dev.off()
  file.copy(output_svg, output_svg_overall, overwrite = TRUE)
  quit(save = "no", status = 0)
}

# subset data is applied after sorting; per-group if grouping is enabled

### add funkyheatmap data

row_info <- NULL
row_groups <- NULL

# Group rows by the specified column if provided
if ("Overall Score" %in% colnames(metrics_tab)) {
  if (!is.null(group_col)) {
    assertthat::assert_that(group_col %in% colnames(metrics_tab),
      msg = sprintf("group_col '%s' not found in metrics_tab", group_col))

    # Compute per-group top-2 scores for sorting
    group_scores <- metrics_tab[!is.na(`Overall Score`),
      .(best = max(`Overall Score`),
        second = sort(`Overall Score`, decreasing = TRUE)[min(2, .N)]),
      by = c(group_col)]

    metrics_tab <- merge(metrics_tab, group_scores, by = group_col, all.x = TRUE)
    metrics_tab[, .group_label := paste0(group_col, "=", fifelse(is.na(get(group_col)), "NA", as.character(get(group_col))))]

    # Resolve group display order, allowing bare values or "col=val" prefixed ones
    all_labels <- unique(metrics_tab$.group_label)
    if (!is.null(group_col_order)) {
      specified <- ifelse(
        grepl("=", group_col_order, fixed = TRUE),
        group_col_order,
        paste0(group_col, "=", group_col_order)
      )
      all_labels <- c(specified[specified %in% all_labels], setdiff(all_labels, specified))
    }

    metrics_tab[, .group_rank := match(.group_label, all_labels)]
    metrics_tab <- metrics_tab[order(.group_rank, -`Overall Score`, -best, -second, na.last = TRUE)]

    # Apply top-n within each group after group-aware sorting.
    metrics_tab <- metrics_tab[, .SD[seq_len(min(.N, n_top))], by = .group_label]
    metrics_tab[, id := as.character(seq_len(.N))]

    row_info   <- metrics_tab[, .(id, group = .group_label)]
    row_groups <- data.table(group = all_labels, Group = all_labels)

    metrics_tab[, c(group_col, "best", "second", ".group_label", ".group_rank") := NULL]
  } else {
    metrics_tab <- metrics_tab[order(-`Overall Score`, na.last = TRUE)]
    metrics_tab <- metrics_tab[1:min(n_top, nrow(metrics_tab))]
    metrics_tab[, id := as.character(seq_len(.N))]
  }
} else {
  metrics_tab <- metrics_tab[1:min(n_top, nrow(metrics_tab))]
  metrics_tab[, id := as.character(seq_len(.N))]
}

# remove uninformative columns
columns_to_keep <- c(metrics, summary_metrics_overall, 'file_name', 'label', 'batch', 'Bio Conservation', 'Batch Correction', 'Overall Score')
columns <- names(metrics_tab)[
  sapply(metrics_tab, uniqueN) != 1 |
  names(metrics_tab) %in% columns_to_keep
]
metrics_tab <- metrics_tab[, ..columns]

get_col_width <- function(col, dt, factor=0.3) {
  width <- max(nchar(as.character(dt[[col]])), na.rm = TRUE)
  if (is.infinite(width) | width < 1) {
    width <- nchar(col)
  }
  return(factor * width + 1) # +1 for padding
}

#add column info metadata for plotting using funkyheatmap
dt1 <- data.table(
  id=integration_setup,
  group="Integration Setup",
  geom='text',
  palette='setup',
  width=sapply(integration_setup, get_col_width, dt=metrics_tab)
)
dt2 <- data.table(
  id="Overall Score",
  group="Overall",
  geom="bar",
  palette="overall",
  width=4
)
dt3 <- data.table(
  id="Batch Correction",
  group="Batch Correction",
  geom='bar',
  palette='batch',
  width=4
)
dt4 <- data.table(
  id=batch_metrics,
  name=batch_metrics,
  group="Batch Correction",
  geom='circle',
  palette='batch',
  width=1
)
dt5 <- data.table(
  id="Bio Conservation",
  group="Bio Conservation",
  geom='bar',
  palette='bio',
  width=4
)
dt6 <- data.table(
  id=bio_metrics,
  name=bio_metrics,
  group="Bio Conservation",
  geom='funkyrect',
  palette='bio',
  width=1
)
column_info <- rbind(dt1, dt2, dt3, dt4, dt5, dt6, fill=TRUE)
column_info <- column_info[column_info$id %in% colnames(metrics_tab)]
column_info$name <- column_info$id
print("column_info")
print(column_info)

n_top <- min(n_top, nrow(metrics_tab))

# Define shared plot configuration
plot_column_groups <- data.table(
  group = c("Integration Setup", "Overall", "Batch Correction", "Bio Conservation"),
  palette = c('setup', 'overall', 'batch', 'bio'),
  level1 = c("Integration Setup", "Overall", "Batch Correction", "Bio Conservation")
)

plot_palettes <- list(
  setup = rev(brewer.pal(5, "Greys")),
  overall = rev(brewer.pal(5, "Greens")),
  bio = rev(brewer.pal(5, "YlOrBr")),
  batch = rev(brewer.pal(5, "Blues"))
)

# Filter palettes to only include those referenced in column_info
used_palettes <- unique(column_info$palette[!is.na(column_info$palette)])
plot_palettes <- plot_palettes[names(plot_palettes) %in% used_palettes]

plot_legends <- list(
  list(title = "Batch Correction", palette = "batch", geom = "circle"),
  list(title = "Bio Conservation", palette = "bio", geom = "funkyrect")
)

# Adapt label spacing to annotation length to avoid overlaps.
max_col_name_len <- max(nchar(as.character(column_info$name)), na.rm = TRUE)
has_row_groups <- !is.null(row_groups) && nrow(row_groups) > 0
plot_position_args <- position_arguments(
  col_annot_offset = max(3, min(10, ceiling(max_col_name_len * 0.4))),
  col_annot_angle = 90
)

# Extract dataset name if available
title_text <- NULL
if ('dataset' %in% colnames(metrics_tab)) {
  datasets <- unique(na.omit(metrics_tab$dataset))
  if (length(datasets) > 0) {
    title_text <- paste("Dataset:", paste(datasets, collapse=", "))
    cat('Title text:', title_text, '\n')
  }
}

order_metric_block <- function(cols, dt) {
  cols <- cols[cols %in% colnames(dt)]
  if (length(cols) <= 1) {
    return(cols)
  }

  col_stats <- data.table(
    id = cols,
    na_count = sapply(cols, function(col) sum(is.na(dt[[col]]))),
    rank_score = sapply(cols, function(col) {
      non_na_rows <- which(!is.na(dt[[col]]))
      if (length(non_na_rows) == 0) {
        return(Inf)
      }
      mean(non_na_rows)
    }),
    name_len = nchar(cols)
  )
  setorder(col_stats, na_count, rank_score, name_len, id)
  col_stats$id
}

# Sort full-plot columns by completeness, then by availability in higher-ranked rows.
## Consolidate user ordering with automatic ordering (allow partial ordering).
ordered_batch_metrics <- NULL
ordered_bio_metrics <- NULL
if (!is.null(metric_order_batch)) {
  specified <- metric_order_batch[metric_order_batch %in% batch_metrics]
  remaining <- setdiff(batch_metrics, specified)
  remaining_ordered <- if (length(remaining) > 0) order_metric_block(remaining, metrics_tab) else character(0)
  ordered_batch_metrics <- c(specified, remaining_ordered)
} else {
  ordered_batch_metrics <- order_metric_block(batch_metrics, metrics_tab)
}

if (!is.null(metric_order_bio)) {
  specified <- metric_order_bio[metric_order_bio %in% bio_metrics]
  remaining <- setdiff(bio_metrics, specified)
  remaining_ordered <- if (length(remaining) > 0) order_metric_block(remaining, metrics_tab) else character(0)
  ordered_bio_metrics <- c(specified, remaining_ordered)
} else {
  ordered_bio_metrics <- order_metric_block(bio_metrics, metrics_tab)
}

full_cols <- c(
  "id",
  intersect(integration_setup, colnames(metrics_tab)),
  intersect("Overall Score", colnames(metrics_tab)),
  intersect("Batch Correction", colnames(metrics_tab)),
  ordered_batch_metrics,
  intersect("Bio Conservation", colnames(metrics_tab)),
  ordered_bio_metrics
)
full_cols <- unique(full_cols[full_cols %in% c("id", column_info$id)])
column_info <- column_info[match(full_cols[full_cols != "id"], id)]
metrics_tab_full <- metrics_tab[, ..full_cols]

# Create full heatmap
g <- funky_heatmap(
  metrics_tab_full,
  row_info = row_info,
  row_groups = row_groups,
  column_info = column_info,
  column_groups = plot_column_groups,
  palettes = plot_palettes,
  legends = plot_legends,
  position_args = plot_position_args,
  scale_column = FALSE
)

if (!is.null(title_text)) g <- g + ggtitle(title_text)
ggsave(
  output_svg,
  g,
  width = g$width,
  height = g$height,
  device = svglite::svglite,
  dpi=dpi
)

# Create overall metrics heatmap (subset with only overall/bio/batch scores)
overall_cols <- c(
  integration_setup,
  "Overall Score",
  "Batch Correction",
  batch_metrics_overall,
  "Bio Conservation",
  bio_metrics_overall
)
overall_cols <- overall_cols[overall_cols %in% colnames(metrics_tab)]

summary_column_info <- rbind(
  data.table(id=batch_metrics_overall, group="Batch Correction", geom='circle', palette='batch', width=1),
  data.table(id=bio_metrics_overall, group="Bio Conservation", geom='funkyrect', palette='bio', width=1),
  fill=TRUE
)
column_info_overall <- rbind(column_info, summary_column_info, fill=TRUE)
column_info_overall <- unique(column_info_overall, by = "id")
overall_cols <- c(
  intersect(integration_setup, overall_cols),
  intersect("Overall Score", overall_cols),
  intersect("Batch Correction", overall_cols),
  order_metric_block(batch_metrics_overall, metrics_tab),
  intersect("Bio Conservation", overall_cols),
  order_metric_block(bio_metrics_overall, metrics_tab)
)
overall_cols <- overall_cols[overall_cols %in% colnames(metrics_tab)]
column_info_overall <- column_info_overall[id %in% overall_cols]
column_info_overall <- column_info_overall[match(overall_cols, id)]
column_info_overall <- column_info_overall[!is.na(id)]
column_info_overall$group <- ""
column_info_overall$name <- column_info_overall$id

g <- funky_heatmap(
  metrics_tab[, ..overall_cols],
  row_info = row_info,
  row_groups = row_groups,
  column_info = column_info_overall,
  palettes = plot_palettes,
  legends = plot_legends,
  position_args = plot_position_args,
  scale_column = FALSE
)

if (!is.null(title_text)) g <- g + ggtitle(title_text)
ggsave(
  output_svg_overall,
  g,
  width = g$width,
  height = g$height,
  device = svglite::svglite,
  dpi=dpi
)
