suppressMessages({
  library(tidyverse)
  library(rjson)      # tree JSON
  library(scales)
  library(ggpubr)
  library(argparse)
  library(ggbeeswarm)
  library(yaml)
  library(corrplot)
})

# -------------------------
# helpers
# -------------------------
relevel_if_present <- function(f, level, after = Inf) {
  # Only relevel if `level` exists, otherwise return unchanged factor
  f <- as.factor(f)
  if (level %in% levels(f)) {
    return(forcats::fct_relevel(f, level, after = after))
  }
  f
}

stop_if_missing <- function(df, cols) {
  missing <- setdiff(cols, colnames(df))
  if (length(missing) > 0) {
    stop("Missing required columns in clades TSV: ", paste(missing, collapse = ", "))
  }
}

# -------------------------
# args
# -------------------------
parser <- ArgumentParser()
parser$add_argument("--tree", required = TRUE)
parser$add_argument("--clades", required = TRUE)
parser$add_argument("--out", required = TRUE)
parser$add_argument("--grid", required = TRUE)
parser$add_argument("--config", required = TRUE, help = "Path to Snakemake config.yaml")
parser$add_argument("--dir-out", required = TRUE)
parser$add_argument("--color", default = "True")
parser$add_argument("--virus", required = TRUE, help = "Virus key in config.yaml (e.g. ev-d68, eva71, cva16, cva10)")
xargs <- parser$parse_args()

output_dir <- xargs$dir_out

# -------------------------
# input
# -------------------------
clades <- tibble(read.csv(xargs$clades, sep = "\t", header = TRUE, check.names = FALSE))
new_tree <- fromJSON(file = xargs$tree)

# YAML: incomplete final line warning is from the file itself (no trailing newline).
# Best fix: add newline at end of config.yaml.
# If you want to avoid the warning in code, you can read file content first:
cfg <- yaml::read_yaml(xargs$config)

# Your YAML structure: top-level keys are viruses
sugg_params <- cfg$viruses[[xargs$virus]]

if (is.null(sugg_params)) {
  stop(sprintf("Virus '%s' not found at top-level of %s", xargs$virus, xargs$config))
}

# -------------------------
# configurable plot behavior (defaults)
# -------------------------
target_subtrees <- sugg_params$target_subtrees %||% NA_real_  # optional in config.yaml
y_limits <- sugg_params$y_limits %||% c(-15, 100)            # optional in config.yaml

# parameter columns expected in clades TSV
param_cols <- c(
  "cutoff",
  "divergence_addition",
  "divergence_scale",
  "min_size",
  "bushiness_branch_scale",
  "branch_length_scale"
)

stop_if_missing(clades, c("new_clade", "old_clade", param_cols))

# Convert params to factors for violin plots
cd <- clades %>%
  mutate(across(all_of(param_cols), as.factor))

# nicer x-axis labels
pretty_param <- function(x) str_to_title(str_replace_all(x, "_", " "))

# -------------------------
# Violin plots (filled by param level)
# -------------------------
make_violin_plot <- function(df, param) {
  p <- ggplot(df, aes(x = .data[[param]], y = .data[["new_clade"]], fill = .data[[param]])) +
    geom_violin(trim = FALSE, color = "black") +
    geom_boxplot(width = 0.06, color = "black", fill = "white", alpha = 0.5, linewidth = 0.6) +
    scale_fill_brewer(palette = "Paired") +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16)
    ) +
    labs(x = pretty_param(param), y = NULL) +
    coord_cartesian(ylim = y_limits)

  if (!is.na(target_subtrees)) {
    p <- p + geom_hline(yintercept = target_subtrees, linetype = "dashed", color = "black")
  }
  p
}

params_order <- c(
  "branch_length_scale", "bushiness_branch_scale",
  "divergence_scale", "cutoff",
  "divergence_addition", "min_size"
)

violin_plots <- setNames(lapply(params_order, \(p) make_violin_plot(cd, p)), params_order)

gga1 <- ggarrange(plotlist = violin_plots, ncol = 3, nrow = 2, align = "hv")
combined1 <- annotate_figure(
  gga1,
  left = grid::textGrob("Number of subtrees\nmedian [IQR]", rot = 90, gp = grid::gpar(fontsize = 18))
)

ggsave(xargs$out, plot = combined1, width = 14, height = 8, dpi = 300, bg = "white")

# -------------------------
# Violin plots + quasirandom points colored by cutoff
# (only makes sense if cutoff exists)
# -------------------------
make_violin_colored <- function(df, param) {
  ggplot(df, aes(x = .data[[param]], y = .data[["new_clade"]])) +
    geom_violin(trim = FALSE, color = "black") +
    geom_quasirandom(aes(color = .data[["cutoff"]]), size = 0.6, groupOnX = TRUE, alpha = 0.8) +
    geom_boxplot(width = 0.06, color = "black", fill = "white", alpha = 0.8, linewidth = 0.6) +
    scale_color_viridis_d(option = "C", begin = 0.2, end = 0.9, direction = -1, name = "Cutoff") +
    theme_classic() +
    theme(
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16)
    ) +
    labs(x = pretty_param(param), y = NULL) +
    coord_cartesian(ylim = y_limits) +
    guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
}

violin_plots_colored <- setNames(lapply(params_order, \(p) make_violin_colored(cd, p)), params_order)

gga2 <- ggarrange(
  plotlist = violin_plots_colored,
  common.legend = TRUE, legend = "right",
  ncol = 3, nrow = 2,
  align = "hv"
)
combined2 <- annotate_figure(
  gga2,
  left = grid::textGrob("Number of subtrees\nmedian [IQR]", rot = 90, gp = grid::gpar(fontsize = 18))
)

ggsave(file.path(output_dir, "violin_plots_cutoff.tiff"),
       plot = combined2, width = 12, height = 8, dpi = 300, bg = "white")

# -------------------------
# Cutoff vs median new_clade by each other parameter (step plots)
# -------------------------
serotype_count <- dplyr::n_distinct(clades$old_clade)

step_params <- c("branch_length_scale", "bushiness_branch_scale", "divergence_addition", "divergence_scale", "min_size")
plot_list <- list()

for (param in step_params) {
  clades_summarised <- clades %>%
    group_by(.data[["cutoff"]], .data[[param]]) %>%
    summarise(new_clade = median(.data[["new_clade"]], na.rm = TRUE), .groups = "drop")

  p <- ggplot(clades_summarised, aes(x = .data[["cutoff"]], color = as.factor(.data[[param]]))) +
    geom_step(aes(y = .data[["new_clade"]]), linewidth = 0.8) +
    geom_hline(yintercept = serotype_count, linetype = "dashed", color = "royalblue", linewidth = 0.6) +
    annotate("text",
             x = max(clades$cutoff, na.rm = TRUE) * 1.02, y = serotype_count,
             label = serotype_count, color = "royalblue", size = 5) +
    scale_color_brewer(palette = "Blues") +
    scale_y_continuous(name = "Number of subtrees [median]") +
    scale_x_continuous(name = "Cutoff value") +
    theme_classic() +
    guides(color = guide_legend(title = pretty_param(param))) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y = element_text(color = "royalblue"),
      axis.ticks.y = element_line(color = "royalblue"),
      axis.line.y = element_line(color = "royalblue"),
      axis.line = element_line(linewidth = 1),
      axis.ticks = element_line(linewidth = 0.9),
      axis.ticks.length = unit(2, "mm"),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    coord_cartesian(ylim = c(-1, max(60, max(clades$new_clade, na.rm = TRUE) + 1)))

  plot_list[[param]] <- p
}

p_all <- ggarrange(plotlist = plot_list, common.legend = FALSE, ncol = 3, nrow = 2, align = "hv")
p_all <- annotate_figure(
  p_all,
  left = text_grob("Number of subtrees [median]", color = "royalblue", rot = 90, size = 18)
)

ggsave(filename = file.path(output_dir, "cutoff_param_grid.tiff"),
       plot = p_all, width = 15, height = 5, dpi = 300, bg = "white")

# -------------------------
# Clade table plot (labels from tree JSON)
# -------------------------
extract_leaf_clades <- function(tree, clade_key = "clade_membership") {
  clade_map <- list()

  walk <- function(node) {
    is_leaf <- !is.null(node$node_attrs$strain)
    if (is_leaf) {
      strain <- node$name
      clade <- node$node_attrs[[clade_key]]$value
      clade_map[[strain]] <<- clade
    }
    if (!is.null(node$children)) {
      lapply(node$children, walk)
    }
  }

  walk(tree$tree)
  clade_map
}

old_map <- extract_leaf_clades(new_tree, "clade_membership")
new_map <- extract_leaf_clades(new_tree, "new-clade")

labels <- tibble(sequence = names(old_map), old_clade = unlist(old_map)) %>%
  inner_join(tibble(sequence = names(new_map), new_clade = unlist(new_map)), by = "sequence") %>%
  transmute(
    tip.label = sequence,
    serotype = old_clade,
    subtree = new_clade
  ) %>%
  mutate(
    # only relevel if present (fixes your warnings)
    serotype = relevel_if_present(serotype, "pre-ABC", after = Inf),
    subtree   = relevel_if_present(subtree,   "pre-ABC", after = Inf)
  ) %>%
  arrange(subtree)

tab <- table(labels$subtree, labels$serotype)
x <- tab / apply(tab, 1, sum)

# order rows/cols by weighted average (your original logic)
xval <- seq_len(ncol(x))
xsum <- apply(x, 1, function(row) sum(xval * row))
io <- order(xsum)
jo <- order(xval)

# column colors (recycle if more columns than palette)
col_scheme <- c("#3f4dcb", "#4681c9", "#5aa4a8", "#78b67e", "#9ebe5a",
                "#c5b945", "#e0a23a", "#e67231", "#dc2f24")
col_scheme <- rep(col_scheme, length.out = ncol(x))

out_grid <- if (tolower(xargs$color) == "false") {
  xargs$grid
} else {
  paste0(str_remove(xargs$grid, "\\.tiff$"), "_colored.tiff")
}

tiff(out_grid, width = 1400, height = 1400, res = 300, bg = "white")
par(mar = c(5, 5, 1, 2), font.lab = 2)

xlab_default <- "Existing clades"
ylab_default <- "Suggested subtrees"
xlab <- sugg_params$table_xlab %||% xlab_default
ylab <- sugg_params$table_ylab %||% ylab_default

shim <- 0.4
plot(NA,
     xlim = c(shim, ncol(x) - shim),
     ylim = c(shim, nrow(x) - shim),
     xaxt = "n", yaxt = "n",
     xlab = xlab, ylab = ylab, bty = "n")

if (tolower(xargs$color) == "false") {
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(ncol(x))) {
      xx <- 1 - x[io[i], jo[j]]
      if (xx < 1) xx <- min(0.9, xx)
      rect(j - 1, i - 1, j, i, col = rgb(xx, xx, xx), border = NA)
      count <- tab[io[i], jo[j]]
      if (count > 0 && xx > 0.1) text(j - 0.5, i - 0.5, label = count, cex = 0.5)
    }
  }
} else {
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(ncol(x))) {
      xx <- 1 - x[io[i], jo[j]]
      if (xx < 1) {
        rect(j - 1, i - 1, j, i, col = col_scheme[jo[j]], border = NA)
      }
    }
  }
}

for (i in 0:(nrow(x) - 1)) {
  abline(v = i, col = "grey85")
  abline(h = i, col = "grey85")
}

axis(side = 1, at = seq_len(ncol(x)) - 0.5, label = colnames(x), cex.axis = 0.8, las = 2)
axis(side = 2, at = seq_len(nrow(x)) - 0.5, label = rownames(x), cex.axis = 0.6, las = 2)
axis(side = 4, at = seq_len(nrow(x)) - 0.5, label = apply(tab[io, , drop = FALSE], 1, sum),
     cex.axis = 0.6, las = 2, lwd = 0, line = -0.8)

# annotate chosen parameters (works for any virus)
x_pos <- grconvertX(0.0, from = "nfc", to = "user")
y_pos <- grconvertY(0.01, from = "nfc", to = "user")
param_lines <- c(
  paste0("virus: ", xargs$virus),
  vapply(params_order, function(k) {
    v <- sugg_params[[k]]
    if (is.null(v)) return(NA_character_)
    paste0(k, ": ", v)
  }, character(1))
)
param_lines <- param_lines[!is.na(param_lines)]

par(mar = c(9, 5, 1, 2), font.lab = 2)

par(xpd = NA)  # allow drawing in the margin
mtext(
  text = paste(param_lines, collapse = "\n"),
  side = 1,        # bottom margin
  line = 6.5,      # vertical position within bottom margin (tune this)
  at = 0,     # anchor to left edge of plot region
  adj = 0,         # left align
  cex = 0.55,
  col = "grey40"
)
par(xpd = FALSE)

invisible(dev.off())

# -------------------------
# Correlation matrix (parameter effects)
# -------------------------
corr_cols <- intersect(
  c("new_clade", "cutoff", "divergence_addition", "divergence_scale",
    "min_size", "bushiness_branch_scale", "branch_length_scale", "score"),
  colnames(clades)
)

corr_df <- clades %>%
  select(all_of(corr_cols)) %>%
  mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.x)))))

# Drop columns with zero variance (or all NA)
sd_vec <- sapply(corr_df, sd, na.rm = TRUE)
keep <- is.finite(sd_vec) & sd_vec > 0
corr_df2 <- corr_df[, keep, drop = FALSE]

if (ncol(corr_df2) >= 2) {
  corr_mat <- cor(corr_df2, use = "pairwise.complete.obs")

  # If any NA remain, clustering will fail -> fall back to original order
  use_hclust <- !anyNA(corr_mat)

  tiff(file.path(output_dir, "parameter_effects_correlation_matrix.tiff"),
       width = 1600, height = 1400, res = 300, bg = "white")
  corrplot(
    corr_mat,
    method = "color",
    type = "upper",
    order = if (use_hclust) "hclust" else "original",
    tl.cex = 0.8,
    tl.col = "black"
  )
  invisible(dev.off())
} else {
  message("Skipping correlation plot: fewer than 2 non-constant numeric columns available.")
}