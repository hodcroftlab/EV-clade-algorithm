suppressMessages({
  library(tidyverse)
  library(rjson)
  library(scales)
  library(ggpubr)
  library(car)
  library(plotly)
  library(argparse)
  library(ggbeeswarm)
})

parser <- ArgumentParser()

parser$add_argument('--tree')
parser$add_argument('--clades')
parser$add_argument('--out')
parser$add_argument('--grid')
parser$add_argument('--params', default = "config/suggestion_params.json")
parser$add_argument('--dir-out')
parser$add_argument('--color')

xargs<- parser$parse_args()


output_dir <- xargs$dir_out
clades  <- tibble(read.csv(xargs$clades, sep="\t", header=T))
new_tree <- fromJSON(file=xargs$tree)
sugg_params <- fromJSON(file=xargs$params)

## melt into long format and use geom_smooth to compare effects
clades_long <- clades %>%
  pivot_longer(cols = c(cutoff, divergence_addition, divergence_scale,
                        min_size, bushiness_branch_scale, branch_length_scale),
               names_to = "parameter", values_to = "value")

# ggplot(clades_long, aes(x = value, y = new_clade)) +
#   geom_point(alpha = 0.3) +
#   geom_smooth(method = "loess", se = FALSE, color = "blue") +
#   facet_wrap(~parameter, scales = "free_x") +
#   ylim(min(clades$new_clade)-1,max(clades$new_clade)+1)+
#   labs(x = "Parameter value", y = "Number of subtrees",
#        title = "Parameter influence on new clade count") +
#   theme_minimal()

####### Violin plots #########
# Violin plots summarizing the number of subtrees per clade for each parameter
cd <- clades %>% dplyr::mutate(
  cutoff = as.factor(cutoff),
  divergence_addition = as.factor(divergence_addition),
  divergence_scale = as.factor(divergence_scale),
  min_size = as.factor(min_size),
  bushiness_branch_scale = as.factor(bushiness_branch_scale),
  branch_length_scale = as.factor(branch_length_scale)
)

violin_plots <- list()
params <- c("branch_length_scale","bushiness_branch_scale", "divergence_scale","cutoff", "divergence_addition", "min_size")

for (param in params) {
  violin_plots[[param]] <- (cd %>%
                              ggplot(aes(x = .data[[param]], y = new_clade, fill = factor(.data[[param]]))) +
                              geom_violin(trim = FALSE,color = "black") +
                              geom_boxplot(width=0.04, color = "black", fill="white", alpha=0.5, size=0.8)+
                              scale_fill_brewer(palette = "Paired") +
                              theme_classic() +
                              geom_hline(yintercept = 9, linetype = "dashed", color = "black")+
                              theme(
                                legend.position = "none",
                                axis.text = element_text(size = 14),
                                axis.title = element_text(size = 16)) +
                              labs(x = str_to_title(str_replace_all(param,"_", " ")), y = NULL)+
                              coord_cartesian(ylim = c(-15, 100)))
}

# combine plots
gga <- ggarrange(
  plotlist = violin_plots,
  ncol = 3, nrow = 2,
  align = "hv",
  font.label = list(size = 14, color = "black", face = "bold")
)
  
combined_plot <- annotate_figure(gga,
                        left = grid::textGrob("Number of subtrees\nmedian [IQR]", rot = 90, gp = grid::gpar(fontsize = 18)))

combined_plot
                        
ggsave(xargs$out,
       plot = combined_plot, width = 14, height = 8, dpi = 300,bg = "white")


violin_plots_colored <- list()
params <- c("branch_length_scale","bushiness_branch_scale","cutoff", "divergence_addition", "divergence_scale", "min_size")

for (param in params) {
  violin_plots_colored[[param]] <- (cd %>%
                              ggplot(aes(x = .data[[param]], y = new_clade)) +
                              geom_violin(trim = FALSE,color = "black") +
                              geom_quasirandom(aes(color = factor(cutoff)), size = 0.5, groupOnX = TRUE) +
                              geom_boxplot(width=0.04, color = "black", fill="white", alpha=0.8)+
                              
                              scale_fill_brewer(palette = "Paired") +
                              scale_color_viridis_d(option = "C", begin = 0.2, end = 0.9, direction = -1)+
                              theme_classic() +
                              geom_hline(yintercept = 9, linetype = "dashed", color = "black")+
                              theme(
                                # legend.position = "none",
                                legend.text = element_text(size = 12),
                                legend.title = element_text(size = 14, face = "bold"),
                                axis.text = element_text(size = 14),
                                axis.title = element_text(size = 16)) +
                              guides(color = guide_legend(override.aes = list(size = 2, alpha = 1), 
                                                           title = "Cutoff")) +
                              labs(x = str_to_title(str_replace_all(param,"_", " ")), y = NULL)+
                              coord_cartesian(ylim = c(-15, 100)))
}

# combine plots
gga <- ggarrange(
  plotlist = violin_plots_colored,
  common.legend = TRUE, legend = "right",
  ncol = 3, nrow = 2,
  align = "hv",
  font.label = list(size = 14, color = "black", face = "bold")
)

combined_plot <- annotate_figure(gga,
                                 left = grid::textGrob("Number of subtrees\nmedian [IQR]", rot = 90, gp = grid::gpar(fontsize = 18)))

ggsave(paste0(output_dir,"/violin_plots_cutoff.tiff"), 
       plot = combined_plot, width = 12, height = 8, dpi = 300,bg = "white")


#### Scaling Parameters vs. score ####



# Poon's internal branch length plot
serotype_count = unique(clades$old_clade)
params <-  c("branch_length_scale","bushiness_branch_scale", "divergence_addition", "divergence_scale", "min_size")

# Store plots
plot_list <- list()

# Loop through parameters
for (param in params) {
  
  clades_summarised <- clades %>%
    group_by(cutoff, .data[[param]]) %>%
    summarise(new_clade = median(new_clade), .groups = "keep")
  
  p <- ggplot(clades_summarised, aes(x = cutoff, color = factor(.data[[param]]))) +
    geom_step(aes(y = new_clade), linewidth = 0.8) +
    geom_linerange(aes(xmin = min(clades$cutoff), xmax = max(clades$cutoff), y = serotype_count),
                   linetype = "dashed", color = "royalblue", linewidth = 0.6) +
    annotate("text", x = max(clades$cutoff) * 1.05, y = serotype_count,
             label = serotype_count, color = "royalblue", size = 5) +
    scale_color_brewer(palette = "Blues", aesthetics = c("color", "fill")) +
    scale_y_continuous(name = "Number of subtrees [median]") +
    scale_x_continuous(name = "Cutoff value", breaks = seq(0.4, 1.5, 0.2)) +
    theme_classic() +
    guides(color = guide_legend(title = str_to_title(str_replace_all(param,"_", " ")))) +
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
    ) + coord_cartesian(ylim = c(-1, 60))
  
  plot_list[[param]] <- p
}

# Arrange all plots in a grid
p_all <- ggarrange(plotlist = plot_list, common.legend = FALSE,
                   ncol = 3, nrow = 2, align = "hv",
                   font.label = list(size = 14, color = "black", face = "bold")
)
p_all <- annotate_figure(p_all,
                         left = text_grob("Number of subtrees [median]", color = "royalblue", rot = 90, size = 18))


# Save final plot
ggsave(filename = paste0(output_dir, "/cutoff_param_grid.tiff"),
       plot = p_all, width = 15, height = 5, dpi = 300, bg = "white")


## Table plot
# -- Recursive function to extract name → clade mapping --
extract_leaf_clades <- function(tree, clade_key = "clade_membership") {
  clade_map <- list()
  
  walk <- function(node) {
    is_leaf <- !is.null(node$node_attrs$strain)
    if (is_leaf) {
      strain <- node$name  # strain ID (e.g. "KX255397")
      clade <- node$node_attrs[[clade_key]]$value
      clade_map[[strain]] <<- clade
    }
    if (!is.null(node$children)) {
      lapply(node$children, walk)
    }
  }
  
  walk(tree$tree)
  return(clade_map)
}

# -- Extract old and new mappings --
old_map <- extract_leaf_clades(new_tree)
new_map <- extract_leaf_clades(new_tree, "new-clade") 

df_old <- tibble(sequence = names(old_map), old_clade = unlist(old_map))
df_new <- tibble(sequence = names(new_map), new_clade = unlist(new_map))

labels <- inner_join(df_old, df_new, by = "sequence")

# labels <- df_merged %>%
#   count(old_clade, new_clade, name = "count")

labels <- labels %>% rename(
  tip.label = sequence,
  serotype = old_clade,
  subtree = new_clade
) %>% 
  # mutate(serotype = factor(serotype, levels = c("A", "A1", "A2", "B", "B1", "B2", "B3", "C", "pre-ABC"))) %>% 
  mutate(subtree =  fct_relevel(as.factor(subtree), "C", after = Inf)) %>%
  mutate(serotype = fct_relevel(as.factor(serotype), "pre-ABC", after = Inf)) %>%
  mutate(subtree =  fct_relevel(subtree, "pre-ABC", after = Inf)) %>%
  arrange(subtree)

tab <- table(labels$subtree, labels$serotype)
x <- tab / apply(tab, 1, sum)
xval <- 1:length(colnames(x))  # not sure if it would be better to keep original names
xsum <- apply(x, 1, function(row) sum(xval*row))
io <- order(xsum)
jo <- order(xval)  # column index


col_scheme <- c("#3f4dcb", "#4681c9", "#5aa4a8", "#78b67e", "#9ebe5a",
                "#c5b945", "#e0a23a", "#e67231", "#dc2f24")


tiff(ifelse(xargs$color == "False", xargs$grid, paste0(str_remove(xargs$grid, ".tiff"), "_colored.tiff")), 
width=1400, height=1400, res=300, bg="white")
if (xargs$color == "False"){
  par(mar=c(5,5,1,2), font.lab=2)
  
  shim <- 0.4 #0.65 or 0.4
  plot(NA, xlim=c(shim, ncol(x)-shim), ylim=c(shim, nrow(x)-shim), 
       xaxt='n', yaxt='n',
       xlab="RIVM labels", ylab="Subtree", bty='n')
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      xx <- 1-x[io[i], jo[j]]
      if (xx < 1) xx <- min(0.9, xx) 
      rect(j-1, i-1, j, i, col=rgb(xx, xx, xx), border=NA)
      count <- tab[io[i], jo[j]]
      if (count>0 && xx>0.1) text(j-0.5, i-0.5, adj=0.5, label=count, cex=0.5)
    }
  }
  for (i in 0:nrow(x)-1) {
    abline(v=i, col='grey80')
    abline(h=i, col='grey80')
  }
  axis(side=1, at=1:ncol(x)-0.5, label=colnames(x), cex.axis=0.8, las=2)
  axis(side=2, at=1:nrow(x)-0.5, label=rownames(x), #label=paste("s", io-1, sep=""), 
       las=2, cex.axis=0.7)
  axis(side=4, at=1:nrow(x)-0.5, label=apply(tab[io,], 1, sum),
       cex.axis=0.6, las=2, lwd=0, line=-0.9)
  
  # Get the bottom-left of the full figure region in user coordinates
  x_pos <- grconvertX(0.01, from = "nfc", to = "user")
  y_pos <- grconvertY(0.01, from = "nfc", to = "user")
  
  # Add your label
  text(x = x_pos, y = y_pos,
      labels = paste('
       "cutoff":', sugg_params$cutoff, ',
       "divergence_addition":', sugg_params$divergence_addition, ', 
       "divergence_scale":', sugg_params$divergence_scale, ',
       "min_size":', sugg_params$min_size, ',
       "bushiness_branch_scale":', sugg_params$bushiness_branch_scale, ',
       "branch_length_scale":', sugg_params$branch_length_scale
      ),
       adj = c(0.0, 0), xpd = NA, cex = 0.6, col = "grey50")
  
} else {
  par(mar=c(5,5,1,2), font.lab=2)
  shim <- 0.4
  plot(NA, xlim=c(shim, ncol(x)-shim), ylim=c(shim, nrow(x)-shim), 
       xaxt='n', yaxt='n',
       xlab="EV-D68 clades (current)", ylab="Subtrees proposed by algorithm", bty='n')
  
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      xx <- 1 - x[io[i], jo[j]]
      
      if (xx < 1) {  # only draw colored box if match
        base_col <- col_scheme[j]
        rect(j-1, i-1, j, i, col=base_col, border=NA)
        
        count <- tab[io[i], jo[j]]
        if (count > 0) {
          # text(j-0.5, i-0.5, adj=0.5, label=count, cex=0.5)
        }
      }
    }
  }
  
  for (i in 0:nrow(x)-1) {
    abline(v=i, col='grey85')
    abline(h=i, col='grey85')
  }
  
  axis(side=1, at=1:ncol(x)-0.5, label=colnames(x), cex.axis=0.8, las=2)
  axis(side=2, at=1:nrow(x)-0.5, label=rownames(x), las=2, cex.axis=0.5)
  axis(side=4, at=1:nrow(x)-0.5, label=apply(tab[io,], 1, sum),
       cex.axis=0.6, las=2, lwd=0, line=-0.7)
  }
invisible(dev.off())


## Bushiness is weird - no effect at all?
# clades %>%
#   group_by(cutoff, bushiness_branch_scale, branch_length_scale) %>%
#   summarise(new_clade = median(new_clade), .groups = "keep")
# 
# clades_long %>% filter(
#   parameter == "bushiness_branch_scale" | parameter == "branch_length_scale" | parameter == "cutoff"
# ) %>% 
#   ggplot(aes(x=value, y=new_clade, color=parameter)) +
#   geom_point()

# interaction linear model
# lm_model <- lm(new_clade ~ cutoff * divergence_addition * divergence_scale * min_size * bushiness_branch_scale * branch_length_scale, data = clades)
# summary(lm_model)

# Check for multicollinearity

# vif_values <- vif(lm_model)
# 
# # Check residuals
# par(mfrow=c(2,2))
# plot(lm_model)


#### effect of divergence on leaf nodes ####
# Recursive function to collect divergence from leaves
# Recursive traversal to collect (node_name, delta_div)
get_delta_divs <- function(tree) {
  delta_list <- list()
  
  walk <- function(node, parent_div = 0) {
    node_div <- node$node_attrs$div %||% NA
    if (!is.null(node_div)) {
      delta <- node_div - parent_div
      delta_list[[node$name]] <<- list(
        name = node$name,
        div = node_div,
        delta_div = delta
      )
    }
    
    if (!is.null(node$children)) {
      for (child in node$children) {
        walk(child, node_div)
      }
    }
  }
  
  walk(tree$tree)
  bind_rows(delta_list)
}

# Run it
delta_df <- get_delta_divs(new_tree)


# Summary
mean_div <- mean(delta_df$delta_div, na.rm = TRUE)
median_div <- median(delta_df$delta_div, na.rm = TRUE)
summary(delta_df)

df <- expand.grid(
  # delta_d = delta_df$delta_div,
  delta_d = seq(0.1, 10, length.out = 100),
  d_add = c(0.55, 1),
  d_scale = c(8,9)
)

df$score <- with(df, (d_add / (delta_d + d_scale)) * delta_d)

ggplot(df, aes(x = delta_d, y = score, color = factor(d_scale), linetype = factor(d_add))) +
  geom_line(size = 1) +
  labs(
    x = expression(Delta*d),
    y = "Score contribution",
    color = "d_scale",
    linetype = "d_add"
  ) +
  theme_minimal(base_size = 14)
ggsave(paste0(output_dir,"/score_contribution_plot.png"), width = 10, height = 6, dpi = 300, bg = "white")


#### Scale vs score with nsubtrees colored
# clades %>%
#   ggplot(aes(x= cutoff , y = score, color = new_clade)) +
#   geom_jitter() +
#   theme_classic()
# 
# # 3D plot: cutoff, new_clade, score
# 
# p <- clades %>%
#   plot_ly(x = ~cutoff, y = ~new_clade, z = ~score, color = ~new_clade, colors = "Set1") %>%
#   add_markers() %>%
#   layout(scene = list(
#     xaxis = list(title = "Cutoff"),
#     yaxis = list(title = "# Subtrees"),
#     zaxis = list(title = "Score")
#   )) %>%
#   config(displayModeBar = FALSE)
# # Save the 3D plot as an HTML file
# htmlwidgets::saveWidget(as_widget(p), file = paste0(output_dir, "/3D_plot.html"), selfcontained = TRUE)



# Create a correlation matrix showing how parameter changes affect outcomes
param_effects <- clades %>%
  select(new_clade, cutoff, divergence_addition, min_size, 
         score) %>%
  cor(use = "complete.obs")

# Heatmap of correlations
library(corrplot)
tiff(paste0(output_dir, "/parameter_effects_correlation_matrix.tiff"), width=1600, height=1400, res=300, bg="white")
corrplot(param_effects, method = "color", type = "upper", 
         order = "hclust", tl.cex = 0.8, tl.col = "black")
invisible(dev.off())


## scales
score_components <- clades %>%
  # Assuming you can calculate these separately:
  mutate(
    bushiness_contrib = bushiness_branch_scale / (bushiness_branch_scale + branch_length_scale + divergence_scale ),
    branch_contrib = branch_length_scale / (bushiness_branch_scale + branch_length_scale + divergence_scale ),
    div_contrib = divergence_scale  / (bushiness_branch_scale + branch_length_scale + divergence_scale )
  ) %>%
  pivot_longer(cols = ends_with("_contrib"), 
               names_to = "component", values_to = "contribution")

ggplot(score_components, aes(x = cutoff, y = contribution, fill = component)) +
  geom_col(position = "stack") +
  facet_wrap(~new_clade, scales = "free_x") +
  labs(title = "Score Component Contributions by Subtree Count")



## clade count sensitivity to cutoff
ggplot(clades, aes(x = cutoff, y = new_clade, color = factor(min_size))) +
  geom_line(aes(group = interaction(min_size, bushiness_branch_scale, branch_length_scale)), alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, size = 1.5) +
  facet_wrap(~bushiness_branch_scale, labeller = label_both) +
  labs(title = "Clade Count Sensitivity to Cutoff",
       x = "Cutoff Threshold", y = "Number of New Clades",
       color = "Min Size") +
  theme_minimal()

