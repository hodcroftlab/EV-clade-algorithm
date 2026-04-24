library(jsonlite)
library(tidyverse)
library(corrplot)

# setwd("//wsl.localhost/SwissTPH-Ubuntu-24.04/home/neunna/workspace/influenza-clade-nomenclature/clade-suggestion-algorithm")

# Recursively count tips and mutations and get median phylo_score (bushiness)
count_tree_stats <- function(node) {
  nt_muts <- length(node$branch_attrs$mutations$nuc %||% list())
  aa_muts <- 0
  if (!is.null(node$branch_attrs$mutations)) {
    # Sum all non-nuc mutation lists
    aa_keys <- setdiff(names(node$branch_attrs$mutations), "nuc")
    for (key in aa_keys) {
      aa_muts <- aa_muts + length(node$branch_attrs$mutations[[key]])
    }
  }
  # get median phylo_score (bushiness)
  phylo_score <- node$branch_attrs$phylo_score %||% NA
  
  n_tips <- 0
  nt_total <- nt_muts
  aa_total <- aa_muts
  phylo_median <- phylo_score
  
  if (!is.null(node$children)) {
    for (child in node$children) {
      stats <- count_tree_stats(child)
      n_tips <- n_tips + stats$n_tips
      nt_total <- nt_total + stats$nt_muts
      aa_total <- aa_total + stats$aa_muts
      phylo_median <- median(phylo_score, stats$phylo_score, na.rm = TRUE)
    }
  } else {
    n_tips <- 1
  }
  
  
  list(n_tips=n_tips, nt_muts=nt_total, aa_muts=aa_total, phylo_score=phylo_median)
}

# Helper for missing fields
`%||%` <- function(x, y) if (is.null(x)) y else x

# ---- MAIN ----
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) < 1) stop("usage: Rscript count_mutations_tree.R tree.json")
# tree_json <- fromJSON(args[1], simplifyVector = FALSE)


trees_paths <-  c(
  "C:/Users/neunna/Downloads/seasonal-flu_h1n1pdm_ha_2y.json",
  "C:/Users/neunna/Downloads/seasonal-flu_h3n2_ha_2y.json",
  "C:/Users/neunna/Downloads/seasonal-flu_h1n1pdm_na_2y.json",
  "C:/Users/neunna/Downloads/seasonal-flu_h3n2_na_2y.json",
  "C:/Users/neunna/Downloads/seasonal-flu_vic_ha_2y.json",
  "C:/Users/neunna/Downloads/seasonal-flu_vic_na_2y.json",
  "C:/Users/neunna/Downloads/enterovirus_d68_genome (2).json"
)

df= data.frame(
  path = NA,
  build = NA,
  n_tips = NA,
  nt_muts = NA,
  aa_muts = NA,
  phylo_score = NA,
  stringsAsFactors = FALSE
)

for (path in trees_paths) {
  tree_json <- fromJSON(path, simplifyVector = FALSE)
  root <- tree_json$tree
  stats <- count_tree_stats(root)
  cat(sprintf("Number of tips: %d\nTotal nucleotide mutations: %d\nTotal amino acid mutations: %d\n",
              stats$n_tips, stats$nt_muts, stats$aa_muts))
  
  name <- gsub(".*Downloads/(.*)\\.json", "\\1", path)
  
  df = add_row(df,
               path = path,
               build = name,
               n_tips = stats$n_tips,
               nt_muts = stats$nt_muts,
               aa_muts = stats$aa_muts,
               phylo_score = stats$phylo_score)
  
}



# params matrix
`seasonal_A-H3N2_HA` <- c("cutoff"= 1.0,
                          "bushiness_branch_scale"= 1,
                          "divergence_scale"= 8,
                          "branch_length_scale"= 4,
                          "divergence_addition"= 1.0,
                          "min_size"= 30)

`seasonal_A-H1N1pdm_HA` <- c("cutoff"= 1.0,
                             "bushiness_branch_scale"= 1,
                             "divergence_scale"= 8,
                             "branch_length_scale"= 4,
                             "divergence_addition"= 1.0,
                             "min_size"= 20)

`seasonal_A-H1N1pdm_NA` <-  c("cutoff"= 1.0,
                              "bushiness_branch_scale"= 1,
                              "divergence_scale"= 8,
                              "branch_length_scale"= 4,
                              "divergence_addition"= 1.0,
                              "min_size"= 30)

`seasonal_A-H3N2_NA` <- c("cutoff"= 1.2,
                          "bushiness_branch_scale"= 1,
                          "divergence_scale"= 10,
                          "branch_length_scale"= 6,
                          "divergence_addition"= 1.0,
                          "min_size"= 30)

`seasonal_B-Vic_HA` <-  c("cutoff"= 1.0,
                          "bushiness_branch_scale"= 1,
                          "divergence_scale"= 8,
                          "branch_length_scale"= 4,
                          "divergence_addition"= 1.0,
                          "min_size"= 20)

`seasonal_B-Vic_NA` <- c("cutoff"= 1.0,
                         "bushiness_branch_scale"= 1,
                         "divergence_scale"= 8,
                         "branch_length_scale"= 4,
                         "divergence_addition"= 1.0,
                         "min_size"= 20)

`enterovirus_d68_genome` <- c("cutoff"= 1.2,
                              "bushiness_branch_scale"= 9.0,
                              "divergence_scale"= 5.0,
                              "branch_length_scale"= 3,
                              "divergence_addition"= 0.5, 
                              "min_size"= 15)



df_merged = df %>% drop_na(path) %>%
  mutate(params = case_when(
    build == "seasonal-flu_h3n2_ha_2y" ~ list(`seasonal_A-H3N2_HA`),
    build == "seasonal-flu_h1n1pdm_ha_2y" ~ list(`seasonal_A-H1N1pdm_HA`),
    build == "seasonal-flu_h1n1pdm_na_2y" ~ list(`seasonal_A-H1N1pdm_NA`),
    build == "seasonal-flu_h3n2_na_2y" ~ list(`seasonal_A-H3N2_NA`),
    build == "seasonal-flu_vic_ha_2y" ~ list(`seasonal_B-Vic_HA`),
    build == "seasonal-flu_vic_na_2y" ~ list(`seasonal_B-Vic_NA`),
    build == "enterovirus_d68_genome (2)" ~ list(`enterovirus_d68_genome`)
  )) %>%
  unnest_wider(params)

df_merged = df_merged %>% 
  mutate(build = str_replace_all(pattern = "seasonal-flu_|_2y| (2)", replacement = "", build))


# plot the dependence between n_tips and bushiness_branch_scale¨
df_merged %>% 
  ggplot(aes(x = bushiness_branch_scale, y = n_tips)) +
  geom_point() +
  scale_y_log10()+
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Dependence of Number of Tips on Bushiness Branch Scale",
       x = "Bushiness Branch Scale",
       y = "Number of Tips")


# plot the dependence between nt_muts & aa_muts and divergence_scale branch_length_scale
df_merged %>% 
  ggplot(aes(x = divergence_scale, y = nt_muts)) +
  geom_point() +
  scale_y_log10()+
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Dependence of Nucleotide Mutations on Divergence Scale",
       x = "Divergence Scale",
       y = "Nucleotide Mutations")

# plot the dependence between nt_muts & aa_muts and branch_length_scale
df_merged %>% 
  ggplot(aes(x = branch_length_scale, y = nt_muts)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_y_log10()+
  labs(title = "Dependence of Nucleotide Mutations on Branch Length Scale",
       x = "Branch Length Scale",
       y = "Nucleotide Mutations")

# Heatmap of correlations
df_merged %>%
  select(cutoff, divergence_addition, divergence_scale, min_size, 
         bushiness_branch_scale, branch_length_scale, n_tips, nt_muts,aa_muts) %>%
  cor(use = "complete.obs") %>% 
  replace(is.na(.), 0) %>%
  corrplot(., method = "color", type = "upper", 
           order = "hclust", tl.cex = 0.8, tl.col = "black")

