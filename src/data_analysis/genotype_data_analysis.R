# script meant to analyse scab genotypic data
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(devtools)
install_other_requirements <- F
if (install_other_requirements) {
  install.packages("BiocManager")
  library(BiocManager)
  BiocManager::install("snpStats")
  BiocManager::install("mixOmicsTeam/mixOmics")
  py_install("umap-learn", pip = T, pip_ignore_installed = T)
}
library(mixOmics)
library(data.table)
library(plotly)
library(ggplot2)
library(umap)
library(dplyr)
library(htmlwidgets)
library(rstudioapi)
library(stringr)
# detect and set script path automatically, and source functions
setwd(dirname(getActiveDocumentContext()$path))
source("../functions.R")
# set options to increase memory
options(expressions = 5e5)
options(warn = -1)

# set paths and booleans
geno_dir_path <- "../../data/genotype_data/"
# output result path for genotype graphics
output_geno_graphics_path <- "../../results/graphics/genomic_graphics/"

# define color function for genotypes of phenotypes
use_rainbow_colors_ <- F
if (use_rainbow_colors_) {
  col_func_genotype <- colorRampPalette(c(
    "black", "blue", "red", "orange",
    "yellow", "green"
  ))
} else {
  col_func_genotype <- colorRampPalette(c(
    "black"
  ))
}

# perform umap and/or pca ?
perform_umap_ <- T
perform_pca_ <- T

# define umap training and plot parameters

# umap parameters, most sensitive ones
random_state_umap_ <- 15
n_neighbors_umap_ <- 15
min_dist_ <- 0.1

# read genotype data
geno_df <- t(as.data.frame(fread(paste0(
  geno_dir_path,
  "genotype_data_family_original.csv"
))))
geno_df <- as.data.frame(geno_df)

# get snp index and name slice from df
idx_snp_index_names <- match(
  c("snp_index", "snp_name"),
  rownames(geno_df)
)
snp_index_name_df <- remove_na_columns(
  geno_df[idx_snp_index_names, ]
)

# drop slice and assign new colnames
geno_df <- geno_df[-idx_snp_index_names, ]
new_rownames_ <- paste0(
  rownames(geno_df), "/",
  geno_df[, 1]
)
rownames(geno_df) <- new_rownames_
geno_df <- geno_df[, -1]
colnames(geno_df) <- snp_index_name_df["snp_name", ]

# replace patterns with values in data frame, note equal length is assumed
# first replacement
first_vect_pattern_ <- c(
  "AA", "GG", "AG", "GA", "CC", "AC", "CA", "GC",
  "CG", "TT", "AT", "TA", "CT", "TC", "GT", "TG"
)
first_vect_replace_val_ <- c(
  "AA", "BB", "AB", "BA", "BB", "AB", "BA", "BB",
  "BB", "AA", "AA", "AA", "BA", "AB", "BA", "AB"
)
geno_df <- replace_pattern(
  geno_df,
  first_vect_pattern_,
  first_vect_replace_val_
)
# second replacement
second_vect_pattern_ <- c(
  "AA", "AB", "BA", "BB"
)
second_vect_replace_val_ <- c(
  "0", "1", "1", "2"
)
geno_df <- replace_pattern(
  geno_df,
  second_vect_pattern_,
  second_vect_replace_val_
)

# convert to numeric
geno_df <- apply(geno_df, 2, as.numeric)
rownames(geno_df) <- new_rownames_

# get an appreciation of missing value rate
100 * sum(is.na(geno_df)) / (nrow(geno_df) * ncol(geno_df))

# impute with column mean for a low missing rate, i.e. don't bother with a more
# sophisticated method such a missForest (which is an amazing method !)
geno_df <- apply(geno_df, 2, impute_mean)

# remove monomorphic markers and convert to data frame
geno_df <- as.data.frame(
  remove_monomorphic_markers(geno_df)$filtered_df
)
dim(geno_df)

# save geno_df to data directory
fwrite(geno_df,
  file = paste0(
    geno_dir_path,
    "genotype_data_family.csv"
  ),
  row.names = T
)

# should umap and/or pca be performed ?
if (perform_umap_) {
  geno_umap_2d <- data.frame(
    umap(
      scale(geno_df,
        center = T,
        scale = T
      ),
      n_components = 2,
      random_state = random_state_umap_,
      n_neighbors = n_neighbors_umap_,
      min_dist = min_dist_
    )[["layout"]]
  )

  # plot umap with genotype as label
  geno_umap_2d$label <- rownames(geno_df)
  labels_ <- unique(geno_umap_2d$label)
  n_genotypes <- length(labels_)

  # define colors for labels
  color_labels_ <- col_func_genotype(
    length(labels_)
  )[1:n_genotypes]
  names(color_labels_) <- labels_

  # define titles and make umap plot
  fig_title_ <-
    "UMAP 2D plot for genotype data"
  label_title_ <- "Genotype"
  fig_x_y <- plot_ly(
    type = "scatter", mode = "markers"
  ) %>%
    layout(
      plot_bgcolor = "#e5ecf6",
      title = fig_title_,
      xaxis = list(title = "first component"),
      yaxis = list(title = "second component")
    )
  # regroup by label
  for (label_ in unique(geno_umap_2d$label)) {
    data_subset <- geno_umap_2d[geno_umap_2d$label == label_, ]
    fig_x_y <- fig_x_y %>%
      add_trace(
        data = data_subset,
        x = ~X1, y = ~X2,
        type = "scatter", mode = "markers",
        marker = list(color = color_labels_[label_]),
        name = label_
      )
  }
  fig_x_y <- fig_x_y %>% layout(
    legend = list(title = list(text = label_title_))
  )
  saveWidget(fig_x_y,
    file = paste0(
      output_geno_graphics_path,
      "genotype_umap_2d.html"
    )
  )
}
if (perform_pca_) {
  geno_pca_obj_ <- pca(geno_df,
    ncomp = 200,
    center = T, scale = T
  )
  geno_pca_mat_ <- as.data.frame(geno_pca_obj_$variates$X)
  geno_pca_exp_var_ <- geno_pca_obj_$prop_expl_var$X
  geno_pca_cum_exp_var_ <- 100 * geno_pca_obj_$cum.var
  # 100*geno_pca_exp_var_[1:2]
  # plot(geno_pca_cum_exp_var_)
  # plotIndiv(geno_pca_obj_)

  # plot pca with genotype as label
  n_comp_plot_ <- 2
  geno_pca_mat_ <- geno_pca_mat_[, 1:n_comp_plot_]
  geno_pca_mat_$label <- rownames(geno_df)
  labels_ <- unique(geno_pca_mat_$label)
  n_genotypes <- length(labels_)

  # define colors for labels
  color_labels_ <- col_func_genotype(
    length(labels_)
  )[1:n_genotypes]
  names(color_labels_) <- labels_

  # define titles and make umap plot
  fig_title_ <-
    "PCA 2D plot for genotype data"
  label_title_ <- "Genotype"
  fig_x_y <- plot_ly(
    type = "scatter", mode = "markers"
  ) %>%
    layout(
      plot_bgcolor = "#e5ecf6",
      title = fig_title_,
      xaxis = list(title = paste0(
        names(geno_pca_exp_var_)[1], ": ",
        signif(100 * as.numeric(geno_pca_exp_var_)[1], 2), "%"
      )),
      yaxis = list(title = paste0(
        names(geno_pca_exp_var_)[2], ": ",
        signif(100 * as.numeric(geno_pca_exp_var_)[2], 2), "%"
      ))
    )
  # regroup by label
  for (label_ in unique(geno_pca_mat_$label)) {
    data_subset <- geno_pca_mat_[geno_pca_mat_$label == label_, ]
    fig_x_y <- fig_x_y %>%
      add_trace(
        data = data_subset,
        x = ~PC1, y = ~PC2,
        type = "scatter", mode = "markers",
        marker = list(color = color_labels_[label_]),
        name = label_
      )
  }
  fig_x_y <- fig_x_y %>% layout(
    legend = list(title = list(text = label_title_))
  )
  saveWidget(fig_x_y,
    file = paste0(
      output_geno_graphics_path,
      "genotype_pca_2d.html"
    )
  )
}
