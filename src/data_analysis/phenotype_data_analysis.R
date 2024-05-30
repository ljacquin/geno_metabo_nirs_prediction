# script meant to analyse scab phenotypic data
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
pheno_dir_path <- "../../data/phenotype_data/"
# output result path for phenotype graphics
output_pheno_graphics_path <- "../../results/graphics/phenotype_graphics/"

# define color function for genotypes of phenotypes
use_rainbow_colors_ <- F
if (use_rainbow_colors_) {
  col_func_phenotype <- colorRampPalette(c(
    "black", "blue", "red", "orange",
    "yellow", "green"
  ))
} else {
  col_func_phenotype <- colorRampPalette(c(
    "black"
  ))
}

# perform umap and/or pca ?
perform_umap_ <- T
perform_pca_ <- T

# define umap training and plot parameters

# umap parameters, most sensitive ones
random_state_umap_ <- 15
n_neighbors_umap_ <- 30
min_dist_ <- 0.1

# read phenotype data
pheno_df <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "phenotype_data_family_original.csv"
)))
pheno_df <- as.data.frame(pheno_df)
pheno_df[pheno_df == "*"] <- NA

# convert to numeric
new_rownames_ <- pheno_df$id
pheno_df <- apply(
  pheno_df[, -1],
  2, as.numeric
)
pheno_df <- as.data.frame(apply(
  pheno_df,
  2, impute_mean
))
rownames(pheno_df) <- new_rownames_

# should umap and/or pca be performed ?
if (perform_umap_) {
  pheno_umap_df <- data.frame(
    umap(
      scale(pheno_df,
        center = T,
        scale = T
      ),
      n_components = 3,
      random_state = random_state_umap_,
      n_neighbors = n_neighbors_umap_,
      min_dist = min_dist_
    )[["layout"]]
  )

  # plot umap with phenotype as label
  pheno_umap_2d <- pheno_umap_df[, 1:2]
  pheno_umap_2d$label <- rownames(pheno_df)
  labels_ <- unique(pheno_umap_2d$label)
  n_genotypes <- length(labels_)

  # define colors for labels
  color_labels_ <- col_func_phenotype(
    length(labels_)
  )[1:n_genotypes]
  names(color_labels_) <- labels_

  # define titles and make umap plot
  fig_title_ <-
    "UMAP 2D plot for phenotype data"
  label_title_ <- "phenotype"
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
  for (label_ in unique(pheno_umap_2d$label)) {
    data_subset <- pheno_umap_2d[pheno_umap_2d$label == label_, ]
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
      output_pheno_graphics_path,
      "phenotype_umap_2d.html"
    )
  )

  # plot umap in 3d with phenotype as label
  pheno_umap_3d <- pheno_umap_df[, 1:3]
  pheno_umap_3d$label <- rownames(pheno_df)
  labels_ <- unique(pheno_umap_3d$label)
  n_genotypes <- length(labels_)

  # define colors for labels
  color_labels_ <- col_func_phenotype(
    length(labels_)
  )[1:n_genotypes]
  names(color_labels_) <- labels_

  # define titles and make umap plot
  fig_title_ <-
    "UMAP 3D plot for phenotype data"
  label_title_ <- "phenotype"
  fig_x_y_z <- plot_ly(
    type = "scatter3d", mode = "markers"
  ) %>%
    layout(
      plot_bgcolor = "#e5ecf6",
      title = fig_title_,
      xaxis = list(title = "first component"),
      yaxis = list(title = "second component"),
      zaxis = list(title = "third component")
    )
  # regroup by label
  for (label_ in unique(pheno_umap_3d$label)) {
    data_subset <- pheno_umap_3d[pheno_umap_3d$label == label_, ]
    fig_x_y_z <- fig_x_y_z %>%
      add_trace(
        data = data_subset,
        x = ~X1, y = ~X2, z = ~X3,
        type = "scatter3d", mode = "markers",
        marker = list(color = color_labels_[label_]),
        name = label_
      )
  }
  fig_x_y_z <- fig_x_y_z %>% layout(
    legend = list(title = list(text = label_title_))
  )
  saveWidget(fig_x_y_z,
    file = paste0(
      output_pheno_graphics_path,
      "phenotype_umap_3d.html"
    )
  )
}
if (perform_pca_) {
  pheno_pca_obj_ <- pca(pheno_df,
    ncomp = 10,
    center = T, scale = T
  )
  pheno_pca_mat_ <- as.data.frame(pheno_pca_obj_$variates$X)
  pheno_pca_exp_var_ <- pheno_pca_obj_$prop_expl_var$X
  pheno_pca_cum_exp_var_ <- 100 * pheno_pca_obj_$cum.var
  # 100*pheno_pca_exp_var_[1:2]
  # plot(pheno_pca_cum_exp_var_)
  # plotIndiv(pheno_pca_obj_)

  # plot pca with phenotype as label
  n_comp_plot_ <- 2
  pheno_pca_mat_ <- pheno_pca_mat_[, 1:n_comp_plot_]
  pheno_pca_mat_$label <- rownames(pheno_df)
  labels_ <- unique(pheno_pca_mat_$label)
  n_phenotypes <- length(labels_)

  # define colors for labels
  color_labels_ <- col_func_phenotype(
    length(labels_)
  )[1:n_phenotypes]
  names(color_labels_) <- labels_

  # define titles and make umap plot
  fig_title_ <-
    "PCA 2D plot for phenotype data"
  label_title_ <- "phenotype"
  fig_x_y <- plot_ly(
    type = "scatter", mode = "markers"
  ) %>%
    layout(
      plot_bgcolor = "#e5ecf6",
      title = fig_title_,
      xaxis = list(title = paste0(
        names(pheno_pca_exp_var_)[1], ": ",
        signif(100 * as.numeric(pheno_pca_exp_var_)[1], 2), "%"
      )),
      yaxis = list(title = paste0(
        names(pheno_pca_exp_var_)[2], ": ",
        signif(100 * as.numeric(pheno_pca_exp_var_)[2], 2), "%"
      ))
    )
  # regroup by label
  for (label_ in unique(pheno_pca_mat_$label)) {
    data_subset <- pheno_pca_mat_[pheno_pca_mat_$label == label_, ]
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
      output_pheno_graphics_path,
      "phenotype_pca_2d.html"
    )
  )
}
