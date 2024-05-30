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

# set nirs file for analysis
phenom_type_ <- "NIRJ5" # NIRJ5  NIRJ-1
phenom_file <- "NIRJ5_merged_wave_ls_mean_df.csv" # NIRJ5_merged_wave_ls_mean_df.csv

# set paths and booleans
phenom_dir_path <- "../../data/phenomic_data/"
# output result path for phenom graphics
output_phenom_graphics_path <- paste0(
  "../../results/graphics/phenomic_graphics/",
  phenom_type_
)

# define color function for nirs (phenom) of phenotypes
use_rainbow_colors_ <- F
if (use_rainbow_colors_) {
  col_func_phenom <- colorRampPalette(c(
    "black", "blue", "red", "orange",
    "yellow", "green"
  ))
} else {
  col_func_phenom <- colorRampPalette(c(
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

# read phenom data
phenom_df <- t(as.data.frame(fread(paste0(
  phenom_dir_path,
  phenom_file
))))[-1, ]

id_ <- rownames(phenom_df)
phenom_df <- apply(phenom_df, 2, as.numeric)

# should umap and/or pca be performed ?
if (perform_umap_) {
  phenom_umap_ <- data.frame(
    umap(
      scale(phenom_df,
        center = T,
        scale = T
      ),
      n_components = 3,
      random_state = random_state_umap_,
      n_neighbors = n_neighbors_umap_,
      min_dist = min_dist_
    )[["layout"]]
  )

  # plot umap 2d with phenom as label
  phenom_umap_2d <- phenom_umap_[, 1:2]
  phenom_umap_2d$label <- id_
  labels_ <- unique(phenom_umap_2d$label)
  n_phenom <- length(labels_)

  # define colors for labels
  color_labels_ <- col_func_phenom(
    length(labels_)
  )[1:n_phenom]
  names(color_labels_) <- labels_

  # define titles and make umap plot
  fig_title_ <- paste0(
    "UMAP 2D plot for ",
    phenom_type_, " phenomic data"
  )
  label_title_ <- paste0(
    phenom_type_,
    " phenomic"
  )
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
  for (label_ in unique(phenom_umap_2d$label)) {
    data_subset <- phenom_umap_2d[phenom_umap_2d$label == label_, ]
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
  # fig_x_y
  saveWidget(fig_x_y,
    file = paste0(
      output_phenom_graphics_path,
      paste0(
        "/",
        phenom_type_,
        "_phenom_umap_2d.html"
      )
    )
  )

  # plot umap 3d with phenom as label
  phenom_umap_3d <- phenom_umap_[, 1:3]
  phenom_umap_3d$label <- id_
  labels_ <- unique(phenom_umap_3d$label)
  n_phenom <- length(labels_)

  # define colors for labels
  color_labels_ <- col_func_phenom(
    length(labels_)
  )[1:n_phenom]
  names(color_labels_) <- labels_

  # define titles and make umap plot
  fig_title_ <- paste0(
    "UMAP 3D plot for ",
    phenom_type_, " phenomic data"
  )
  label_title_ <- paste0(
    phenom_type_,
    " phenomic"
  )
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
  for (label_ in unique(phenom_umap_3d$label)) {
    data_subset <- phenom_umap_3d[phenom_umap_3d$label == label_, ]
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
  fig_x_y_z
  # fig_x_y_z
  saveWidget(fig_x_y_z,
    file = paste0(
      output_phenom_graphics_path,
      paste0(
        "/",
        phenom_type_,
        "_phenom_umap_3d.html"
      )
    )
  )
}
if (perform_pca_) {
  phenom_pca_obj_ <- pca(phenom_df,
    ncomp = 2,
    center = T, scale = T
  )
  phenom_pca_mat_ <- as.data.frame(phenom_pca_obj_$variates$X)
  phenom_pca_exp_var_ <- phenom_pca_obj_$prop_expl_var$X
  phenom_pca_cum_exp_var_ <- 100 * phenom_pca_obj_$cum.var
  # 100*phenom_pca_exp_var_[1:2]
  # plot(phenom_pca_cum_exp_var_)
  # plotIndiv(phenom_pca_obj_)

  # plot pca with phenom as label
  n_comp_plot_ <- 2
  phenom_pca_mat_ <- phenom_pca_mat_[, 1:n_comp_plot_]
  phenom_pca_mat_$label <- id_
  labels_ <- unique(phenom_pca_mat_$label)
  n_phenoms <- length(labels_)

  # define colors for labels
  color_labels_ <- col_func_phenom(
    length(labels_)
  )[1:n_phenom]
  names(color_labels_) <- labels_

  # define titles and make umap plot
  fig_title_ <-
    paste0(
      "PCA 2D plot for ",
      phenom_type_, " phenomic data"
    )
  label_title_ <- paste0(
    phenom_type_,
    " phenomic"
  )
  fig_x_y <- plot_ly(
    type = "scatter", mode = "markers"
  ) %>%
    layout(
      plot_bgcolor = "#e5ecf6",
      title = fig_title_,
      xaxis = list(title = paste0(
        names(phenom_pca_exp_var_)[1], ": ",
        signif(100 * as.numeric(phenom_pca_exp_var_)[1], 2), "%"
      )),
      yaxis = list(title = paste0(
        names(phenom_pca_exp_var_)[2], ": ",
        signif(100 * as.numeric(phenom_pca_exp_var_)[2], 2), "%"
      ))
    )
  # regroup by label
  for (label_ in unique(phenom_pca_mat_$label)) {
    data_subset <- phenom_pca_mat_[phenom_pca_mat_$label == label_, ]
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
      output_phenom_graphics_path,
      paste0(
        "/",
        phenom_type_,
        "_phenom_pca_2d.html"
      )
    )
  )
}
