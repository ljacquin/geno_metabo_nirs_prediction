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

# set metabolome file for analysis
metabo_type_ <- "consti"            # consti, b04, BCZ14_phefileB2, BCZ14_phefileB3
metabo_file <- "CONSTI_phefile.csv"  # CONSTI_phefile.csv, B04_LsMeans.csv, BCZ14_phefileB2.csv, BCZ14_phefileB3.csv

# set paths and booleans
metabo_dir_path <- "../../data/metabolome_data/"
# output result path for metabo graphics
output_metabo_graphics_path <- paste0(
  "../../results/graphics/metabolome_graphics/",
  metabo_type_
)

# define color function for metabos of phenotypes
use_rainbow_colors_ <- F
if (use_rainbow_colors_) {
  col_func_metabo <- colorRampPalette(c(
    "black", "blue", "red", "orange",
    "yellow", "green"
  ))
} else {
  col_func_metabo <- colorRampPalette(c(
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

# read metabo data
metabo_df <- as.data.frame(fread(paste0(
  metabo_dir_path,
  metabo_file
)))

metabo_df <- remove_na_columns(metabo_df)
id_ <- metabo_df$id
metabo_df <- apply(metabo_df[, -1], 2, as.numeric)

# should umap and/or pca be performed ?
if (perform_umap_) {
  metabo_umap_2d <- data.frame(
    umap(
      scale(metabo_df,
        center = T,
        scale = T
      ),
      n_components = 2,
      random_state = random_state_umap_,
      n_neighbors = n_neighbors_umap_,
      min_dist = min_dist_
    )[["layout"]]
  )

  # plot umap with metabo as label
  metabo_umap_2d$label <- id_
  labels_ <- unique(metabo_umap_2d$label)
  n_metabos <- length(labels_)

  # define colors for labels
  color_labels_ <- col_func_metabo(
    length(labels_)
  )[1:n_metabos]
  names(color_labels_) <- labels_

  # define titles and make umap plot
  fig_title_ <- paste0(
    "UMAP 2D plot for ",
    metabo_type_, " metabolome data"
  )
  label_title_ <- paste0(
    metabo_type_,
    " metabolome"
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
  for (label_ in unique(metabo_umap_2d$label)) {
    data_subset <- metabo_umap_2d[metabo_umap_2d$label == label_, ]
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
      output_metabo_graphics_path,
      paste0("/",
        metabo_type_,
        "_metabo_umap_2d.html"
      )
    )
  )
}
if (perform_pca_) {
  metabo_pca_obj_ <- pca(metabo_df,
    ncomp = 2,
    center = T, scale = T
  )
  metabo_pca_mat_ <- as.data.frame(metabo_pca_obj_$variates$X)
  metabo_pca_exp_var_ <- metabo_pca_obj_$prop_expl_var$X
  metabo_pca_cum_exp_var_ <- 100 * metabo_pca_obj_$cum.var
  # 100*metabo_pca_exp_var_[1:2]
  # plot(metabo_pca_cum_exp_var_)
  # plotIndiv(metabo_pca_obj_)

  # plot pca with metabo as label
  n_comp_plot_ <- 2
  metabo_pca_mat_ <- metabo_pca_mat_[, 1:n_comp_plot_]
  metabo_pca_mat_$label <- id_
  labels_ <- unique(metabo_pca_mat_$label)
  n_metabos <- length(labels_)

  # define colors for labels
  color_labels_ <- col_func_metabo(
    length(labels_)
  )[1:n_metabos]
  names(color_labels_) <- labels_

  # define titles and make umap plot
  fig_title_ <-
    paste0(
      "PCA 2D plot for ",
      metabo_type_, " metabolome data"
    )
  label_title_ <- paste0(
    metabo_type_,
    " metabolome"
  )
  fig_x_y <- plot_ly(
    type = "scatter", mode = "markers"
  ) %>%
    layout(
      plot_bgcolor = "#e5ecf6",
      title = fig_title_,
      xaxis = list(title = paste0(
        names(metabo_pca_exp_var_)[1], ": ",
        signif(100 * as.numeric(metabo_pca_exp_var_)[1], 2), "%"
      )),
      yaxis = list(title = paste0(
        names(metabo_pca_exp_var_)[2], ": ",
        signif(100 * as.numeric(metabo_pca_exp_var_)[2], 2), "%"
      ))
    )
  # regroup by label
  for (label_ in unique(metabo_pca_mat_$label)) {
    data_subset <- metabo_pca_mat_[metabo_pca_mat_$label == label_, ]
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
      output_metabo_graphics_path,
      paste0("/",
        metabo_type_,
        "_metabo_pca_2d.html"
      )
    )
  )
}
