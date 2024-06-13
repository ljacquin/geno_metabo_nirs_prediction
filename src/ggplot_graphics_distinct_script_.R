library(data.table)
library(plotly)
library(htmlwidgets)
library(ggplot2)
library(tidyr)
library(dplyr)
library(plotly)
library(rstudioapi)
# install.packages("remotes")
# remotes::install_github("trevorld/ggpattern")
library(ggpattern)
setwd(dirname(getActiveDocumentContext()$path))

# set path for reading results and saving graphics
genomic_prediction_result_path <- "../results/genomic_prediction/"
metabolomic_prediction_result_path <- "../results/metabolomic_prediction/"
phenomic_prediction_result_path <- "../results/phenomic_prediction/"
output_pred_graphics_path <- "../results/manuscript_graphics/"

# file suffixe for bloc and position
file_suffix <- ""

# use spectra data corrected for bloc and position ?
use_bloc_pos_correction_ <- F
if (use_bloc_pos_correction_) {
  file_suffix <- "_bloc_pos_corrected_for_nirs"
}

# get genomic prediction results
res_b04_18 <- as.data.frame(fread(paste0(
  genomic_prediction_result_path,
  "genomic_pred_results_8347_SNP_BLUPaudpc.eau.EUB04.18.csv"
)))[, "GBLUP"]

res_b04_19 <- as.data.frame(fread(paste0(
  genomic_prediction_result_path,
  "genomic_pred_results_8347_SNP_BLUPaudpc.eau.EUB04.19.csv"
)))[, "GBLUP"]

res_bcz14_19 <- as.data.frame(fread(paste0(
  genomic_prediction_result_path,
  "genomic_pred_results_8347_SNP_BLUPaudpc.eau.09BCZ14.prin19.csv"
)))[, "GBLUP"]

df_res_gblup_ <- as.data.frame(cbind(
  res_b04_18,
  res_b04_19,
  res_bcz14_19
))
n_row_ <- nrow(df_res_gblup_)
df_res_gblup_ <- df_res_gblup_[-c(n_row_ - 1, n_row_), ]
colnames(df_res_gblup_) <- c(
  "AUDPC B04 2018 / GBLUP",
  "AUDPC B04 2019 / GBLUP",
  "AUDPC BCZ14 / GBLUP"
)

# get metabolomic prediction results

# consti
res_b04_18 <- fread(paste0(
  metabolomic_prediction_result_path, "consti/",
  "metabolomic_pred_results_1576_metabolites_BLUPaudpc.eau.EUB04.18.csv"
))[, "MBLUP"]

res_b04_19 <- fread(paste0(
  metabolomic_prediction_result_path, "consti/",
  "metabolomic_pred_results_1576_metabolites_BLUPaudpc.eau.EUB04.19.csv"
))[, "MBLUP"]

res_bcz14_19 <- fread(paste0(
  metabolomic_prediction_result_path, "consti/",
  "metabolomic_pred_results_1576_metabolites_BLUPaudpc.eau.09BCZ14.prin19.csv"
))[, "MBLUP"]

df_res_mblup_consti_ <- as.data.frame(cbind(
  res_b04_18,
  res_b04_19,
  res_bcz14_19
))
n_row_ <- nrow(df_res_mblup_consti_)
df_res_mblup_consti_ <- df_res_mblup_consti_[-c(n_row_ - 1, n_row_), ]
colnames(df_res_mblup_consti_) <- c(
  "AUDPC B04 2018 / MBLUP based on CONSTI",
  "AUDPC B04 2019 / MBLUP based on CONSTI",
  "AUDPC BCZ14 / MBLUP based on CONSTI"
)

# b04
res_b04_18 <- fread(paste0(
  metabolomic_prediction_result_path, "b04/",
  "metabolomic_pred_results_2102_metabolites_BLUPaudpc.eau.EUB04.18.csv"
))[, "MBLUP"]

res_b04_19 <- fread(paste0(
  metabolomic_prediction_result_path, "b04/",
  "metabolomic_pred_results_2102_metabolites_BLUPaudpc.eau.EUB04.19.csv"
))[, "MBLUP"]

# bcz14_b2
res_bcz14_b2 <- as.data.frame(fread(paste0(
  metabolomic_prediction_result_path, "bcz14_b2/",
  "metabolomic_pred_results_2955_metabolites_BLUPaudpc.eau.09BCZ14.prin19.csv"
)))[, "MBLUP"]

# bcz14_b3
res_bcz14_b3 <- as.data.frame(fread(paste0(
  metabolomic_prediction_result_path, "bcz14_b3/",
  "metabolomic_pred_results_2955_metabolites_BLUPaudpc.eau.09BCZ14.prin19.csv"
)))[, "MBLUP"]

df_res_mblup_ <- as.data.frame(cbind(
  res_b04_18,
  res_b04_19,
  res_bcz14_b2,
  res_bcz14_b3
))
n_row_ <- nrow(df_res_mblup_)
df_res_mblup_ <- df_res_mblup_[-c(n_row_ - 1, n_row_), ]
colnames(df_res_mblup_) <- c(
  "AUDPC B04 2018 / MBLUP based on B04 ls-means",
  "AUDPC B04 2019 / MBLUP based on B04 ls-means",
  "AUDPC BCZ14 / MBLUP based on BCZ14 B2",
  "AUDPC BCZ14 / MBLUP based on BCZ14 B3"
)

if (use_bloc_pos_correction_) {
  # get phenomic prediction results
  nirs_type_ <- "NIRJ-1"

  res_b04_18 <- as.data.frame(fread(paste0(
    paste0(phenomic_prediction_result_path, nirs_type_, "/"),
    "phenomic_pred_results_2151_wave_lengths_bloc_pos_corrected_BLUPaudpc.eau.EUB04.18.csv"
  )))[, "PBLUP"]

  res_b04_19 <- as.data.frame(fread(paste0(
    paste0(phenomic_prediction_result_path, nirs_type_, "/"),
    "phenomic_pred_results_2151_wave_lengths_bloc_pos_corrected_BLUPaudpc.eau.EUB04.19.csv"
  )))[, "PBLUP"]

  res_bcz14_19 <- as.data.frame(fread(paste0(
    paste0(phenomic_prediction_result_path, nirs_type_, "/"),
    "phenomic_pred_results_2151_wave_lengths_bloc_pos_corrected_BLUPaudpc.eau.09BCZ14.prin19.csv"
  )))[, "PBLUP"]

  df_res_pblup_nirj_1 <- as.data.frame(cbind(
    res_b04_18,
    res_b04_19,
    res_bcz14_19
  ))
  n_row_ <- nrow(df_res_pblup_nirj_1)
  df_res_pblup_nirj_1 <- df_res_pblup_nirj_1[-c(n_row_ - 1, n_row_), ]
  colnames(df_res_pblup_nirj_1) <- c(
    "AUDPC B04 2018 / PBLUP based on NIRJ-1 wave length ls-means",
    "AUDPC B04 2019 / PBLUP based on NIRJ-1 wave length ls-means",
    "AUDPC BCZ14 / PBLUP based on NIRJ-1 wave length ls-means"
  )

  nirs_type_ <- "NIRJ5"

  res_b04_18 <- as.data.frame(fread(paste0(
    paste0(phenomic_prediction_result_path, nirs_type_, "/"),
    "phenomic_pred_results_2151_wave_lengths_bloc_pos_corrected_BLUPaudpc.eau.EUB04.18.csv"
  )))[, "PBLUP"]

  res_b04_19 <- as.data.frame(fread(paste0(
    paste0(phenomic_prediction_result_path, nirs_type_, "/"),
    "phenomic_pred_results_2151_wave_lengths_bloc_pos_corrected_BLUPaudpc.eau.EUB04.19.csv"
  )))[, "PBLUP"]

  res_bcz14_19 <- as.data.frame(fread(paste0(
    paste0(phenomic_prediction_result_path, nirs_type_, "/"),
    "phenomic_pred_results_2151_wave_lengths_bloc_pos_corrected_BLUPaudpc.eau.09BCZ14.prin19.csv"
  )))[, "PBLUP"]

  df_res_pblup_nirj5 <- as.data.frame(cbind(
    res_b04_18,
    res_b04_19,
    res_bcz14_19
  ))
  n_row_ <- nrow(df_res_pblup_nirj5)
  df_res_pblup_nirj5 <- df_res_pblup_nirj5[-c(n_row_ - 1, n_row_), ]
  colnames(df_res_pblup_nirj5) <- c(
    "AUDPC B04 2018 / PBLUP based on NIRJ5 wave length ls-means",
    "AUDPC B04 2019 / PBLUP based on NIRJ5 wave length ls-means",
    "AUDPC BCZ14 / PBLUP based on NIRJ5 wave length ls-means"
  )
} else {
  # get phenomic prediction results
  nirs_type_ <- "NIRJ-1"

  res_b04_18 <- as.data.frame(fread(paste0(
    paste0(phenomic_prediction_result_path, nirs_type_, "/"),
    "phenomic_pred_results_2151_wave_lengths_BLUPaudpc.eau.EUB04.18.csv"
  )))[, "PBLUP"]

  res_b04_19 <- as.data.frame(fread(paste0(
    paste0(phenomic_prediction_result_path, nirs_type_, "/"),
    "phenomic_pred_results_2151_wave_lengths_BLUPaudpc.eau.EUB04.19.csv"
  )))[, "PBLUP"]

  res_bcz14_19 <- as.data.frame(fread(paste0(
    paste0(phenomic_prediction_result_path, nirs_type_, "/"),
    "phenomic_pred_results_2151_wave_lengths_BLUPaudpc.eau.09BCZ14.prin19.csv"
  )))[, "PBLUP"]

  df_res_pblup_nirj_1 <- as.data.frame(cbind(
    res_b04_18,
    res_b04_19,
    res_bcz14_19
  ))
  n_row_ <- nrow(df_res_pblup_nirj_1)
  df_res_pblup_nirj_1 <- df_res_pblup_nirj_1[-c(n_row_ - 1, n_row_), ]
  colnames(df_res_pblup_nirj_1) <- c(
    "AUDPC B04 2018 / PBLUP based on NIRJ-1 wave length ls-means",
    "AUDPC B04 2019 / PBLUP based on NIRJ-1 wave length ls-means",
    "AUDPC BCZ14 / PBLUP based on NIRJ-1 wave length ls-means"
  )

  nirs_type_ <- "NIRJ5"

  res_b04_18 <- as.data.frame(fread(paste0(
    paste0(phenomic_prediction_result_path, nirs_type_, "/"),
    "phenomic_pred_results_2151_wave_lengths_BLUPaudpc.eau.EUB04.18.csv"
  )))[, "PBLUP"]

  res_b04_19 <- as.data.frame(fread(paste0(
    paste0(phenomic_prediction_result_path, nirs_type_, "/"),
    "phenomic_pred_results_2151_wave_lengths_BLUPaudpc.eau.EUB04.19.csv"
  )))[, "PBLUP"]

  res_bcz14_19 <- as.data.frame(fread(paste0(
    paste0(phenomic_prediction_result_path, nirs_type_, "/"),
    "phenomic_pred_results_2151_wave_lengths_BLUPaudpc.eau.09BCZ14.prin19.csv"
  )))[, "PBLUP"]

  df_res_pblup_nirj5 <- as.data.frame(cbind(
    res_b04_18,
    res_b04_19,
    res_bcz14_19
  ))
  n_row_ <- nrow(df_res_pblup_nirj5)
  df_res_pblup_nirj5 <- df_res_pblup_nirj5[-c(n_row_ - 1, n_row_), ]
  colnames(df_res_pblup_nirj5) <- c(
    "AUDPC B04 2018 / PBLUP based on NIRJ5 wave length ls-means",
    "AUDPC B04 2019 / PBLUP based on NIRJ5 wave length ls-means",
    "AUDPC BCZ14 / PBLUP based on NIRJ5 wave length ls-means"
  )
}

# add a column method to the data frame
df_res_gblup_$Method <- "GBLUP"
df_res_mblup_$Method <- "MBLUP"
df_res_mblup_consti_$Method <- "MBLUP_CONSTI"
df_res_pblup_nirj_1$Method <- "PBLUP_NIRJ1"
df_res_pblup_nirj5$Method <- "PBLUP_NIRJ5"

# combine data frame into a long format
df_combined <- bind_rows(
  df_res_gblup_ %>% pivot_longer(-Method,
    names_to = "Variable",
    values_to = "Value"
  ),
  df_res_mblup_ %>% pivot_longer(-Method,
    names_to = "Variable",
    values_to = "Value"
  ),
  df_res_mblup_consti_ %>% pivot_longer(-Method,
    names_to = "Variable",
    values_to = "Value"
  ),
  df_res_pblup_nirj_1 %>% pivot_longer(-Method,
    names_to = "Variable",
    values_to = "Value"
  ),
  df_res_pblup_nirj5 %>% pivot_longer(-Method,
    names_to = "Variable",
    values_to = "Value"
  )
)

# define colors for each method
colors <- c(
  "GBLUP" = "green", "MBLUP" = "blue", "MBLUP_CONSTI" = "red",
  "PBLUP_NIRJ1" = "orange", "PBLUP_NIRJ5" = "purple"
)

# create boxplots with ggplot2
ggplot_obj_ <-
  ggplot(df_combined, aes(x = Variable, y = Value, fill = Method)) +
  geom_boxplot() +
  labs(
    title = "Boxplots of PA for differents methods based on different datasets (i.e. genomic, metabolomic and phenomic) for different AUDPC",
    x = "",
    y = "Predictive ability (PA)",
    fill = "Method"
  ) +
  scale_fill_manual(values = colors) + # set colors manually
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme_bw() +
  ylim(-1, 1) +
  facet_grid(. ~ Method, scales = "free_x", space = "free_x") + # add facet_grid
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(t = 20)),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )
ggsave(
  filename = paste0(
    output_pred_graphics_path, "genomic_metabolomic_phenomic_prediction",
    file_suffix,
    ".png"
  ),
  width = 20,
  height = 8,
  limitsize = F,
  plot = ggplot_obj_
)

# convert ggplot to plotly object and save html
ggplotly_obj_ <- ggplotly(ggplot_obj_)

saveWidget(ggplotly_obj_, file = paste0(
  output_pred_graphics_path,
  "genomic_metabolomic_phenomic_prediction",
  file_suffix,
  ".html"
))
