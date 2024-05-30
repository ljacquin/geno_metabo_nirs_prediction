library(data.table)
library(plotly)
library(htmlwidgets)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

# set path for reading results and saving graphics
genomic_prediction_result_path <- "../results/genomic_prediction/"
metabolomic_prediction_result_path <- "../results/metabolomic_prediction/"
phenomic_prediction_result_path <- "../results/phenomic_prediction/"
output_pred_graphics_path <- "../results/manuscript_graphics/"

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

df_res_ <- as.data.frame(cbind(
  res_b04_18,
  res_b04_19,
  res_bcz14_19
))
n_row_ <- nrow(df_res_)
df_res_ <- df_res_[-c(n_row_ - 1, n_row_), ]
colnames(df_res_) <- c(
  "AUDPC B04 2018",
  "AUDPC B04 2019",
  "AUDPC BCZ14"
)

# get conditions names
condition_names <- colnames(df_res_)
df_res_ <- as.data.frame(apply(
  df_res_,
  2, as.numeric
))

# initialize plot_ly boxplot graphic
boxplots_pa_ <- plot_ly()

# add boxplots
for (condition_ in condition_names) {
  boxplots_pa_ <- add_boxplot(
    boxplots_pa_,
    y = df_res_[[condition_]],
    name = condition_,
    boxpoints = "all",
    jitter = 0.3,
    pointpos = -1.8
  )
}
# add layout
boxplots_pa_ <- boxplots_pa_ %>%
  layout(
    title = paste0(
      "GBLUP PA distributions of conditions for scenarios based on 5-folds CV and 20 shuffling, using 8347 SNPs"
    ),
    yaxis = list(title = "Predictive ability (PA)", range = c(-1, 1)),
    legend = list(title = list(text = "Condition"))
  )

# save boxplots_pa_ graphics
saveWidget(boxplots_pa_, file = paste0(
  output_pred_graphics_path,
  "genomic_predictive_ability_for_conditions.html"
))


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

df_res_ <- cbind(
  res_b04_18,
  res_b04_19,
  res_bcz14_19
)
n_row_ <- nrow(df_res_)
df_res_ <- df_res_[-c(n_row_ - 1, n_row_), ]
colnames(df_res_) <- c(
  "AUDPC B04 2018",
  "AUDPC B04 2019",
  "AUDPC BCZ14"
)

# get conditions names
condition_names <- colnames(df_res_)
df_res_ <- as.data.frame(apply(
  df_res_,
  2, as.numeric
))

# initialize plot_ly boxplot graphic
boxplots_pa_ <- plot_ly()

# add boxplots
for (condition_ in condition_names) {
  boxplots_pa_ <- add_boxplot(
    boxplots_pa_,
    y = df_res_[[condition_]],
    name = condition_,
    boxpoints = "all",
    jitter = 0.3,
    pointpos = -1.8
  )
}
# add layout
boxplots_pa_ <- boxplots_pa_ %>%
  layout(
    title = paste0(
      "MBLUP PA distributions of conditions for scenarios based on 5-folds CV and 20 shuffling, based on 1576 metabolites from CONSTI"
    ),
    yaxis = list(title = "Predictive ability (PA)", range = c(-1, 1)),
    legend = list(title = list(text = "Condition"))
  )

# save boxplots_pa_ graphics
saveWidget(boxplots_pa_, file = paste0(
  output_pred_graphics_path,
  "consti_metabolic_predictive_ability_for_conditions.html"
))


# b04

res_b04_18 <- fread(paste0(
  metabolomic_prediction_result_path, "b04/",
  "metabolomic_pred_results_2102_metabolites_BLUPaudpc.eau.EUB04.18.csv"
))[, "MBLUP"]

res_b04_19 <- fread(paste0(
  metabolomic_prediction_result_path, "b04/",
  "metabolomic_pred_results_2102_metabolites_BLUPaudpc.eau.EUB04.19.csv"
))[, "MBLUP"]

df_res_ <- cbind(
  res_b04_18,
  res_b04_19
)
n_row_ <- nrow(df_res_)
df_res_ <- df_res_[-c(n_row_ - 1, n_row_), ]
colnames(df_res_) <- c(
  "AUDPC B04 2018",
  "AUDPC B04 2019"
)

# get conditions names
condition_names <- colnames(df_res_)
df_res_ <- as.data.frame(apply(
  df_res_,
  2, as.numeric
))

# initialize plot_ly boxplot graphic
boxplots_pa_ <- plot_ly()

# add boxplots
for (condition_ in condition_names) {
  boxplots_pa_ <- add_boxplot(
    boxplots_pa_,
    y = df_res_[[condition_]],
    name = condition_,
    boxpoints = "all",
    jitter = 0.3,
    pointpos = -1.8
  )
}
# add layout
boxplots_pa_ <- boxplots_pa_ %>%
  layout(
    title = paste0(
      "MBLUP PA distributions of conditions for scenarios based on 5-folds CV and 20 shuffling, using 2102 metabolites from B04"
    ),
    yaxis = list(title = "Predictive ability (PA)", range = c(-1, 1)),
    legend = list(title = list(text = "Condition"))
  )

# save boxplots_pa_ graphics
saveWidget(boxplots_pa_, file = paste0(
  output_pred_graphics_path,
  "b04_metabolic_predictive_ability_for_conditions.html"
))

# bcz14_b2

res_bcz14_b2 <- as.data.frame(fread(paste0(
  metabolomic_prediction_result_path, "bcz14_b2/",
  "metabolomic_pred_results_2955_metabolites_BLUPaudpc.eau.09BCZ14.prin19.csv"
)))[, "MBLUP"]
res_bcz14_b2 <- res_bcz14_b2[-c(
  length(res_bcz14_b2) - 1,
  length(res_bcz14_b2)
)]

# initialize plot_ly boxplot graphic
boxplots_pa_ <- plot_ly()

# add boxplots
boxplots_pa_ <- add_boxplot(
  boxplots_pa_,
  y = res_bcz14_b2,
  name = "AUDPC BCZ14",
  boxpoints = "all",
  jitter = 0.3,
  pointpos = -1.8
)

# add layout
boxplots_pa_ <- boxplots_pa_ %>%
  layout(
    title = paste0(
      "MBLUP PA distributions of condition for scenarios based on 5-folds CV and 20 shuffling, using 2955 metabolites from BCZ14_B2"
    ),
    yaxis = list(title = "Predictive ability (PA)", range = c(-1, 1)),
    legend = list(title = list(text = "Condition"))
  )

# save boxplots_pa_ graphics
saveWidget(boxplots_pa_, file = paste0(
  output_pred_graphics_path,
  "bcz14_b2_metabolic_predictive_ability_for_conditions.html"
))


# bcz14_b3

res_bcz14_b3 <- as.data.frame(fread(paste0(
  metabolomic_prediction_result_path, "bcz14_b3/",
  "metabolomic_pred_results_2955_metabolites_BLUPaudpc.eau.09BCZ14.prin19.csv"
)))[, "MBLUP"]
res_bcz14_b3 <- res_bcz14_b3[-c(
  length(res_bcz14_b3) - 1,
  length(res_bcz14_b3)
)]

# initialize plot_ly boxplot graphic
boxplots_pa_ <- plot_ly()

# add boxplots
boxplots_pa_ <- add_boxplot(
  boxplots_pa_,
  y = res_bcz14_b3,
  name = "AUDPC BCZ14",
  boxpoints = "all",
  jitter = 0.3,
  pointpos = -1.8
)

# add layout
boxplots_pa_ <- boxplots_pa_ %>%
  layout(
    title = paste0(
      "MBLUP PA distributions of condition for scenarios based on 5-folds CV and 20 shuffling, using 2955 metabolites from BCZ14_B3"
    ),
    yaxis = list(title = "Predictive ability (PA)", range = c(-1, 1)),
    legend = list(title = list(text = "Condition"))
  )

# save boxplots_pa_ graphics
saveWidget(boxplots_pa_, file = paste0(
  output_pred_graphics_path,
  "bcz14_b3_metabolic_predictive_ability_for_conditions.html"
))


# get phenomic prediction results
vect_nirs_type_ <- c("NIRJ-1", "NIRJ5")

for (nirs_type_ in vect_nirs_type_) {
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

  df_res_ <- as.data.frame(cbind(
    res_b04_18,
    res_b04_19,
    res_bcz14_19
  ))
  n_row_ <- nrow(df_res_)
  df_res_ <- df_res_[-c(n_row_ - 1, n_row_), ]
  colnames(df_res_) <- c(
    "AUDPC B04 2018",
    "AUDPC B04 2019",
    "AUDPC BCZ14"
  )

  # get conditions names
  condition_names <- colnames(df_res_)
  df_res_ <- as.data.frame(apply(
    df_res_,
    2, as.numeric
  ))

  # initialize plot_ly boxplot graphic
  boxplots_pa_ <- plot_ly()

  # add boxplots
  for (condition_ in condition_names) {
    boxplots_pa_ <- add_boxplot(
      boxplots_pa_,
      y = df_res_[[condition_]],
      name = condition_,
      boxpoints = "all",
      jitter = 0.3,
      pointpos = -1.8
    )
  }
  # add layout
  boxplots_pa_ <- boxplots_pa_ %>%
    layout(
      title = paste0(
        "PBLUP PA distributions of conditions for scenarios based on 5-folds CV and 20 shuffling, using ",
        nirs_type_, " 2151 wave lengths"
      ),
      yaxis = list(title = "Predictive ability (PA)", range = c(-1, 1)),
      legend = list(title = list(text = "Condition"))
    )

  # save boxplots_pa_ graphics
  saveWidget(boxplots_pa_, file = paste0(
    output_pred_graphics_path,
    paste0(
      nirs_type_,
      "_phenomic_predictive_ability_for_conditions.html"
    )
  ))
}
