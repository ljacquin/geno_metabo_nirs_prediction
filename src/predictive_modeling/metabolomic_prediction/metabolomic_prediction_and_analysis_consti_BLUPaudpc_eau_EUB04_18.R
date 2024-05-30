# script meant to perform metabolomic prediction and analyses for apple scab resistance
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(MASS)
library(data.table)
library(stringr)
library(doParallel)
library(doRNG)
library(robustbase)
library(foreach)
library(parallel)
library(Matrix)
library(rgl)
library(cvTools)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(dplyr)
library(KRMM)
library(kernlab)
library(glmnet)
library(ranger)
library(stringr)
library(devtools)
install_other_requirements <- F
if (install_other_requirements) {
  devtools::install_github("ljacquin/KRMM")
}
local_computation_ <- F
if (local_computation_) {
  library(rstudioapi)
  setwd(dirname(getActiveDocumentContext()$path))
}
source("../../functions.R")
options(warn = -1)

# define function(s) and package(s) to export for parallelization
pkgs_to_export_ <- c(
  "ranger",
  "kernlab",
  "KRMM",
  "glmnet",
  "foreach",
  "cvTools"
)
# set input paths
metabo_dir_path <- "../../../data/metabolome_data/"
pheno_dir_path <- "../../../data/phenotype_data/"

# output result path for metabolome graphics
output_pred_results_path <- "../../../results/metabolomic_prediction/"
output_pred_graphics_path <- "../../../results/graphics/metabolomic_prediction_graphics/"

# define selected traits and the analyzed one defined by trait_
selected_traits_ <- c(
  "BLUPaudpc.eau.EUB04.18",
  "BLUPaudpc.eau.EUB04.19",
  "BLUPaudpc.eau.09BCZ14.prin19"
)
trait_ <- "BLUPaudpc.eau.EUB04.18"

# define metabo_type_ for analysis
metabo_type_ <- "consti" # consti, b04, BCZ14_phefileB2, BCZ14_phefileB3
metabo_file <- "CONSTI_phefile.csv" # CONSTI_phefile.csv, B04_LsMeans.csv, BCZ14_phefileB2.csv, BCZ14_phefileB3.csv

# define shift seed value by
mult_seed_by_ <- 100

# set k for K-folds cv
k_folds_ <- 5

# define number of shuffles
n_shuff_ <- 20

# define color function for metabo
col_func_metabo <- colorRampPalette(c(
  "black", "blue", "red", "orange",
  "yellow", "green"
))

# get phenotype and metabolome data
pheno_df <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "phenotype_data_family.csv"
)))
# read metabo data
metabo_df <- as.data.frame(fread(paste0(
  metabo_dir_path,
  metabo_file
)))
metabo_df <- remove_na_columns(metabo_df)

# merge pheno_df and metabo_df for integrity of analyses and slice the merged df
merged_df <- merge(pheno_df, metabo_df, by = "id")
pheno_df <- merged_df[, c("id", selected_traits_)]
metabo_df <- merged_df[, -match(selected_traits_, colnames(merged_df))]
metabo_df <- apply(metabo_df[, -1], 2, as.numeric)

# remove na for analyzed trait_ and corresponding rows for marker data
idx_na_trait_ <- which(is.na(pheno_df[, trait_]))
if (length(idx_na_trait_) > 0) {
  pheno_df <- pheno_df[-idx_na_trait_, ]
  metabo_df <- metabo_df[-idx_na_trait_, ]
}
metabo_df <- scale(metabo_df)

# get number of phenotypes
n <- nrow(pheno_df)

# register parallel backend
cl <- makeCluster(detectCores())
registerDoParallel(cl)

df_result_ <- foreach(
  shuff_ = 1:n_shuff_,
  .packages = pkgs_to_export_,
  .combine = rbind
) %dopar% {
  # set seed, define a new set of indices, and create folds for k-fold cross-validation
  set.seed(shuff_ * mult_seed_by_)
  print(paste0("shuff_ : ", shuff_))
  idx_ <- sample(1:n, size = n, replace = FALSE)
  Folds <- cvFolds(n, K = k_folds_, type = "consecutive")

  fold_results <- foreach(
    fold_ = 1:k_folds_,
    .packages = pkgs_to_export_,
    .combine = rbind
  ) %dopar% {
    idx_val_fold <- which(Folds$which == fold_)
    idx_val <- idx_[idx_val_fold]
    idx_train <- idx_[-idx_val_fold]

    # initialize vector of results for the current fold
    fold_result <- c(
      "RF" = NA, "SVR" = NA, "RKHS" = NA, "GBLUP" = NA, "LASSO" = NA,
      "SVR_support_vectors" = NA
    )

    # train and predict with Random Forest
    rf_model <- ranger(
      y = pheno_df[idx_train, trait_],
      x = metabo_df[idx_train, ],
      mtry = ncol(metabo_df) / 3,
      num.trees = 1000
    )
    f_hat_val_rf <- predict(
      rf_model,
      metabo_df[idx_val, ]
    )
    fold_result["RF"] <- cor(
      f_hat_val_rf$predictions,
      pheno_df[idx_val, trait_]
    )

    # train and predict with SVR
    c_par <- max(
      abs(mean(pheno_df[idx_train, trait_])
      + 3 * sd(pheno_df[idx_train, trait_])),
      abs(mean(pheno_df[idx_train, trait_])
      - 3 * sd(pheno_df[idx_train, trait_]))
    )
    gaussian_svr_model <- ksvm(
      x = as.matrix(metabo_df[idx_train, ]),
      y = pheno_df[idx_train, trait_],
      scaled = FALSE, type = "eps-svr",
      kernel = "rbfdot",
      kpar = "automatic", C = c_par, epsilon = 0.1
    )
    idx_sv_ <- SVindex(gaussian_svr_model)
    sv_ <- pheno_df[idx_train, "id"][idx_sv_]
    f_hat_val_gaussian_svr <- predict(
      gaussian_svr_model,
      as.matrix(metabo_df[idx_val, ])
    )
    fold_result["SVR"] <- cor(
      f_hat_val_gaussian_svr,
      pheno_df[idx_val, trait_]
    )
    fold_result["SVR_support_vectors"] <- paste0(sv_, collapse = ", ")

    # train and predict with GBLUP (linear kernel krmm)
    linear_krmm_model <- krmm(
      Y = pheno_df[idx_train, trait_],
      Matrix_covariates = metabo_df[idx_train, ],
      method = "GBLUP"
    )
    f_hat_val_linear_krmm <- predict_krmm(linear_krmm_model,
      Matrix_covariates = metabo_df[idx_val, ],
      add_flxed_effects = T
    )
    fold_result["GBLUP"] <- cor(
      f_hat_val_linear_krmm,
      pheno_df[idx_val, trait_]
    )

    # train and predict with RKHS (non-linear Gaussian kernel krmm)
    gaussian_krmm_model <- krmm(
      Y = pheno_df[idx_train, trait_],
      Matrix_covariates = metabo_df[idx_train, ],
      method = "RKHS", kernel = "Gaussian",
      rate_decay_kernel = 0.1
    )
    f_hat_val_gaussian_krmm <- predict_krmm(gaussian_krmm_model,
      Matrix_covariates = metabo_df[idx_val, ],
      add_flxed_effects = T
    )
    fold_result["RKHS"] <- cor(
      f_hat_val_gaussian_krmm,
      pheno_df[idx_val, trait_]
    )

    # train and predict with LASSO
    cv_fit_lasso_model <- cv.glmnet(
      intercept = TRUE, y = pheno_df[idx_train, trait_],
      x = as.matrix(metabo_df[idx_train, ]),
      type.measure = "mse", alpha = 1.0, nfold = 10,
      parallel = TRUE
    )
    f_hat_val_lasso <- predict(cv_fit_lasso_model,
      newx = as.matrix(metabo_df[idx_val, ]),
      s = "lambda.min"
    )
    fold_result["LASSO"] <- cor(
      f_hat_val_lasso,
      pheno_df[idx_val, trait_]
    )

    fold_result
  }

  fold_results
}

# stop the parallel backend
stopCluster(cl)
registerDoSEQ()

# change column names
colnames(df_result_) <- c(
  "RF", "SVR", "RKHS",
  "MBLUP", "LASSO", "SVR_support_vectors"
)

# create directory for trait_ graphics if it does not exist
if (!dir.exists(paste0(
  output_pred_graphics_path, metabo_type_, "/",
  trait_, "/"
))) {
  dir.create(paste0(
    output_pred_graphics_path, metabo_type_, "/",
    trait_, "/"
  ))
}

# convert to data frame format
df_result_ <- as.data.frame(df_result_)
df_ <- df_result_[, -match(
  "SVR_support_vectors",
  colnames(df_result_)
)]

# get methods names
df_ <- as.data.frame(apply(df_, 2, as.numeric))
method_names <- colnames(df_)

# initialize plot_ly boxplot graphic
boxplots_pa_ <- plot_ly()

# add boxplots
for (method_ in method_names) {
  boxplots_pa_ <- add_boxplot(
    boxplots_pa_,
    y = df_[[method_]],
    name = method_,
    boxpoints = "all",
    jitter = 0.3,
    pointpos = -1.8
  )
}
# add layout
boxplots_pa_ <- boxplots_pa_ %>%
  layout(
    title = paste0(
      "Metabolomic prediction PA distributions of methods for ",
      trait_, ", based on ", ncol(metabo_df), " metabolites across ",
      n_shuff_, " shuffling and ", k_folds_, "-folds scenarios"
    ),
    yaxis = list(title = "Predictive ability (PA)", range = c(-1, 1)),
    legend = list(title = list(text = "Prediction method"))
  )

# save boxplots_pa_ graphics
saveWidget(boxplots_pa_, file = paste0(
  output_pred_graphics_path, metabo_type_, "/",
  trait_, "/pa_",
  trait_, "_", ncol(metabo_df), "_metabolites", ".html"
))

# save predictive ability results
rownames(df_) <- paste0("pa_shuff_", 1:nrow(df_))
df_stat <- as.data.frame(rbind(apply(df_, 2, mean), apply(df_, 2, sd)))
rownames(df_stat) <- c("pa_mean", "pa_sd")
df_ <- rbind(df_, df_stat)
fwrite(df_,
  file = paste0(
    output_pred_results_path, metabo_type_,
    "/metabolomic_pred_results_",
    ncol(metabo_df), "_metabolites_",
    trait_, ".csv"
  ), row.names = T
)

# get genotypes used as support vectors
sv_metabo <- unlist(str_split(df_result_$SVR_support_vectors,
  pattern = ", "
))

# create a frequency table
freq_table <- table(sv_metabo)

# convert the frequency table to a dataframe
df_sv_ <- data.frame(
  Genotype = names(freq_table),
  Count = as.numeric(freq_table)
)

# create barplot_sv_orig with colored bars based on origins
barplot_sv_orig <- plot_ly(df_sv_,
  x = ~Genotype, y = ~Count, type = "bar",
  color = ~Genotype,
  colors = col_func_metabo(nrow(df_sv_))
) %>%
  layout(
    title = paste0(
      "Counts of genotypes identified as support vectors across ",
      n_shuff_, " Gaussian SVR models, \n used for metabolomic prediction of ",
      trait_, " across ", n_shuff_, " shuffling and ",
      k_folds_, "-folds scenarios"
    ),
    xaxis = list(
      categoryorder = "total descending", title = "Genotype",
      tickangle = 300
    ),
    yaxis = list(title = "Count"),
    legend = list(title = list(text = "Origin"))
  )
# save barplot_sv_orig graphics
saveWidget(barplot_sv_orig, file = paste0(
  output_pred_graphics_path, metabo_type_, "/",
  trait_, "/sv_",
  trait_, "_", ncol(metabo_df), "_metabolites", ".html"
))
