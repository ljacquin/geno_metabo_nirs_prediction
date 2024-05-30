# script meant to perform phenomic prediction and analyses for apple scab resistance
# note: text is formatted from Addins using Style active file from styler package

# clear memory and source libraries
rm(list = ls())
library(MASS)
library(data.table)
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
library(KRMM)
library(kernlab)
library(glmnet)
library(ranger)
library(geomtextpath)
library(dplyr)
library(tidyr)
library(viridis)
library(stringr)
library(prospectr)
library(signal)
library(FactoMineR)
library(factoextra)
library(lme4)
library(lsmeans)
library(png)
library(devtools)
install_other_requirements <- F
if (install_other_requirements) {
  devtools::install_github("ljacquin/KRMM")
}
local_computation_ <- T
if (local_computation_) {
  library(rstudioapi)
  setwd(dirname(getActiveDocumentContext()$path))
}
source("../../functions.R")

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
nirs_dir_path <- "../../../data/phenomic_data/"
pheno_dir_path <- "../../../data/phenotype_data/"

# output result path for phenomic graphics
output_pred_results_path <- "../../../results/phenomic_prediction/"
output_graphics_path <- "../../../results/graphics/phenomic_graphics/"
output_pred_graphics_path <- "../../../results/graphics/phenomic_prediction_graphics/"

# define nirs type for analysis
nirs_type_ <- "NIRJ-1" #"NIRJ5" # "NIRJ-1"

# define selected blocs for experiment
selected_blocs_ <- c("B2", "B8") #c("B4", "B6") # c("B2", "B8")

# should de-trending be applied to spectra ?
apply_detrend_ <- T

# should Savitzky-Golay filtering be applied to spectra ?
apply_sg_filter_ <- T

# define derivative order, polynome order and number of support points for
# Savitzky-Golay filtering
deriv_order_ <- 1 # 2
polynome_order_ <- 2 # 3
num_supp_ <- 37 # 61
vect_order_char_ <- c("first", "second")

# define selected traits, and the analyzed one defined by trait_
selected_traits_ <- c(
  "BLUPaudpc.eau.EUB04.18",
  "BLUPaudpc.eau.EUB04.19",
  "BLUPaudpc.eau.09BCZ14.prin19"
)
trait_ <- "BLUPaudpc.eau.EUB04.19" # "BLUPaudpc.eau.09BCZ14.prin19"

# define experimental vars
experim_vars_ <- c(
  "geno", "bloc", "kinetic",
  "enveloppe", "position"
)

# define shift seed value by
mult_seed_by_ <- 100

# set k for K-folds cv
k_folds_ <- 5

# define number of shuffles
n_shuff_ <- 20

# define color function for genotypes
col_func_genotype <- colorRampPalette(c(
  "black", "blue", "red", "orange",
  "yellow", "green"
))

# get phenotype, nirs and bloc data

# pheno data
pheno_df <- as.data.frame(fread(paste0(
  pheno_dir_path,
  "phenotype_data_family.csv"
)))

# nirs data
nirs_df <- as.data.frame(fread(
  paste0(
    nirs_dir_path,
    paste0(nirs_type_, "_FAMF.txt")
  ),
  header = T
))
nirs_df$Wavelength <- paste0(
  nirs_type_, "_", unlist(str_extract_all(
    nirs_df$Wavelength,
    pattern = "[0-9_]{3,4}_[0-9_]{1}"
  ))
)

# geno/bloc data
geno_bloc_df <- as.data.frame(fread(paste0(
  nirs_dir_path,
  "labels_alu_spring_2024.csv"
), header = T))
colnames(geno_bloc_df) <- c("bloc", "geno", "enveloppe")

# split wave length column into "kinetic", "enveloppe" and "position"
nirs_df <- separate(nirs_df, Wavelength, c(
  "kinetic",
  "enveloppe",
  "position"
), sep = "_")

# select genotypes according to defined blocs if they are defined
if (length(selected_blocs_) > 0) {
  geno_bloc_df <- geno_bloc_df[
    geno_bloc_df$bloc %in% selected_blocs_,
  ]
}

# filter nirs data according to selected blocs and merge by enveloppe
nirs_filtered_df <- merge(nirs_df, geno_bloc_df, by = "enveloppe")
nirs_filtered_df <- nirs_filtered_df %>% dplyr::select(
  experim_vars_,
  everything()
)

# clean filtered nirs data for columns having na on all rows
nirs_filtered_df <- nirs_filtered_df[
  , colSums(is.na(nirs_filtered_df)) != nrow(nirs_filtered_df)
]

# correct for geno, bloc and  position effect
fitted_col_names_ <- colnames(nirs_filtered_df)[
  (length(experim_vars_)
  + 1):ncol(nirs_filtered_df)
]
waves_fitted_values_list <- lapply(
  fitted_col_names_,
  fit_geno_bloc_pos_and_extract,
  df = nirs_filtered_df
)
waves_fitted_values_df <- do.call(
  cbind,
  waves_fitted_values_list
)
nirs_filtered_df <- cbind(
  nirs_filtered_df[, experim_vars_],
  waves_fitted_values_df
)
colnames(nirs_filtered_df) <- c(
  experim_vars_,
  fitted_col_names_
)

# rename wave length vars associated to reflectance
wave_len_orig_vars_ <- colnames(nirs_filtered_df)[-match(
  experim_vars_,
  colnames(nirs_filtered_df)
)]
wave_len_vars_ <- paste0("x_", wave_len_orig_vars_)
colnames(nirs_filtered_df) <- c(
  experim_vars_,
  wave_len_vars_
)

# compute adjusted ls-means, for each wavelength, for each genotype over blocs
# note : this is a long computation, hence the results should be saved and reused later
if (!file.exists(paste0(
  nirs_dir_path,
  nirs_type_, "_", "merged_wave_ls_mean_df.csv"
))) {
  # define a list to save all data frames and perform ls-means for each wavelength,
  # by genotype, over all blocs
  list_ls_mean_df_ <- vector("list", length(wave_len_vars_))
  names(list_ls_mean_df_) <- wave_len_vars_

  for (wave_var_ in wave_len_vars_) {
    lm_model <- lm(
      formula(paste0(
        wave_var_,
        "~ geno + bloc"
      )),
      data = nirs_filtered_df
    )
    wave_var_ls_means_df <- as.data.frame(
      lsmeans(
        lm_model, ~geno
      )
    )[, c("geno", "lsmean")]
    colnames(wave_var_ls_means_df) <- c("geno", wave_var_)
    list_ls_mean_df_[[wave_var_]] <- wave_var_ls_means_df
  }

  # merge list of data frames
  merged_wave_ls_mean_df <- Reduce(merge_by_geno, list_ls_mean_df_)
  colnames(merged_wave_ls_mean_df) <- c("geno", wave_len_orig_vars_)

  # transpose merged_wave_ls_mean_df data frame
  merged_wave_ls_mean_df <- as.data.frame(t(merged_wave_ls_mean_df))
  colnames(merged_wave_ls_mean_df) <- merged_wave_ls_mean_df[1, ]
  merged_wave_ls_mean_df <- merged_wave_ls_mean_df[-1, ]
  merged_wave_ls_mean_df <- cbind(
    as.numeric(rownames(merged_wave_ls_mean_df)),
    merged_wave_ls_mean_df
  )
  colnames(merged_wave_ls_mean_df)[1] <- "lambda"
  fwrite(merged_wave_ls_mean_df, file = paste0(
    nirs_dir_path,
    nirs_type_, "_",
    "merged_wave_ls_mean_df.csv"
  ))
} else {
  merged_wave_ls_mean_df <- as.data.frame(fread(paste0(
    nirs_dir_path,
    nirs_type_, "_",
    "merged_wave_ls_mean_df.csv"
  )))
}

# save the graphical ls means of the spectra for each genotype accross blocs
png(
  paste0(
    output_graphics_path,
    nirs_type_, "/",
    nirs_type_, "_spectra_ls_means.PNG"
  ),
  width = 1200, height = 600
)
matplot(
  merged_wave_ls_mean_df$lambda,
  merged_wave_ls_mean_df[, -1],
  type = "l", lty = 1, pch = 0,
  xlab = "Lambda (nm)", ylab = "Reflectance",
  main = paste0(
    nirs_type_,
    " reflectance ls-means for each genotype accross blocs"
  ),
  xlim = c(
    min(merged_wave_ls_mean_df$lambda),
    max(merged_wave_ls_mean_df$lambda)
  )
)
dev.off()

# get spectra resolution in nm
spec_res_ <- (max(merged_wave_ls_mean_df$lambda)
- min(merged_wave_ls_mean_df$lambda)) /
  (length(merged_wave_ls_mean_df$lambda) - 1)

# apply different transformations to the ls-means spectra
wave_trans_df <- merged_wave_ls_mean_df
rm(merged_wave_ls_mean_df)
wave_trans_df[, 2:ncol(wave_trans_df)] <- scale(
  apply(
    wave_trans_df[, 2:ncol(wave_trans_df)], 2,
    as.numeric
  )
)

# apply detrending ?
if (apply_detrend_) {
  lambda_ <- wave_trans_df$lambda
  # apply detrending
  wave_trans_df <- as.data.frame(t(detrend(
    X = t(apply(wave_trans_df[, -1], 2, as.numeric)),
    wav = lambda_
  )))
  wave_trans_df <- cbind(
    lambda_,
    wave_trans_df
  )
  colnames(wave_trans_df)[1] <- "lambda"
  # save de-trended spectra
  png(
    paste0(
      output_graphics_path,
      nirs_type_, "/",
      nirs_type_, "_spectra_ls_means_detrended.PNG"
    ),
    width = 1200, height = 600
  )
  matplot(
    wave_trans_df$lambda,
    wave_trans_df[, -1],
    type = "l", lty = 1, pch = 0,
    xlab = "Lambda (nm)", ylab = "Reflectance",
    main = paste0(
      nirs_type_,
      " detrended reflectance ls-means for each genotype accross blocs"
    ),
    xlim = c(
      min(wave_trans_df$lambda),
      max(wave_trans_df$lambda)
    )
  )
  dev.off()
}

# apply Savitzky-Golay filtering ?
if (apply_sg_filter_) {
  # apply derivative, for the chosen order, to the ls means of the spectra
  wave_trans_df <- t(as.matrix(
    apply(
      scale(apply(wave_trans_df[, -1], 2, as.numeric)), 2, function(x) {
        sgolayfilt(x,
          p = polynome_order_,
          n = num_supp_,
          m = deriv_order_,
          ts = spec_res_
        )
      }
    )
  )) # p is the polynome order, m-th derivative, n support points
  wave_trans_df <- as.data.frame(cbind(
    rownames(wave_trans_df),
    wave_trans_df
  ))
  colnames(wave_trans_df)[1] <- "id"
  # save filtered spectra
  if (apply_detrend_) {
    out_path_ <- paste0(
      output_graphics_path,
      nirs_type_, "/",
      nirs_type_,
      "_spectra_detrended_ls_means_and_sgolay_",
      vect_order_char_[deriv_order_],
      "_order_filtered.PNG"
    )
    graphic_title_ <- paste0(
      nirs_type_,
      " detrended and Savitzky–Golay ",
      vect_order_char_[deriv_order_],
      " order filtered reflectance ls-means for each genotype accross blocs"
    )
  } else {
    out_path_ <- paste0(
      output_graphics_path,
      nirs_type_, "/",
      nirs_type_,
      "_spectra_ls_means_sgolay_",
      vect_order_char_[deriv_order_],
      "_order_filtered.PNG"
    )
    graphic_title_ <- paste0(
      nirs_type_,
      " Savitzky–Golay ",
      vect_order_char_[deriv_order_],
      " order filtered reflectance ls-means for each genotype accross blocs"
    )
  }
  wave_trans_sg_plot_df <- as.data.frame(t(wave_trans_df))[-1, ]
  wave_trans_sg_plot_df <- cbind(
    as.numeric(wave_len_orig_vars_),
    wave_trans_sg_plot_df
  )
  colnames(wave_trans_sg_plot_df)[1] <- "lambda"
  png(out_path_, width = 1200, height = 600)
  matplot(
    wave_trans_sg_plot_df$lambda,
    wave_trans_sg_plot_df[, -1],
    type = "l", lty = 1, pch = 0,
    xlab = "Lambda (nm)", ylab = "Reflectance",
    main = graphic_title_,
    xlim = c(
      min(wave_trans_sg_plot_df$lambda),
      max(wave_trans_sg_plot_df$lambda)
    )
  )
  dev.off()
} else {
  wave_trans_df <- t(as.matrix(wave_trans_df[, -1]))
  wave_trans_df <- as.data.frame(cbind(
    rownames(wave_trans_df),
    wave_trans_df
  ))
  colnames(wave_trans_df)[1] <- "id"
}

# merge pheno_df and wave_trans_df for integrity of analyses
# and slice the merged df
merged_df <- merge(pheno_df, wave_trans_df, by = "id")
pheno_df <- merged_df[, c("id", selected_traits_)]
wave_trans_vars_df <- merged_df[, -match(selected_traits_, colnames(merged_df))]
wave_trans_vars_df <- apply(wave_trans_vars_df[, -1], 2, as.numeric)

# remove na for analyzed trait_ and corresponding rows for marker data
idx_na_trait_ <- which(is.na(pheno_df[, trait_]))
if (length(idx_na_trait_) > 0) {
  pheno_df <- pheno_df[-idx_na_trait_, ]
  wave_trans_vars_df <- wave_trans_vars_df[-idx_na_trait_, ]
}

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
      "RF" = NA, "SVR" = NA, "RKHS" = NA, "PBLUP" = NA, "LASSO" = NA,
      "SVR_support_vectors" = NA
    )

    # train and predict with Random Forest
    rf_model <- ranger(
      y = pheno_df[idx_train, trait_],
      x = wave_trans_vars_df[idx_train, ],
      mtry = ncol(wave_trans_vars_df) / 3,
      num.trees = 1000
    )
    f_hat_val_rf <- predict(
      rf_model,
      wave_trans_vars_df[idx_val, ]
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
      x = as.matrix(wave_trans_vars_df[idx_train, ]),
      y = pheno_df[idx_train, trait_],
      scaled = FALSE, type = "eps-svr",
      kernel = "rbfdot",
      kpar = "automatic", C = c_par, epsilon = 0.01
    )
    idx_sv_ <- SVindex(gaussian_svr_model)
    sv_ <- pheno_df[idx_train, "id"][idx_sv_]
    f_hat_val_gaussian_svr <- predict(
      gaussian_svr_model,
      as.matrix(wave_trans_vars_df[idx_val, ])
    )
    fold_result["SVR"] <- cor(
      f_hat_val_gaussian_svr,
      pheno_df[idx_val, trait_]
    )
    fold_result["SVR_support_vectors"] <- paste0(sv_, collapse = ", ")

    # train and predict with GBLUP (linear kernel krmm),
    # note PBLUB is equivalent to GBLUP method here
    linear_krmm_model <- krmm(
      Y = pheno_df[idx_train, trait_],
      Matrix_covariates = wave_trans_vars_df[idx_train, ],
      method = "GBLUP"
    )
    f_hat_val_linear_krmm <- predict_krmm(linear_krmm_model,
      Matrix_covariates = wave_trans_vars_df[idx_val, ],
      add_flxed_effects = T
    )
    fold_result["PBLUP"] <- cor(
      f_hat_val_linear_krmm,
      pheno_df[idx_val, trait_]
    )

    # train and predict with RKHS (non-linear Gaussian kernel krmm)
    gaussian_krmm_model <- krmm(
      Y = pheno_df[idx_train, trait_],
      Matrix_covariates = wave_trans_vars_df[idx_train, ],
      method = "RKHS", kernel = "Gaussian",
      rate_decay_kernel = 1,
      init_sigma2K = 1e9, init_sigma2E = 1
    )
    f_hat_val_gaussian_krmm <- predict_krmm(gaussian_krmm_model,
      Matrix_covariates = wave_trans_vars_df[idx_val, ],
      add_flxed_effects = T
    )
    fold_result["RKHS"] <- cor(
      f_hat_val_gaussian_krmm,
      pheno_df[idx_val, trait_]
    )

    # train and predict with LASSO
    cv_fit_lasso_model <- cv.glmnet(
      intercept = TRUE, y = pheno_df[idx_train, trait_],
      x = as.matrix(wave_trans_vars_df[idx_train, ]),
      type.measure = "mse", alpha = 1.0, nfold = 10,
      parallel = TRUE
    )
    f_hat_val_lasso <- predict(cv_fit_lasso_model,
      newx = as.matrix(wave_trans_vars_df[idx_val, ]),
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

# create directory for trait_ graphics if it does not exist
if (!dir.exists(paste0(
  output_pred_graphics_path,
  nirs_type_, "/", trait_, "/"
))) {
  dir.create(paste0(
    output_pred_graphics_path,
    nirs_type_, "/", trait_, "/"
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
      "Phenomic prediction PA distributions of methods for ",
      trait_, ",\n based on ", nirs_type_, " ",
      ncol(wave_trans_vars_df), " wave lengths across ",
      n_shuff_, " shuffling and ", k_folds_, "-folds scenarios"
    ),
    yaxis = list(title = "Predictive ability (PA)", range = c(-1, 1)),
    legend = list(title = list(text = "Prediction method"))
  )

# save boxplots_pa_ graphics
saveWidget(boxplots_pa_, file = paste0(
  output_pred_graphics_path,
  nirs_type_, "/",
  trait_, "/pa_",
  trait_, "_", ncol(wave_trans_vars_df), "_wave_lengths", ".html"
))
boxplots_pa_

# save predictive ability results
rownames(df_) <- paste0("pa_shuff_", 1:nrow(df_))
df_stat <- as.data.frame(rbind(apply(df_, 2, mean), apply(df_, 2, sd)))
rownames(df_stat) <- c("pa_mean", "pa_sd")
df_ <- rbind(df_, df_stat)
fwrite(df_,
  file = paste0(
    output_pred_results_path,
    nirs_type_, "/",
    "phenomic_pred_results_", ncol(wave_trans_vars_df), "_wave_lengths_",
    trait_, ".csv"
  ), row.names = T
)

# get genotypes used as support vectors
sv_geno <- unlist(str_split(df_result_$SVR_support_vectors,
  pattern = ", "
))

# create a frequency table
freq_table <- table(sv_geno)

# convert the frequency table to a dataframe
df_sv_ <- data.frame(
  Genotype = names(freq_table),
  Count = as.numeric(freq_table)
)

# create barplot_sv_orig with colored bars based on origins
barplot_sv_orig <- plot_ly(df_sv_,
  x = ~Genotype, y = ~Count, type = "bar",
  color = ~Genotype,
  colors = col_func_genotype(nrow(df_sv_))
) %>%
  layout(
    title = paste0(
      "Counts of genotypes identified as support vectors across ",
      n_shuff_, " Gaussian SVR models, \n used for phenomic prediction of ",
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
  output_pred_graphics_path,
  nirs_type_, "/",
  trait_, "/sv_for_",
  trait_, "_", ncol(wave_trans_vars_df),
  "_wave_lengths", ".html"
))
