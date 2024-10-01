
suppressPackageStartupMessages(library(arm))
suppressPackageStartupMessages(library(ranger))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(rootSolve))
suppressPackageStartupMessages(library(cubature))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(glm2))
suppressPackageStartupMessages(library(simcausal))
suppressPackageStartupMessages(library(hal9001))
suppressPackageStartupMessages(library(extraDistr))
suppressPackageStartupMessages(library(mc2d))
suppressPackageStartupMessages(library(SuperLearner))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(randomForest))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(xgboost))
suppressPackageStartupMessages(library(mlr3superlearner))
suppressPackageStartupMessages(library(nnet))
suppressPackageStartupMessages(library(glm2))
suppressPackageStartupMessages(library(simcausal))
suppressPackageStartupMessages(library(mlr3extralearners))
suppressPackageStartupMessages(library(lightgbm))
suppressPackageStartupMessages(library(kknn))
suppressPackageStartupMessages(library(earth))
suppressPackageStartupMessages(library(Formula))
suppressPackageStartupMessages(library(plotmo))
suppressPackageStartupMessages(library(plotrix))
suppressPackageStartupMessages(library(TeachingDemos))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(np))
library(here)



##############################################################################################################################


rm(list=ls())
source(here("R/utils.R"))
load(here("Case_Study/data","paper_data_full.RData"))
data <- data %>%
  mutate(across(everything(), as.numeric))

qs <- as.numeric(quantile(unlist(data[,grep("dex_a_",colnames(data))]), c(0.2, 0.4, 0.6, 0.8)))
## scale treatment to 0-1 across all time points ##
data <- data %>% mutate(across(starts_with("dex_a_"), ~ case_when(
  . < qs[1] ~ 0,
  . >= qs[1] & . < qs[2] ~ 1,
  . >= qs[2] & . < qs[3] ~ 2,
  . >= qs[3] & . < qs[4] ~ 3,
  . >= qs[4] ~ 4
) ))


# types #
Y_type = "cum_cases_percap_"
prev_Y = "cum_cases_percap_"

# # range of data #
start_time = 2
end_time = 25
# stacks #
stackr_temp_c = list("mean", "lightgbm", "multinom","xgboost", "nnet",
                     "knn", "rpart", "svm", "lda","naivebayes","glmnet", "cv_glmnet",
                     list("randomforest", ntree = 1000, id = "randomforest"),
                     list("ranger", num.trees = 1000, id = "ranger")
)

learners_ranger <- create.Learner('SL.ranger', params = list(num.trees = 100))
learners_rf <- create.Learner('SL.randomForest', params = list(ntree = 100))
stackm_temp_c = c('SL.mean','SL.glm','SL.glmnet',
                  'SL.bayesglm','SL.earth','SL.rpart',learners_ranger$names,learners_rf$names)

stackr_temp <- list(stackr_temp_c,stackr_temp_c,stackr_temp_c,stackr_temp_c,stackr_temp_c,stackr_temp_c,stackr_temp_c,stackr_temp_c)
stackm_temp <- list(stackm_temp_c,stackm_temp_c,stackm_temp_c,stackm_temp_c,stackm_temp_c,stackm_temp_c,stackm_temp_c,stackm_temp_c)

# seed #
s = 123

A_ex = paste0("dex_a_",seq(start_time, end_time, by = 3))
end_time = max(seq(start_time, end_time, by = 3))
tau = length(A_ex)
# covariates #
covs_base <- c(c("log_pop_dens_", "log_pop_", "pct_rural_2010_",
                 "ppl_below_poverty_", "hh_below_poverty_", "pct_black_",
                 "pct_hispanic_", "pct_asian_", "pct_under_18_", "pct_65plus_",
                 "some_college_", "pct_female_",
                 "inc_ineq_", "pct_dem_2016_",
                 "meatpack_cat4_", "meatpack_cat5_", "meatpack_4or5_",
                 "drive_alone_to_work_", "prison_",
                 "pct_crowded_housing_", "svi_ses_", "avg_hh_size_", "unemployment_",
                 "uninsured_adults_", "adult_smoke_rate_", "log_med_hh_inc_",
                 "obesity_pct_", "diabetes_pct_"),
               paste0("state_abb_",unique(gsub("_\\d+","",gsub("state_abb_","",colnames(data)[grep("state", colnames(data))]))),"_")
)
cov_var <- c("mask_level_")
L_ex <- list(c(paste0(covs_base, start_time-1), paste0(prev_Y,cov_var, start_time-1)),
             paste0(c(prev_Y,cov_var), start_time+2),
             paste0(c(prev_Y,cov_var), start_time+5),
             paste0(c(prev_Y,cov_var), start_time+8),
             paste0(c(prev_Y,cov_var), start_time+11),
             paste0(c(prev_Y,cov_var), start_time+14),
             paste0(c(prev_Y,cov_var), start_time+17),
             paste0(c(prev_Y,cov_var), start_time+20)
)

# outcome #
Y_ex = paste0(Y_type,"25")

### RUN ###
run_tmp <-
  estimators(data = data,
             tau = tau,
             A = A_ex,
             L = L_ex,
             Y = Y_ex,
             est = c("sdr","tml"),
             A_type = "categorical",
             outcome = "gaussian",
             phi = function(a) cbind(1, rowMeans(a)),
             lambda = T,
             nfolds = 1,
             valSets = NULL,
             stackr = stackr_temp,
             stackm = stackm_temp,
             k = Inf,
             trace = T,
             cores = 1,
             ken.dens = T,
             approx_int = T,
             trunc = F,
             trim = F,
             U2_method = "Nelder-Mead",
             seed = s, A.min = 0, A.max = 4)
# run_tmp$out %>% select(sdr, varsdr, sdrlo, sdrhi, tml, vartml, tmllo, tmlhi, ipw)


save(run_tmp, file = here("Case_Study/results/NM_categorical_notrim_6months.Rdata"))

