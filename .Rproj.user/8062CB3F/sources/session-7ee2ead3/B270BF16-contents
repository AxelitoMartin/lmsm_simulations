rm(list=ls())
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
suppressPackageStartupMessages(library(e1071))
library(here)

# load #
load(here('data/lambda_fits_simple_mult.Rdata'))
load(here('data/true_discrete_cat5_simple_mult.Rdata'))
source(here('R/utils.R'))

stackr = list("mean", "lightgbm", "multinom","xgboost", "nnet",
              "knn", "rpart", "naivebayes","glmnet",
              list("randomforest", ntree = 250, id = "randomforest"),
              list("ranger", num.trees = 250, id = "ranger")
)

stackm = c('SL.mean','SL.glmnet','SL.earth','SL.glm.interaction')

tau = 4
A.min <- 0
A.max <- 4
increment <- 1
quants <- sapply(1:tau,
                 function(j)
                   seq(A.min, A.max, by = increment))
A_span <- as.matrix(expand.grid(
  lapply(1:tau,function(j) quants[,j]
  )))
binsizes <- 1
bindim <- 1


temp <- as.data.frame(A_span)
colnames(temp) <- c("A_1","A_2","A_3","A_4")
lambda = lambda_fits
la <- sapply(1:tau,function(j)lambda_cust(j, A_name = c("A_1","A_2","A_3","A_4"), temp, k = Inf, lambda, trunc = F, A_type = "categorical"))
la <- apply(la, 1, prod)

# function #
simula <- function(j){
  print(j)
  seed <- task[j, 'seed']
  set.seed(seed)
  n <- task[j, 'n']
  data <- hush(datagen_mult(n, size_A = 5, tau = tau, mu = mu_mult, pl = pl_mult, m4 = m4_mult, L1_size = 5))
  
  # with stab weights #
  results_la <- hush(simul_mult(data = data$data, stackr = stackr, stackm = stackm, lambda = lambda, seed = seed, la = la,
                                A_span = A_span, binsizes = binsizes, bindim = bindim, 
                                get_la = T, trueW = NULL))
  sink()
  
  results <- list(results_la)
  return(results)
}

set.seed(1)
task <- expand.grid(n = c(500), seed = sample(10^6, 250))

system.time(
  res <- mclapply(1:nrow(task), simula, mc.cores = 5)
)

save(res, file = here("./results/Mult_n500.Rdata"))