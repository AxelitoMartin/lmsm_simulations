
########### Message suppression #######

#' hush
#'
#' Simple function to suppress spurious output
#' @param code Chunck of code with output to console to be suppressed
#' @return output of funciton without console output
#' @export
hush=function(code){
  sink("/dev/null") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}


########## Bounds ############
#' bound
#'
#' Values to be bounded from above and below
#' @param x Values to be bounded
#' @p boundary
#' @return Bounded values
#' @export
bound <- function(x, p = 1e-5) pmax(pmin(x, 1 - p), p)

#' truncate
#'
#' Values to be truncated from below
#' @param x Values to be bounded
#' @p boundary
#' @return Truncated values
#' @export
truncate <- function(x, p = 1e-2) pmax(x, p)


############### ESTIMATING FUNCTIONS ##################
#' lambda_cust_int
#'
#' Function to perform sum/integration of U1 including the stabilized weights provided by the user
#' @param j time point
#' @param data data to be used to perform prediction of lambda values
#' @param k Last time point to consider when getting prior history
#' @param fits user provided stabilizing lambda value fits (e.g.: multinomial fits at each time point considering prior treatment history)
#' @return predicted lambdas for vlaues of treatment
#' @export
lambda_cust_int <- function(j, data, A_name, k, fits, trunc = F) {
  
  nameA <- "A"
  a <- data[, nameA]
  n <- length(a)
  
  if(length(unique(fits[[j]])) == 1 || length(unique(a)) == 1)
    return(rep(1, n))
  
  if(j == 1){
    temp_fit <- fits[[j]]
    pihat <- c()
    for(r in 1:n){
      pihat[r] <- temp_fit[(a[r] + 1)]
    }
    return(pihat)
  }
  
  # select all prior trt history
  if(j > 1) namesA <- namesA <- A_name[max(1, (j-k)):(j-1)] else namesA <- NULL
  # stack in a single dataset
  l <- data %>% select(one_of(c(namesA)))
  temp_fit <- fits[[j]]
  preds <- predict(temp_fit, type = "probs", newdata = l)
  
  pihat <- c()
  for(r in 1:n){
    if(!is.null(dim((preds)) )){
      temp_preds <- preds[r,]
    } else temp_preds <- preds
    pihat[r] <- as.numeric(temp_preds[(a[r] + 1)])
  }
  
  if(trunc) pihat <- bound(pihat, p = 0.01)
  return(pihat)
  
}


#' SL_ladata_int
#'
#' Function to perform sum/integration of U1 including the stabilized weights provided by the user
#' @param j time point
#' @param data data to be used to perform prediction of lambda values
#' @param k Last time point to consider when getting prior history
#' @param fits user provided stabilizing lambda value fits (e.g.: multinomial fits at each time point considering prior treatment history)
#' @return predicted lambdas for vlaues of treatment
#' @export
SL_ladata_int <- function(j, data, A_name, lambda , k = Inf, trunc) {
  
  
  # nameA <- "A"
  # a <- data[, nameA]
  # n <- length(a)
  #
  # # select all prior trt/covariate history
  # namesA <- A_name[max(1, (j-k)):(j-1)]
  #
  # approx.fn <- function(x,y,z){ predict(smooth.spline(x,y),x=z)$y }
  #
  #
  #
  # l <- data %>% select(one_of(c(namesA)))
  #
  # if(ncol(l) == 1) {
  #   temp <- colnames(l)
  #   l <- cbind(0,l)
  #   colnames(l) <- c("Intercept",temp)
  # }
  #
  #
  # # rm(lambda)
  # # lambda <- lambda_fits
  # # lambda[[j]]$pimod$data <- lambda[[j]]$pimod$xlevels <- lambda[[j]]$pimod$R <- lambda[[j]]$pimod$residuals <-
  # #   lambda[[j]]$pimod$fitted.values <- lambda[[j]]$pimod$weights <- lambda[[j]]$pimod$var <-
  # #   lambda[[j]]$pimod$additive.predictors <- lambda[[j]]$pimod$prior.weights <- lambda[[j]]$pimod$y <- NULL
  #
  #
  # pimod.vals <- predict(lambda[[j]]$pimod,newdata = l)
  # # pimod.vals <- predict(lambda[[j]]$pimod,newdata = l)$pred
  # sq.res <- (a-pimod.vals[1:n])^2
  #
  # pi2mod.vals <-predict(lambda[[j]]$pi2mod,newdata = l)
  # # pi2mod.vals <-predict(lambda[[j]]$pi2mod,newdata = l)$pred
  # # construct estimated pi/varpi and mu/m values
  # a.std         <- (a-pimod.vals)/sqrt(abs(pi2mod.vals))
  # pihat.vals    <- approx.fn(density(na.omit(a.std[1:n]))$x, density(na.omit(a.std)[1:n])$y,na.omit(a.std[1:n]))
  # dens_a         <- (pihat.vals/sqrt(abs(pi2mod.vals)))[1:n]
  
  nameA <- "A"
  a <- data[, nameA]
  n <- length(a)
  
  # make dat with only current and prior treatments #
  namesA <- A_name[max(1, (j-k)):(j-1)]
  l <- data %>% select(one_of(c(namesA,nameA)))
  
  if(j == 1){
    dens_a <- approx(lambda[[j]]$density_model$x, lambda[[j]]$density_model$y, xout = a)$y
    # dens_a <- predict(lambda[[j]]$density_model, newdata = l)
  } else {
    l <- data %>% select(one_of(c(namesA)))
    
    if(ncol(l) == 1) {
      temp <- colnames(l)
      l <- cbind(0,l)
      colnames(l) <- c("Intercept",temp)
    }
    pimod.vals <- predict(lambda[[j]]$pimod,newdata = l)
    # pimod.vals <- predict(lambda[[j]]$pimod,newdata = l)$pred
    sq.res <- (a-pimod.vals[1:n])^2
    
    pi2mod.vals <-predict(lambda[[j]]$pi2mod,newdata = l)
    # pi2mod.vals <-predict(lambda[[j]]$pi2mod,newdata = l)$pred
    
    
    approx.fn     <- function(x,y,z){ predict(smooth.spline(x,y),x=z)$y }
    # construct estimated pi/varpi and mu/m values
    a.std         <- (a-pimod.vals)/sqrt(abs(pi2mod.vals))
    pihat.vals    <- approx.fn(density(na.omit(a.std[1:n]))$x, density(na.omit(a.std)[1:n])$y,na.omit(a.std[1:n]))
    dens_a         <- (pihat.vals/sqrt(abs(pi2mod.vals)))[1:n]
  }
  
  return(dens_a)
  
}


#' integral
#'
#' Function to get sum/integrated values of treatment over given prior history and predicted outcome including stabilized weights
#' @param fit a SuperLearner fit of the predicted outcome at time point t
#' @param H matrix to be used for prediction
#' @return integrated/sum values of H over fit multiplied by stabilizing lambda
#' @export
integral <- function(fit, H, A, lambda, A_type, A.min, A.max, trunc, size, A_span = NULL, binsizes = NULL, k = Inf, approx_int, get_la){
  
  
  if(A_type == "categorical"){
    # Generate the sequence for A_span outside the loop
    # A_span <- t(seq(A.min, A.max, by = increment))
    # A_span <- as.matrix(A_span)
    # create all counterfactual data #
    newXX <- do.call('rbind',lapply(seq(nrow(H)), function(H_row){
      data.frame(A = A_span[1,], H[H_row, , drop = FALSE])
    }))
    
    # get all predictions for counterfactuals #
    if(get_la){
      t_temp <- sum(A %in% colnames(newXX))
      la <- lambda_cust_int((t_temp+1), newXX, A_name = A, k, fits = lambda, trunc)
    } else la <- 1
    # tmp <- predict(fit, newdata = newXX)$pred[,1] * la
    
    # weird potential issue #
    tmp <- try(predict(fit, newdata = newXX)$pred[,1] * la, silent = T)
    if(is.character(tmp)) as.numeric(predict(fit)$pred[1,1]) * la
    
    # sum counterfactuals within each respective group #
    seq_sum <- c(seq(1,length(tmp), by = size),length(tmp) + 1)
    preds <- sapply(1:(length(seq_sum) - 1), function(i) {
      sum(tmp[seq_sum[i]:(seq_sum[i + 1] - 1)])
    })
    
  } else if(A_type == "continuous"){
    
    if(approx_int){
      
      # Generate the sequence for A_span outside the loop
      # A_span <- t(seq(A.min, A.max, by = increment))
      
      # create all counterfactual data #
      newXX <- do.call('rbind',lapply(seq(nrow(H)), function(H_row){
        data.frame(A = A_span[1,], H[H_row, , drop = FALSE])
      }))
      
      # get all predictions for counterfactuals #
      if(get_la){
        t_temp <- sum(A %in% colnames(newXX))
        if(is.list(lambda)){
          la <- SL_ladata_int((t_temp+1), newXX, A, lambda, k, trunc)
        } else if(is.function(lambda)){
          la <- lambda(data[,A[t_temp+1]])
        } else la <- 1
      } else la <- 1
      tmp <- predict(fit, newdata = newXX)$pred[,1] * la
      
      # sum counterfactuals within each respective group #
      seq_sum <- c(seq(1,length(tmp), by = length(A_span)),length(tmp) + 1)
      preds <- sapply(1:(length(seq_sum) - 1), function(i) {
        sum(tmp[seq_sum[i]:(seq_sum[i + 1] - 1)] * binsizes)
      })
      
    } else{
      preds <- sapply(seq(nrow(H)),
                      function(i) {
                        integrand <- function(a){
                          # print(a)
                          newXX <- data.frame(A = a[1,], H[i, , drop = FALSE])
                          if(get_la) {
                            t_temp <- sum(A %in% colnames(newXX))
                            la <- SL_ladata_int((t_temp+1), newXX, A, lambda, k, trunc)
                          } else la <- 1
                          tmp <- matrix(predict(fit, newdata = newXX)$pred[,1] * la, nrow = 1)
                          # print(dim(tmp))
                          return(tmp)
                        }
                        
                        out <-  cubature::hcubature(integrand,
                                                    lower = A.min,
                                                    upper = A.max,
                                                    vectorInterface = TRUE)
                        return(out$integral)
                      })
    }
  }
  
  
  ##############
  return(preds)
}



#' Jac
#'
#' Function to estimate Jacobian matrix of estimators using user provided stabilizing lambda
#' @param b value/vector of estimated values of U1+U2
#' @return Jacobian matrix
#' @export
Jac <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int){
  
  if(outcome == "binomial"){
    
    if(approx_int){
      fun12 <- function(a, A, lambda, phi, m, trunc, tau, la){
        pp <- phi(a)[,1] * phi(a)[,2]
        mm <- m(a, b, phi)
        return(as.matrix(pp * mm * (1 - mm) * la, nrow = 1))
      }
      fun11 <- function(a, A, lambda, phi, m, trunc, tau, la){
        pp <- phi(a)[,1]^2
        mm <- m(a, b, phi)
        return(as.matrix(pp * mm * (1 - mm) * la, nrow = 1))
      }
      fun22 <- function(a, A, lambda, phi, m, trunc, tau, la){
        pp <- phi(a)[,2]^2
        mm <- m(a, b, phi)
        return(as.matrix(pp * mm * (1 - mm) * la, nrow = 1))
      }
      
      # CHANGE INTEGRALS TO SUMS OVER DOMAIN OF A #
      U12 <- sum(fun12(A_span, A = A, lambda = lambda, phi = phi, m = m, trunc = trunc, tau = tau, la = la)*
                   bindim)
      U11 <- sum(fun11(A_span, A = A, lambda = lambda, phi = phi, m = m, trunc = trunc, tau = tau, la = la)*
                   bindim)
      U22 <- sum(fun22(A_span, A = A, lambda = lambda, phi = phi, m = m, trunc = trunc, tau = tau, la = la)*
                   bindim)
      
      U2h  <- matrix(c(U11, U12, U12, U22), ncol = 2)
      return(- U2h)
    } else{
      fun12 <- function(a){
        pp <- phi(t(a))[,1] * phi(t(a))[,2]
        mm <- m(t(a), b, phi)
        if(la == 1) {
          la_temp <- 1 # do.call('cbind',mclapply(1:tau,function(j)SL_ladata(j, temp, A_name = A, lambda, k, trunc),
          # mc.cores = cores))
        } else la_temp <- 1
        return(as.matrix(pp * mm * (1 - mm) * la_temp, nrow = 1))
      }
      fun11 <- function(a){
        pp <- phi(t(a))[,1]^2
        mm <- m(t(a), b, phi)
        if(la) {
          la_temp <- 1 # do.call('cbind',mclapply(1:tau,function(j)SL_ladata(j, temp, A_name = A, lambda, k, trunc),
          # mc.cores = cores))
        } else la_temp <- 1
        return(as.matrix(pp * mm * (1 - mm) * la_temp, nrow = 1))
      }
      fun22 <- function(a){
        pp <- phi(t(a))[,2]^2
        mm <- m(t(a), b, phi)
        if(la) {
          la_temp <- 1 # do.call('cbind',mclapply(1:tau,function(j)SL_ladata(j, temp, A_name = A, lambda, k, trunc),
          # mc.cores = cores))
        } else la_temp <- 1
        return(as.matrix(pp * mm * (1 - mm) * la_temp, nrow = 1))
      }
      U12 <- cubature::hcubature(fun12,
                                 lower = rep(A.min, tau),
                                 upper = rep(A.max, tau),
                                 vectorInterface = TRUE)
      U11 <- cubature::hcubature(fun11,
                                 lower = rep(A.min, tau),
                                 upper = rep(A.max, tau),
                                 vectorInterface = TRUE)
      U22 <- cubature::hcubature(fun22,
                                 lower = rep(A.min, tau),
                                 upper = rep(A.max, tau),
                                 vectorInterface = TRUE)
      
      U2h  <- matrix(c(U11$integral, U12$integral, U12$integral, U22$integral), ncol = 2)
      return(- U2h)
    }
  } else if(outcome == "gaussian"){
    
    if(approx_int){
      fun12 <- function(a, A, lambda, phi, m, trunc, tau, la){
        pp <- phi(a)[,1] * phi(a)[,2]
        return(as.matrix(pp * la, nrow = 1))
      }
      fun11 <- function(a, A, lambda, phi, m, trunc, tau, la){
        pp <- phi(a)[,1]^2
        return(as.matrix(pp * la, nrow = 1))
      }
      fun22 <- function(a, A, lambda, phi, m, trunc, tau, la){
        pp <- phi(a)[,2]^2
        return(as.matrix(pp * la, nrow = 1))
      }
      
      # CHANGE INTEGRALS TO SUMS OVER DOMAIN OF A #
      U12 <- sum(fun12(A_span, A = A, lambda = lambda, phi = phi, m = m, trunc = trunc, tau = tau, la = la)*
                   bindim)
      U11 <- sum(fun11(A_span, A = A, lambda = lambda, phi = phi, m = m, trunc = trunc, tau = tau, la = la)*
                   bindim)
      U22 <- sum(fun22(A_span, A = A, lambda = lambda, phi = phi, m = m, trunc = trunc, tau = tau, la = la)*
                   bindim)
      
      U2h  <- matrix(c(U11, U12, U12, U22), ncol = 2)
      return(- U2h)
    } else{
      fun12 <- function(a){
        pp <- phi(t(a))[,1] * phi(t(a))[,2]
        mm <- m(t(a), b, phi)
        if(la == 1) {
          la_temp <- 1 # do.call('cbind',mclapply(1:tau,function(j)SL_ladata(j, temp, A_name = A, lambda, k, trunc),
          # mc.cores = cores))
        } else la_temp <- 1
        return(as.matrix(pp * mm * (1 - mm) * la_temp, nrow = 1))
      }
      fun11 <- function(a){
        pp <- phi(t(a))[,1]^2
        mm <- m(t(a), b, phi)
        if(la) {
          la_temp <- 1 # do.call('cbind',mclapply(1:tau,function(j)SL_ladata(j, temp, A_name = A, lambda, k, trunc),
          # mc.cores = cores))
        } else la_temp <- 1
        return(as.matrix(pp * mm * (1 - mm) * la_temp, nrow = 1))
      }
      fun22 <- function(a){
        pp <- phi(t(a))[,2]^2
        mm <- m(t(a), b, phi)
        if(la) {
          la_temp <- 1 # do.call('cbind',mclapply(1:tau,function(j)SL_ladata(j, temp, A_name = A, lambda, k, trunc),
          # mc.cores = cores))
        } else la_temp <- 1
        return(as.matrix(pp * mm * (1 - mm) * la_temp, nrow = 1))
      }
      U12 <- cubature::hcubature(fun12,
                                 lower = rep(A.min, tau),
                                 upper = rep(A.max, tau),
                                 vectorInterface = TRUE)
      U11 <- cubature::hcubature(fun11,
                                 lower = rep(A.min, tau),
                                 upper = rep(A.max, tau),
                                 vectorInterface = TRUE)
      U22 <- cubature::hcubature(fun22,
                                 lower = rep(A.min, tau),
                                 upper = rep(A.max, tau),
                                 vectorInterface = TRUE)
      
      U2h  <- matrix(c(U11$integral, U12$integral, U12$integral, U22$integral), ncol = 2)
      return(- U2h)
    }
  }
  
}


#' U2fun
#'
#' Function to estimate U1 + U2 using stabilizing lambda and vectorizing process (for quicker processing)
#' @param b value/vector of starting values (from IPW fit)
#' @return Vector of length b as solutions of U1 + U2
#' @export
U2fun <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int) {
  
  if(A_type == "categorical"){
    fU2 <- function(a, A, lambda_int, phi_int, m_int, trunc, tau, outcome, la) {
      pp <- phi_int(a)
      mm <- m_int(a, b, phi_int)
      return( - mm * pp * la)
    }
    
    out <- apply(fU2(a = A_span, A = A, lambda_int = lambda, phi_int = phi, m_int = m, trunc = trunc, tau = tau, outcome = outcome, la = la),2,sum)
    return(out)
    
  } else if(A_type == "continuous"){
    
    if(approx_int){
      
      fU2 <- function(a, A, lambda_int, phi_int, m_int, trunc, tau, outcome, la) {
        pp <- phi_int(a)
        mm <- m_int(a, b, phi_int)
        return( - mm * pp * la)
      }
      
      out <- apply(fU2(a = A_span, A = A, lambda_int = lambda, phi_int = phi,
                       m_int = m, trunc = trunc, tau = tau, outcome = outcome,  la = la) *
                     bindim,2,sum)
      return(out)
      
    } else{
      
      fU2 <- function(a) {
        pp <- c(1, sum(a))
        mm <- plogis(pp %*% b)[,1]
        return( - mm * pp * la)
      }
      
      fint1 <- function(a) fU2(a)[1]
      fint2 <- function(a) fU2(a)[2]
      # CAN VECTORIZE THE FUNCTIONS FOR FASTER COMPUTATION
      out1 <- cubintegrate(fint1, lower = rep(A.min, tau), upper = rep(A.max, tau), method = "cuhre")
      out2 <- cubintegrate(fint2, lower = rep(A.min, tau), upper = rep(A.max, tau), method = "cuhre")
      return(c(out1$integral, out2$integral))
      
    }
    
  }
  
}




#' get_int_info
#'
#' Generate information for reimann integration
#' @param data data used
#' @param A vector of variable names
#' @param tau number of variables
#' @param quant_bin bin size for quantiles
#' @return output of funciton without console output
#' @export

get_int_info <- function(data, A, tau, quant_bin = 0.01, A.min = 0, A.max = 1){
  
  # get quantiles #
  quants <- sapply(1:tau,
                   function(j)
                     as.numeric(quantile(data[,A[j]], probs = seq(0, 1, by = quant_bin))))
  
  
  
  # Compute middle points and differences (bin sizes)
  A_diffs <- lapply(1:tau, function(j) {
    # Middle points
    mids <- sapply(2:length(quants[, j]), function(i) {
      (quants[i - 1, j] + quants[i, j]) / 2
    })
    
    # Differences (bin sizes)
    diffs <- diff(quants[, j])
    
    list(middle = mids, diff = diffs)
  })
  
  # Extract only middle points to create the grid
  A_span <- as.matrix(expand.grid(
    lapply(A_diffs, function(item) item$middle)
  )) %>% as.data.frame()
  
  # Create a similar grid for bin sizes
  binsizes <- expand.grid(
    lapply(A_diffs, function(item) item$diff)
  ) %>% as.data.frame()
  
  # dimensions to adjust for #
  bindim <- apply(binsizes,1,prod)
  
  
  return(list(A_span = A_span, binsizes = binsizes, bindim = bindim, quants = quants))
  
  
  # # old generation #
  # # generate A_span from quantiles #
  # A_span <- as.matrix(expand.grid(
  #   lapply(1:tau,function(j) quants[,j]
  #   ))) %>% as.data.frame()
  #
  # # create corresponding bins #
  # # Calculate midpoints between quantiles for each variable and get size of each bin #
  # binsizes <- as.matrix(expand.grid(
  #   lapply(1:tau,function(j) {
  #     tmp <- as.numeric(c(quantile(data[,A[j]],0),
  #                         sapply(2:length(quants[,j]), function(i) {
  #                           (quants[i-1,j] + quants[i,j]) / 2
  #                         }),
  #                         quantile(data[,A[j]],1)))
  #     sapply(2:length(tmp), function(i) {
  #       tmp[i] - tmp[i-1]
  #     })
  #   }
  #   ))) %>% as.data.frame()
  #
  # # dimensions to adjust for #
  # bindim <- apply(binsizes,1,prod)
  
}



#' #' perform_checks
#' #'
#' #' Checks the arguments of the function
#' perform_checks <- function(args_list){
#'   if(sum(c("sdr","tml") %in% est) == 0)
#'     stop('Please select at least one of the available estimators: c("sdr","tml").')
#'   if(!is.function(phi)) stop("Argument phi must be a function.")
#'   # if(!is.function(m)) stop("Argument m must be a function.")
#'   if(!(isTRUE(lambda) || isFALSE(lambda) || is.function(lambda)))
#'     if(!is.list(lambda)){ stop("Argument lambda must be a boolean, a function or a list.")
#'     } else if(length(lambda) != tau)
#'       warning("Length of lambda fits different than tau, default multinomial method will be used instead.")
#'   if(nfolds < 1 || round(nfolds) != nfolds) stop("nfolds must be a positive integer.")
#'   if(nfolds < 2) warning("In order to meet the Donsker class requirements cross-fitting must be performed. Note that your estimate may not have a causal interpretation.")
#'   if(length(stackr) < 1){
#'     warning("Argument stackr was empty, the default set of learners will be used instead.")
#'     stackr = c("mean","randomforest","glmnet","knn","lightgbm","rpart","xgboost",
#'                "naivebayes","svm","lda","nnet","ranger","multinom")
#'   }
#'   if(length(stackm) < 1){
#'     warning("Argument stackm was empty, the default set of learners will be used instead.")
#'     stackm = c('SL.mean','SL.glm','SL.caretEarth','SL.glmnet',
#'                'SL.glm.interaction','SL.gam','SL.randomForest',
#'                'SL.rpart','SL.bayesglm',"SL.svm")
#'   }
#'
#'   if(outcome == "gaussian"){
#'     if(approx_int) warning("Approximate integration will be used, results may not be accurate. Ensure that increment covers important regions of your exposure.")
#'   }
#'   if(cores < 1 || round(cores) != cores) stop("cores must be a positive integer.")
#'
#'   if(!(U2_method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent")))
#'     warning("U2_method for `optim()` not recognized, defaulting to Nelder-Mead.")
#'
#'   # add warning if A names different than those provided in fits lambda set lambda to TRUE instead #
#'
#'   ### data checks and fix ###
#'   # Checks #
#'   if(length(A) != tau) stop("The number of treatments did not equal the number of time points.")
#'   if(length(L) != tau) stop("The number of covariates did not equal the number of time points.")
#'   if(is.na(match(Y, colnames(data)))) stop("The outcome variable specified was not found in the data.")
#'   # renaming #
#'   colnames(data)[match(Y, colnames(data))] <- "Y"
#'   if(rescale_outcome)
#'     data$Y <- (data$Y - min(data$Y)) / (max(data$Y) - min(data$Y))
#'   data <- as.data.frame(data)
#'
#'   # creating folds #
#'   if(is.null(valSets))
#'     valSets <- split(sample(seq_len(nrow(data))), rep(seq_len(nfolds), length = nrow(data)))
#'
#'   if(!(is.logical(trim) || is.numeric(trim))) stop("trim must be a boolean or a numeric value.")
#'   if(!is.logical(trim) && (trim < 0 || trim > 1)) stop("trim must be a positive value between 0 and 1.")
#'
#'   # update values #
#'   # set function m as function of the outcome #
#'   if(outcome == "binomial"){
#'     m = function(a, b, phi) plogis(phi(a) %*% b)[,1]
#'   } else if(outcome == "gaussian") {
#'     m = function(a, b, phi) (phi(a) %*% b)[,1]
#'   }
#'
#'   args_list_out <- list(data = data, tau = tau, A = A, L = L, Y = Y, est = est,
#'                         outcome = outcome, A_type = A_type,
#'                         phi = phi,
#'                         lambda = lambda, nfolds = nfolds,
#'                         valSets = valSets,
#'                         stackr = stackr,
#'                         stackm = stackm, k = k,
#'                         trunc = trunc, trim = trim,
#'                         scale_weights = scale_weights, rescale_outcome = rescale_outcome,
#'                         merge_results = merge_results, trace = trace,
#'                         int_bounds = int_bounds, increment = increment,
#'                         cores = cores, U2_method = U2_method,
#'                         approx_int = approx_int, ken.dens = ken.dens, force_la_1 = force_la_1, m = m)
#'
#'   return(args_list_out)
#' }


##################################






################ MULTICLASS ESTIMATION ###################

#' mult_densities
#'
#' Function to estimate densities of categorical treatment values accross time points as a function of previous
#' history using mlr3superlearner package.
#' @param j time point
#' @param data data used
#' @param stack SuperLeaner libraries to be used for fitting
#' @param valSets sets for cross-fitting
#' @param k Last time point to consider when getting prior history
#' @param trunc boolean specifying if estimates should be truncated.
#' @return estimate of probability mass for observed treatment
#' @export
#' @import
#' mlr3superlearner
#' mlr3extralearners
#' SuperLearner
#' dplyr

mult_densities <- function(j, A_name, L_name, data,
                           stack , valSets, k = Inf, trunc, seed = 123) {
  
  
  set.seed(seed)
  
  # a <- unlist(data[, A_name[j]])
  a <- data[, A_name[j]]
  n <- length(a)
  
  # select all prior trt/covariate history
  if(j > 1) namesA <- A_name[max(1, (j-k)):(j-1)] else namesA <- NULL
  namesL <- unlist(L_name[max(1, j-k+1):j])
  # stack in a single dataset
  l <- data %>% select(one_of(c(namesA, namesL)))
  temp_dat <- as.data.frame(cbind(a,l)) %>%
    mutate(a = as.character(a))
  
  # case with single treatment at given time point #
  if(length(unique(temp_dat$a)) == 1)
    return(rep(1, n))
  
  
  pihat <- c()
  for(i in 1:length(valSets)){
    if(length(valSets) > 1){
      temp_dat_train <- temp_dat[-valSets[[i]],]
      temp_dat_test <- temp_dat[valSets[[i]],]
    } else{
      temp_dat_train <- temp_dat_test <- temp_dat[valSets[[i]],]
    }
    if(j > 1)
      fit <- hush(mlr3superlearner(data = temp_dat_train, target = "a",outcome_type = "multiclass",
                                   library = stack,
                                   folds = 5,discrete = T,info = F))
    
    if(j == 1){
      
      # if(is.list(stack)){
      #
      # } else{
      #   if(sum(stack != "SL.mean") > 1 && length(grep("glmnet",stack)) > 0)
      #     stack <- stack[-grep("glmnet",stack)]
      # }
      if(sum(stack != "SL.mean") > 1 && length(grep("glmnet",stack)) > 0)
        stack <- stack[-grep("glmnet",stack)]
      fit <- hush(mlr3superlearner(data = temp_dat_train, target = "a",outcome_type = "multiclass",
                                   library = stack,
                                   folds = 5,discrete = T,info = F))
    }
    
    preds <- predict(fit, newdata = temp_dat_test)
    counter <- 1
    for(r in valSets[[i]]){
      temp_preds <- preds[counter,]
      if(trunc) pihat[r] <- bound(as.numeric(temp_preds[(a[r] + 1)]), p = 0.01) else pihat[r] <- as.numeric(temp_preds[(a[r] + 1)])
      counter <- counter + 1
    }
  }
  return(pihat)
  
}



################ STABILIZING LAMBDA ################

#' lambda_cust
#'
#' Function to estimate lambdas the values to be used to stabilize weights, requires use to provide a list of
#' fits of length equal to the number of time points of interest. This can be done by providing a multinomial fit
#' or a mlr3superlearner fit. At the time point of the first treatment empirical means of each treamtent will be used.
#' @param j time point
#' @param data data used
#' @param k Last time point to consider when getting prior history
#' @param fits list of length equal to number of time points of interest containing fits to be used to predict lambdas for dataset
#' @return estimate of probability mass for observed treatment based on prior treatment history (lambda)
#' @export
#' @import dplyr
lambda_cust <- function(j, A_name, data, k = Inf, fits, trunc = F, A_type) {
  
  # nameA <- paste0('A_', j, sep = '')
  a <- data[, A_name[j]]
  n <- length(a)
  temp_fit <- fits[[j]]
  
  if(A_type == "categorical"){
    if((length(unique(fits[[j]]))  == 1) || (length(unique(a))  == 1))
      return(rep(1, n))
    
    
    if(j == 1){
      pihat <- c()
      for(r in 1:n){
        pihat[r] <- temp_fit[(a[r] + 1)]
      }
      return(pihat)
    }
    
    # select all prior trt/covariate history
    namesA <- A_name[max(1, (j-k)):(j-1)] # paste('A_', max(1, (j-k)):(j-1), sep = '') else namesA <- NULL
    
    l <- data %>% select(one_of(c(namesA)))
    preds <- predict(temp_fit, type = "probs", newdata = l)
    
    pihat <- c()
    for(r in 1:n){
      if(!is.null(dim((preds)) )){
        temp_preds <- preds[r,]
      } else temp_preds <- preds
      pihat[r] <- as.numeric(temp_preds[(a[r] + 1)])
    }
    
    if(trunc) pihat <- bound(pihat, p = 0.01)
    return(pihat)
  }
  
  
  if(A_type == "continuous"){
    # select all prior trt/covariate history
    namesA <- A_name[max(1, (j-k)):(j-1)]
    
    if(j == 1){
      # make dat with only current and prior treatments #
      l <- data %>% select(one_of(c(namesA,A_name[j]))) %>%
        rename(A = A_name[j])
      # Predict using the model
      dens_a <- approx(temp_fit$density_model$x, temp_fit$density_model$y, xout = a)$y
      # dens_a <- predict(temp_fit$density_model, newdata = l)
    } else {
      l <- data %>% select(one_of(c(namesA)))
      
      if(ncol(l) == 1) {
        temp <- colnames(l)
        l <- cbind(0,l)
        colnames(l) <- c("Intercept",temp)
      }
      # pimod.vals <- predict(temp_fit$pimod,newdata = l)$pred
      pimod.vals <- predict(temp_fit$pimod,newdata = l)
      sq.res <- (a-pimod.vals[1:n])^2
      
      # pi2mod.vals <-predict(temp_fit$pi2mod,newdata = l)$pred
      pi2mod.vals <-predict(temp_fit$pi2mod,newdata = l)
      
      approx.fn     <- function(x,y,z){ predict(smooth.spline(x,y),x=z)$y }
      # construct estimated pi/varpi and mu/m values
      a.std         <- (a-pimod.vals)/sqrt(abs(pi2mod.vals))
      pihat.vals    <- approx.fn(density(na.omit(a.std[1:n]))$x, density(na.omit(a.std)[1:n])$y,na.omit(a.std[1:n]))
      dens_a         <- (pihat.vals/sqrt(abs(pi2mod.vals)))[1:n]
      
    }
    
    return(dens_a)
  }
  
}


### Add function for default stabilizing weights with multinomial ###
#' lambda_est
#'
#' Function to estimate lambdas the values to be used to stabilize weights, requires use to provide a list of
#' fits of length equal to the number of time points of interest. This can be done by providing a multinomial fit
#' or a mlr3superlearner fit. At the time point of the first treatment empirical means of each treamtent will be used.
#' @param j time point
#' @param data data used
#' @param k Last time point to consider when getting prior history
#' @param fits list of length equal to number of time points of interest containing fits to be used to predict lambdas for dataset
#' @return estimate of probability mass for observed treatment based on prior treatment history (lambda)
#' @export
#' @import dplyr nnet
lambda_est <- function(j,A_name , data, k){
  
  print(j)
  a <- data[, A_name[j]]
  n <- length(a)
  
  if(length(unique(a)) == 1)
    return(rep(1, n))
  
  if(j == 1){
    preds <- summary(as.factor(a))/n
    fit <- preds
  }
  
  else{
    namesA <- A_name[max(1, (j-k)):(j-1)]# paste('A_', max(1, (j-k)):(j-1), sep = '')
    l <- data %>% select(one_of(c(namesA)))
    
    temp_dat <- as.data.frame(cbind(a, l))
    fit <- multinom(a ~ ., data = temp_dat)
    fit$residuals <- NULL
    fit$weights <- NULL
  }
  return(fit)
}




################ CONTINUOUS TREATMENT #############
### Not currently available ###

#' OR_densities
#'
#' Function to estimate densities of continuous treatment values accross time points as a function of previous
#' history using data duplication and odds ratios.
#' @param j time point
#' @param data data used
#' @param stack SuperLeaner libraries to be used for fitting
#' @param valSets sets for cross-fitting
#' @param k Last time point to consider when getting prior history
#' @return Odds ratio as estimate of density
OR_densities <- function(j, data, A_name, L_name, stack, valSets, k = Inf) {
  
  A <- data[, A_name[j]]
  
  # select all prior trt/covariate history
  namesL <- unlist(L_name[max(1, j-k+1):j])
  if(j > 1) {
    namesA <- A_name[max(1, (j-k)):(j-1)]
    H <- data %>% select(one_of(c(namesA, namesL)))# %>% mutate_all(as.factor)
  } else H <- data %>% select(one_of(c(namesL)))# %>% mutate_all(as.factor)
  
  # stack in a single dataset
  # H <- data %>% select(one_of(c(namesA, namesL)))# %>% mutate_all(as.factor)
  
  n <- length(A)
  
  # create random set of treatment values
  Ad <- runif(n)
  id <- rep(1:n, 2)
  Delta <- c(rep(0, n), rep(1, n))
  # duplicate data #
  X <- data.frame(A = c(Ad, A), rbind(H, H))
  
  # fit model with outcome being true observed vs random generated #
  fit <- SuperLearner(Y = Delta, X = X, newX = X,family = binomial(), SL.library = stack)
  # fit <- SuperLearner(Y = Delta, X = X, newX = X,family = gaussian(), SL.library = stack)
  # fit <- mySL(Y = Delta, X = X, id = id, SL.library = stack, family = binomial(), validRows = valSets)
  
  # get prediction for outcome (in the group that had the observed treatment)
  u <- bound(fit$SL.predict[Delta == 1, 1], p = 0.01)
  
  # this is the density? #
  rat <- u  / (1 - u)
  
  return(rat)
  
}



#' SL_densities
#'
#' Function to estimate densities of continuous treatment values accross time points as a function of previous
#' history using the method described by Kennedy et al.
#' @param j time point
#' @param data data used
#' @param stack SuperLeaner libraries to be used for fitting
#' @param valSets sets for cross-fitting
#' @param k Last time point to consider when getting prior history
#' @return estimate of density
SL_densities <- function(j, data, A_name, L_name, stack, valSets, k = Inf, trunc, seed = 1, length.out = 100) {
  
  set.seed(seed)
  # print(paste("Dens: ",j))
  a <- data[, A_name[j]]
  n <- length(a)
  sl.lib_gps <- stack
  
  # select all prior trt/covariate history
  namesL <- unlist(L_name[max(1, j-k+1):j])
  if(j > 1) {
    namesA <- A_name[max(1, (j-k)):(j-1)]
    l <- data %>% select(one_of(c(namesA, namesL)))# %>% mutate_all(as.factor)
  } else {
    l <- data %>% select(one_of(c(namesL)))# %>% mutate_all(as.factor)
    if('SL.glmnet' %in% sl.lib_gps) {
      sl.lib_gps <- sl.lib_gps[!sl.lib_gps %in% 'SL.glmnet']
    }
  }
  
  # set up evaluation points & matrices for predictions
  n <- length(a)
  # set up evaluation points & matrices for predictions
  a.min <- min(a); a.max <- max(a)
  # generate a set of potential values for treatment falling between the bounds of the observed treatment
  a.vals <- seq(a.min,a.max,length.out=length.out)
  temp <- as.data.frame(cbind( l[rep(1:n,length(a.vals)),],a=rep(a.vals,rep(n,length(a.vals)))))
  colnames(temp) <- colnames(cbind(l,a))
  la.new            <- data.frame(rbind(cbind(l,a),  temp))
  res_GP            <- 1
  l.new             <- as.data.frame(la.new[,-dim(la.new)[2]])
  colnames(l.new) <- colnames(l)
  
  ###### GPS
  #print("SuperLearner for pi")
  pimod <- SuperLearner(Y=a, X=l, SL.library=sl.lib_gps, newX = l.new, family=gaussian)
  pimod.vals <- as.numeric(pimod$SL.predict[,1])
  
  
  #print("SuperLearner for pi2")
  sq.res <- (a-pimod.vals[1:n])^2
  pi2mod <- SuperLearner(Y=sq.res,X=l, SL.library=sl.lib_gps, newX=l.new,family=gaussian)
  pi2mod.vals <- as.numeric(pi2mod$SL.predict[,1])
  
  approx.fn     <- function(x,y,z){ predict(smooth.spline(x,y),x=z)$y }
  
  ## Current ##
  a.std         <- (a-pimod.vals)/sqrt(pi2mod.vals)
  pihat.vals    <- approx.fn(density(na.omit(a.std[1:n]))$x, density(na.omit(a.std)[1:n])$y,na.omit(a.std))
  pihat         <- (pihat.vals/sqrt(pi2mod.vals))[1:n]
  
  if(trunc) pihat <- truncate(pihat, p = 0.001)
  
  return(pihat)
  
}



#' SL_densities_cf
#'
#' Function to estimate densities of continuous treatment values accross time points as a function of previous
#' history using the method described by Kennedy et al.
#' @param j time point
#' @param data data used
#' @param stack SuperLeaner libraries to be used for fitting
#' @param valSets sets for cross-fitting
#' @param k Last time point to consider when getting prior history
#' @return estimate of density
SL_densities_cf <- function(j, data, A_name, L_name, stack, valSets, k = Inf, trunc, seed = 1, length.out = 100) {
  
  set.seed(seed)
  sl.lib_gps <- stack
  
  pihat_out <- c()
  for(i in 1:length(valSets)){
    if(length(valSets) > 1){
      temp_dat_train <- data[-valSets[[i]],]
      temp_dat_test <- data[valSets[[i]],]
    } else{
      n = nrow(data)
      valSets <- list(1:n)
      temp_dat_train <- temp_dat_test <- data[valSets[[i]],]
    }
    
    a <- temp_dat_train[, A_name[j]]
    a.test <- temp_dat_test[, A_name[j]]
    
    # select all prior trt/covariate history
    namesL <- unlist(L_name[max(1, j-k+1):j])
    if(j > 1) {
      namesA <- A_name[max(1, (j-k)):(j-1)]
      l <- temp_dat_train %>% select(one_of(c(namesA, namesL)))
      l.test <- temp_dat_test %>% select(one_of(c(namesA, namesL)))
    } else {
      l <- temp_dat_train %>% select(one_of(c(namesL)))
      l.test <- temp_dat_test %>% select(one_of(c(namesL)))
      if('SL.glmnet' %in% sl.lib_gps) {
        sl.lib_gps <- sl.lib_gps[!sl.lib_gps %in% 'SL.glmnet']
      }
    }
    
    # set up evaluation points & matrices for predictions
    n <- length(a)
    n.test <- length(a.test)
    a.min <- min(data[, A_name[j]]); a.max <- max(data[, A_name[j]])
    # generate a set of potential values for treatment falling between the bounds of the observed treatment
    a.vals <- seq(a.min,a.max,length.out=length.out) # 500
    
    temp <- as.data.frame(cbind( l.test[rep(1:n.test,length(a.vals)),],a=rep(a.vals,rep(n.test,length(a.vals)))))
    colnames(temp) <- colnames(cbind(l.test,a.test))
    la.new            <- data.frame(rbind(cbind(l.test,a.test),  temp))
    res_GP            <- 1
    l.new             <- as.data.frame(la.new[,-dim(la.new)[2]])
    colnames(l.new) <- colnames(l)
    
    ###### GPS
    #print("SuperLearner for pi")
    pimod <- SuperLearner(Y=a, X=l, SL.library=sl.lib_gps, newX = l.new, family=gaussian)
    pimod.vals.test <- as.numeric(pimod$SL.predict[,1])
    pimod.vals <- predict(pimod, newdata=l)$pred[,1]
    
    
    #print("SuperLearner for pi2")
    sq.res <- (a-pimod.vals[1:n])^2
    pi2mod <- SuperLearner(Y=sq.res,X=l, SL.library=sl.lib_gps, newX=l.new,family=gaussian)
    pi2mod.vals.test <- as.numeric(pi2mod$SL.predict[,1])
    pi2mod.vals <- predict(pi2mod, newdata=l)$pred[,1]
    
    approx.fn     <- function(x,y,z){ predict(smooth.spline(x,y),x=z)$y }
    
    if(sum(pi2mod.vals.test[1:n.test] < 0) > 0) warning("Negative variance estimate. ")
    
    ## Current ##
    a.std.test         <- (a.test-pimod.vals.test)/sqrt(abs(pi2mod.vals.test))
    a.std         <- (a-pimod.vals)/sqrt(abs(pi2mod.vals))
    pihat.vals    <- approx.fn(density(na.omit(a.std[1:n]))$x, density(na.omit(a.std[1:n]))$y,na.omit(a.std.test))
    pihat         <- (pihat.vals/sqrt(abs(pi2mod.vals.test)))[1:n.test]
    
    pihat_out[valSets[[i]]] <- pihat
    
  }
  
  if(trunc) pihat_out <- truncate(pihat_out, p = 0.001)
  
  return(pihat_out)
  
}









#' SL_densities_lambda
#'
#' Function to estimate densities of continuous treatment values accross time points as a function of previous
#' history using the method described by Kennedy et al.
#' @param j time point
#' @param data data used
#' @param stack SuperLeaner libraries to be used for fitting
#' @param valSets sets for cross-fitting
#' @param k Last time point to consider when getting prior history
#' @return estimate of density
SL_densities_lambda <- function(j, data, A_name, stack, k = Inf, length.out = 500) {
  
  sl.lib_gps <- stack
  a <- data[, A_name[j]]
  n <- length(a)
  
  # select all prior trt/covariate history
  namesA <- A_name[max(1, (j-k)):(j-1)]
  
  
  if(j == 1){
    l <- data %>% select(one_of(c(namesA,A_name[j]))) %>%
      rename(A = A_name[j])
    
    
    density_model <- density(l$A)
    return(list("density_model" = density_model))
    
    # a.min <- min(a); a.max <- max(a)
    # # generate a set of potential values for treatment falling between the bounds of the observed treatment
    # a.vals <- seq(a.min,a.max,length.out=length.out)
    # temp <- as.data.frame(cbind( l[rep(1:n,length(a.vals)),],a=rep(a.vals,rep(n,length(a.vals)))))
    # colnames(temp) <- colnames(cbind(l,a))
    # la.new            <- data.frame(rbind(cbind(l,a),  temp))
    # res_GP            <- 1
    # l.new             <- as.data.frame(la.new[,-dim(la.new)[2]])
    # colnames(l.new) <- colnames(l)
    #
    #
    # pimod <- gam(A ~ 1, data = l)
    # pimod.vals <- predict(pimod, newdata = l.new)
    #
    # sq.res <- (a-pimod.vals[1:n])^2
    # l$sq.res <- sq.res
    # pi2mod <- gam(sq.res ~ 1, data = l)
    #
    # return(list("pimod" = pimod,
    #             "pi2mod"= pi2mod))
    # bw <- npudensbw(~A, data = l)
    # density_model <- npudens(bw)
    
    # density_model <- density(l$A)
    # return(list("density_model" = density_model))
  } else{
    l <- data %>% select(one_of(c(namesA)))
    
    # set up evaluation points & matrices for predictions
    a.min <- min(a); a.max <- max(a)
    # generate a set of potential values for treatment falling between the bounds of the observed treatment
    a.vals <- seq(a.min,a.max,length.out=length.out)
    temp <- as.data.frame(cbind( l[rep(1:n,length(a.vals)),],a=rep(a.vals,rep(n,length(a.vals)))))
    colnames(temp) <- colnames(cbind(l,a))
    la.new            <- data.frame(rbind(cbind(l,a),  temp))
    res_GP            <- 1
    l.new             <- as.data.frame(la.new[,-dim(la.new)[2]])
    colnames(l.new) <- colnames(l)
    
    # make formula #
    formula_pimod <- as.formula(paste0("a~",paste0("s(",colnames(l),")",collapse = "+")))
    pimod <- gam(formula_pimod, data = l)
    pimod.vals <- predict(pimod, newdata = l.new)
    
    sq.res <- (a-pimod.vals[1:n])^2
    formula_pi2mod <- as.formula(paste0("sq.res~",paste0("s(",colnames(l),")",collapse = "+")))
    l$sq.res <- sq.res
    pi2mod <- gam(formula_pi2mod, data = l)
    
    
    return(list("pimod" = pimod,
                "pi2mod"= pi2mod))
    
    # ###### GPS
    # #print("SuperLearner for pi")
    # if(ncol(l) == 1){
    #   X <- cbind(0,l)
    #   newX <- cbind(0,l.new)
    #   colnames(X) <- colnames(newX) <- c("Intercept",colnames(l))
    #   pimod <- SuperLearner(Y=a, X=X, SL.library=sl.lib_gps, newX = newX, family=gaussian)
    # } else {
    #   pimod <- SuperLearner(Y=a, X=l, SL.library=sl.lib_gps, newX = l.new, family=gaussian)
    # }
    # pimod.vals <- pimod$SL.predict
    #
    # #print("SuperLearner for pi2")
    # sq.res <- (a-pimod.vals[1:n])^2
    # if(ncol(l) == 1){
    #   X <- cbind(0,l)
    #   newX <- cbind(0,l.new)
    #   colnames(X) <- colnames(newX) <- c("Intercept",colnames(l))
    #   pi2mod <- SuperLearner(Y=sq.res,X=X, SL.library=sl.lib_gps, newX=newX,family=gaussian)
    # } else {
    #   pi2mod <- SuperLearner(Y=sq.res,X=l, SL.library=sl.lib_gps, newX=l.new,family=gaussian)
    # }
    #
    # return(list("pimod" = pimod,
    #             "pi2mod"= pi2mod))
  }
  
}

#' SL_densities_ladata
#' Function to estimate densities of continuous treatment values accross time points as a function of previous
SL_densities_ladata <- function(j, data, A_name, lambda , k = Inf, trunc) {
  
  # a <- data[, A_name[j]]
  # n <- length(a)
  #
  # # select all prior trt/covariate history
  # namesA <- A_name[max(1, (j-k)):(j-1)]
  #
  # approx.fn <- function(x,y,z){ predict(smooth.spline(x,y),x=z)$y }
  #
  #
  #
  # l <- data %>% select(one_of(c(namesA)))
  #
  # if(ncol(l) == 1) {
  #   temp <- colnames(l)
  #   l <- cbind(0,l)
  #   colnames(l) <- c("Intercept",temp)
  # }
  #
  #
  # # rm(lambda)
  # # lambda <- lambda_fits
  # # lambda[[j]]$pimod$data <- lambda[[j]]$pimod$xlevels <- lambda[[j]]$pimod$R <- lambda[[j]]$pimod$residuals <-
  # #   lambda[[j]]$pimod$fitted.values <- lambda[[j]]$pimod$weights <- lambda[[j]]$pimod$var <-
  # #   lambda[[j]]$pimod$additive.predictors <- lambda[[j]]$pimod$prior.weights <- lambda[[j]]$pimod$y <- NULL
  #
  #
  # pimod.vals <- predict(lambda[[j]]$pimod,newdata = l)
  # # pimod.vals <- predict(lambda[[j]]$pimod,newdata = l)$pred
  # sq.res <- (a-pimod.vals[1:n])^2
  #
  # pi2mod.vals <-predict(lambda[[j]]$pi2mod,newdata = l)
  # # pi2mod.vals <-predict(lambda[[j]]$pi2mod,newdata = l)$pred
  # # construct estimated pi/varpi and mu/m values
  # a.std         <- (a-pimod.vals)/sqrt(abs(pi2mod.vals))
  # pihat.vals    <- approx.fn(density(na.omit(a.std[1:n]))$x, density(na.omit(a.std)[1:n])$y,na.omit(a.std[1:n]))
  # dens_a         <- (pihat.vals/sqrt(abs(pi2mod.vals)))[1:n]
  
  
  
  a <- data[, A_name[j]]
  n <- length(a)
  
  # select all prior trt/covariate history
  namesA <- A_name[max(1, (j-k)):(j-1)]
  
  approx.fn <- function(x,y,z){ predict(smooth.spline(x,y),x=z)$y }
  if(j == 1){
    # make dat with only current and prior treatments #
    l <- data %>% select(one_of(c(namesA,A_name[j]))) %>%
      rename(A = A_name[j])
    
    # dens_a <- approx.fn(lambda[[j]]$density_model$x, lambda[[j]]$density_model$y, a)
    dens_a <- approx(lambda[[j]]$density_model$x, lambda[[j]]$density_model$y, xout = a)$y
    # Predict using the model
    # dens_a <- predict(lambda[[j]]$density_model, newdata = l)
  } else {
    l <- data %>% select(one_of(c(namesA)))
    
    if(ncol(l) == 1) {
      temp <- colnames(l)
      l <- cbind(0,l)
      colnames(l) <- c("Intercept",temp)
    }
    
    
    # rm(lambda)
    # lambda <- lambda_fits
    # lambda[[j]]$pimod$data <- lambda[[j]]$pimod$xlevels <- lambda[[j]]$pimod$R <- lambda[[j]]$pimod$residuals <-
    #   lambda[[j]]$pimod$fitted.values <- lambda[[j]]$pimod$weights <- lambda[[j]]$pimod$var <-
    #   lambda[[j]]$pimod$additive.predictors <- lambda[[j]]$pimod$prior.weights <- lambda[[j]]$pimod$y <- NULL
    
    
    pimod.vals <- predict(lambda[[j]]$pimod,newdata = l)
    # pimod.vals <- predict(lambda[[j]]$pimod,newdata = l)$pred
    sq.res <- (a-pimod.vals[1:n])^2
    
    pi2mod.vals <-predict(lambda[[j]]$pi2mod,newdata = l)
    # pi2mod.vals <-predict(lambda[[j]]$pi2mod,newdata = l)$pred
    # construct estimated pi/varpi and mu/m values
    a.std         <- (a-pimod.vals)/sqrt(abs(pi2mod.vals))
    pihat.vals    <- approx.fn(density(na.omit(a.std[1:n]))$x, density(na.omit(a.std)[1:n])$y,na.omit(a.std[1:n]))
    dens_a         <- (pihat.vals/sqrt(abs(pi2mod.vals)))[1:n]
    
  }
  
  return(as.numeric(dens_a))
  
}


#' SL_lambda
#'
#' Function to estimate densities of continuous treatment values accross time points as a function of previous
#' history using the method described by Kennedy et al.
#' @param j time point
#' @param data data used
#' @param stack SuperLeaner libraries to be used for fitting
#' @param valSets sets for cross-fitting
#' @param k Last time point to consider when getting prior history
#' @return estimate of density
SL_lambda <- function(j, data, A_name, k = Inf) {
  
  # print(j)
  # make dat with only current and prior treatments #
  namesA <- A_name[max(1, (j-k)):(j-1)]
  l <- data %>% select(one_of(c(namesA,A_name[j]))) %>%
    rename(A = A_name[j])
  
  # l <- as.data.frame(l)
  if(j ==1) {
    bw <- npudensbw(~A, data = l)
    density_model <- npudens(bw)
  } else {
    formula_str <- paste("A", "~", paste(namesA, collapse = " + "))
    model_formula <- as.formula(formula_str)
    model_formula <<- as.formula(formula_str)
    density_model <- npcdens(model_formula, data = l)
  }
  
  
  return(list("density_model" = density_model))
  
}


#' SL_ladata
#' Function to estimate densities of continuous treatment values accross time points as a function of previous

SL_ladata <- function(j, data, A_name, lambda , k = Inf, trunc) {
  
  # make dat with only current and prior treatments #
  namesA <- A_name[max(1, (j-k)):(j-1)]
  l <- data %>% select(one_of(c(namesA,A_name[j]))) %>%
    rename(A = A_name[j])
  
  
  # Predict using the model
  dens_a <- predict(lambda[[j]]$density_model, newdata = l)
  
  # if(trunc) dens_a <- truncate(dens_a, p = 0.01)
  
  # # Estimate variance
  # residuals <- (l$A - predictions)^2
  # l$residuals <- residuals
  # mse <- predict(lambda[[j]]$density_model, newdata = l)
  #
  # approx.fn     <- function(x,y,z){ predict(smooth.spline(x,y),x=z)$y }
  # # construct estimated pi/varpi and mu/m values
  # a.std         <- (residuals)/sqrt(mse)
  # pihat.vals    <- approx.fn(density(na.omit(a.std))$x, density(na.omit(a.std))$y,na.omit(a.std))
  # dens_a         <- (pihat.vals/sqrt(mse))
  
  return(dens_a)
  
}



#' estimators
#'
#' Function to estimate user specified transformation phi functional. Estimators currently avalaible
#' include sequentially doubly robust estimator, targeted maximum likelihood estimator and substituion estimator.
#' @name estimators
#' @param data dataset of interest
#' @param tau number of timepoints of interest
#' @param A vector character specifying the treatment to be used in analysis at each time point (must be of length equal to tau)
#' @param L list of character vectors specifying covariates to be included in analysis at each time point (must be of length equal to tau)
#' @param Y character specifying name of the outcome
#' @param est estimates to be included in analysis, options are "sdr" (sequentially doubly robust), "tml" (targeted maximum likelihood).
#' @param outcome character specifying outcome type. Currently available "binomial".
#' @param A_type character specifying if treatment is categorical or continuous. Only categorical available currently.
#' @param phi use define transformation functional of \theta(\bar{a},v) that will be the target of the causal estimator
#' @param m parametric model for \theta(\bar{a},v)
#' @param lambda can be either a list (of length tau) of models with covariates prior treatment history to get stabilizng weights. Or can be boolean, with TRUE, internal estimation of lambda using multinomial fits, or FALSE, to use no stabilizing weights. Default is TRUE.
#' @param nfolds number of folds to be used in cross-fitting
#' @param valSets sets for cross-fitting (default `split(sample(seq_len(nrow(data))), rep(seq_len(nfolds), length = nrow(data)))`)
#' @param stackr mlr3superlearner libraries used for probability mass estimation in the denominator of the weights.
#' @param stackm SuperLearner libraries used for outcome model prediction
#' @param k number of prior time points to consider when using prior history in fits. Default is Inf which is all prior time points.
#' @param trunc boolean specifying if density or probability estimates should be truncated. Useful when weights are large.
#' @param scale_weights boolean specifying if weights estimated should be scaled by their mean. Useful when weights are large.
#' @param rescale_outcome boolean specifying if outcome should be rescaled.
#' @param merge_results should results all be merged in a single dataframe (boolean).
#' @param trace should print statements be made to show progress (boolean).
#' @return list of the estimated values, for U1, and U2, of the estimators chosen. Along with variance and confidence intervals.
#' @export

estimators  <- function(data, tau, A, L, Y, est = c("sdr","tml"),
                        outcome = "binomial", A_type = "categorical", A.min = 0, A.max = 1,
                        phi = function(a) cbind(1, rowSums(a)),
                        lambda = TRUE, nfolds = 5,
                        valSets = NULL,
                        stackr = c("mean","randomforest","glmnet","knn","lightgbm","rpart","xgboost",
                                   "naivebayes","svm","lda","nnet","ranger","multinom"),
                        stackm = c('SL.mean','SL.glm','SL.caretEarth','SL.glmnet',
                                   'SL.glm.interaction','SL.gam','SL.randomForest',
                                   'SL.rpart','SL.bayesglm',"SL.svm"), k = Inf,
                        trunc = F, trim = F, scale_weights = F, rescale_outcome = F,
                        merge_results = T, trace = T,
                        int_bounds = NULL, increment = 1,
                        cores = 1, U2_method = "CG", known_weights = NULL,
                        approx_int = F, ken.dens = T, force_la_1 = F, seed = 123) {
  
  
  # replication #
  set.seed(seed)
  
  # arguments checks #
  if(trace) print("Performing checks")
  
  if(anyNA(data)) stop("Data contains missing values. Not supported by SuperLearner, please impute or perform complete case.")
  if(sum(c("sdr","tml") %in% est) == 0)
    stop('Please select at least one of the available estimators: c("sdr","tml").')
  if(!is.function(phi)) stop("Argument phi must be a function.")
  # if(!is.function(m)) stop("Argument m must be a function.")
  if(!(isTRUE(lambda) || isFALSE(lambda) || is.function(lambda)))
    if(!is.list(lambda)){ stop("Argument lambda must be a boolean, a function or a list.")
    } else if(length(lambda) != tau)
      warning("Length of lambda fits different than tau, default multinomial method will be used instead.")
  if(nfolds < 1 || round(nfolds) != nfolds) stop("nfolds must be a positive integer.")
  if(nfolds < 2) warning("In order to meet the Donsker class requirements cross-fitting must be performed. Note that your estimate may not have a causal interpretation.")
  if(length(stackr) < 1){
    warning("Argument stackr was empty, the default set of learners will be used instead.")
    stackr = c("mean","randomforest","glmnet","knn","lightgbm","rpart","xgboost",
               "naivebayes","svm","lda","nnet","ranger","multinom")
  }
  if(length(stackm) < 1){
    warning("Argument stackm was empty, the default set of learners will be used instead.")
    stackm = c('SL.mean','SL.glm','SL.caretEarth','SL.glmnet',
               'SL.glm.interaction','SL.gam','SL.randomForest',
               'SL.rpart','SL.bayesglm',"SL.svm")
  }
  
  if(cores < 1 || round(cores) != cores) stop("cores must be a positive integer.")
  
  if(!(U2_method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN","Brent")))
    warning("U2_method for `optim()` not recognized, defaulting to Nelder-Mead.")
  
  # add warning if A names different than those provided in fits lambda set lambda to TRUE instead #
  
  ### data checks and fix ###
  # Checks #
  if(length(A) != tau) stop("The number of treatments did not equal the number of time points.")
  if(length(L) != tau) stop("The number of covariates did not equal the number of time points.")
  if(is.na(match(Y, colnames(data)))) stop("The outcome variable specified was not found in the data.")
  # renaming #
  colnames(data)[match(Y, colnames(data))] <- "Y"
  # if(rescale_outcome)
  #   data$Y <- (data$Y - min(data$Y)) / (max(data$Y) - min(data$Y))
  data <- as.data.frame(data)
  
  
  
  # approx_int stuff #
  if(A_type == "continuous"){
    if(approx_int) {
      warning("Approximate integration will be used, results may not be accurate. Ensure that increment covers important regions of your exposure.")
      if(is.null(increment)) {
        warning("No increment specified for Riemann integration. Default will generate segments based on density of exposures.")
        
        
        int_info <- get_int_info(data, A, tau, quant_bin)
        A_span <- int_info$A_span
        binsizes <- int_info$binsizes
        bindim <- int_info$bindim
      }
    }
  }
  
  # creating folds #
  if(is.null(valSets))
    valSets <- split(sample(seq_len(nrow(data))), rep(seq_len(nfolds), length = nrow(data)))
  
  if(!(is.logical(trim) || is.numeric(trim))) stop("trim must be a boolean or a numeric value.")
  if(!is.logical(trim) && (trim < 0 || trim > 1)) stop("trim must be a positive value between 0 and 1.")
  
  # update values #
  # set function m as function of the outcome #
  if(outcome == "binomial"){
    m = function(a, b, phi) plogis(phi(a) %*% b)[,1]
  } else if(outcome == "gaussian") {
    m = function(a, b, phi) (phi(a) %*% b)[,1]
  }
  
  
  #### Generating weights ####
  if(trace) print("Step 1: Estimating weights.")
  if(is.null(known_weights)){
    # get lambda / stabilizing weights #
    if(is.list(lambda)){
      get_la <- T
      ladata <- sapply(1:tau, function(j)lambda_cust(j, A_name = A, data, k, fits = lambda, trunc, A_type))
    } else if(is.function(lambda)){
      ladata <- sapply(A, function(j)lambda(data[,j]))
      get_la <- T
      # lambda <- c() make lambda the bounds of integration here
    } else if(!lambda){
      ladata <- 1
      get_la <- F
    } else if(lambda){
      get_la <- T
      if(A_type == "categorical"){
        lambda <- hush(sapply(1:tau, function(j)lambda_est(j, A_name = A , data, k)))
        ladata <- hush(sapply(1:tau, function(j)lambda_cust(j, A_name = A, data, k, fits = lambda, trunc, A_type)))
      } else if(A_type == "continuous"){
        
        lambda <- hush(mclapply(1:tau, function(j)
          SL_densities_lambda(j, data, A_name = A, stack = stackr[[j]], k),
          mc.cores = cores
        ))
        
        ladata <- hush(do.call('cbind',mclapply(1:tau, function(j)
          SL_densities_ladata(j, data, A_name = A, lambda = lambda, k, trunc),
          mc.cores = cores)))
        
      }
    }
    
    
    
    # get g / density / probability mass #
    if(A_type == "categorical") {
      set.seed(seed)
      dens <- do.call('cbind',mclapply(1:tau, function(j)
        mult_densities(j, A_name = A, L_name = L, data,
                       stack = stackr[[j]], valSets, k, trunc), mc.cores = cores))
    } else if(A_type == "continuous"){
      if(!ken.dens){
        dens <- sapply(1:tau, function(j)OR_densities(j, data, A_name = A, L_name = L, stack = stackr[[j]], valSets, k))
      } else if(ken.dens){
        set.seed(seed)
        dens <- do.call('cbind',mclapply(1:tau, function(j)
          SL_densities_cf(j, data, A_name = A, L_name = L, stack = stackr[[j]], valSets, k, trunc),
          mc.cores = cores))
        
      }
    }
    
    # if(get_la) Zn <- t(apply(ladata / dens, 1, cumprod)) else Zn <- t(apply(ladata / dens, 1, cumprod))
    Zn <- t(apply(ladata / dens, 1, cumprod))
    # trim weights if needed #
    if(is.logical(trim)){
      if (isTRUE(trim)) { Zn <- apply(Zn, 2, function(x)pmin(x, quantile(x, 0.99)))
      } else Zn <- apply(Zn, 2, function(x)pmin(x, quantile(x, 1)))
    } else if(is.numeric(trim)){ Zn <- apply(Zn, 2, function(x)pmin(x, quantile(x, trim))) }
  } else{
    if(!(nrow(known_weights) == nrow(data) && ncol(known_weights) == tau))
      stop("Known weights must be of size n by tau.")
    dens <- known_weights
    Zn <- t(apply(1 / dens, 1, cumprod))
    if(is.list(lambda)) {
      get_la <- T
      ladata <- sapply(1:tau, function(j)lambda_cust(j, A_name = A, data, k, fits = lambda, trunc, A_type))
      Zn <- t(apply(ladata / dens, 1, cumprod))
    } else get_la <- F
  }
  if(tau == 1) Zn <- t(t(Zn))
  n <- nrow(data)
  
  # apply user defined transfomation #
  phidata <- as.matrix(data %>% select(all_of(A))) %>% phi()
  # Then to outcome #
  dataphi <- do.call('cbind',lapply(1:ncol(phidata), function(j)phidata[,j] * data[, 'Y']))
  
  
  eifsdr <- eifsub <- eiftml <- matrix(NA, ncol = ncol(dataphi), nrow = n)
  sdr <- sdr_bis <- sub <- tml <- rep(NA, ncol(dataphi))
  fitsub <- fitsdr <- fittml <- list()
  
  #### Generating outcome model ####
  if(trace) print("Step 2: Estimating outcome model:")
  
  if(force_la_1) get_la <- F
  ### Generate A_span if doesn't exist and continuous exposure with Riemann integration ###
  if(inherits(try(A_span, silent = T), 'try-error') && !is.null(increment)){
    
    
    quants <- sapply(1:tau,
                     function(j)
                       seq(A.min, A.max, by = increment))
    
    
    # generate A_span from quantiles #
    A_span <- as.matrix(expand.grid(
      lapply(1:tau,function(j) quants[,j]
      ))) #%>% as.data.frame()
    
    # create corresponding bins #
    binsizes <- 1
    
    # dimensions to adjust for #
    bindim <- 1
  }
  
  
  
  ###################################################
  
  for(kk in 1:ncol(dataphi)) {
    
    if(trace) print(paste0("Column phi: ",kk,"..."))
    if(rescale_outcome){
      mmax <- max(dataphi[, kk])
      mmin <- min(dataphi[, kk])
    } else{
      mmax <- 1
      mmin <- 0
    }
    
    Y <- (dataphi[, kk] - mmin) / (mmax - mmin)
    
    msubn <- msubd <- mtmln <- mtmld <-
      cbind(matrix(NA, nrow = n, ncol = tau), Y)
    
    
    for(t in tau:1){
      
      if(trace) print(paste0("Time point: ",t,"..."))
      A_temp <- data %>% pull(.data[[A[t]]])
      
      namesA <- paste('A_', max(1, t-k+1):(t-1), sep = '')
      namesA <- A[max(1, t-k+1):(t-1)]
      if(t == 1) namesA <- NULL
      namesL <- unlist(L[max(1, t-k+1):t])
      
      H <- data %>% select(one_of(c(namesA, namesL)))
      fam <- gaussian()
      XX  <- data.frame(A = A_temp, H)
      
      # check if variance exists - if not switch stack to mean #
      if(sd(msubd[, t + 1]) > .Machine$double.eps) stacksub <- stackm[[t]] else stacksub <- c('SL.mean')
      
      if(A_type == "categorical") size <- length(unique(A_temp))
      
      sub_fit <- mclapply(1:length(valSets), function(i){
        if(length(valSets) > 1){
          temp_dat_train <- XX[-valSets[[i]],]
          temp_dat_test <- XX[valSets[[i]],]
          Y_fit <- msubd[-valSets[[i]],t + 1]
        } else{
          temp_dat_train <- temp_dat_test <- XX[valSets[[i]],]
          Y_fit <- msubd[valSets[[i]], t + 1]
        }
        set.seed(seed)
        fitsub <- SuperLearner(Y = Y_fit, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stacksub)
        msubn_temp <- predict(fitsub)$pred[, 1]
        
        msubd_temp <- integral(fitsub, H[valSets[[i]],,drop = FALSE], A,
                               lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                               k, approx_int, get_la)
        return(list("msubn_temp"=msubn_temp, "msubd_temp"=msubd_temp))
      }, mc.cores = cores)
      
      for(i in 1:length(valSets)){
        msubn[valSets[[i]], t] <- sub_fit[[i]]$msubn_temp
        msubd[valSets[[i]], t] <- sub_fit[[i]]$msubd_temp
      }
      
      
      if(t == tau) {
        
        msdrn <- cbind(matrix(NA, nrow = n, ncol = tau - 1), msubn[, tau], Y)
        msdrd <- cbind(matrix(NA, nrow = n, ncol = tau - 1), msubd[, tau], Y)
        
        mtmlnini <- msubn[, tau]
        mtmldini <- msubd[, tau]
        
      }
      
      if(t < tau) {
        
        
        ### TMLE ###
        if("tml" %in% est){
          if(sd(mtmld[, t + 1]) > .Machine$double.eps) {
            stacktml <- stackm[[t]]
          } else {
            stacktml <- c('SL.mean')
          }
          
          
          tml_fit <- mclapply(1:length(valSets), function(i){
            if(length(valSets) > 1){
              temp_dat_train <- XX[-valSets[[i]],]
              temp_dat_test <- XX[valSets[[i]],]
              Y_fit <- mtmld[-valSets[[i]],t + 1]
            } else{
              temp_dat_train <- temp_dat_test <- XX[valSets[[i]],]
              Y_fit <- mtmld[valSets[[i]], t + 1]
            }
            set.seed(seed)
            fittml <- SuperLearner(Y = Y_fit, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stacktml)
            mtmlnini_temp <- predict(fittml)$pred[, 1]
            mtmldini_temp <- integral(fittml, H[valSets[[i]],,drop = FALSE], A,
                                      lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                      k, approx_int, get_la)
            return(list("mtmlnini_temp"=mtmlnini_temp, "mtmldini_temp"=mtmldini_temp))
            
          }, mc.cores = cores)
          
          
          for(i in 1:length(valSets)){
            mtmlnini[valSets[[i]]] <- tml_fit[[i]]$mtmlnini_temp
            mtmldini[valSets[[i]]] <- tml_fit[[i]]$mtmldini_temp
          }
          
        }
        
        
        ### Sequential regression estimator ###
        if("sdr" %in% est){
          
          Zs <- t(apply((ladata / dens)[, (t+1):tau, drop = FALSE], 1, cumprod))
          if(t == tau - 1) Zs <- t(Zs)
          
          outsdr <- rowSums(Zs * (msdrd[, (t + 2):(tau + 1), drop = FALSE] -
                                    msdrn[, (t + 1):tau, drop = FALSE])) + msdrd[, t + 1]
          if(sd(outsdr) > .Machine$double.eps) stacksdr <- stackm[[t]] else stacksdr <- c('SL.mean')
          
          
          sdr_fit <- mclapply(1:length(valSets), function(i){
            if(length(valSets) > 1){
              temp_dat_train <- XX[-valSets[[i]],]
              temp_dat_test <- XX[valSets[[i]],]
              # Y_fit <- msdrn[-valSets[[i]],t + 1]
              Y_fit <- outsdr[-valSets[[i]]]
            } else{
              temp_dat_train <- temp_dat_test <- XX[valSets[[i]],]
              Y_fit <- outsdr[valSets[[i]]]
            }
            set.seed(seed)
            fitsdr <- SuperLearner(Y = Y_fit, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stacksdr)
            msdrn_temp <- predict(fitsdr)$pred[, 1]
            msdrd_temp <- integral(fitsdr, H[valSets[[i]],,drop = FALSE], A,
                                   lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                   k, approx_int, get_la)
            return(list("msdrn_temp"=msdrn_temp, "msdrd_temp"=msdrd_temp))
            
          }, mc.cores = cores)
          
          
          for(i in 1:length(valSets)){
            msdrn[valSets[[i]], t] <- sdr_fit[[i]]$msdrn_temp
            msdrd[valSets[[i]], t] <- sdr_fit[[i]]$msdrd_temp
          }
          
        }
      }
      
      tilt <- lm(mtmld[, t + 1] ~ offset(mtmlnini),
                 weights = Zn[, t])
      
      mtmld[, t] <- mtmldini + coef(tilt)
      mtmln[, t] <- mtmlnini + coef(tilt)
      
    }
    
    
    if("sdr" %in% est){
      eifsdr[, kk] <- rowSums(Zn * (msdrd[, 2:(tau + 1), drop = FALSE] -
                                      msdrn[, 1:tau, drop = FALSE])) + msdrd[, 1]
      eifsdr[, kk] <- eifsdr[, kk] * (mmax - mmin) + mmin
      sdr[kk] <- mean(eifsdr[, kk])
    }
    if("tml" %in% est){
      eiftml[, kk] <- rowSums(Zn * (mtmld[, 2:(tau + 1), drop = FALSE] -
                                      mtmln[, 1:tau, drop = FALSE])) + mtmld[, 1]
      eiftml[, kk] <- eiftml[, kk] * (mmax - mmin) + mmin
      tml[kk] <- mean(eiftml[, kk])
    }
    
    sub[kk] <- mean(msubd[, 1]) * (mmax - mmin) + mmin
  }
  
  IPW <- glm2(data[, 'Y'] ~ 0 + phidata, weights = Zn[, tau], family = outcome)
  naive_fit <- glm(data[, 'Y'] ~ 0 + phidata, family = outcome)
  
  if(get_la){
    
    if(approx_int && A_type == "continuous"){
      if(trace) print("Step 2.5: Getting stabilizing weights...")
      # A_span <- expand.grid(lapply(1:tau, function(j) seq(A.min, A.max, by = increment)))
      
      temp <- as.data.frame(A_span)
      colnames(temp) <- A
      
      if(is.list(lambda)){
        la <- hush(do.call('cbind',mclapply(1:tau, function(j)
          SL_densities_ladata(j, temp, A_name = A, lambda = lambda, k, trunc),
          mc.cores = cores)))
        # la <- do.call('cbind',mclapply(1:tau,function(j)SL_ladata(j, temp, A_name = A, lambda, k, trunc),
        #                                mc.cores = cores))
        la <- apply(la, 1, prod)
      } else if(is.function(lambda)){
        la <- sapply(1:tau, function(j)lambda(A_span[,j]))
        la <- apply(la, 1, prod)
      }
    } else if(A_type == "categorical"){
      # A_span <- expand.grid(lapply(1:tau, function(j) seq(A.min, A.max, by = increment)))
      
      temp <- as.data.frame(A_span)
      colnames(temp) <- A
      la <- do.call('cbind',mclapply(1:tau,function(j)lambda_cust(j, A_name = A, temp, k, lambda, trunc, A_type),
                                     mc.cores = cores))
      la <- apply(la, 1, prod)
    } else if(!approx_int && A_type == "continuous"){
      la <- T
    }
  } else {la <- 1}
  
  if(force_la_1) la <- 1
  
  fsub1 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int){
    sum(abs(sub +  U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)))}
  fsdr1 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int) {
    tt <- sdr + U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)
    sqrt((t(tt) %*% solve(var(eifsdr)) %*% tt)[, 1])
  }
  ftml1 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int) {
    tt <- tml + U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)
    sqrt((t(tt) %*% solve(var(eiftml)) %*% tt)[, 1])
  }
  
  
  
  if(trace) print("Step 3: Solving...")
  ssub <- optim(par = coef(IPW), A = A, lambda = lambda, phi = phi, m = m,
                A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                fn = fsub1, method = U2_method,
                control = list(abstol = 0.01 * log(n)/n, maxit = 5000))
  if("sdr" %in% est)
    ssdr <- optim(par = coef(IPW), A = A, lambda = lambda, phi = phi, m = m,
                  A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                  la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                  fn = fsdr1, method = U2_method,
                  control = list(abstol = 0.01 * log(n)/n, maxit = 10^6))
  if("tml" %in% est)
    stml <- optim(par = coef(IPW), A = A, lambda = lambda, phi = phi, m = m,
                  A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                  la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                  fn = ftml1, method = U2_method,
                  control = list(abstol = 0.01 * log(n)/n, maxit = 10^6))
  
  
  if(trace) print("Step 4: Estimating variances")
  
  
  # variance and CI of SDR #
  sdrU1lo <- sdr - qnorm(0.975) * sqrt(diag(var(eifsdr))/n)
  sdrU1hi <- sdr + qnorm(0.975) * sqrt(diag(var(eifsdr))/n)
  
  Jsdr <- Jac(ssdr$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
              trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jisdr <- try(solve(Jsdr), silent = T)
  
  if(inherits(Jisdr, 'try-error')){
    Jisdr <- ginv(Jsdr)
    if(inherits(Jisdr, 'try-error')){
      sdrlo <- sdrhi <- rep(NA, ncol(dataphi))
      varsdr <- rep(NA, ncol(dataphi))
    } else {
      varsdr <- Jisdr %*% var(eifsdr) %*% t(Jisdr)
      sdrlo <- ssdr$par - qnorm(0.975) * sqrt(diag(varsdr)/n)
      sdrhi <- ssdr$par + qnorm(0.975) * sqrt(diag(varsdr)/n)
      varsdr = diag(varsdr) / n
    }
  } else {
    varsdr <- Jisdr %*% var(eifsdr) %*% t(Jisdr)
    sdrlo <- ssdr$par - qnorm(0.975) * sqrt(diag(varsdr)/n)
    sdrhi <- ssdr$par + qnorm(0.975) * sqrt(diag(varsdr)/n)
    varsdr = diag(varsdr) / n
  }
  
  
  # variance and CI of TMLE #
  tmlU1lo <- tml - qnorm(0.975) * sqrt(diag(var(eiftml))/n)
  tmlU1hi <- tml + qnorm(0.975) * sqrt(diag(var(eiftml))/n)
  
  Jtml <- Jac(stml$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
              trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jitml <- try(solve(Jtml), silent = T)
  
  if(inherits(Jitml, 'try-error')){
    Jitml <- ginv(Jtml)
    if(inherits(Jitml, 'try-error')){
      tmllo <- tmlhi <- rep(NA, ncol(dataphi))
      vartml <- rep(NA, ncol(dataphi))
    } else {
      vartml <- Jitml %*% var(eiftml) %*% t(Jitml)
      tmllo <- stml$par - qnorm(0.975) * sqrt(diag(vartml)/n)
      tmlhi <- stml$par + qnorm(0.975) * sqrt(diag(vartml)/n)
      vartml = diag(vartml) / n
    }
  } else {
    vartml <- Jitml %*% var(eiftml) %*% t(Jitml)
    tmllo <- stml$par - qnorm(0.975) * sqrt(diag(vartml)/n)
    tmlhi <- stml$par + qnorm(0.975) * sqrt(diag(vartml)/n)
    vartml = diag(vartml) / n
  }
  
  type <- paste0("Column phi:",1:ncol(dataphi))
  out <- list()
  if("sdr" %in% est){
    sdr_out <- tibble(
      type = type,
      sdr = ssdr$par,
      varsdr = varsdr,
      sdrlo = sdrlo,
      sdrhi = sdrhi,
      sdrU1 = sdr,
      sdrU1lo = sdrU1lo,
      sdrU1hi = sdrU1hi
    )
    out$SDR <- sdr_out
  }
  if("tml" %in% est){
    tml_out <- tibble(
      type = type,
      tml = stml$par,
      vartml = vartml,
      tmllo = tmllo,
      tmlhi = tmlhi,
      tmlU1 = tml,
      tmlU1lo = tmlU1lo,
      tmlU1hi = tmlU1hi
    )
    out$TML <- tml_out
  }
  
  ### futher values ###
  out$IPW <- tibble(
    type = type,
    ipw = coef(IPW),
    varipw = as.numeric(summary(IPW)$coefficients[,2])^2,
    naive = coef(naive_fit),
    varnaive = as.numeric(summary(naive_fit)$coefficients[,2])^2,
    sub   = ssub$par,
  )
  
  out$weights <- Zn
  
  if(get_la)
    out$ladata <- ladata
  if(rescale_outcome)
    out$Y <- data$Y
  
  # add arguments used #
  # out$args <-
  
  if(trace) print("Complete.")
  
  
  if(merge_results){
    out <- tibble(n = n,
                  type = type,
                  subU1 = sub,
                  sub   = ssub$par,
                  varipw = as.numeric(summary(IPW)$coefficients[,2])^2,
                  ipw   = coef(IPW),
                  sdrU1 = sdr,
                  sdrU1lo = sdrU1lo,
                  sdrU1hi = sdrU1hi,
                  sdr     = ssdr$par,
                  varsdr  = varsdr, #diag(varsdr) / n,
                  sdrlo = sdrlo,
                  sdrhi = sdrhi,
                  tmlU1 = tml,
                  tmlU1lo = tmlU1lo,
                  tmlU1hi = tmlU1hi,
                  tml     = stml$par,
                  vartml  = vartml,# diag(vartml) / n,
                  tmllo = tmllo,
                  tmlhi = tmlhi)
    return(list(out= out,
                eifsdr = eifsdr,
                eiftml = eiftml,
                # eifsub = eifsub,
                Zn = Zn,
                Zs = Zs,
                dens = dens,
                Jtml = Jtml,
                Jsdr = Jsdr,
                ssdr = ssdr,
                stml = stml,
                ssub = ssub,
                msdrd = msdrd,
                msdrn = msdrn,
                mtmld = mtmln,
                mtmld = mtmld
    ))
    
  }
  
  
  return(out)
}





# rm(list=ls())
#
# library(dplyr)
# library(parallel)
# library(SuperLearner)
# suppressPackageStartupMessages(library(mlr3superlearner))
# suppressPackageStartupMessages(library(mlr3extralearners))

################################################################################
#################### 1. Make functions to gen data #############################
################################################################################


# Covariates generation #

#' pl_mult
#'
#' Generate binary covariate.
#' @param a treatment at previous time point
#' @param l covariate at current time point
#' @param lp covariate at previous time point
#' @export
pl_mult  <- function(a, l, lp, size_A) {
  a <- a / size_A
  pp <- plogis((-0.5 * a + lp +  2 * l - a * l + lp * l))
  browser(expr = {any(pp < 0 | pp > 1)})
  return(pp)
}


# Treatment generation #

#' mu_mult
#'
#' Generate a treatment probabilities.
#' @param covs Covariate set
#' @param coefs Coefficients for each of the different possible treatments
#' @param size_A Number of possible treatments
#' @export
mu_mult <- function(covs,coefs,size_A){
  probs <- matrix(NA, nrow = nrow(covs), ncol(coefs))
  for(i in 1:ncol(coefs)){
    probs[,i] <- coefs[1,i] + coefs[2,i]*covs[,1]/size_A + coefs[3,i]*covs[,2] + coefs[4,i]*covs[,3]/size_A + coefs[5,i] * covs[,1]/size_A * covs[,2] + coefs[6,i] * covs[,4]
  }
  total_prob <- apply(probs, 1, function(x){x[1] + sum(exp(x[-1]))})
  for(i in 1:nrow(probs)){
    probs[i,1] <- probs[i,1]/total_prob[i]
    probs[i,2:ncol(probs)] <- exp(probs[i,2:ncol(probs)])/total_prob[i]
  }
  
  return(probs)
}



# Outcome generation #

#' m4_mult
#'
#' Time point 3 model.
#' @export
m4_mult <- function(a4, a3, a2, a1, l4, l3, l2, l1,size_A) {
  a4 <- a4 / size_A
  a3 <- a3 / size_A
  a2 <- a2 / size_A
  a1 <- a1 / size_A
  
  pp <- plogis( - 2.5 +
                  5*(l4 + l3 + l2 + l1/5) + 5 * l4 * l3 * l2 * l1/5 -
                  a4 * (5 - 2 * a3 * l4 * l3) -
                  a3 * (4 - 1.5 * a2 * l3 * l2) -
                  a2 * (3 - 1 * a1 * l2 * l1/5) -
                  a1 * (2 - 1 * l1/5)
  )
  # pp <- plogis(
  #   -1 + 2 * (l4 + l3 + l2 + l1) + (l4 * l3 + l2 * l1) -
  #     3 * (a4 + a3 + a2 + a1)
  # )
  return(pp)
}


#' m3_mult
#'
#' Time point 3 model.
#' @export
m3_mult <- function(a4, a3, a2, a1, l3, l2, l1,size_A) {
  m4_mult(a4, a3, a2, a1, 1, l3, l2, l1,size_A) * pl_mult(a3, l3, l2,size_A) +
    m4_mult(a4, a3, a2, a1, 0, l3, l2, l1,size_A) * (1 - pl_mult(a3, l3, l2,size_A))
}

#' m2_mult
#'
#' Time point 2 model.
#' @export
m2_mult <- function(a4, a3, a2, a1, l2, l1,size_A) {
  m3_mult(a4, a3,  a2, a1, 1, l2, l1,size_A) * pl_mult(a2, l2, l1,size_A) +
    m3_mult(a4, a3, a2, a1, 0, l2, l1,size_A) * (1 - pl_mult(a2, l2, l1,size_A))
}

#' m1_mult
#'
#' Time point 1 model.
#' @export
m1_mult <- function(a4, a3, a2, a1, l1,size_A) {
  m2_mult(a4, a3, a2, a1, 1, l1,size_A) * pl_mult(a1, l1, 0,size_A) +
    m2_mult(a4, a3, a2, a1, 0, l1,size_A) * (1 - pl_mult(a1, l1, 0,size_A))
}

#' m0_mult
#'
#' Time point 0 model.
#' @export
m0_mult <- function(a4, a3, a2, a1, size_A, L1_size = 5) {
  m1_mult(a4, a3, a2, a1, 0/size_A, size_A) * 0.2 + m1_mult(a4, a3, a2, a1, 1/size_A, size_A) * 0.2 +  m1_mult(a4, a3, a2, a1, 2/size_A, size_A) * 0.2 +
    m1_mult(a4, a3, a2, a1, 3/size_A, size_A) * 0.2 +  m1_mult(a4, a3, a2, a1, 4/size_A, size_A) * 0.2
}

#' g_mult
#'
#' Get true probability of observed treatment
#' @param data data observed
#' @param coefs true coefficients of data generation
#' @param size size of the treatment
#' @param j time point of interest
#' @export

g_mult <- function(data, coefs, size, j){
  
  if(j == 1){
    data$`A_-1` <- 0
    data$L_0 <- 0
    data$A_0 <- 0
    data$L_1 <- data$L_1/5
  } else if(j == 2){
    data$A_0 <- 0
    data$L_0 <- 0
    data$L_1 <- data$L_1/5
  }
  
  covs_keep <- c( paste0("A_", (j-1)), paste0("L_", j), paste0("A_", (j-2)), paste0("L_", (j-1)) )
  data_temp <- data %>% select(one_of(covs_keep))
  
  probs <- matrix(NA, nrow = nrow(data), ncol(coefs))
  for(i in 1:ncol(coefs)){
    probs[,i] <- coefs[1,i] + coefs[2,i]*data_temp[,1]/size +
      coefs[3,i]*data_temp[,2] + coefs[4,i]*data_temp[,3]/size +
      coefs[5,i] * data_temp[,1]/size * data_temp[,2] + coefs[6,i] * data_temp[,4]
  }
  
  total_prob <- apply(probs, 1, function(x){x[1] + sum(exp(x[-1]))})
  for(i in 1:nrow(probs)){
    probs[i,1] <- probs[i,1]/total_prob[i]
    probs[i,2:ncol(probs)] <- exp(probs[i,2:ncol(probs)])/total_prob[i]
  }
  
  true_prob <- c()
  for(r in 1:nrow(data_temp)){
    temp_preds <- probs[r,]
    true_prob[r] <- as.numeric(temp_preds[(data[r,paste0("A_",j)] + 1)])
  }
  # match observed treatment with it's prob
  return(true_prob)
}


#' datagen_mult
#'
#' Generate a dataset with treatment from a mulitnomial distribution for a specified number of time points with a binary outcome.
#' @param n Sample size of the dataset to be generated
#' @param size Size of the categorical treatment to be generated
#' @param tau Number of time points to be generated
#' @param mu Functional specifying the form of the probabilities to be generated for each treatment
#' @param pl Functional specifying the form of the probabilities to be generated for each covariate
#' @param m4 Functional specifying the form of the outcome
#' @export
datagen_mult <- function(n = 250, size_A = 5, tau = 4, mu, pl, m4, L1_size = 5){
  
  
  # T = 1 #
  L_1 <- sample(0:(L1_size-1), n, prob = rep(1/L1_size, L1_size), replace = T)/(L1_size-1)
  coefs_1 = matrix(c(c(1,0,0,0,0, 0),
                     c(0,0.1,0.25,-0.1,0.1, 1),
                     c(0,-0.1,0.5,0.1,-0.1, 0.75),
                     c(0,0.2,0.75,0.2,-0.2, 0/5),
                     c(0,-0.2,1,-0.2,0.2, 0.25)), ncol = size_A)
  p_1 <- mu(cbind(0.5,L_1,0.5,0.5), coefs_1,size_A = size_A - 1)
  A_1 <- as.double(apply(p_1, MARGIN = 1, function(x)sample(x = 0:(size_A - 1), size = 1, prob = x)))
  # get true p1 #
  true_p1 <- c()
  for(r in 1:n){
    temp_preds <- p_1[r,]
    true_p1[r] <- as.numeric(temp_preds[(A_1[r] + 1)])
  }
  
  # T = 2 #
  L_2 <- rbinom(n, 1, pl(A_1, L_1, 0, size_A = size_A - 1))
  # coefs_2 = matrix(c(c(1,0,0,0,0, 0),
  #                    c(0,0.1,0.25,-0.1,0.1, 1),
  #                    c(0,-0.1,0.5,0.1,-0.1, 0.75),
  #                    c(0,0.2,0.75,0.2,-0.2, 0.5),
  #                    c(0,-0.2,1,-0.2,0.2, 0.25)), ncol = size_A)
  coefs_2 = matrix(c(c(1,0,0,0,0, 0),
                     c(0,0.1,0.25,-0.1,0.1, 1),
                     c(0,-0.1,0.5,0.1,-0.1, 0.75),
                     c(0,0.2,0.75,0.2,-0.2, 0.5),
                     c(0,-0.2,1,-0.2,0.2, 0.25)), ncol = size_A)*1.5
  p_2 <- mu(cbind(A_1,L_2,0, L_1), coefs_2,size_A = size_A - 1)
  A_2 <- as.double(apply(p_2, MARGIN = 1, function(x)sample(x = 0:(size_A - 1), size = 1, prob = x)))
  true_p2 <- c()
  for(r in 1:n){
    temp_preds <- p_2[r,]
    true_p2[r] <- as.numeric(temp_preds[(A_2[r] + 1)])
  }
  
  # T = 3 #
  L_3 <- rbinom(n, 1, pl(A_2, L_2, L_1, size_A = size_A - 1))
  
  coefs_3 = matrix(c(c(1,0,0,0,0, 0),
                     c(0,0.1,0.25,-0.1,0.1, 1),
                     c(0,-0.1,0.5,0.1,-0.1, 0.75),
                     c(0,0.2,0.75,0.2,-0.2, 0.5),
                     c(0,-0.2,1,-0.2,0.2, 0.25)), ncol = size_A)*1.5
  p_3 <- mu(cbind(A_2,L_3,A_1, L_2), coefs_3, size_A = size_A - 1)
  A_3 <- as.double(apply(p_3, MARGIN = 1, function(x)sample(x = 0:(size_A - 1), size = 1, prob = x)))
  true_p3 <- c()
  for(r in 1:n){
    temp_preds <- p_3[r,]
    true_p3[r] <- as.numeric(temp_preds[(A_3[r] + 1)])
  }
  
  # T = 4 #
  L_4 <- rbinom(n, 1, pl(A_3, L_3, L_2, size_A = size_A - 1))
  
  coefs_4 = matrix(c(c(1,0,0,0,0, 0),
                     c(0,0.1,0.25,-0.1,0.1, 1),
                     c(0,-0.1,0.5,0.1,-0.1, 0.75),
                     c(0,0.2,0.75,0.2,-0.2, 0.5),
                     c(0,-0.2,1,-0.2,0.2, 0.25)), ncol = size_A)*1.5
  p_4 <- mu(cbind(A_3,L_4,A_2, L_3), coefs_4,size_A = size_A - 1)
  A_4 <- as.double(apply(p_4, MARGIN = 1, function(x)sample(x = 0:(size_A - 1), size = 1, prob = x)))
  true_p4 <- c()
  for(r in 1:n){
    temp_preds <- p_4[r,]
    true_p4[r] <- as.numeric(temp_preds[(A_4[r] + 1)])
  }
  
  
  data <- as.data.frame(cbind(L_1,A_1,L_2,A_2,L_3,A_3,L_4,A_4))
  
  Y <- extraDistr::rbern(n = n, prob = with(data,m4(A_4, A_3, A_2, A_1, L_4, L_3, L_2, L_1,size_A = size_A - 1)))
  data$Y = Y
  
  # true values #
  true_ps <- cbind(true_p1,true_p2,true_p3,true_p4)
  return(list(data=data,
              true_ps = true_ps,
              coefs_4 = coefs_4,
              coefs_3 = coefs_3,
              coefs_2 = coefs_2,
              coefs_1 = coefs_1
  ))
  
}

############## GENERATE TRUE VALUES ##################


#' true_mult
#'
#' Generate a large dataset with treatment from a multinomial distribution for a specified number of time points with a binary outcome to get the true underlying value of the treatment effect for a specified phi function.
#' @param n Sample size of the dataset to be generated
#' @param size Size of the categorical treatment to be generated
#' @param tau Number of time points to be generated
#' @param s seed
#' @param phi functional
#' @param m functional
#' @param save_models should models generated be saved
#' @param mu Functional specifying the form of the probabilities to be generated for each treatment
#' @param pl Functional specifying the form of the probabilities to be generated for each covariate
#' @param m4 Functional specifying the form of the outcome
#' @export
true_mult <- function(n, size, tau = 4, s, phi, m, save_models = T, mu = mu_mult, pl = pl_mult, m4 = m4_mult, cores = 4) {
  
  
  s = 123
  n = 10^6
  size_A = 5
  tau = 4
  phi = function(a) cbind(1, rowSums(a))
  m = function(a, b, phi) plogis(phi(a) %*% b)[,1]
  save_models = T
  mu = mu_mult
  pl = pl_mult
  m4 = m4_mult
  increment <- 1
  A.min <- 0
  A.max <- size_A - 1
  tau = 4
  k = Inf
  trunc = F
  cores = 4
  A <- c("A_1","A_2","A_3","A_4")
  L = list(c("L_1"),c("L_2"),c("L_3"),c("L_4"))
  data <- datagen_mult(n = n, size_A = size_A, tau = tau, mu = mu, pl = pl, m4 = m4)
  dat <- data$data
  
  library(dplyr)
  library(parallel)
  library(SuperLearner)
  suppressPackageStartupMessages(library(mlr3superlearner))
  suppressPackageStartupMessages(library(mlr3extralearners))
  
  # Discrete #
  set.seed(s)
  data <- datagen_mult(n = n, size_A = size_A, tau = tau, mu = mu, pl = pl, m4 = m4)
  dat <- data$data
  dd <- data$true_ps
  
  # for lambda fit in the population here and keep it --> feed it to the actual estimators
  k = Inf
  lambda_fits <- list()
  weightsnum <- matrix(NA, nrow = n, ncol = tau)
  for(j in 1:tau){
    
    nameA <- paste0('A_', j, sep = '')
    a <- dat[, nameA]
    
    if(j == 1){
      preds <- summary(as.factor(a))/n
      lambda_fits[[j]] <- preds
      for(r in 1:n){
        weightsnum[r,j] <- as.numeric(preds[(a[r] + 1)])
      }
    }
    
    else{
      if(j > 1) namesA <- paste('A_', max(1, (j-k)):(j-1), sep = '') else namesA <- NULL
      namesL <- paste('L_', max(1, j-k+1):j, sep = '')
      # stack in a single dataset
      l <- dat %>% select(one_of(c(namesA)))
      
      temp_dat <- as.data.frame(cbind(a, l))
      fit <- nnet::multinom(a ~ ., data = temp_dat)
      fit$residuals <- NULL
      fit$weights <- NULL
      lambda_fits[[j]] <- fit
      preds <- predict(fit, newdata = temp_dat, type = "probs")
      for(r in 1:n){
        temp_preds <- preds[r,]
        weightsnum[r,j] <- as.numeric(temp_preds[(a[r] + 1)])
      }
    }
  }
  
  if(save_models)
    save(lambda_fits, file = here::here("data/lambda_fits_simple_mult.Rdata"))
  
  
  # True IPW #
  Zn <- t(apply(1/dd, 1, cumprod))
  summary(Zn)
  weights <- Zn[,4]
  
  model <- glm2::glm2(dat$Y ~ 0 + phi(with(dat, data.frame(A_1, A_2, A_3, A_4))),
                      family = binomial, weights = weights)
  
  # True IPW la #
  Zn_la <- t(apply(weightsnum/dd, 1, cumprod))
  summary(Zn_la)
  weights_la <- Zn_la[,4]
  
  model_la <- glm2::glm2(dat$Y ~ 0 + phi(with(dat, data.frame(A_1, A_2, A_3, A_4))),
                         family = binomial, weights = weights_la)
  
  
  # Bias IPW #
  valSets <- list(1:n)
  dens_w <- do.call('cbind',mclapply(1:tau, function(j)
    mult_densities(j, A_name = A, L_name = L, dat,
                   stack = "mean", valSets, k, trunc), mc.cores = cores))
  bias_Zn <- t(apply(1/dens_w, 1, cumprod))
  bias_weights <- bias_Zn[,4]
  model1 <- glm2::glm2(dat$Y ~ 0 + phi(with(dat, data.frame(A_1, A_2, A_3, A_4))),
                       family = binomial, weights = bias_Zn[,4])
  
  # Bias IPW la #
  bias_Zn_la <- t(apply(weightsnum/dens_w, 1, cumprod))
  bias_weights_la <- bias_Zn_la[,4]
  model1_la <- glm2::glm2(dat$Y ~ 0 + phi(with(dat, data.frame(A_1, A_2, A_3, A_4))),
                          family = binomial, weights = bias_Zn_la[,4])
  
  
  
  ######## TRUE VALUE ########
  fU1 <- function(a) {
    
    mm0 <- m0_mult(a[4,], a[3,], a[2,], a[1,], size_A = size_A-1)
    pp <- phi(t(a))
    
    temp <- as.data.frame(t(a))
    colnames(temp) <- paste0("A_",1:tau)
    la <- 1
    
    return(mm0 * pp * la)
    
  }
  
  A_span <- t(expand.grid(a1 = seq(A.min,A.max,by = increment),
                          a2 = seq(A.min,A.max,by = increment),
                          a3 = seq(A.min,A.max,by = increment),
                          a4 = seq(A.min,A.max,by = increment)))
  U1 <- apply(fU1(A_span), 2, sum)
  
  
  # i <- 1
  f <- function(b) {
    # print(i)
    # i <<- i + 1
    fU2 <- function(a) {
      
      temp <- as.data.frame(t(a))
      colnames(temp) <- paste0("A_",1:tau)
      la <- 1
      
      pp <- phi(t(a))
      mm <- m(t(a), b, phi)
      
      out <- - mm * pp * la
      return(out)
    }
    
    
    U2 <- apply(fU2(A_span), 2, sum)
    
    
    return(sum(abs(U1 + U2)))
    
  }
  
  tt_CG <- optim(par = coef(model), fn = f, method = "CG",
                 control = list(abstol = 0.01 * 1/n, maxit = 10^6, reltol = .Machine$double.eps))
  tt_NM <- optim(par = coef(model), fn = f, method = "Nelder-Mead",
                 control = list(abstol = 0.01 * 1/n, maxit = 10^6, reltol = .Machine$double.eps))
  
  
  
  ######## TRUE VALUE LA ########
  fU1_la <- function(a) {
    
    mm0 <- m0_mult(a[4,], a[3,], a[2,], a[1,], size_A = size_A-1)
    pp <- phi(t(a))
    
    temp <- as.data.frame(t(a))
    colnames(temp) <- paste0("A_",1:tau)
    la <- sapply(1:tau, function(j)lambda_cust(j, A_name = colnames(temp), data = temp, k = Inf, fits = lambda_fits, A_type = "categorical"))
    la <- apply(la, 1, prod)
    
    return(mm0 * pp * la)
    
  }
  
  A_span <- t(expand.grid(a1 = seq(A.min,A.max,by = increment),
                          a2 = seq(A.min,A.max,by = increment),
                          a3 = seq(A.min,A.max,by = increment),
                          a4 = seq(A.min,A.max,by = increment)))
  U1_la <- apply(fU1_la(A_span), 2, sum)
  
  
  # i <- 1
  f <- function(b) {
    # print(i)
    # i <<- i + 1
    fU2_la <- function(a) {
      
      temp <- as.data.frame(t(a))
      colnames(temp) <- paste0("A_",1:tau)
      la <- sapply(1:tau, function(j)lambda_cust(j, A_name = colnames(temp), data = temp, k = Inf, fits = lambda_fits, A_type = "categorical"))
      la <- apply(la, 1, prod)
      
      pp <- phi(t(a))
      mm <- m(t(a), b, phi)
      
      out <- - mm * pp * la
      return(out)
    }
    
    
    U2_la <- apply(fU2_la(A_span), 2, sum)
    
    
    return(sum(abs(U1_la + U2_la)))
    
  }
  
  tt_CG_la <- optim(par = coef(model), fn = f, method = "CG",
                    control = list(abstol = 0.01 * 1/n, maxit = 10^6, reltol = .Machine$double.eps))
  tt_NM_la <- optim(par = coef(model), fn = f, method = "Nelder-Mead",
                    control = list(abstol = 0.01 * 1/n, maxit = 10^6, reltol = .Machine$double.eps))
  
  
  ############################################################################################################
  
  
  
  
  
  truevalue <- list(
    U1 = U1, par_CG = tt_CG$par, par_NM = tt_NM$par,
    ipw = coef(model),
    U1_la = U1_la, par_CG_la = tt_CG_la$par, par_NM_la = tt_NM_la$par,
    ipw = coef(model_la),
    trueweights = summary(Zn),
    trueweights_la = summary(Zn_la),
    ipw_bias = coef(model1),
    ipw_bias_la = coef(model1_la)
  )
  
  
  
  
  # checks #
  # truevalue
  # sum(dat$Y)/n
  #
  # summary(data$true_ps)
  #
  # # fix data #
  # dat$L_1 <- dat$L_1/5
  # dat$A_1_tmp <- dat$A_1/size
  # dat$A_2_tmp <- dat$A_2/size
  # dat$A_3_tmp <- dat$A_3/size
  # dat$A_4_tmp <- dat$A_4/size
  # # make models #
  # m1 <- multinom(A_1 ~ L_1, data = dat, trace = FALSE)
  # coef(m1)
  # m2 <- multinom(A_2 ~ A_1_tmp + L_2 + A_1_tmp * L_2 + L_1 , data = dat, trace = FALSE)
  # coef(m2)
  # m3 <- multinom(A_3_tmp ~ A_2_tmp + L_3 + A_1_tmp + A_2_tmp * L_3 + L_2 , data = dat, trace = FALSE)
  # coef(m3)
  # m4 <- multinom(A_4_tmp ~ A_3_tmp + L_4 + A_2_tmp + A_3_tmp * L_4 + L_3 , data = dat, trace = FALSE)
  # coef(m4)
  # # mY <- glm(Y ~ A_4_tmp * A_3_tmp * A_2_tmp * A_1_tmp * L_4 * L_3 * L_2 * L_1, data = dat, family = binomial, trace = FALSE)
  # # coef(mY)
  # # looks fine #
  
  
  
  if(save_models)
    save(truevalue, file = here::here('./data/true_discrete_cat5_simple_mult.Rdata'))
  return(truevalue)
}


# true_mult(n = 10^6, size = 5, tau = 4, s = 210793,
#           phi = function(a) cbind(1, rowSums(a)),
#           m = function(a, b, phi) plogis(phi(a) %*% b)[,1],
#           save_models = T, mu = mu_mult, pl = pl_mult, m4 = m4_mult)





#' simul_mult
#' Simulation function to get results from scenarios 1-5
#' @param data
#' @param stackr
#' @param stackm
#' @param seed
#' @export

simul_mult <- function(data, stackr, stackm, lambda, seed, la,
                       A_span, binsizes, bindim,
                       get_la = T, trueW = NULL){
  
  
  # make the different scenarios #
  stackr_c <- stackr
  stackm_c <- stackm
  stackm_w <- stackr_w <- c("SL.mean")
  # stacks scenarios #
  # 1 all correct #
  stackr_1 <- list(stackr_c,stackr_c,stackr_c,stackr_c)
  stackm_1 <- list(stackm_c,stackm_c,stackm_c,stackm_c)
  # 2 correct outcome #
  stackr_2 <- list(stackr_w,stackr_w,stackr_w,stackr_w)
  stackm_2 <- list(stackm_c,stackm_c,stackm_c,stackm_c)
  # 3 correct weights #
  stackr_3 <- list(stackr_c,stackr_c,stackr_c,stackr_c)
  stackm_3 <- list(stackm_w,stackm_w,stackm_w,stackm_w)
  # 4
  stackr_4 <- list(stackr_w, stackr_w, stackr_c, stackr_c)
  stackm_4 <- list(stackm_c, stackm_c, stackm_w, stackm_w)
  # 5
  stackr_5 <- list(stackr_c, stackr_c, stackr_w, stackr_w)
  stackm_5 <- list(stackm_w, stackm_w, stackm_c, stackm_c)
  
  
  
  ### estimator function ###
  # set.seed(seed)
  
  # arguments needed #
  outcome = "binomial"
  A_type = "categorical"
  data <- as.data.frame(data)
  approx_int = T
  valSets = NULL
  k = Inf
  tau <- 4
  A = c("A_1","A_2","A_3","A_4")
  L = list(c("L_1"),c("L_2"),c("L_3"),c("L_4"))
  Y = "Y"
  cores = 1
  nfolds = 1
  trunc = F
  phi = function(a) cbind(1, rowSums(a))
  rescale_outcome = F
  trace = T
  trim = F
  
  # creating folds #
  if(is.null(valSets))
    valSets <- list(1:nrow(data))
  
  # update values #
  # set function m as function of the outcome #
  if(outcome == "binomial"){
    m = function(a, b, phi) plogis(phi(a) %*% b)[,1]
  } else if(outcome == "gaussian") {
    m = function(a, b, phi) (phi(a) %*% b)[,1]
  }
  
  #### Generating weights ####
  # stab weights #
  if(get_la){
    ladata <- hush(do.call('cbind',mclapply(1:tau, function(j)
      lambda_cust(j, A_name = A, data, k, fits = lambda, trunc, A_type),
      mc.cores = cores))) } else {
        la <- rep(1,n)
        ladata <- 1
      }
  
  if(is.null(trueW)){
    dens_c <- do.call('cbind',mclapply(1:tau, function(j)
      mult_densities(j, data, A_name = A, L_name = L, stack = stackr_c, valSets, k, trunc),
      mc.cores = cores)) } else {
        dens_c <- trueW
      }
  
  # incorrect #
  # set.seed(seed)
  dens_w <- do.call('cbind',mclapply(1:tau, function(j)
    mult_densities(j, A_name = A, L_name = L, data,
                   stack = "mean", valSets, k, trunc),
    mc.cores = cores))
  
  # weights for each scenarios #
  Zn_1 <- t(apply(ladata / dens_c, 1, cumprod))
  Zn_2 <- t(apply(ladata / dens_w, 1, cumprod))
  Zn_3 <- t(apply(ladata / dens_c, 1, cumprod))
  Zn_4 <- t(apply(ladata / cbind(dens_w[,1:2], dens_c[,3:4]), 1, cumprod))
  Zn_5 <- t(apply(ladata / cbind(dens_c[,1:2], dens_w[,3:4]), 1, cumprod))
  
  
  # Outcome models #
  n <- nrow(data)
  
  # apply user defined transformation #
  phidata <- data %>% select(all_of(A)) %>% phi()
  
  # Then to outcome #
  dataphi <- sapply(1:ncol(phidata), function(j)phidata[,j] * data[, 'Y'])
  
  eifsdr_1 <- eifsub_1 <- eiftml_1 <-
    eifsdr_2 <- eifsub_2 <- eiftml_2 <-
    eifsdr_3 <- eifsub_3 <- eiftml_3 <-
    eifsdr_4 <- eifsub_4 <- eiftml_4 <-
    eifsdr_5 <- eifsub_5 <- eiftml_5 <- matrix(NA, ncol = ncol(dataphi), nrow = n)
  sdr_1 <- sub_1 <- tml_1 <-
    sdr_2 <- sub_2 <- tml_2 <-
    sdr_3 <- sub_3 <- tml_3 <-
    sdr_4 <- sub_4 <- tml_4 <-
    sdr_5 <- sub_5 <- tml_5 <- rep(NA, ncol(dataphi))
  fitsub <- fitsdr <- fittml <- list()
  
  
  for(kk in 1:ncol(dataphi)) {
    
    if(trace) print(paste0("Column phi: ",kk,"..."))
    if(rescale_outcome){
      mmax <- max(dataphi[, kk])
      mmin <- min(dataphi[, kk])
    } else{
      mmax <- 1
      mmin <- 0
    }
    
    Y <- (dataphi[, kk] - mmin) / (mmax - mmin)
    
    
    # correct/wrong #
    msubn_c <- msubd_c <- mtmln_c <- mtmld_c <-
      msubn_w <- msubd_w <- mtmln_w <- mtmld_w <-
      cbind(matrix(NA, nrow = n, ncol = tau), Y)
    
    # scenarios #
    msubn_1 <- msubd_1 <- mtmln_1 <- mtmld_1 <-
      msubn_2 <- msubd_2 <- mtmln_2 <- mtmld_2 <-
      msubn_3 <- msubd_3 <- mtmln_3 <- mtmld_3 <-
      msubn_4 <- msubd_4 <- mtmln_4 <- mtmld_4 <-
      msubn_5 <- msubd_5 <- mtmln_5 <- mtmld_5 <-
      cbind(matrix(NA, nrow = n, ncol = tau), Y)
    
    
    
    for(t in tau:1){
      
      if(trace) print(paste0("Time point: ",t,"..."))
      A_temp <- data %>% pull(.data[[A[t]]])
      if(A_type == "categorical") size <- length(unique(A_temp))
      
      namesA <- paste('A_', max(1, t-k+1):(t-1), sep = '')
      namesA <- A[max(1, t-k+1):(t-1)]
      if(t == 1) namesA <- NULL
      namesL <- unlist(L[max(1, t-k+1):t])
      
      H <- data %>% select(one_of(c(namesA, namesL)))
      # if(outcome == "binomial" && t == tau) fam <- binomial() else fam <- gaussian()
      fam <- gaussian()
      XX  <- data.frame(A = A_temp, H)
      
      ########################################################################
      
      if(t == tau) {
        if(A_type == "categorical") size <- length(unique(A_temp))
        
        
        ########################################################################
        sub_fit <- lapply(1:length(valSets), function(i){
          
          temp_dat_train <- temp_dat_test <- XX[valSets[[i]],]
          # scenario 1 #
          Y_fit_1 <- msubd_1[valSets[[i]], t + 1]
          if(sd(Y_fit_1) > .Machine$double.eps) stackm_1[[t]] <- stackm_1[[t]] else stackm_1[[t]] <- c('SL.mean')
          # set.seed(seed)
          fitsub_1 <- hush(SuperLearner(Y = Y_fit_1, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_1[[t]]))
          if(sum(fitsub_1$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fitsub_1$errorsInCVLibrary)
            fitsub_1 <- SuperLearner(Y = Y_fit_1, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_1[[t]][valid_coefs])
          }
          msubn_temp_1 <- predict(fitsub_1)$pred[, 1]
          msubd_temp_1 <- integral(fitsub_1, H[valSets[[i]],,drop = FALSE], A,
                                   lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                   k, approx_int, get_la)
          
          # scenario 2 #
          Y_fit_2 <- msubd_2[valSets[[i]], t + 1]
          if(sd(Y_fit_2) > .Machine$double.eps) stackm_2[[t]] <- stackm_2[[t]] else stackm_2[[t]] <- c('SL.mean')
          # set.seed(seed)
          fitsub_2 <- hush(SuperLearner(Y = Y_fit_2, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_2[[t]]))
          if(sum(fitsub_2$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fitsub_2$errorsInCVLibrary)
            fitsub_2 <- SuperLearner(Y = Y_fit_2, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_2[[t]][valid_coefs])
          }
          msubn_temp_2 <- predict(fitsub_2)$pred[, 1]
          msubd_temp_2 <- integral(fitsub_2, H[valSets[[i]],,drop = FALSE], A,
                                   lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                   k, approx_int, get_la)
          
          # scenario 3 #
          Y_fit_3 <- msubd_3[valSets[[i]], t + 1]
          if(sd(Y_fit_3) > .Machine$double.eps) stackm_3[[t]] <- stackm_3[[t]] else stackm_3[[t]] <- c('SL.mean')
          # set.seed(seed)
          fitsub_3 <- hush(SuperLearner(Y = Y_fit_3, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_3[[t]]))
          if(sum(fitsub_3$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fitsub_3$errorsInCVLibrary)
            fitsub_3 <- SuperLearner(Y = Y_fit_3, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_3[[t]][valid_coefs])
          }
          msubn_temp_3 <- predict(fitsub_3)$pred[, 1]
          msubd_temp_3 <- integral(fitsub_3, H[valSets[[i]],,drop = FALSE], A,
                                   lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                   k, approx_int, get_la)
          
          # scenario 4 #
          Y_fit_4 <- msubd_4[valSets[[i]], t + 1]
          if(sd(Y_fit_4) > .Machine$double.eps) stackm_4[[t]] <- stackm_4[[t]] else stackm_4[[t]] <- c('SL.mean')
          # set.seed(seed)
          fitsub_4 <- hush(SuperLearner(Y = Y_fit_4, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_4[[t]]))
          if(sum(fitsub_4$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fitsub_4$errorsInCVLibrary)
            fitsub_4 <- SuperLearner(Y = Y_fit_4, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_4[[t]][valid_coefs])
          }
          msubn_temp_4 <- predict(fitsub_4)$pred[, 1]
          msubd_temp_4 <- integral(fitsub_4, H[valSets[[i]],,drop = FALSE], A,
                                   lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                   k, approx_int, get_la)
          
          
          # scenario 5 #
          Y_fit_5 <- msubd_5[valSets[[i]], t + 1]
          if(sd(Y_fit_5) > .Machine$double.eps) stackm_5[[t]] <- stackm_5[[t]] else stackm_5[[t]] <- c('SL.mean')
          # set.seed(seed)
          fitsub_5 <- hush(SuperLearner(Y = Y_fit_5, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_5[[t]]))
          if(sum(fitsub_5$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fitsub_5$errorsInCVLibrary)
            fitsub_5 <- SuperLearner(Y = Y_fit_5, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_5[[t]][valid_coefs])
          }
          msubn_temp_5 <- predict(fitsub_5)$pred[, 1]
          msubd_temp_5 <- integral(fitsub_5, H[valSets[[i]],,drop = FALSE], A,
                                   lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                   k, approx_int, get_la)
          
          # return all values @#
          return(list("msubn_temp_1"=msubn_temp_1, "msubd_temp_1"=msubd_temp_1,
                      "msubn_temp_2"=msubn_temp_2, "msubd_temp_2"=msubd_temp_2,
                      "msubn_temp_3"=msubn_temp_3, "msubd_temp_3"=msubd_temp_3,
                      "msubn_temp_4"=msubn_temp_4, "msubd_temp_4"=msubd_temp_4,
                      "msubn_temp_5"=msubn_temp_5, "msubd_temp_5"=msubd_temp_5
          ))
        })
        
        for(i in 1:length(valSets)){
          # scanrio 1 #
          msubn_1[valSets[[i]], t] <- sub_fit[[i]]$msubn_temp_1
          msubd_1[valSets[[i]], t] <- sub_fit[[i]]$msubd_temp_1
          # scenario 2 #
          msubn_2[valSets[[i]], t] <- sub_fit[[i]]$msubn_temp_2
          msubd_2[valSets[[i]], t] <- sub_fit[[i]]$msubd_temp_2
          # scenario 3 #
          msubn_3[valSets[[i]], t] <- sub_fit[[i]]$msubn_temp_3
          msubd_3[valSets[[i]], t] <- sub_fit[[i]]$msubd_temp_3
          # scenario 4 #
          msubn_4[valSets[[i]], t] <- sub_fit[[i]]$msubn_temp_4
          msubd_4[valSets[[i]], t] <- sub_fit[[i]]$msubd_temp_4
          # scenario 5 #
          msubn_5[valSets[[i]], t] <- sub_fit[[i]]$msubn_temp_5
          msubd_5[valSets[[i]], t] <- sub_fit[[i]]$msubd_temp_5
        }
        
        
        
        ########################################################################
        # init tml #
        mtmlnini_1 <- msubn_1[, tau]
        mtmldini_1 <- msubd_1[, tau]
        mtmlnini_2 <- msubn_2[, tau]
        mtmldini_2 <- msubd_2[, tau]
        mtmlnini_3 <- msubn_3[, tau]
        mtmldini_3 <- msubd_3[, tau]
        mtmlnini_4 <- msubn_4[, tau]
        mtmldini_4 <- msubd_4[, tau]
        mtmlnini_5 <- msubn_5[, tau]
        mtmldini_5 <- msubd_5[, tau]
        
        # init sdr #
        msdrn_1 <- cbind(matrix(NA, nrow = n, ncol = tau - 1), msubn_1[, tau], Y)
        msdrd_1 <- cbind(matrix(NA, nrow = n, ncol = tau - 1), msubd_1[, tau], Y)
        msdrn_2 <- cbind(matrix(NA, nrow = n, ncol = tau - 1), msubn_2[, tau], Y)
        msdrd_2 <- cbind(matrix(NA, nrow = n, ncol = tau - 1), msubd_2[, tau], Y)
        msdrn_3 <- cbind(matrix(NA, nrow = n, ncol = tau - 1), msubn_3[, tau], Y)
        msdrd_3 <- cbind(matrix(NA, nrow = n, ncol = tau - 1), msubd_3[, tau], Y)
        msdrn_4 <- cbind(matrix(NA, nrow = n, ncol = tau - 1), msubn_4[, tau], Y)
        msdrd_4 <- cbind(matrix(NA, nrow = n, ncol = tau - 1), msubd_4[, tau], Y)
        msdrn_5 <- cbind(matrix(NA, nrow = n, ncol = tau - 1), msubn_5[, tau], Y)
        msdrd_5 <- cbind(matrix(NA, nrow = n, ncol = tau - 1), msubd_5[, tau], Y)
        ########################################################################
      }
      
      if(t < tau) {
        
        
        ########################################################################
        ### TMLE ###
        ########################################################################
        tml_fit <- lapply(1:length(valSets), function(i){
          
          
          temp_dat_train <- temp_dat_test <- XX[valSets[[i]],]
          # scenario 1 #
          Y_fit_1 <- mtmld_1[valSets[[i]], t + 1]
          if(sd(Y_fit_1) > .Machine$double.eps) stackm_1[[t]] <- stackm_1[[t]] else stackm_1[[t]] <- c('SL.mean')
          # set.seed(seed)
          fittml_1 <- hush(SuperLearner(Y = Y_fit_1, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_1[[t]]))
          if(sum(fittml_1$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fittml_1$errorsInCVLibrary)
            fittml_1 <- SuperLearner(Y = Y_fit_1, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_1[[t]][valid_coefs])
          }
          mtmlnini_temp_1 <- predict(fittml_1)$pred[, 1]
          mtmldini_temp_1 <- integral(fittml_1, H[valSets[[i]],,drop = FALSE], A,
                                      lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                      k, approx_int, get_la)
          
          # scenario 2 #
          Y_fit_2 <- mtmld_2[valSets[[i]], t + 1]
          if(sd(Y_fit_2) > .Machine$double.eps) stackm_2[[t]] <- stackm_2[[t]] else stackm_2[[t]] <- c('SL.mean')
          # set.seed(seed)
          fittml_2 <- hush(SuperLearner(Y = Y_fit_2, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_2[[t]]))
          if(sum(fittml_2$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fittml_2$errorsInCVLibrary)
            fittml_2 <- SuperLearner(Y = Y_fit_2, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_2[[t]][valid_coefs])
          }
          mtmlnini_temp_2 <- predict(fittml_2)$pred[, 1]
          mtmldini_temp_2 <- integral(fittml_2, H[valSets[[i]],,drop = FALSE], A,
                                      lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                      k, approx_int, get_la)
          
          # scenario 3 #
          Y_fit_3 <- mtmld_3[valSets[[i]], t + 1]
          if(sd(Y_fit_3) > .Machine$double.eps) stackm_3[[t]] <- stackm_3[[t]] else stackm_3[[t]] <- c('SL.mean')
          # set.seed(seed)
          fittml_3 <- hush(SuperLearner(Y = Y_fit_3, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_3[[t]]))
          if(sum(fittml_3$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fittml_3$errorsInCVLibrary)
            fittml_3 <- SuperLearner(Y = Y_fit_3, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_3[[t]][valid_coefs])
          }
          mtmlnini_temp_3 <- predict(fittml_3)$pred[, 1]
          mtmldini_temp_3 <- integral(fittml_3, H[valSets[[i]],,drop = FALSE], A,
                                      lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                      k, approx_int, get_la)
          
          # scenario 4 #
          Y_fit_4 <- mtmld_4[valSets[[i]], t + 1]
          if(sd(Y_fit_4) > .Machine$double.eps) stackm_4[[t]] <- stackm_4[[t]] else stackm_4[[t]] <- c('SL.mean')
          # set.seed(seed)
          fittml_4 <- hush(SuperLearner(Y = Y_fit_4, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_4[[t]]))
          if(sum(fittml_4$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fittml_4$errorsInCVLibrary)
            fittml_4 <- SuperLearner(Y = Y_fit_4, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_4[[t]][valid_coefs])
          }
          mtmlnini_temp_4 <- predict(fittml_4)$pred[, 1]
          mtmldini_temp_4 <- integral(fittml_4, H[valSets[[i]],,drop = FALSE], A,
                                      lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                      k, approx_int, get_la)
          
          
          # scenario 5 #
          Y_fit_5 <- mtmld_5[valSets[[i]], t + 1]
          if(sd(Y_fit_5) > .Machine$double.eps) stackm_5[[t]] <- stackm_5[[t]] else stackm_5[[t]] <- c('SL.mean')
          # set.seed(seed)
          fittml_5 <- hush(SuperLearner(Y = Y_fit_5, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_5[[t]]))
          if(sum(fittml_5$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fittml_5$errorsInCVLibrary)
            fittml_5 <- SuperLearner(Y = Y_fit_5, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_5[[t]][valid_coefs])
          }
          mtmlnini_temp_5 <- predict(fittml_5)$pred[, 1]
          mtmldini_temp_5 <- integral(fittml_5, H[valSets[[i]],,drop = FALSE], A,
                                      lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                      k, approx_int, get_la)
          
          # return all values @#
          return(list("mtmlnini_temp_1"=mtmlnini_temp_1, "mtmldini_temp_1"=mtmldini_temp_1,
                      "mtmlnini_temp_2"=mtmlnini_temp_2, "mtmldini_temp_2"=mtmldini_temp_2,
                      "mtmlnini_temp_3"=mtmlnini_temp_3, "mtmldini_temp_3"=mtmldini_temp_3,
                      "mtmlnini_temp_4"=mtmlnini_temp_4, "mtmldini_temp_4"=mtmldini_temp_4,
                      "mtmlnini_temp_5"=mtmlnini_temp_5, "mtmldini_temp_5"=mtmldini_temp_5
          ))
          
        })
        
        
        for(i in 1:length(valSets)){
          # scanrio 1 #
          mtmlnini_1[valSets[[i]]] <- tml_fit[[i]]$mtmlnini_temp_1
          mtmldini_1[valSets[[i]]] <- tml_fit[[i]]$mtmldini_temp_1
          # scenario 2 #
          mtmlnini_2[valSets[[i]]] <- tml_fit[[i]]$mtmlnini_temp_2
          mtmldini_2[valSets[[i]]] <- tml_fit[[i]]$mtmldini_temp_2
          # scenario 3 #
          mtmlnini_3[valSets[[i]]] <- tml_fit[[i]]$mtmlnini_temp_3
          mtmldini_3[valSets[[i]]] <- tml_fit[[i]]$mtmldini_temp_3
          # scenario 4 #
          mtmlnini_4[valSets[[i]]] <- tml_fit[[i]]$mtmlnini_temp_4
          mtmldini_4[valSets[[i]]] <- tml_fit[[i]]$mtmldini_temp_4
          # scenario 5 #
          mtmlnini_5[valSets[[i]]] <- tml_fit[[i]]$mtmlnini_temp_5
          mtmldini_5[valSets[[i]]] <- tml_fit[[i]]$mtmldini_temp_5
        }
        
        
        ########################################################################
        
        
        ### Sequential regression estimator ###
        
        # weights for all scenarios #
        Zs_1 <- t(apply((ladata / dens_c)[, (t+1):tau, drop = FALSE], 1, cumprod))
        Zs_2 <- t(apply((ladata / dens_w)[, (t+1):tau, drop = FALSE], 1, cumprod))
        Zs_3 <- t(apply((ladata / dens_c)[, (t+1):tau, drop = FALSE], 1, cumprod))
        Zs_4 <- t(apply((ladata / cbind(dens_w[,1:2], dens_c[,3:4]))[, (t+1):tau, drop = FALSE], 1, cumprod))
        Zs_5 <- t(apply((ladata / cbind(dens_c[,1:2], dens_w[,3:4]))[, (t+1):tau, drop = FALSE], 1, cumprod))
        if(t == tau - 1) {
          Zs_1 <- t(Zs_1)
          Zs_2 <- t(Zs_2)
          Zs_3 <- t(Zs_3)
          Zs_4 <- t(Zs_4)
          Zs_5 <- t(Zs_5)
        }
        
        # outcome for all scenarios #
        outsdr_1 <- rowSums(Zs_1 * (msdrd_1[, (t + 2):(tau + 1), drop = FALSE] -
                                      msdrn_1[, (t + 1):tau, drop = FALSE])) + msdrd_1[, t + 1]
        outsdr_2 <- rowSums(Zs_2 * (msdrd_2[, (t + 2):(tau + 1), drop = FALSE] -
                                      msdrn_2[, (t + 1):tau, drop = FALSE])) + msdrd_2[, t + 1]
        outsdr_3 <- rowSums(Zs_3 * (msdrd_3[, (t + 2):(tau + 1), drop = FALSE] -
                                      msdrn_3[, (t + 1):tau, drop = FALSE])) + msdrd_3[, t + 1]
        outsdr_4 <- rowSums(Zs_4 * (msdrd_4[, (t + 2):(tau + 1), drop = FALSE] -
                                      msdrn_4[, (t + 1):tau, drop = FALSE])) + msdrd_4[, t + 1]
        outsdr_5 <- rowSums(Zs_5 * (msdrd_5[, (t + 2):(tau + 1), drop = FALSE] -
                                      msdrn_5[, (t + 1):tau, drop = FALSE])) + msdrd_5[, t + 1]
        
        
        
        sdr_fit <- lapply(1:length(valSets), function(i){
          
          temp_dat_train <- temp_dat_test <- XX[valSets[[i]],]
          # scenario 1 #
          Y_fit_1 <- outsdr_1[valSets[[i]]]
          if(sd(Y_fit_1) > .Machine$double.eps) stackm_1[[t]] <- stackm_1[[t]] else stackm_1[[t]] <- c('SL.mean')
          # set.seed(seed)
          fitsdr_1 <- hush(SuperLearner(Y = Y_fit_1, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_1[[t]]))
          if(sum(fitsdr_1$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fitsdr_1$errorsInCVLibrary)
            fitsdr_1 <- SuperLearner(Y = Y_fit_1, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_1[[t]][valid_coefs])
          }
          msdrn_temp_1 <- predict(fitsdr_1)$pred[, 1]
          msdrd_temp_1 <- integral(fitsdr_1, H[valSets[[i]],,drop = FALSE], A,
                                   lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                   k, approx_int, get_la)
          
          # scenario 2 #
          Y_fit_2 <- outsdr_2[valSets[[i]]]
          if(sd(Y_fit_2) > .Machine$double.eps) stackm_2[[t]] <- stackm_2[[t]] else stackm_2[[t]] <- c('SL.mean')
          # set.seed(seed)
          fitsdr_2 <- hush(SuperLearner(Y = Y_fit_2, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_2[[t]]))
          if(sum(fitsdr_2$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fitsdr_2$errorsInCVLibrary)
            fitsdr_2 <- SuperLearner(Y = Y_fit_2, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_2[[t]][valid_coefs])
          }
          msdrn_temp_2 <- predict(fitsdr_2)$pred[, 1]
          msdrd_temp_2 <- integral(fitsdr_2, H[valSets[[i]],,drop = FALSE], A,
                                   lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                   k, approx_int, get_la)
          
          # scenario 3 #
          Y_fit_3 <- outsdr_3[valSets[[i]]]
          if(sd(Y_fit_3) > .Machine$double.eps) stackm_3[[t]] <- stackm_3[[t]] else stackm_3[[t]] <- c('SL.mean')
          # set.seed(seed)
          fitsdr_3 <- hush(SuperLearner(Y = Y_fit_3, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_3[[t]]))
          if(sum(fitsdr_3$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fitsdr_3$errorsInCVLibrary)
            fitsdr_3 <- SuperLearner(Y = Y_fit_3, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_3[[t]][valid_coefs])
          }
          msdrn_temp_3 <- predict(fitsdr_3)$pred[, 1]
          msdrd_temp_3 <- integral(fitsdr_3, H[valSets[[i]],,drop = FALSE], A,
                                   lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                   k, approx_int, get_la)
          
          # scenario 4 #
          Y_fit_4 <- outsdr_4[valSets[[i]]]
          if(sd(Y_fit_4) > .Machine$double.eps) stackm_4[[t]] <- stackm_4[[t]] else stackm_4[[t]] <- c('SL.mean')
          # set.seed(seed)
          fitsdr_4 <- hush(SuperLearner(Y = Y_fit_4, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_4[[t]]))
          if(sum(fitsdr_4$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fitsdr_4$errorsInCVLibrary)
            fitsdr_4 <- SuperLearner(Y = Y_fit_4, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_4[[t]][valid_coefs])
          }
          msdrn_temp_4 <- predict(fitsdr_4)$pred[, 1]
          msdrd_temp_4 <- integral(fitsdr_4, H[valSets[[i]],,drop = FALSE], A,
                                   lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                   k, approx_int, get_la)
          
          
          # scenario 5 #
          Y_fit_5 <- outsdr_5[valSets[[i]]]
          if(sd(Y_fit_5) > .Machine$double.eps) stackm_5[[t]] <- stackm_5[[t]] else stackm_5[[t]] <- c('SL.mean')
          # set.seed(seed)
          fitsdr_5 <- hush(SuperLearner(Y = Y_fit_5, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_5[[t]]))
          if(sum(fitsdr_5$errorsInCVLibrary) > 0){
            valid_coefs <- which(!fitsdr_5$errorsInCVLibrary)
            fitsdr_5 <- SuperLearner(Y = Y_fit_5, X = temp_dat_train, newX = temp_dat_test,family = fam, SL.library = stackm_5[[t]][valid_coefs])
          }
          msdrn_temp_5 <- predict(fitsdr_5)$pred[, 1]
          msdrd_temp_5 <- integral(fitsdr_5, H[valSets[[i]],,drop = FALSE], A,
                                   lambda, A_type, A.min, A.max, trunc, size, A_span = t(unique(A_span[,t])), binsizes = unique(binsizes[,t]),
                                   k, approx_int, get_la)
          
          # return all values @#
          return(list("msdrn_temp_1"=msdrn_temp_1, "msdrd_temp_1"=msdrd_temp_1,
                      "msdrn_temp_2"=msdrn_temp_2, "msdrd_temp_2"=msdrd_temp_2,
                      "msdrn_temp_3"=msdrn_temp_3, "msdrd_temp_3"=msdrd_temp_3,
                      "msdrn_temp_4"=msdrn_temp_4, "msdrd_temp_4"=msdrd_temp_4,
                      "msdrn_temp_5"=msdrn_temp_5, "msdrd_temp_5"=msdrd_temp_5
          ))
          
        })
        
        
        for(i in 1:length(valSets)){
          # scanrio 1 #
          msdrn_1[valSets[[i]], t] <- sdr_fit[[i]]$msdrn_temp_1
          msdrd_1[valSets[[i]], t] <- sdr_fit[[i]]$msdrd_temp_1
          # scenario 2 #
          msdrn_2[valSets[[i]], t] <- sdr_fit[[i]]$msdrn_temp_2
          msdrd_2[valSets[[i]], t] <- sdr_fit[[i]]$msdrd_temp_2
          # scenario 3 #
          msdrn_3[valSets[[i]], t] <- sdr_fit[[i]]$msdrn_temp_3
          msdrd_3[valSets[[i]], t] <- sdr_fit[[i]]$msdrd_temp_3
          # scenario 4 #
          msdrn_4[valSets[[i]], t] <- sdr_fit[[i]]$msdrn_temp_4
          msdrd_4[valSets[[i]], t] <- sdr_fit[[i]]$msdrd_temp_4
          # scenario 5 #
          msdrn_5[valSets[[i]], t] <- sdr_fit[[i]]$msdrn_temp_5
          msdrd_5[valSets[[i]], t] <- sdr_fit[[i]]$msdrd_temp_5
        }
        
      }
      
      ########################################################################
      
      # scenario 1 #
      tilt_1 <- lm(mtmld_1[, t + 1] ~ offset(mtmlnini_1),
                   weights = Zn_1[, t])
      mtmld_1[, t] <- mtmldini_1 + coef(tilt_1)
      mtmln_1[, t] <- mtmlnini_1 + coef(tilt_1)
      # scenario 2 #
      tilt_2 <- lm(mtmld_2[, t + 1] ~ offset(mtmlnini_2),
                   weights = Zn_2[, t])
      mtmld_2[, t] <- mtmldini_2 + coef(tilt_2)
      mtmln_2[, t] <- mtmlnini_2 + coef(tilt_2)
      # scenario 3 #
      tilt_3 <- lm(mtmld_3[, t + 1] ~ offset(mtmlnini_3),
                   weights = Zn_3[, t])
      mtmld_3[, t] <- mtmldini_3 + coef(tilt_3)
      mtmln_3[, t] <- mtmlnini_3 + coef(tilt_3)
      # scenario 4 #
      tilt_4 <- lm(mtmld_4[, t + 1] ~ offset(mtmlnini_4),
                   weights = Zn_4[, t])
      mtmld_4[, t] <- mtmldini_4 + coef(tilt_4)
      mtmln_4[, t] <- mtmlnini_4 + coef(tilt_4)
      # scenario 5 #
      tilt_5 <- lm(mtmld_5[, t + 1] ~ offset(mtmlnini_5),
                   weights = Zn_5[, t])
      mtmld_5[, t] <- mtmldini_5 + coef(tilt_5)
      mtmln_5[, t] <- mtmlnini_5 + coef(tilt_5)
      
      
    }
    
    ### SDR scenarios ###
    # scenario 1 #
    eifsdr_1[, kk] <- rowSums(Zn_1 * (msdrd_1[, 2:(tau + 1), drop = FALSE] -
                                        msdrn_1[, 1:tau, drop = FALSE])) + msdrd_1[, 1]
    eifsdr_1[, kk] <- eifsdr_1[, kk] * (mmax - mmin) + mmin
    sdr_1[kk] <- mean(eifsdr_1[, kk])
    # scenario 2 #
    eifsdr_2[, kk] <- rowSums(Zn_2 * (msdrd_2[, 2:(tau + 1), drop = FALSE] -
                                        msdrn_2[, 1:tau, drop = FALSE])) + msdrd_2[, 1]
    eifsdr_2[, kk] <- eifsdr_2[, kk] * (mmax - mmin) + mmin
    sdr_2[kk] <- mean(eifsdr_2[, kk])
    # scenario 3 #
    eifsdr_3[, kk] <- rowSums(Zn_3 * (msdrd_3[, 2:(tau + 1), drop = FALSE] -
                                        msdrn_3[, 1:tau, drop = FALSE])) + msdrd_3[, 1]
    eifsdr_3[, kk] <- eifsdr_3[, kk] * (mmax - mmin) + mmin
    sdr_3[kk] <- mean(eifsdr_3[, kk])
    # scenario 4 #
    eifsdr_4[, kk] <- rowSums(Zn_4 * (msdrd_4[, 2:(tau + 1), drop = FALSE] -
                                        msdrn_4[, 1:tau, drop = FALSE])) + msdrd_4[, 1]
    eifsdr_4[, kk] <- eifsdr_4[, kk] * (mmax - mmin) + mmin
    sdr_4[kk] <- mean(eifsdr_4[, kk])
    # scenario 5 #
    eifsdr_5[, kk] <- rowSums(Zn_5 * (msdrd_5[, 2:(tau + 1), drop = FALSE] -
                                        msdrn_5[, 1:tau, drop = FALSE])) + msdrd_5[, 1]
    eifsdr_5[, kk] <- eifsdr_5[, kk] * (mmax - mmin) + mmin
    sdr_5[kk] <- mean(eifsdr_5[, kk])
    
    ### TML scenarios ###
    # scenario 1 #
    eiftml_1[, kk] <- rowSums(Zn_1 * (mtmld_1[, 2:(tau + 1), drop = FALSE] -
                                        mtmln_1[, 1:tau, drop = FALSE])) + mtmld_1[, 1]
    eiftml_1[, kk] <- eiftml_1[, kk] * (mmax - mmin) + mmin
    tml_1[kk] <- mean(eiftml_1[, kk])
    # scenario 2 #
    eiftml_2[, kk] <- rowSums(Zn_2 * (mtmld_2[, 2:(tau + 1), drop = FALSE] -
                                        mtmln_2[, 1:tau, drop = FALSE])) + mtmld_2[, 1]
    eiftml_2[, kk] <- eiftml_2[, kk] * (mmax - mmin) + mmin
    tml_2[kk] <- mean(eiftml_2[, kk])
    # scenario 3 #
    eiftml_3[, kk] <- rowSums(Zn_3 * (mtmld_3[, 2:(tau + 1), drop = FALSE] -
                                        mtmln_3[, 1:tau, drop = FALSE])) + mtmld_3[, 1]
    eiftml_3[, kk] <- eiftml_3[, kk] * (mmax - mmin) + mmin
    tml_3[kk] <- mean(eiftml_3[, kk])
    # scenario 4 #
    eiftml_4[, kk] <- rowSums(Zn_4 * (mtmld_4[, 2:(tau + 1), drop = FALSE] -
                                        mtmln_4[, 1:tau, drop = FALSE])) + mtmld_4[, 1]
    eiftml_4[, kk] <- eiftml_4[, kk] * (mmax - mmin) + mmin
    tml_4[kk] <- mean(eiftml_4[, kk])
    # scenario 5 #
    eiftml_5[, kk] <- rowSums(Zn_5 * (mtmld_5[, 2:(tau + 1), drop = FALSE] -
                                        mtmln_5[, 1:tau, drop = FALSE])) + mtmld_5[, 1]
    eiftml_5[, kk] <- eiftml_5[, kk] * (mmax - mmin) + mmin
    tml_5[kk] <- mean(eiftml_5[, kk])
    
    ### SUB scenarios ###
    # scenario 1 #
    sub_1[kk] <- mean(msubd_1[, 1]) * (mmax - mmin) + mmin
    # scenario 2 #
    sub_2[kk] <- mean(msubd_2[, 1]) * (mmax - mmin) + mmin
    # scenario 3 #
    sub_3[kk] <- mean(msubd_3[, 1]) * (mmax - mmin) + mmin
    # scenario 4 #
    sub_4[kk] <- mean(msubd_4[, 1]) * (mmax - mmin) + mmin
    # scenario 5 #
    sub_5[kk] <- mean(msubd_5[, 1]) * (mmax - mmin) + mmin
  }
  
  # IPW fits #
  # scenario 1 #
  IPW_1 <- glm2(data[, 'Y'] ~ 0 + phidata, weights = Zn_1[, tau], family = outcome)
  # scenario 2 #
  IPW_2 <- glm2(data[, 'Y'] ~ 0 + phidata, weights = Zn_2[, tau], family = outcome)
  # scenario 3 #
  IPW_3 <- glm2(data[, 'Y'] ~ 0 + phidata, weights = Zn_3[, tau], family = outcome)
  # scenario 4 #
  IPW_4 <- glm2(data[, 'Y'] ~ 0 + phidata, weights = Zn_4[, tau], family = outcome)
  # scenario 5 #
  IPW_5 <- glm2(data[, 'Y'] ~ 0 + phidata, weights = Zn_5[, tau], family = outcome)
  # naive fit #
  naive_fit <- glm(data[, 'Y'] ~ 0 + phidata, family = outcome)
  
  
  
  ##############################################################################
  # make functions for solving U2 #
  # scenario 1 #
  fsub1_1 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int){
    sum(abs(sub_1 +  U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)))}
  fsdr1_1 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int) {
    tt <- sdr_1 + U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)
    sqrt((t(tt) %*% solve(var(eifsdr_1)) %*% tt)[, 1])
  }
  ftml1_1 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int) {
    tt <- tml_1 + U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)
    sqrt((t(tt) %*% solve(var(eiftml_1)) %*% tt)[, 1])
  }
  # scenario 2 #
  fsub1_2 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int){
    sum(abs(sub_2 +  U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)))}
  fsdr1_2 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int) {
    tt <- sdr_2 + U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)
    sqrt((t(tt) %*% solve(var(eifsdr_2)) %*% tt)[, 1])
  }
  ftml1_2 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int) {
    tt <- tml_2 + U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)
    sqrt((t(tt) %*% solve(var(eiftml_2)) %*% tt)[, 1])
  }
  # scenario 3 #
  fsub1_3 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int){
    sum(abs(sub_3 +  U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)))}
  fsdr1_3 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int) {
    tt <- sdr_3 + U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)
    sqrt((t(tt) %*% solve(var(eifsdr_3)) %*% tt)[, 1])
  }
  ftml1_3 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int) {
    tt <- tml_3 + U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)
    sqrt((t(tt) %*% solve(var(eiftml_3)) %*% tt)[, 1])
  }
  # scenario 4 #
  fsub1_4 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int){
    sum(abs(sub_4 +  U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)))}
  fsdr1_4 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int) {
    tt <- sdr_4 + U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)
    sqrt((t(tt) %*% solve(var(eifsdr_4)) %*% tt)[, 1])
  }
  ftml1_4 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int) {
    tt <- tml_4 + U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)
    sqrt((t(tt) %*% solve(var(eiftml_4)) %*% tt)[, 1])
  }
  # scenario 5 #
  fsub1_5 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int){
    sum(abs(sub_5 +  U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)))}
  fsdr1_5 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int) {
    tt <- sdr_5 + U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)
    sqrt((t(tt) %*% solve(var(eifsdr_5)) %*% tt)[, 1])
  }
  ftml1_5 <- function(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int) {
    tt <- tml_5 + U2fun(b, A, lambda, phi, m, A_type, A.min, A.max, trunc, tau, outcome, la, A_span, bindim, approx_int)
    sqrt((t(tt) %*% solve(var(eiftml_5)) %*% tt)[, 1])
  }
  
  
  
  
  
  ##############################################################################
  if(trace) print("Step 3: Solving NM...")
  # scenario 1 #
  ssdr_1 <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                  A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                  la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                  fn = fsdr1_1, method = "Nelder-Mead",
                  control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  stml_1 <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                  A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                  la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                  fn = ftml1_1, method = "Nelder-Mead",
                  control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  
  # scenario 2 #
  ssdr_2 <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                  A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                  la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                  fn = fsdr1_2, method = "Nelder-Mead",
                  control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  stml_2 <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                  A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                  la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                  fn = ftml1_2, method = "Nelder-Mead",
                  control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  
  # scenario 3 #
  ssdr_3 <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                  A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                  la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                  fn = fsdr1_3, method = "Nelder-Mead",
                  control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  stml_3 <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                  A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                  la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                  fn = ftml1_3, method = "Nelder-Mead",
                  control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  
  # scenario 4 #
  ssdr_4 <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                  A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                  la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                  fn = fsdr1_4, method = "Nelder-Mead",
                  control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  stml_4 <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                  A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                  la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                  fn = ftml1_4, method = "Nelder-Mead",
                  control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  # scenario 5 #
  ssdr_5 <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                  A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                  la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                  fn = fsdr1_5, method = "Nelder-Mead",
                  control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  stml_5 <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                  A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                  la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                  fn = ftml1_5, method = "Nelder-Mead",
                  control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  
  ##############################################################################
  
  
  if(trace) print("Step 4: Estimating variances NM")
  # Scenario 1 #
  # SDR #
  sdrU1lo_1 <- sdr_1 - qnorm(0.975) * sqrt(diag(var(eifsdr_1))/n)
  sdrU1hi_1 <- sdr_1 + qnorm(0.975) * sqrt(diag(var(eifsdr_1))/n)
  
  Jsdr_1 <- Jac(ssdr_1$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jisdr_1 <- try(solve(Jsdr_1), silent = T)
  if(inherits(Jisdr_1, 'try-error')){
    Jisdr_1 <- ginv(Jsdr_1)
    if(inherits(Jisdr_1, 'try-error')){
      sdrlo_1 <- sdrhi_1 <- rep(NA, ncol(dataphi))
      varsdr_1 <- rep(NA, ncol(dataphi))
    }
  } else {
    varsdr_1 <- Jisdr_1 %*% var(eifsdr_1) %*% t(Jisdr_1)
    sdrlo_1 <- ssdr_1$par - qnorm(0.975) * sqrt(diag(varsdr_1)/n)
    sdrhi_1 <- ssdr_1$par + qnorm(0.975) * sqrt(diag(varsdr_1)/n)
    varsdr_1 = diag(varsdr_1) / n
  }
  # TMLE #
  tmlU1lo_1 <- tml_1 - qnorm(0.975) * sqrt(diag(var(eiftml_1))/n)
  tmlU1hi_1 <- tml_1 + qnorm(0.975) * sqrt(diag(var(eiftml_1))/n)
  
  Jtml_1 <- Jac(stml_1$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jitml_1 <- try(solve(Jtml_1), silent = T)
  
  if(inherits(Jitml_1, 'try-error')){
    Jitml_1 <- ginv(Jtml_1)
    if(inherits(Jitml_1, 'try-error')){
      tmllo_1 <- tmlhi_1 <- rep(NA, ncol(dataphi))
      vartml_1 <- rep(NA, ncol(dataphi))
    }
  } else {
    vartml_1 <- Jitml_1 %*% var(eiftml_1) %*% t(Jitml_1)
    tmllo_1 <- stml_1$par - qnorm(0.975) * sqrt(diag(vartml_1)/n)
    tmlhi_1 <- stml_1$par + qnorm(0.975) * sqrt(diag(vartml_1)/n)
    vartml_1 = diag(vartml_1) / n
  }
  
  # Scenario 2 #
  # SDR #
  sdrU1lo_2 <- sdr_2 - qnorm(0.975) * sqrt(diag(var(eifsdr_2))/n)
  sdrU1hi_2 <- sdr_2 + qnorm(0.975) * sqrt(diag(var(eifsdr_2))/n)
  
  Jsdr_2 <- Jac(ssdr_2$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jisdr_2 <- try(solve(Jsdr_2), silent = T)
  if(inherits(Jisdr_2, 'try-error')){
    Jisdr_2 <- ginv(Jsdr_2)
    if(inherits(Jisdr_2, 'try-error')){
      sdrlo_2 <- sdrhi_2 <- rep(NA, ncol(dataphi))
      varsdr_2 <- rep(NA, ncol(dataphi))
    }
  } else {
    varsdr_2 <- Jisdr_2 %*% var(eifsdr_2) %*% t(Jisdr_2)
    sdrlo_2 <- ssdr_2$par - qnorm(0.975) * sqrt(diag(varsdr_2)/n)
    sdrhi_2 <- ssdr_2$par + qnorm(0.975) * sqrt(diag(varsdr_2)/n)
    varsdr_2 = diag(varsdr_2) / n
  }
  # TMLE #
  tmlU1lo_2 <- tml_2 - qnorm(0.975) * sqrt(diag(var(eiftml_2))/n)
  tmlU1hi_2 <- tml_2 + qnorm(0.975) * sqrt(diag(var(eiftml_2))/n)
  
  Jtml_2 <- Jac(stml_2$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jitml_2 <- try(solve(Jtml_2), silent = T)
  
  if(inherits(Jitml_2, 'try-error')){
    Jitml_2 <- ginv(Jtml_2)
    if(inherits(Jitml_2, 'try-error')){
      tmllo_2 <- tmlhi_2 <- rep(NA, ncol(dataphi))
      vartml_2 <- rep(NA, ncol(dataphi))
    }
  } else {
    vartml_2 <- Jitml_2 %*% var(eiftml_2) %*% t(Jitml_2)
    tmllo_2 <- stml_2$par - qnorm(0.975) * sqrt(diag(vartml_2)/n)
    tmlhi_2 <- stml_2$par + qnorm(0.975) * sqrt(diag(vartml_2)/n)
    vartml_2 = diag(vartml_2) / n
  }
  
  # Scenario 3 #
  # SDR #
  sdrU1lo_3 <- sdr_3 - qnorm(0.975) * sqrt(diag(var(eifsdr_3))/n)
  sdrU1hi_3 <- sdr_3 + qnorm(0.975) * sqrt(diag(var(eifsdr_3))/n)
  
  Jsdr_3 <- Jac(ssdr_3$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jisdr_3 <- try(solve(Jsdr_3), silent = T)
  if(inherits(Jisdr_3, 'try-error')){
    Jisdr_3 <- ginv(Jsdr_3)
    if(inherits(Jisdr_3, 'try-error')){
      sdrlo_3 <- sdrhi_3 <- rep(NA, ncol(dataphi))
      varsdr_3 <- rep(NA, ncol(dataphi))
    }
  } else {
    varsdr_3 <- Jisdr_3 %*% var(eifsdr_3) %*% t(Jisdr_3)
    sdrlo_3 <- ssdr_3$par - qnorm(0.975) * sqrt(diag(varsdr_3)/n)
    sdrhi_3 <- ssdr_3$par + qnorm(0.975) * sqrt(diag(varsdr_3)/n)
    varsdr_3 = diag(varsdr_3) / n
  }
  # TMLE #
  tmlU1lo_3 <- tml_3 - qnorm(0.975) * sqrt(diag(var(eiftml_3))/n)
  tmlU1hi_3 <- tml_3 + qnorm(0.975) * sqrt(diag(var(eiftml_3))/n)
  
  Jtml_3 <- Jac(stml_3$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jitml_3 <- try(solve(Jtml_3), silent = T)
  
  if(inherits(Jitml_3, 'try-error')){
    Jitml_3 <- ginv(Jtml_3)
    if(inherits(Jitml_3, 'try-error')){
      tmllo_3 <- tmlhi_3 <- rep(NA, ncol(dataphi))
      vartml_3 <- rep(NA, ncol(dataphi))
    }
  } else {
    vartml_3 <- Jitml_3 %*% var(eiftml_3) %*% t(Jitml_3)
    tmllo_3 <- stml_3$par - qnorm(0.975) * sqrt(diag(vartml_3)/n)
    tmlhi_3 <- stml_3$par + qnorm(0.975) * sqrt(diag(vartml_3)/n)
    vartml_3 = diag(vartml_3) / n
  }
  
  # Scenario 4 #
  # SDR #
  sdrU1lo_4 <- sdr_4 - qnorm(0.975) * sqrt(diag(var(eifsdr_4))/n)
  sdrU1hi_4 <- sdr_4 + qnorm(0.975) * sqrt(diag(var(eifsdr_4))/n)
  
  Jsdr_4 <- Jac(ssdr_4$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jisdr_4 <- try(solve(Jsdr_4), silent = T)
  if(inherits(Jisdr_4, 'try-error')){
    Jisdr_4 <- ginv(Jsdr_4)
    if(inherits(Jisdr_4, 'try-error')){
      sdrlo_4 <- sdrhi_4 <- rep(NA, ncol(dataphi))
      varsdr_4 <- rep(NA, ncol(dataphi))
    }
  } else {
    varsdr_4 <- Jisdr_4 %*% var(eifsdr_4) %*% t(Jisdr_4)
    sdrlo_4 <- ssdr_4$par - qnorm(0.975) * sqrt(diag(varsdr_4)/n)
    sdrhi_4 <- ssdr_4$par + qnorm(0.975) * sqrt(diag(varsdr_4)/n)
    varsdr_4 = diag(varsdr_4) / n
  }
  # TMLE #
  tmlU1lo_4 <- tml_4 - qnorm(0.975) * sqrt(diag(var(eiftml_4))/n)
  tmlU1hi_4 <- tml_4 + qnorm(0.975) * sqrt(diag(var(eiftml_4))/n)
  
  Jtml_4 <- Jac(stml_4$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jitml_4 <- try(solve(Jtml_4), silent = T)
  
  if(inherits(Jitml_4, 'try-error')){
    Jitml_4 <- ginv(Jtml_4)
    if(inherits(Jitml_4, 'try-error')){
      tmllo_4 <- tmlhi_4 <- rep(NA, ncol(dataphi))
      vartml_4 <- rep(NA, ncol(dataphi))
    }
  } else {
    vartml_4 <- Jitml_4 %*% var(eiftml_4) %*% t(Jitml_4)
    tmllo_4 <- stml_4$par - qnorm(0.975) * sqrt(diag(vartml_4)/n)
    tmlhi_4 <- stml_4$par + qnorm(0.975) * sqrt(diag(vartml_4)/n)
    vartml_4 = diag(vartml_4) / n
  }
  
  # Scenario 5 #
  # SDR #
  sdrU1lo_5 <- sdr_5 - qnorm(0.975) * sqrt(diag(var(eifsdr_5))/n)
  sdrU1hi_5 <- sdr_5 + qnorm(0.975) * sqrt(diag(var(eifsdr_5))/n)
  
  Jsdr_5 <- Jac(ssdr_5$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jisdr_5 <- try(solve(Jsdr_5), silent = T)
  if(inherits(Jisdr_5, 'try-error')){
    Jisdr_5 <- ginv(Jsdr_5)
    if(inherits(Jisdr_5, 'try-error')){
      sdrlo_5 <- sdrhi_5 <- rep(NA, ncol(dataphi))
      varsdr_5 <- rep(NA, ncol(dataphi))
    }
  } else {
    varsdr_5 <- Jisdr_5 %*% var(eifsdr_5) %*% t(Jisdr_5)
    sdrlo_5 <- ssdr_5$par - qnorm(0.975) * sqrt(diag(varsdr_5)/n)
    sdrhi_5 <- ssdr_5$par + qnorm(0.975) * sqrt(diag(varsdr_5)/n)
    varsdr_5 = diag(varsdr_5) / n
  }
  # TMLE #
  tmlU1lo_5 <- tml_5 - qnorm(0.975) * sqrt(diag(var(eiftml_5))/n)
  tmlU1hi_5 <- tml_5 + qnorm(0.975) * sqrt(diag(var(eiftml_5))/n)
  
  Jtml_5 <- Jac(stml_5$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jitml_5 <- try(solve(Jtml_5), silent = T)
  
  if(inherits(Jitml_5, 'try-error')){
    Jitml_5 <- ginv(Jtml_5)
    if(inherits(Jitml_5, 'try-error')){
      tmllo_5 <- tmlhi_5 <- rep(NA, ncol(dataphi))
      vartml_5 <- rep(NA, ncol(dataphi))
    }
  } else {
    vartml_5 <- Jitml_5 %*% var(eiftml_5) %*% t(Jitml_5)
    tmllo_5 <- stml_5$par - qnorm(0.975) * sqrt(diag(vartml_5)/n)
    tmlhi_5 <- stml_5$par + qnorm(0.975) * sqrt(diag(vartml_5)/n)
    vartml_5 = diag(vartml_5) / n
  }
  
  
  ##############################################################################
  
  
  ##############################################################################
  if(trace) print("Step 5: Solving CG...")
  # scenario 1 #
  ssdr_1_CG <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                     A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                     la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                     fn = fsdr1_1, method = "CG",
                     control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  stml_1_CG <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                     A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                     la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                     fn = ftml1_1, method = "CG",
                     control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  
  # scenario 2 #
  ssdr_2_CG <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                     A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                     la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                     fn = fsdr1_2, method = "CG",
                     control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  stml_2_CG <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                     A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                     la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                     fn = ftml1_2, method = "CG",
                     control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  
  # scenario 3 #
  ssdr_3_CG <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                     A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                     la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                     fn = fsdr1_3, method = "CG",
                     control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  stml_3_CG <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                     A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                     la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                     fn = ftml1_3, method = "CG",
                     control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  
  # scenario 4 #
  ssdr_4_CG <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                     A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                     la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                     fn = fsdr1_4, method = "CG",
                     control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  stml_4_CG <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                     A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                     la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                     fn = ftml1_4, method = "CG",
                     control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  # scenario 5 #
  ssdr_5_CG <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                     A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                     la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                     fn = fsdr1_5, method = "CG",
                     control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  stml_5_CG <- optim(par = coef(IPW_1), A = A, lambda = lambda, phi = phi, m = m,
                     A_type = A_type, A.min = A.min, A.max = A.max, trunc = trunc, tau = tau, outcome = outcome,
                     la = la, A_span = A_span, bindim = bindim, approx_int = approx_int,
                     fn = ftml1_5, method = "CG",
                     control = list(abstol = 0.001 * log(n)/n, maxit = 10^6, reltol = .Machine$double.eps))
  
  ##############################################################################
  
  
  if(trace) print("Step 6: Estimating variances CG")
  # Scenario 1 #
  # SDR #
  Jsdr_1_CG <- Jac(ssdr_1_CG$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                   trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jisdr_1_CG <- try(solve(Jsdr_1_CG), silent = T)
  if(inherits(Jisdr_1_CG, 'try-error')){
    Jisdr_1_CG <- ginv(Jsdr_1_CG)
    if(inherits(Jisdr_1_CG, 'try-error')){
      sdrlo_1_CG <- sdrhi_1_CG <- rep(NA, ncol(dataphi))
      varsdr_1_CG <- rep(NA, ncol(dataphi))
    }
  } else {
    varsdr_1_CG <- Jisdr_1_CG %*% var(eifsdr_1) %*% t(Jisdr_1_CG)
    sdrlo_1_CG <- ssdr_1_CG$par - qnorm(0.975) * sqrt(diag(varsdr_1_CG)/n)
    sdrhi_1_CG <- ssdr_1_CG$par + qnorm(0.975) * sqrt(diag(varsdr_1_CG)/n)
    varsdr_1_CG = diag(varsdr_1_CG) / n
  }
  # TMLE #
  Jtml_1_CG <- Jac(stml_1_CG$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                   trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jitml_1_CG <- try(solve(Jtml_1_CG), silent = T)
  
  if(inherits(Jitml_1_CG, 'try-error')){
    Jitml_1_CG <- ginv(Jtml_1_CG)
    if(inherits(Jitml_1_CG, 'try-error')){
      tmllo_1_CG <- tmlhi_1_CG <- rep(NA, ncol(dataphi))
      vartml_1_CG <- rep(NA, ncol(dataphi))
    }
  } else {
    vartml_1_CG <- Jitml_1_CG %*% var(eiftml_1) %*% t(Jitml_1_CG)
    tmllo_1_CG <- stml_1_CG$par - qnorm(0.975) * sqrt(diag(vartml_1_CG)/n)
    tmlhi_1_CG <- stml_1_CG$par + qnorm(0.975) * sqrt(diag(vartml_1_CG)/n)
    vartml_1_CG = diag(vartml_1_CG) / n
  }
  
  # Scenario 2 #
  # SDR #
  Jsdr_2_CG <- Jac(ssdr_2_CG$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                   trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jisdr_2_CG <- try(solve(Jsdr_2_CG), silent = T)
  if(inherits(Jisdr_2_CG, 'try-error')){
    Jisdr_2_CG <- ginv(Jsdr_2_CG)
    if(inherits(Jisdr_2_CG, 'try-error')){
      sdrlo_2_CG <- sdrhi_2_CG <- rep(NA, ncol(dataphi))
      varsdr_2_CG <- rep(NA, ncol(dataphi))
    }
  } else {
    varsdr_2_CG <- Jisdr_2_CG %*% var(eifsdr_2) %*% t(Jisdr_2_CG)
    sdrlo_2_CG <- ssdr_2_CG$par - qnorm(0.975) * sqrt(diag(varsdr_2_CG)/n)
    sdrhi_2_CG <- ssdr_2_CG$par + qnorm(0.975) * sqrt(diag(varsdr_2_CG)/n)
    varsdr_2_CG = diag(varsdr_2_CG) / n
  }
  # TMLE #
  Jtml_2_CG <- Jac(stml_2_CG$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                   trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jitml_2_CG <- try(solve(Jtml_2_CG), silent = T)
  
  if(inherits(Jitml_2_CG, 'try-error')){
    Jitml_2_CG <- ginv(Jtml_2_CG)
    if(inherits(Jitml_2_CG, 'try-error')){
      tmllo_2_CG <- tmlhi_2_CG <- rep(NA, ncol(dataphi))
      vartml_2_CG <- rep(NA, ncol(dataphi))
    }
  } else {
    vartml_2_CG <- Jitml_2_CG %*% var(eiftml_2) %*% t(Jitml_2_CG)
    tmllo_2_CG <- stml_2_CG$par - qnorm(0.975) * sqrt(diag(vartml_2_CG)/n)
    tmlhi_2_CG <- stml_2_CG$par + qnorm(0.975) * sqrt(diag(vartml_2_CG)/n)
    vartml_2_CG = diag(vartml_2_CG) / n
  }
  
  # Scenario 3 #
  # SDR #
  Jsdr_3_CG <- Jac(ssdr_3_CG$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                   trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jisdr_3_CG <- try(solve(Jsdr_3_CG), silent = T)
  if(inherits(Jisdr_3_CG, 'try-error')){
    Jisdr_3_CG <- ginv(Jsdr_3_CG)
    if(inherits(Jisdr_3_CG, 'try-error')){
      sdrlo_3_CG <- sdrhi_3_CG <- rep(NA, ncol(dataphi))
      varsdr_3_CG <- rep(NA, ncol(dataphi))
    }
  } else {
    varsdr_3_CG <- Jisdr_3_CG %*% var(eifsdr_3) %*% t(Jisdr_3_CG)
    sdrlo_3_CG <- ssdr_3_CG$par - qnorm(0.975) * sqrt(diag(varsdr_3_CG)/n)
    sdrhi_3_CG <- ssdr_3_CG$par + qnorm(0.975) * sqrt(diag(varsdr_3_CG)/n)
    varsdr_3_CG = diag(varsdr_3_CG) / n
  }
  # TMLE #
  Jtml_3_CG <- Jac(stml_3_CG$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                   trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jitml_3_CG <- try(solve(Jtml_3_CG), silent = T)
  
  if(inherits(Jitml_3_CG, 'try-error')){
    Jitml_3_CG <- ginv(Jtml_3_CG)
    if(inherits(Jitml_3_CG, 'try-error')){
      tmllo_3_CG <- tmlhi_3_CG <- rep(NA, ncol(dataphi))
      vartml_3_CG <- rep(NA, ncol(dataphi))
    }
  } else {
    vartml_3_CG <- Jitml_3_CG %*% var(eiftml_3) %*% t(Jitml_3_CG)
    tmllo_3_CG <- stml_3_CG$par - qnorm(0.975) * sqrt(diag(vartml_3_CG)/n)
    tmlhi_3_CG <- stml_3_CG$par + qnorm(0.975) * sqrt(diag(vartml_3_CG)/n)
    vartml_3_CG = diag(vartml_3_CG) / n
  }
  
  # Scenario 4 #
  # SDR #
  Jsdr_4_CG <- Jac(ssdr_4_CG$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                   trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jisdr_4_CG <- try(solve(Jsdr_4_CG), silent = T)
  if(inherits(Jisdr_4_CG, 'try-error')){
    Jisdr_4_CG <- ginv(Jsdr_4_CG)
    if(inherits(Jisdr_4_CG, 'try-error')){
      sdrlo_4_CG <- sdrhi_4_CG <- rep(NA, ncol(dataphi))
      varsdr_4_CG <- rep(NA, ncol(dataphi))
    }
  } else {
    varsdr_4_CG <- Jisdr_4_CG %*% var(eifsdr_4) %*% t(Jisdr_4_CG)
    sdrlo_4_CG <- ssdr_4_CG$par - qnorm(0.975) * sqrt(diag(varsdr_4_CG)/n)
    sdrhi_4_CG <- ssdr_4_CG$par + qnorm(0.975) * sqrt(diag(varsdr_4_CG)/n)
    varsdr_4_CG = diag(varsdr_4_CG) / n
  }
  # TMLE #
  Jtml_4_CG <- Jac(stml_4_CG$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                   trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jitml_4_CG <- try(solve(Jtml_4_CG), silent = T)
  
  if(inherits(Jitml_4_CG, 'try-error')){
    Jitml_4_CG <- ginv(Jtml_4_CG)
    if(inherits(Jitml_4_CG, 'try-error')){
      tmllo_4_CG <- tmlhi_4_CG <- rep(NA, ncol(dataphi))
      vartml_4_CG <- rep(NA, ncol(dataphi))
    }
  } else {
    vartml_4_CG <- Jitml_4_CG %*% var(eiftml_4) %*% t(Jitml_4_CG)
    tmllo_4_CG <- stml_4_CG$par - qnorm(0.975) * sqrt(diag(vartml_4_CG)/n)
    tmlhi_4_CG <- stml_4_CG$par + qnorm(0.975) * sqrt(diag(vartml_4_CG)/n)
    vartml_4_CG = diag(vartml_4_CG) / n
  }
  
  # Scenario 5 #
  # SDR #
  Jsdr_5_CG <- Jac(ssdr_5_CG$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                   trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jisdr_5_CG <- try(solve(Jsdr_5_CG), silent = T)
  if(inherits(Jisdr_5_CG, 'try-error')){
    Jisdr_5_CG <- ginv(Jsdr_5_CG)
    if(inherits(Jisdr_5_CG, 'try-error')){
      sdrlo_5_CG <- sdrhi_5_CG <- rep(NA, ncol(dataphi))
      varsdr_5_CG <- rep(NA, ncol(dataphi))
    }
  } else {
    varsdr_5_CG <- Jisdr_5_CG %*% var(eifsdr_5) %*% t(Jisdr_5_CG)
    sdrlo_5_CG <- ssdr_5_CG$par - qnorm(0.975) * sqrt(diag(varsdr_5_CG)/n)
    sdrhi_5_CG <- ssdr_5_CG$par + qnorm(0.975) * sqrt(diag(varsdr_5_CG)/n)
    varsdr_5_CG = diag(varsdr_5_CG) / n
  }
  # TMLE #
  Jtml_5_CG <- Jac(stml_5_CG$par, A = A, lambda = lambda,  phi = phi, m = m, A_type = A_type, A.min = A.min, A.max = A.max,
                   trunc = trunc, tau = tau, outcome = outcome, la = la, A_span = A_span, bindim = bindim, approx_int = approx_int)
  Jitml_5_CG <- try(solve(Jtml_5_CG), silent = T)
  
  if(inherits(Jitml_5_CG, 'try-error')){
    Jitml_5_CG <- ginv(Jtml_5_CG)
    if(inherits(Jitml_5_CG, 'try-error')){
      tmllo_5_CG <- tmlhi_5_CG <- rep(NA, ncol(dataphi))
      vartml_5_CG <- rep(NA, ncol(dataphi))
    }
  } else {
    vartml_5_CG <- Jitml_5_CG %*% var(eiftml_5) %*% t(Jitml_5_CG)
    tmllo_5_CG <- stml_5_CG$par - qnorm(0.975) * sqrt(diag(vartml_5_CG)/n)
    tmlhi_5_CG <- stml_5_CG$par + qnorm(0.975) * sqrt(diag(vartml_5_CG)/n)
    vartml_5_CG = diag(vartml_5_CG) / n
  }
  
  
  
  
  ##############################################################################
  
  
  # make results #
  type <- c("Intercept","slope")
  
  # scenario 1 #
  out_1 <- data.frame(scen = 1, type = type,
                      # NM #
                      sdr = ssdr_1$par, sdrlo = sdrlo_1, sdrhi = sdrhi_1, sdrU1 = sdr_1, sdrU1lo = sdrU1lo_1, sdrU1hi = sdrU1hi_1, varsdr = varsdr_1,
                      tml = stml_1$par, tmllo = tmllo_1, tmlhi = tmlhi_1, tmlU1 = tml_1, tmlU1lo = tmlU1lo_1, tmlU1hi = tmlU1hi_1, vartml = vartml_1,
                      # CG #
                      sdr_CG = ssdr_1_CG$par, sdrlo_CG = sdrlo_1_CG, sdrhi_CG = sdrhi_1_CG, varsdr_CG = varsdr_1_CG,
                      tml_CG = stml_1_CG$par, tmllo_CG = tmllo_1_CG, tmlhi_CG = tmlhi_1_CG, vartml_CG = vartml_1_CG,
                      n = n,
                      IPW = coef(IPW_1), naive = coef(naive_fit))
  # scenario 2 #
  out_2 <- data.frame(scen = 2, type = type,
                      # NM #
                      sdr = ssdr_2$par, sdrlo = sdrlo_2, sdrhi = sdrhi_2, sdrU1 = sdr_2, sdrU1lo = sdrU1lo_2, sdrU1hi = sdrU1hi_2, varsdr = varsdr_2,
                      tml = stml_2$par, tmllo = tmllo_2, tmlhi = tmlhi_2, tmlU1 = tml_2, tmlU1lo = tmlU1lo_2, tmlU1hi = tmlU1hi_2, vartml = vartml_2,
                      # CG #
                      sdr_CG = ssdr_2_CG$par, sdrlo_CG = sdrlo_2_CG, sdrhi_CG = sdrhi_2_CG, varsdr_CG = varsdr_2_CG,
                      tml_CG = stml_2_CG$par, tmllo_CG = tmllo_2_CG, tmlhi_CG = tmlhi_2_CG, vartml_CG = vartml_2_CG,
                      n = n,
                      IPW = coef(IPW_2), naive = coef(naive_fit))
  # scenario 3 #
  out_3 <- data.frame(scen = 3, type = type,
                      # NM #
                      sdr = ssdr_3$par, sdrlo = sdrlo_3, sdrhi = sdrhi_3, sdrU1 = sdr_3, sdrU1lo = sdrU1lo_3, sdrU1hi = sdrU1hi_3, varsdr = varsdr_3,
                      tml = stml_3$par, tmllo = tmllo_3, tmlhi = tmlhi_3, tmlU1 = tml_3, tmlU1lo = tmlU1lo_3, tmlU1hi = tmlU1hi_3, vartml = vartml_3,
                      # CG #
                      sdr_CG = ssdr_3_CG$par, sdrlo_CG = sdrlo_3_CG, sdrhi_CG = sdrhi_3_CG, varsdr_CG = varsdr_3_CG,
                      tml_CG = stml_3_CG$par, tmllo_CG = tmllo_3_CG, tmlhi_CG = tmlhi_3_CG, vartml_CG = vartml_3_CG,
                      n = n,
                      IPW = coef(IPW_3), naive = coef(naive_fit))
  
  # scenario 4 #
  out_4 <- data.frame(scen = 4, type = type,
                      # NM #
                      sdr = ssdr_4$par, sdrlo = sdrlo_4, sdrhi = sdrhi_4, sdrU1 = sdr_4, sdrU1lo = sdrU1lo_4, sdrU1hi = sdrU1hi_4, varsdr = varsdr_4,
                      tml = stml_4$par, tmllo = tmllo_4, tmlhi = tmlhi_4, tmlU1 = tml_4, tmlU1lo = tmlU1lo_4, tmlU1hi = tmlU1hi_4, vartml = vartml_4,
                      # CG #
                      sdr_CG = ssdr_4_CG$par, sdrlo_CG = sdrlo_4_CG, sdrhi_CG = sdrhi_4_CG, varsdr_CG = varsdr_4_CG,
                      tml_CG = stml_4_CG$par, tmllo_CG = tmllo_4_CG, tmlhi_CG = tmlhi_4_CG, vartml_CG = vartml_4_CG,
                      n = n,
                      IPW = coef(IPW_4), naive = coef(naive_fit))
  
  # scenario 5 #
  out_5 <- data.frame(scen = 5, type = type,
                      # NM #
                      sdr = ssdr_5$par, sdrlo = sdrlo_5, sdrhi = sdrhi_5, sdrU1 = sdr_5, sdrU1lo = sdrU1lo_5, sdrU1hi = sdrU1hi_5, varsdr = varsdr_5,
                      tml = stml_5$par, tmllo = tmllo_5, tmlhi = tmlhi_5, tmlU1 = tml_5, tmlU1lo = tmlU1lo_5, tmlU1hi = tmlU1hi_5, vartml = vartml_5,
                      # CG #
                      sdr_CG = ssdr_5_CG$par, sdrlo_CG = sdrlo_5_CG, sdrhi_CG = sdrhi_5_CG, varsdr_CG = varsdr_5_CG,
                      tml_CG = stml_5_CG$par, tmllo_CG = tmllo_5_CG, tmlhi_CG = tmlhi_5_CG, vartml_CG = vartml_5_CG,
                      n = n,
                      IPW = coef(IPW_5), naive = coef(naive_fit))
  
  # merge #
  out <- rbind(out_1, out_2, out_3, out_4, out_5)
  
  return(list(
    # combined results #
    "out" = out,
    # EIFs #
    # sdr #
    "eifsdr_1" = eifsdr_1, "eifsdr_2" = eifsdr_2, "eifsdr_3" = eifsdr_3, "eifsdr_4" = eifsdr_4, "eifsdr_5" = eifsdr_5,
    # tml #
    "eiftml_1" = eiftml_1, "eiftml_2" = eiftml_2, "eiftml_3" = eiftml_3, "eiftml_4" = eiftml_4, "eiftml_5" = eiftml_5,
    # weights #
    "Zn_1" = Zn_1, "Zn_2" = Zn_2, "Zn_3" = Zn_3, "Zn_4" = Zn_4, "Zn_5" = Zn_5
  ))
  
  
}


