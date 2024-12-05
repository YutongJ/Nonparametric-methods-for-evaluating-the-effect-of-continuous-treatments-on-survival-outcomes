#' Propensity score estimation
#'
#' Estimate the generalized propensity score using kernel smoothing
#' 
#' @param a Predictor of interest, an n-dimensional vector.
#' @param w Covariates for adjustment, an n x p matrix.
#' @param no.folds Number of folds to be used for selection of bandwidth via cross-validation
#' @param no.eval Number of evaluation points (levels of predictor of interest) at which the conditional density will be estimated.
#'                Estimates at intermediate points are estimated using linear interpolation.
#' @param bw Value for bandwidth, if selected a priori. Set to NULL by default.
#'
#'  @return A list containing the following: 
#'         \describe{
#'            \item{cond.dens} Estimate of conditional density of predictor of interest, given covariates.
#'            \item{marg.dens} Estimate of marginal density of predictor of interest.
#'            \item{prop.score} Ratio of marginal density to conditional density.
#'          }             
#' @return An n-dimensional vector containing the conditional density, marginal density, and their ratior.
#'         
#' @export
EstPropScore <- function(a, w, no.folds = 10, no.eval = 20, bw = NULL) {
  n <- length(a)
  a.eval <- seq(min(a), max(a), length.out = no.eval)
  
  # use bandwidth that is optimal for kernel estimation of margianl density
  # is sensible choice if con'l density is about as smooth in a as marginal density
  # could try other choices via cross-validation, but this is comp'l cheaper
  # con'l density estimation is generally hard, and needs further study
  
  if(is.null(bw)) {
    # select bandwidth via cross-validation
    folds <- sample(cut(1:n, breaks = no.folds, labels = FALSE))
    bw.seq <- exp(seq(log(.01 * sd(a)), log(25 * sd(a)), length.out = 100))
    
    pred.error <- matrix(NA, nrow = 100, ncol = no.folds) # prediction error
    for(i in 1:no.folds) {
      for(j in 1:100) {
        a.train <- a[folds != i]
        a.test <- a[folds == i]
        
        n.test <- length(a.test)
        fit.test <- numeric(n.test)
        
        for(k in 1:n.test) {
          kern <- dnorm(a[folds != i], mean = a.test[k], sd = bw.seq[j])
          fit.test[k] <- mean(kern)
        }
        
        pred.error[j,i] <- mean(-log(fit.test))
      }
    }
    avg.pred.error <- apply(pred.error, 1, mean)
    se.pred.error <- apply(pred.error, 1, sd)/sqrt(no.folds)
    bw.opt <- max(bw.seq[avg.pred.error < min(avg.pred.error) + se.pred.error[which.min(avg.pred.error)]])
    # bw.opt <- bw.seq[which.min(apply(pred.error, 1, mean))]
  } else {
    # set bandwidth to a pre-specified value, if one is provided
    bw.opt <- bw
  }
  
  marg.dens <- numeric(n)
  for(i in 1:n) {
    marg.dens[i] <- mean(dnorm(a, mean = a[i], sd = 2 * bw.opt))
  }
  
  f.eval <- matrix(NA, nrow = n, ncol = no.eval)
  for(j in 1:no.eval) {
    kern.aj <- dnorm(a, mean = a.eval[j], sd = bw.opt)
    # hal.fit <- fit_hal(X = w, Y = kern.aj)
    hal.fit <- fit_hal(X = w, Y = kern.aj,
                       family = "poisson", return_x_basis = TRUE)
    f.eval[,j] <- predict(hal.fit, new_data = w)
    # f.eval[,j] <- exp(hal.fit$X.basis %*% hal.fit$coefs)
  }
  
  # set cond.dens[i,j] = cond.dens(a[i], w[j])
  cond.dens <- matrix(NA, n, n)
  for(i in 1:n) {
    k <- min(which(a.eval - a[i] >= 0))
    for(j in 1:n) {
      if(k == 1) {
        cond.dens[i,j] <- f.eval[j,1]
      } else if(k == no.eval) {
        cond.dens[i,j] <- f.eval[j,no.eval]
      } else {
        t <- (a.eval[k] - a[i])/(a.eval[k] - a.eval[k-1])
        cond.dens[i,j] <- f.eval[j,k] * (1-t) + f.eval[j,k-1] * t
      }
    }
  }
  
  marg.dens <- apply(cond.dens, 1, mean)
  prop.score <- diag(marg.dens) %*% (cond.dens^(-1))
  
  out <- list(cond.dens = cond.dens, marg.dens = marg.dens,
              prop.score = prop.score)
  return(out)
}
