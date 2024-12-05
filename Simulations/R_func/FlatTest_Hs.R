#' Test flat null using one-step estimator
#'
#' @param res_OneStep Output from `OneStep_H.R`
#' @param sign We can use sign function for L1-norm + monotone. Default is FALSE.
#' @param structure structure constraint: "flexible" or "monotone" omega's.
#' @param lambda scale parameter required for structure constraint on omega's.
#' @param norm scale constraint: L1-norm or L2-norm.
#' @param null_size The size of null sets for getting empirical p value.
#'                  
#' @return A list containing the following: 
#'         \describe{
#'            \item{nullset} A random set of $\Psi$ under null distribution.
#'            \item{t_stat} Test statistic for $H_0$: 
#'            \item{eif} Estimate of efficient influence function for expectation of drf times any given basis function.
#'            \item{test.stat} Test statistic for flat null.
#'            \item{bs.samples} Bootstrap samples from null limiting distribution of test statistic.
#'            \item{p.val} p-value for test of flat null hypothesis.
#'          }
#'
#' @export
FlatTest <- function(
    res_OneStep,
    sign = FALSE,
    structure = NULL,
    lambda = NULL,
    norm  = NULL,
    # structure = "flexible",
    # lambda = 2,
    # norm  = "L1",
    null_size = 5000
  ){
  # sample size
  n <- dim(res_OneStep$G)[1]
  # size of H sets
  d <- dim(res_OneStep$G)[2]
  
  
  # (1) monotone + L1: sign function
  if (sign){
    #--------------
    # empirical p
    #--------------
    nullset <- apply(MASS::mvrnorm(null_size, mu = rep(0, d), Sigma = res_OneStep$cov),
                     1, function(x){max(abs(x))})
    # test statistic
    t_stat <- max(abs(res_OneStep$onestep))
    # p-value
    p <- mean(nullset >= t_stat)
  }else{
    # weights (from CVXR) for one-step estimator
    w_onestep <- weight_optimize(psi_G = res_OneStep$onestep,
                                 n = n,
                                 G = res_OneStep$G,
                                 G_c = res_OneStep$G_c,
                                 structure = structure,
                                 lambda = lambda,
                                 norm  = norm)
    
    #--------------
    # empirical p
    #--------------
    # null sets of one-step estimator of indicator_knots function
    null_psiG <- MASS::mvrnorm(null_size, mu = rep(0, d), Sigma = res_OneStep$cov)
    # null weights for each observed null null set
    null_weights <- apply(null_psiG, 1, weight_optimize, 
                          n = n,
                          G = res_OneStep$G,
                          G_c = res_OneStep$G_c,
                          structure = structure,
                          lambda = lambda,
                          norm  = norm)
    
    # alternative for-loop [equivalent, a little bit slower than apply function]
    # null_weights <- vector(mode = "list", length = 100)
    # for ( i in 1:30){
    #   null_weights[[i]] <- weight_optimize(psi_G = null_psiG[i,],
    #                                        n = length(A),
    #                                        G = res_OneStep$G,
    #                                        G_c = res_OneStep$G_c,
    #                                        structure = structure,
    #                                        lambda = lambda,
    #                                        norm  = norm,
    #                                        ClosedForm = FALSE)
    #   print(i)
    # }
    
    # null set of h: abs(t(weights) %*% G)
    # weights are different for each observed set of psi_G
    nullset <- abs(do.call(c, lapply(null_weights, function(x){x[[2]]})))
    # test statistic: sup |psi(h)|, h \in H
    t_stat <- abs(w_onestep$value)
    # empirical p-value
    p <- mean(nullset >= t_stat) # 0.07507508
  }
  
  
  # list of results
  res <- list(nullset = nullset,
              t_stat = t_stat, 
              p = p)
  
  return(res)
}
