#' Test flat null using one-step estimator
#'
#' @param ftime Survival time, an n-dimesional vector
#' @param ftype Survival censoring indicator. 1 for event.
#' @param A treatment, an n-dimensional vector.
#' @param W Covariates for adjustment, an n x p matrix.

#' @param S Estimate of conditional survival summary given covariates and treatment
#' @param int_vals S*Integral to a pre-specified time point tau. Used in calculating EIF.
#' @param ps_ratio Ratio of marginal density of exposure to conditional density, given covariates
#' @param As All posible treatment levels.

#' @param G_func A set of transformations $g(A)$.
#' @param G_type Type of basis function.
#'                  
#' @return A list containing the following: 
#'         \describe{
#'            \item{drf.est} Estimate of the dose response function, evaluated at a.
#'            \item{tmle} Estimate of expectation of drf times any given basis function.
#'            \item{eif} Estimate of efficient influence function for expectation of drf times any given basis function.
#'            \item{test.stat} Test statistic for flat null.
#'            \item{bs.samples} Bootstrap samples from null limiting distribution of test statistic.
#'            \item{p.val} p-value for test of flat null hypothesis.
#'          }
#'
#' @export
OneStep <- function(
  ftime,
  ftype, 
  A, 
  W, 
  S = pred_SS[[t0_index+1]], 
  int_vals = IFs_SS[[t0_index]],
  ps_ratio,
  As,
  G_func,
  G_type
  # continuous
  # no.bs = 1000,
  # no.folds = 10L,
  # gamma.hat = NULL
  ){
  # sample size
  n <- length(A)
  d <- length(G_func)
  
  # Estimate of counterfactual mean at each observed trt level - vector:3*1
  # theta(a)
  drf_est_a <- apply(S, 1, mean)
  # theta(A_i)
  drf_est <- drf_est_a[match(A, As)]
  # E{ theta_bar(A) }
  drf_est_c <- drf_est - mean(drf_est)
  
  # collection of transformations (binary trt - n*3)
  G <- do.call(cbind,
               lapply(G_func, function(x){x(A)}))
  
  # center columns of G: h(a) − E[h(A)] - dim: n*3
  # G_c <- G - matrix(colMeans(G), nrow = n, ncol = length(G_func), byrow = TRUE)
  G_c <- (diag(n) - matrix(1/n, n, n)) %*% G
  
  # Compute plug-in & one-step estimator
  # plug-in: E_P{ ( theta_bar(a) ) * ( h(a) − E[h(A)] ) }  - dim: 5*1
  plugin <- t(G) %*% (drf_est_c) / n
  
  # theta_bar(A, tau) + integral part adapted from Ted 
  part1 <- drf_est_c + ps_ratio[cbind(match(A, As), 1:n)] * int_vals[cbind(match(A, As), 1:n)]

  # E_{P,A} [ S(t|A,w) * ( h(A) − E_P[h(A)] ) ]
  part2 <- t(S[match(A,As),]) %*% G_c / n
  
  # 2 * E_{P,A} [ h(a) * theta_bar(a) ]
  part3 <- 2 * matrix(plugin, nrow = n, ncol = d, byrow = TRUE)
  
  # Efficient influence function (EIF) - (should be 5*2)
  eif <- diag(part1) %*% G_c + part2 - part3
  
  
  # one-step - dim: 5*1
  onestep <- plugin + colMeans(eif)
  # estimated variance-covariance matrix of proposed estimands
  cov <- var(eif) / n
  

  # list of results
  res <- list(drf_est = drf_est,
              plugin = plugin,
              G = G,
              G_c = G_c,
              eif = eif,
              onestep = onestep,
              cov = cov)
  
  return(res)
}