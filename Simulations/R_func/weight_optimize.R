#' Test flat null using one-step estimator
#'
#' @param psi_G one-step estimators for psi.
#' @param n Sample size.
#' @param G n by d matrix. For column i is applying function g_i(.) to a_j, where j is 1,...,n.
#' @param G_c Mean-centered G matrix, an n x p matrix.
#' @param structure structure constraint: flexible omega's or monotone increasing omega's.
#' @param lambda scale parameter required for structure constraint on omega's.
#' @param norm scale constraint: L1-norm or L2-norm.
#'
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
weight_optimize <- function(
    psi_G,
    n,
    G,
    G_c,
    # res_OneStep,
    structure = "flexible",
    # structure = "monotone",
    lambda = 2,
    norm  = "L1"
    ){
  # size of G sets
  d <- length(psi_G)
  
  # contrast for monotone c1: beta_{i-1} - beta_i <= 0
  C <- diag(1, nrow = d-1, ncol = d)
  C[cbind(1:(d-1), 2:d)] <- -1
  
  # pre-calculation for L2-norm c2
  V <- t(G_c) %*% G_c
  
  
  # for infinite lambda, we have closed form for both
  if (is.infinite(lambda) & norm =="L1"){
    # if ( norm == "L1"){
      # weights
      omega <- sign(psi_G)
      # h = sum of (weight * psi_G)
      value <- t(omega) %*% psi_G
      
      res <- list(omega_ClosedForm = c(omega),
                  value = c(value))
    # }else{ # norm == "L2"
    #   omega <- sqrt(n) * MASS::ginv(V) %*% psi_G/c(sqrt(t(psi_G) %*% MASS::ginv(V) %*% psi_G))
    #   value <- t(omega) %*% psi_G
    # 
    #   res <- list(omega_ClosedForm = c(omega),
    #               value = c(value))
    # }
  }else{
    # define the weights of basis functions to be estimated
    omega <- CVXR::Variable(d)
    
    #------------------------------
    # constraints for optimization
    #------------------------------
    # structure constraint: the complexity of H is bounded
    # (i) h is not too complex (differentiate the most influencer)
    # (i) h is not too complex (differentiate the most influencer)
    if (structure  == "flexible"){
      # sum |w_i+1 - w_i| <= lambda
      c1 <- tv(omega) <= lambda
      # sum(abs(G %*% omega)) <= lambda
    }
    # (ii) assume a monotone increasing effect trend: w1 < w2 < ... < wd
    if (structure == "monotone"){
      c1 <- C %*% omega <= 0
      # c1 <- diff(omega) <= 0
      # c12 <- C %*% omega >= 0
    } 
    
    
    #scale constraint: supremum is well-defined; define the norm
    # L1-norm: 
    if (norm  == "L1"){
      # max (|omega|)<=1
      c2 <- cvxr_norm(omega, "inf") <= 1
      # c2 <- max(abs(omega)) <= 1
    }
    # L2-norm:  var(h) = sum(h^2) = 1
    if (norm == "L2"){
      c2 <-  quad_form(omega, V) <= n
    }
    
    #------------------------------
    #    optimization
    #------------------------------
    # constraints
    if (is.infinite(lambda)){
      cs <- list(c2)
    }else{
      cs <- list(c1, c2)
    }
    
    
    # supreme of |\psi(h)|, h in H
    objective1 <- CVXR::Maximize( t(omega) %*% psi_G )
    objective2 <- CVXR::Minimize( t(omega) %*% psi_G )
    
    
    # the optimization results for both direction
    suppressWarnings(result_max <- solve(Problem(objective1,
                                                 constraints = cs),
                                         num_iter = 1000, reltol = 1e-02, feastol = 1e-04, abstol = 1e-04)
    )
    # final result list
    result <- result_max
    
    if (structure == "monotone"){
      suppressWarnings(result_min <- solve(Problem(objective2,
                                                   constraints = cs),
                                           num_iter = 1000, reltol = 1e-02, feastol = 1e-04, abstol = 1e-04))
      
      # refit for another direction if assuming monotone effect
      if (result_min$status != "solver_error"){
        if (result_max$status != "solver_error"){ # min exists, max exists
          if (abs(result_max$value) < abs(result_min$value)){
            result <- result_min
          }
        }else{ # min exists, max doesn't exists
          result <- result_min
        }
      }
    }
    
    # output list
    if (result$status != "solver_error" ){
      res <- list(omega_cvxr = as.vector(result$getValue(omega)),
                  value = result$value)
    }else{
      res <- NULL
    }
    
  }
  
  # if (ClosedForm){
  #   # it works for small finite set of H
  #   coef_G_closeform <- solve(t(G) %*% G) %*% psi_G
  #   # return both weights
  #   res <- list(omega_cvxr = as.vector(result$getValue(omega)),
  #               value = result$value,
  #               omega_ClosedForm = coef_G_closeform)
  # }
  
  
  return(res)
}
