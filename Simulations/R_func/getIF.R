#-----------------------
# get influence function (IF)
#-----------------------
#' Data Generate in Simulation
#'
#' Creates a dataset including Y, Delta, S, Trt, W1 and W2.
#'
#' @param Y observed time in the study.
#' @param Delta survival censoring indicator.
#' @param A treatment. Make sure input treatment is binary - only target counterfactual is encoded as 1 (all other treatments are coded as 0).
#' @param a0 prespecified treatment level to examine.
#' @param fit_ftime time grid for which you want to estimate counterfactuals.
#' @param eval_ftime observed time points within the period (0, max of fit_ftime). Both event & censoring time points count.
#' @param gn_Ts n by t matrix of predictions of survival(events) for all individuals at each evaluated time point.
#' @param gn_Cs n by t matrix of predictions of survival(censor)  for all individuals at each evaluated time point.
#' @param gn_As n by t matrix of predictions of survival(censor)  for all individuals at each evaluated time point.
#' 
#'
#' @return A data.frame consists of Y, A, W1 and W2.
#'
#' @examples
#'
#' @export



# Y = ftime
# gn_Ts=pred_taw_ls[[1]]
# gn_Cs=pred_caw_ls[[1]]
# gn_A=gn_As[[1]]
# fit_ftime=event_ftime
# rm(Y,Delta, fit_ftime, gn_Ts, gn_Cs, gn_A)
getIF <- function(Y, 
                  Delta, 
                  A, 
                  a0,
                  fit_ftime, 
                  eval_ftime, 
                  gn_Ts, 
                  gn_Cs, 
                  gn_A){
  
  # sample size
  n <- length(Y)
  
  # sorted time points to be evaluated
  ord <- order(eval_ftime)
  eval_ftime <- eval_ftime[ord]
  # reorder predictions of T's, C's, and propensity scores accordingly
  gn_Ts <- gn_Ts[,ord]
  gn_Cs <- gn_Cs[,ord]
  
  # calculate integral of H for each individual i
  H_integral <- function(i) {
    # mark: why?
    vals <- diff(1/gn_Ts[i, ]) * 1/gn_Cs[i, -ncol(gn_Cs)]
    if (any(eval_ftime[-1] > Y[i])) 
      vals[eval_ftime[-1] > Y[i]] <- 0
    c(0, cumsum(vals))
  }
  # the integral value of H for each individual up to observed time points
  H_vals <- t(sapply(1:n, H_integral))
  
  # stepfun creates a function that can be applied to each observed Y_i
  # the intervals are closed on the left (and open on the right)
  gn_Ts_Y <- sapply(1:n, 
                    function(i){stepfun(eval_ftime, 
                                        c(1, gn_Ts[i, ]), 
                                        right = FALSE)(Y[i])})
  # the intervals are closed on the right (and open on the left)
  gn_Cs_Y <- sapply(1:n, 
                    function(i){stepfun(eval_ftime, 
                                        c(1, gn_Cs[i, ]), 
                                        right = TRUE)(Y[i])})
  
  # matrix for IF at each timepoint for each individual
  IFs <- matrix(NA, nrow = n, ncol = length(fit_ftime))
  IF_parts <- matrix(NA, nrow = n, ncol = length(fit_ftime))
  # IF_ted <- matrix(NA, nrow = n, ncol = length(fit_ftime))
  # target estimator
  surv <- rep(NA, length(fit_ftime))

  for (t0 in fit_ftime) {
    # the smallest time points closing but larger than the evaluated timepoint
    k <- min(which(eval_ftime >= t0))
    # estimated S for all individuals at t0 
    gn_Ts_t0 <- gn_Ts[, k]
    
    # H(t^y, A, W)
    part1 <- H_vals[, k]

    # I(Y<=t, Delta=1) / [ S_0(y|A,W) * G_0(y-|A,W)]
    part2 <- as.numeric(Y <= t0 & Delta == 1)/(gn_Ts_Y * gn_Cs_Y)
    # Influnce Function
    # here we recode target A as 1 outside of this function
    IF <- gn_Ts_t0 + gn_Ts_t0 * as.numeric(A == a0)/gn_A * (part1 - part2)
    IF_part <- gn_Ts_t0 * (part1 - part2)
    
    # locate certain time point to store the results
    t0_ind <- which(fit_ftime == t0)
    # point estimates for survival probablity at different time points
    surv[t0_ind] <- mean(IF)
    # sets of influence function
    # IF_ted[, t0_ind] <- IF - surv[t0_ind]
    IFs[, t0_ind] <- IF
    IF_parts[, t0_ind] <- IF_part
  }
  
  res <- list(times = fit_ftime, 
              # bound predictions of survival outcome between [0, 1]
              surv = pmin(1, pmax(0, surv)), 
              # surv = surv, 
              IFs = IFs,
              IF_parts = IF_parts
              # IF_ted = IF_ted
              )
  
  return(res)
}




