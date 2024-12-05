#' Test flat null using one-step estimator
#'
#' @param ftime Survival time - time-to-event outcome.
#' @param ftype Censoring indicator - 1 for event.
#' @param A Treatment of interest, an n-dimensional vector.
#' @param W Covariates for adjustment, an n x p matrix.
#' @param tau Time point of interest to detect the difference among treatment levels.
#' @param pred_ftime A vector of time points we want to have estimation and inference on.
#' @param As A vector of treatments considered in the dataset. For discrete treatment, include each treatment value once. For continuous treatments, include all treatment values.
#' @param continuous Whether if the treatment is continuous. It affects the estimation of P(A|W).
#' @param gn_As The estimated densities. Default is NULL.
#' @param event.SL.library Library of candidate learners to use to estimate the conditional survival of the event.
#' @param cens.SL.library Library of candidate learners to use to estimate the conditional survival of the censoring.
#' @param obsWeights Optional `n x 1` vector of observation weights. Default is NULL.
#' @param G_type The type of h used for constructing H set. We have several options: sign, indicator, indicator_knots, basis.
#' @param knots If G_type is indicator_knots, need to specify this value. Default is 10 and NULL for others.
#' @param structure structure constraint: "flexible" or "monotone" omega's.
#' @param lambda scale parameter required for structure constraint on omega's.
#' @param norm scale constraint: L1-norm or L2-norm.
#' @param null_size The size of null sets for getting empirical p value.
#' @param verbose Include procedure information or not.
#'                  
#' @return A list containing the following: 
#'         \describe{
#'            \item{drf.est} Estimate of the dose response function, evaluated at a.
#'            \item{tmle} Estimate of expectation of drf times any given basis function.
#'            \item{eif} Estimate of efficient influence function for expectation of drf times any given basis function.
#'            \item{test.stat} Test statistic for flat null.
#'          }
#'
#' @export


# Code development -- variables of interest
# ftime = dat$Y
# ftype = dat$delta
# A = dat$A
# # W = data.frame(dat[,c("W1")])
# W = data.frame(dat[,c("W1", "W2")])
# # pred_ftime = NULL
# pred_ftime = times
# verbose = TRUE
# As = unique(sort(dat$A))
# continuous = TRUE
# G_type = "basis"
# rm(ftime, ftype, A, W, verbose, As, pred_ftime, continuous)
causalSS <- function(
  ftime, 
  ftype,
  A,
  W,  
  pred_ftime = NULL,
  As = c(0,1),
  continuous = FALSE,
  gn_As = NULL,
  event.SL.library,
  cens.SL.library,
  obsWeights = NULL,
  G_type = "indicator_knots",
  knots = 10,
  structure = "flexible",
  lambda = 2,
  norm  = "L1",
  null_size = 5000,
  verbose = TRUE){
  
  # number of times on file
  n <- length(ftime)
  
  
  # print error when all events are in the same censoring category
  if(any(table(ftype) == n)) {
    message("Error: all events have the same censoring type.")
  }
  
  # we may want to include a sanity check function to check the input dataset ahead of fitting
  # check.input <- function(){}
  
  
  # sorted unique event times
  if (is.null(pred_ftime)){
    pred_ftime <- sort(unique(ftime[ftime > 0 & ftime < max(ftime[ftype == 1])]))
  }
  
  # all evaluated sorted unique times (including event and censoring)
  eval_ftime <- sort(unique(c(0, ftime[ftime > 0 & ftime <= max(pred_ftime)], max(pred_ftime))))
  
  
  #----------------------
  # Estimate PS: A ~ W
  #----------------------
  if(verbose) message("Estimating propensity scores")
  
  if(is.null(gn_As)){
    if (continuous){
      fit_aw <- EstPropScore(a = A, w = W, no.eval = 10, no.folds = 5)
      
      # generate prediction of E(A=a | W)
      gn_As <- fit_aw$cond.dens
    }else{
      if (length(As) == 2){
        # logistic regression for binary
        fit_aw <- glm(A ~ ., family = "binomial", data = W)
        
        # generate prediction of E(A=a | W)
        gn_a1 <- stats::predict(fit_aw,
                                newdata = data.frame(A, W),
                                type = "response")
        # list of predictions
        gn_As <- list(gn_a0 = 1-gn_a1,
                      gn_a1 = gn_a1)
      }else{
        # multinomial regression for discrete treatments
        fit_aw <- nnet::multinom(A ~ ., data = W)
        
        # generate prediction of E(A=a | W)
        gn_ <- stats::predict(fit_aw, 
                              newdata = data.frame(A, W), 
                              type = "prob")
        # list of predictions
        gn_As <- asplit(gn_[,match(colnames(gn_), As)],
                        2)
      }
      names(gn_As) <- As
    }
  }
  
  
  
  
  # # using random forests for multiple level treatment A
  # fit_aw <- randomForest::randomForest(factor(A) ~ ., data=W, maxnodes=10,ntree=100)
  # # get P(A = xxx | W) 
  # gn_W <- predict(fit_aw, newdata = W, type = "prob")
  # # generate gn list
  # gn_As <- split(gn_W, col(gn_W))
  
  
  #------------------------------
  # Estimate conditional survival: T ~ A, W   &   C ~ A, W
  #------------------------------
  if(verbose) message("Estimating conditional survival function of the event/censoring time")
  
  # fit survSuperLearner model for conditional survival
  newAW <- cbind(A = rep(As, each = nrow(W)), W)
  # a=proc.time()
  suppressWarnings(
    fit_surv <- survSuperLearner::survSuperLearner(
      time = ftime,
      event = ftype,
      X = cbind(A, W),
      newX = newAW,
      new.times = eval_ftime,
      event.SL.library = event.SL.library,
      cens.SL.library = cens.SL.library,
      obsWeights = obsWeights,
      # event.SL.library = c("survSL.coxph"),
      # cens.SL.library = c("survSL.coxph"),
      verbose = verbose,
      # control = list(initWeightAlg = "survSL.rfsrc", verbose=FALSE), 
      cvControl = list(V = 2)))
  # proc.time() - a
  
  # predictions for event
  pred_taw_ls <- lapply(seq(1, n^2, by=n), function(i){fit_surv$event.SL.predict[i:(i + n - 1), ]})
  names(pred_taw_ls) <- As
  
  # predictions for censoring
  pred_caw_ls <- lapply(seq(1, n^2, by=n), function(i){fit_surv$cens.SL.predict[i:(i + n - 1), ]})
  names(pred_caw_ls) <- As
  
  
  
  #------------------------------
  # Estimate counterfactual survival
  #------------------------------
  if(verbose) message("Estimating counterfactual survival")
  # i : an index of treatments
  getIF_by_trt <- function(i){
    res <- getIF(Y = ftime,
                 Delta = ftype,
                 A = A,
                 a0 = As[i],
                 fit_ftime = pred_ftime, # time points which we want to estimate at
                 eval_ftime = eval_ftime, # need to estimated all time points
                 gn_Ts = pred_taw_ls[[i]],
                 gn_Cs = pred_caw_ls[[i]],
                 gn_A = gn_As[[i]])
    return(res)
  }
  
  # estimation by treatment
  res_IFs <- lapply(seq_along(As), getIF_by_trt)
  names(res_IFs) <- As
  
  

  #------------------------------
  # Hypothesis testing - reformating
  #------------------------------
  if(verbose) message("Conducting hypothesis testing")
  # calculate tatio of marginal density to conditional density -- I(A=a)/p(A|W)
  # ratio_ps <- diag(sapply(gn_As, mean, simplify = TRUE)) %*% (do.call(rbind, gn_As)^(-1))
  ratio_ps <- diag(rowMeans(gn_As)) %*% (gn_As^(-1))
  # ratio_ps <- t(fit_aw$prop.score)
  
  # survival counterfactual prediction for time-points of interest
  # 3 dimension: timepoints * (trt * Ws)
  pred_SS_ <- lapply(seq_along(eval_ftime), 
                 function(i){do.call(rbind, 
                                     lapply(pred_taw_ls, function(x){x[,i]}))})
  # subset timepoints we want to predict (specified in pred_ftime)
  pred_SS <- vector(mode = "list", length = length(pred_ftime))
  for (i in seq_along(pred_ftime)) {
    # the smallest time points closing but larger than the evaluated timepoint
    k <- min(which(eval_ftime >= pred_ftime[i]))
    # estimated S for all individuals at t0 
    pred_SS[[i]] <- pred_SS_[[k]]
  }
  
  IFs_SS <- lapply(seq_along(pred_ftime), 
                 function(i){do.call(rbind, 
                                     lapply(res_IFs, function(x){x$IF_parts[,i]}))})
  
  
  
  #-----------------------------------------------
  # H set: transformation of estimands, e.g. h(A)'s
  #-----------------------------------------------
  # we provide three types of transformation set
  # (1) For L1 + monotone, we can use sign function.
  #  <i>  sign at each observed treatment level
  if (G_type == "sign_As"){
    G_sets <- lapply(As, 
                     function(x){function(a){ifelse(a>= x, 1, -1)}})
  }
  #  <ii> sign at pre-specified knots that are evenly distributed over the treatment range
  if (G_type == "sign_knots"){
    G_sets <- lapply(seq(knots),
                     function(x){function(a){
                       (-1)^(findInterval(a,
                                          seq(min(A),max(A),length.out = knots+1),
                                          left.open=FALSE) >= x)}})
  }
  
  # (2) For all others, we can use indicator function.
  #  <i>  the separating knots are exactly the observed treatment levels
  if (G_type=="indicator_As"){
    G_sets <- lapply(As, 
                     function(x){function(a){as.numeric(a == x)}})
  }
  #  <ii> the separating knots are evenly distributed over the treatment range
  if (G_type=="indicator_knots"){
    G_sets <- lapply(seq(knots), 
                     function(x){function(a){
                       as.numeric(findInterval(a, 
                                               seq(min(A),max(A),length.out = knots+1),
                                               left.open=FALSE) == x)}})
  }
  
  
  
  #------------------------------
  # Hypothesis testing - for different time points of interest
  #------------------------------
  # test for each time point individually (or only for time point of interest)
  res_onesteps <- vector(mode = "list", length = length(pred_ftime))
  res_tests <- vector(mode = "list", length = length(pred_ftime))
  for (t0_index in seq_along(pred_ftime)){
    # One-step estimator and eif
    res_OneStep <- OneStep(ftime = ftime,
                           ftype = ftype, 
                           A = A, 
                           W = W, 
                           S = pred_SS[[t0_index]], 
                           int_vals = IFs_SS[[t0_index]],
                           ps_ratio = ratio_ps,
                           As = As,
                           G_func = G_sets,
                           G_type = G_type)
    
    # hypothesis testing
    if (G_type %in% c("sign_As", "sign_knots")){
      res_test <- FlatTest(res_OneStep = res_OneStep,
                           sign = TRUE,
                           structure = structure,
                           lambda = lambda,
                           norm  = norm,
                           null_size = null_size)
    }else{
      res_test <- FlatTest(res_OneStep = res_OneStep,
                           sign = FALSE,
                           structure = structure,
                           lambda = lambda,
                           norm  = norm,
                           null_size = null_size)
    }
    
    # store results
    res_onesteps[[t0_index]] <- res_OneStep
    res_tests[[t0_index]] <- res_test
    
  }
  
  out <- list(time = pred_ftime,
              pred_SS = pred_SS,
              res_onesteps = res_onesteps,
              res_tests = res_tests)
  
  return(out)
}

