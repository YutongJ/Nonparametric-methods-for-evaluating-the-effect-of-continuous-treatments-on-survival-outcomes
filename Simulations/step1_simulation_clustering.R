#!/bin/bash
options(echo=FALSE) # if you want to see commands in output file
args=(commandArgs(TRUE))
print(args)
arguments = matrix(unlist(strsplit(args,"=")),ncol=2,byrow = T)

for (args_i in 1:length(args)) {
  if(nchar(arguments[args_i,2])>5){
    assign(arguments[args_i,1],arguments[args_i,2])
  }
  else{
    assign(arguments[args_i,1],as.numeric(arguments[args_i,2]))
  }
}

# ind = 10001


# --------------
#   Parameters
# --------------
params1 <- expand.grid(sim = 1:1000,
                       n_samples = c(100, 300, 500, 1000),   # sample sizes
                       G_type = c("sign_As"),
                       knots = NA,
                       structure = NA,
                       lambda = NA,
                       norm  = NA,
                       null_size = 1000)

params2 <- expand.grid(sim = 1:1000,
                       n_samples = c(100, 300, 500, 1000),   # sample sizes
                       G_type = c("sign_knots"),
                       knots = c(10, 20, 50),
                       structure = NA,
                       lambda = NA,
                       norm  = NA,
                       null_size = 1000)

params3 <- expand.grid(sim = 1:1000,
                       n_samples = c(100, 300, 500, 1000),   # sample sizes
                       G_type = c("indicator_knots"),
                       knots = c(10, 20, 50),
                       structure = c("flexible"),
                       lambda = c(4,6,Inf),
                       norm  = c("L1"),
                       null_size = 1000)

params4 <- expand.grid(sim = 1:1000,
                       n_samples = c(100, 300, 500, 1000),   # sample sizes
                       G_type = c("indicator_knots"),
                       knots = c(10, 20, 50),
                       structure = c("flexible"),
                       lambda = c(4,6,Inf),
                       norm  = c("L2"),
                       null_size = 1000)

params5 <- expand.grid(sim = 1:1000,
                       n_samples = c(100, 300, 500, 1000),   # sample sizes
                       G_type = c("indicator_knots"),
                       knots = c(10, 20, 50),
                       structure = c("monotone"),
                       lambda = NA,
                       norm  = c("L2"),
                       null_size = 1000)

# paramss <- rbind(params1, params2, params3, params4, params5)
params <- rbind(params4[-c(1:24000),])

# specific parameter needed for simulation
(sim <- params$sim[ind])
(n <- params$n_samples[ind])
(g <- params$G_type[ind])
(k <- params$knots[ind])
(s <- params$structure[ind])
(l <- params$lambda[ind])
(m <- params$norm[ind])
(B <- params$null_size[ind])
  
# --------------
#   additional parameters
# --------------
# number of discrete treatment groups
# n_trt <- 5
# time points that we want to predict, used for "pred_ftime"
# times <- seq(from = 10, to = 30, by = 15) # 10, 25
times <- 25 # 10, 25


# --------------
#   library
# --------------
library(mgcv, warn.conflicts = FALSE)
library(hal9001)
library(CVXR, warn.conflicts=FALSE)
library(survSuperLearner, warn.conflicts=FALSE)



# rm(list = ls())
# --------------
#   Path
# --------------
dir <- path_read
dir_fun <- paste0(dir, "R_func/")
dir_save <- path_save
dir_data <- "path/dat_and_PS/"


# --------------
#   source
# --------------
# data generating process
# source(paste0(dir_fun, "DGP_complex_ceiling.R"))
# source(paste0(dir_fun, "EstPropScore.R"))
source(paste0(dir_fun, "getIF.R"))

source(paste0(dir_fun, "causalSS_SL.R"))
source(paste0(dir_fun, "OneStep_H.R"))
source(paste0(dir_fun, "weight_optimize.R"))
source(paste0(dir_fun, "FlatTest_Hs.R"))


# ----------------------
#   Generating Dataset
# ----------------------

# assign variables of interest
# ftime = dat$Y
# ftype = dat$delta
# A = dat$A
# # W = data.frame(dat[,c("W1")])
# W = data.frame(dat[,c("W1", "W2")])
# # pred_ftime = NULL
# pred_ftime = times
# As = unique(sort(dat$A))
# continuous = TRUE
# event.SL.library = c("survSL.km", "survSL.coxph","survSL.gam", "survSL.rfsrc")
# cens.SL.library = c("survSL.km", "survSL.coxph","survSL.gam", "survSL.rfsrc")
# # G_type = "indicator_knots"
# G_type = g
# knots = k
# structure = s
# lambda = l
# norm  = m
# null_size = B
# verbose = TRUE
# rm(ftime, ftype, A, W, pred_ftime,As, continuous,
#    event.SL.library, cens.SL.library,
#    G_type, knots, structure, lambda, norm, null_size,verbose)

one_sim <- function(i){
  # load data
  dat <- readRDS(file = paste0(dir_data, "dat_n_",n,"_seed",i,".rds"))
  gn_As <- readRDS(file = paste0(dir_data, "gn_As_n_",n,"_seed",i,".rds"))
  
  # debug(causalSS)
  res <- causalSS(ftime = dat$Y, 
                  ftype = dat$delta,
                  A = dat$A,
                  W = data.frame(dat[,c("W1", "W2")]),
                  pred_ftime = times, # seq(1, 100, by=2),
                  As = unique(sort(dat$A)),
                  continuous = TRUE,
                  gn_As = gn_As,
                  # gn_As = NULL,
                  # event.SL.library = c("survSL.coxph"),
                  # cens.SL.library = c("survSL.coxph"),
                  event.SL.library = c("survSL.km", "survSL.coxph","survSL.gam", "survSL.rfsrc"),
                  cens.SL.library = c("survSL.km", "survSL.coxph","survSL.gam", "survSL.rfsrc"),
                  G_type = g,
                  knots = k,
                  structure = s,
                  lambda = l,
                  norm  = m,
                  null_size = B,
                  verbose = FALSE)
  return(res)
}



# run single simulation
X = proc.time()
out <- one_sim(sim)
proc.time()-X
# save results
save(out, file = paste0(dir_save, "res/res_",g,"_knots",k,"_structure_",s,"_lambda_",l,"_norm_",m,"_n",n,"_seed", sim, ".RData"))
