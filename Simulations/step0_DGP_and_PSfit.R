#!/bin/bash
options(echo=TRUE) # if you want to see commands in output file
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

# ind = 1


# --------------
#   Parameters
# --------------
params <- expand.grid(sim = 1:1000,
                      n_samples = c(100, 300, 500, 1000))

# specific parameter needed for simulation
(sim <- params$sim[ind])
(n <- params$n_samples[ind])


# --------------
#   library
# --------------
library(hal9001)



# cluster
dir <- "path/Codes/"
dir_fun <- paste0(dir, "R_func/")
dir_save <- "path/Results/"



# --------------
#   source
# --------------
source(paste0(dir_fun, "DGP_complex_ceiling.R"))
source(paste0(dir_fun, "EstPropScore.R"))



# simple DGP
# null
dat <- data_generate(seed = sim, n = n, contA = TRUE, effect = "null")
# monotone
dat <- data_generate(seed = sim, n = n, contA = TRUE, effect = "monotone")
# concave
dat <- data_generate(seed = sim, n = n, contA = TRUE, effect = "concave")
saveRDS(dat, file = paste0(dir_save, "dat_n_",n,"_seed",sim,".rds"))

# PS fit
fit_aw <- EstPropScore(a = dat$A, w = data.frame(W1=dat$W1,W2=dat$W2),
                       no.folds = 10,
                       no.eval = 10)
# generate prediction of E(A=a | W)
gn_As <- fit_aw$cond.dens
saveRDS(gn_As, file = paste0(dir_save, "gn_As_n_",n,"_seed",sim,".rds"))


