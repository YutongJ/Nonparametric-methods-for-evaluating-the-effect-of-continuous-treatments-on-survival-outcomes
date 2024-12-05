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


# --------------
#  library
# --------------
library(survSuperLearner)
library(CVXR)
library(mgcv, warn.conflicts = FALSE)

# --------------
#   Path
# --------------
dir <- "path/Codes/"
dir_fun <- paste0(dir, "R_func/")
dir_save <- paste0("path/Results/")


# --------------
#   Parameters
# --------------
params <- expand.grid(k = c(10, 20, 50, 100),
                      w = c(4),
                      lambda = c(4, 6))

(k <- params$k[ind])
(w <- params$w[ind])
(lambda <- params$lambda[ind])


# # --------------
# #   function
# # --------------
source(paste0(dir_fun, "EstPropScore.R"))
source(paste0(dir_fun, "getIF.R"))
source(paste0(dir_fun, "causalSS_SL.R"))
source(paste0(dir_fun, "OneStep_H.R"))
source(paste0(dir_fun, "weight_optimize.R"))
source(paste0(dir_fun, "FlatTest_Hs.R"))

# # read data
dat <- readRDS(file = paste0(dir_save, "dat_pooled_onehot.RDS"))



# if (w == 1){
#   Ws <- data.frame(dat[,c("geo")])
#   Ws_g <- data.frame(dat[,c(paste0("geo_",1:4))])
# }
# if (w == 2){
#   Ws <- data.frame(dat[,c("agegrp")])
#   Ws_g <- data.frame(dat[,c(paste0("agegrp_",1:4))])
# }
# if (w == 3){
#   Ws <- data.frame(dat[,c("weight")])
#   Ws_g <- data.frame(dat[,c("weight")])
# }
if (w == 4){
  Ws <- data.frame(dat[,c("weight","geo")])
  Ws_g <- data.frame(dat[,c("weight",
                            paste0("geo_",1:4))])
}
# if (w == 5){
#   Ws <- data.frame(dat[,c("weight","agegrp")])
#   Ws_g <- data.frame(dat[,c("weight",
#                             paste0("agegrp_",1:4))])
# }
# if (w == 6){
#   Ws <- data.frame(dat[,c("weight","geo","agegrp")])
#   Ws_g <- data.frame(dat[,c("weight",
#                             paste0("geo_",1:4),
#                             paste0("agegrp_",1:4))])
# }

# if (w == 4){
#   Ws <- data.frame(dat[,c("weight","geo","agegrp","risk_score")])
#   Ws_g <- data.frame(dat[,c(paste0("geo_",1:4))])
# }



# PS fit
# set.seed(241121)
# fit_aw <- EstPropScore(a = dat$D61,
#                        w = Ws_g,
#                        no.folds = 10,
#                        no.eval = 20)
# # generate prediction of E(A=a | W)
# gn_As <- fit_aw$cond.dens
# saveRDS(gn_As, file = paste0(dir_save, "res/gn_As.rds"))

# PS read
gn_As <- readRDS(file = paste0(dir_save, "res/gn_As.rds")) # adjust for weight, geo


set.seed(241121)
# proposed method
res <- causalSS(ftime = dat$time, 
                ftype = dat$status,
                A = dat$D61,
                W = Ws,
                pred_ftime = 595, # seq(1, 100, by=2),
                As = unique(sort(dat$D61)),
                continuous = TRUE,
                gn_As = gn_As,
                # gn_As = NULL,
                # event.SL.library = c("survSL.km", "survSL.coxph"),
                # cens.SL.library = c("survSL.km", "survSL.coxph"),
                event.SL.library = c("survSL.coxph", "survSL.gam", "survSL.rfsrc"),
                cens.SL.library = c("survSL.coxph", "survSL.gam", "survSL.rfsrc"),
                # obsWeights = NULL,
                obsWeights = dat$wt,
                # G_type = "sign_As",
                # G_type = "sign_knots",
                G_type = "indicator_knots",
                # knots = 10,
                knots = k,
                structure = "flexible",
                lambda = lambda,
                norm  = "L1",
                null_size = 1000,
                verbose = TRUE)

saveRDS(res, file = paste0(dir_save, "res/res_knots",k, "_lambda", lambda,".rds"))




# calculate weighted survival probability
survprob_sort <- res$pred_SS[[1]] %*% dat$wt
survprob_sort <- survprob_sort / sum(dat$wt)
survprob_wt <- survprob_sort[match(dat$D61, unique(sort(dat$D61)))]

# plotting counterfactual survival probability vs D61 dose
dat_plot <- data.frame(concentration = dat$D61,
                       survprob_wt = survprob_wt, # weighted
                       survprob = res$res_onesteps[[1]]$drf_est,# equally weighted
                       survprob_rowmean = rowMeans(res$pred_SS[[1]])[match(dat$D61, unique(sort(dat$D61)))],# equally weighted
                       p = res$res_tests[[1]]$p)
dat_plot <- dat_plot[order(dat_plot$concentration),]
rownames(dat_plot) <- NULL

saveRDS(dat_plot, file = paste0(dir_save, "res/dat_plot_knots",k, "_lambda", lambda,".rds"))

