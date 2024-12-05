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

# ii = 1

# --------------
#   additional parameters
# --------------
# sample sizes
n <- c(100, 300, 500, 1000)
# if (ii < 14){
#   n <- c(100, 300, 500, 1000)
# }else{
#   n <- c(1000)
# }


# simulations
B <- 1000
# time points that we want to predict, used for "pred_ftime"
times <- 25
# times <- seq(from = 10, to = 30, by = 15) # 10, 25

# type of results
pre_res <- c(
  "res_sign_As_knotsNA_structure_NA_lambda_NA_norm_NA", #1
  "res_sign_knots_knots10_structure_NA_lambda_NA_norm_NA",
  "res_sign_knots_knots20_structure_NA_lambda_NA_norm_NA",
  "res_sign_knots_knots50_structure_NA_lambda_NA_norm_NA",
  "res_indicator_knots_knots10_structure_flexible_lambda_4_norm_L1",#5
  "res_indicator_knots_knots20_structure_flexible_lambda_4_norm_L1",
  "res_indicator_knots_knots50_structure_flexible_lambda_4_norm_L1",
  "res_indicator_knots_knots10_structure_flexible_lambda_6_norm_L1",#8
  "res_indicator_knots_knots20_structure_flexible_lambda_6_norm_L1",
  "res_indicator_knots_knots50_structure_flexible_lambda_6_norm_L1",
  "res_indicator_knots_knots10_structure_flexible_lambda_Inf_norm_L1",#11
  "res_indicator_knots_knots20_structure_flexible_lambda_Inf_norm_L1",
  "res_indicator_knots_knots50_structure_flexible_lambda_Inf_norm_L1",
  "res_indicator_knots_knots10_structure_flexible_lambda_4_norm_L2",#14
  "res_indicator_knots_knots20_structure_flexible_lambda_4_norm_L2",
  "res_indicator_knots_knots50_structure_flexible_lambda_4_norm_L2",
  "res_indicator_knots_knots10_structure_flexible_lambda_6_norm_L2",#17
  "res_indicator_knots_knots20_structure_flexible_lambda_6_norm_L2",
  "res_indicator_knots_knots50_structure_flexible_lambda_6_norm_L2",
  "res_indicator_knots_knots10_structure_flexible_lambda_Inf_norm_L2",#20
  "res_indicator_knots_knots20_structure_flexible_lambda_Inf_norm_L2",
  "res_indicator_knots_knots50_structure_flexible_lambda_Inf_norm_L2",
  "res_indicator_knots_knots10_structure_monotone_lambda_NA_norm_L2",#23
  "res_indicator_knots_knots20_structure_monotone_lambda_NA_norm_L2",
  "res_indicator_knots_knots50_structure_monotone_lambda_NA_norm_L2"
  )


set <- pre_res[ii]
# --------------
#   Path
# --------------
dir <- path_read
dir_fun <- paste0(dir, "R_func/")
dir_data <- "path/dat_and_PS/"
setwd(dir_data)


for (j in seq_along(n)){
  # make a blank list to store results
  lout <- vector(mode = "list", length = B)
  
  # load 1000 simulation results
  for (b in seq(B)){
    if(file.exists(paste0("res/", set, "_n",n[j],"_seed", b, ".RData"))){
      lout[[b]] <- get(load(paste0("res/", set, "_n",n[j],"_seed", b, ".RData")))
    }
  }
  
  p_timepoints <- lapply(seq_along(times),
                         function(i){sapply(lout, function(x){ifelse(is.null(x), NA, x$res_tests[[i]]$p)}, simplify = TRUE)})
  # p_timepoints <- sapply(lout, function(x){ifelse(is.null(x), NA, x$res_tests[[1]]$p)}, simplify = TRUE)
  names(p_timepoints) <- paste0("t_", times)
  
  # lst_p[[j]] <- p_timepoints
  
  save(p_timepoints, file = paste0(dir_data,"pval/", set, "_n",n[j],"_lst_pval.RData"))
}
