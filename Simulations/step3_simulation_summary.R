rm(list = ls())


# --------------
#  library
# --------------
library(ggplot2)


# --------------
#   additional parameters
# --------------
# sample sizes
n <- c(100, 300, 500, 1000)

# simulations
B <- 1000
# time points that we want to predict, used for "pred_ftime"
times <- 25
# times <- seq(from = 10, to = 30, by = 15) # 10, 25


# type of results
pre_res <- c("res_sign_As_knotsNA_structure_NA_lambda_NA_norm_NA", #1
             "res_sign_knots_knots10_structure_NA_lambda_NA_norm_NA",
             "res_sign_knots_knots20_structure_NA_lambda_NA_norm_NA",
             "res_sign_knots_knots50_structure_NA_lambda_NA_norm_NA",
             "res_indicator_knots_knots10_structure_flexible_lambda_4_norm_L1", #5
             "res_indicator_knots_knots20_structure_flexible_lambda_4_norm_L1",
             "res_indicator_knots_knots50_structure_flexible_lambda_4_norm_L1",
             "res_indicator_knots_knots10_structure_flexible_lambda_6_norm_L1", #8
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
             "res_indicator_knots_knots50_structure_monotone_lambda_NA_norm_L2")


# --------------
#   library
# --------------


# --------------
#   Path
# --------------
dir <- path_read
dir_fun <- paste0(dir, "R_func/")
dir_save <- path_save
dir_data <- "path/dat_and_PS/"


# #-------------------------------------------------------------------------------
#----------------------------
# simulation results check
#----------------------------
# dir_data <- paste0("X:/fast/hudson_a/user/yjin2/Causal_Survival/Results/20240724_DGPmono/")
# dir_data <- paste0("X:/fast/hudson_a/user/yjin2/Causal_Survival/Results/20240725_DGPcc/")
# 
# all_files <- list.files(paste0(dir_data, "res/"), full.names = FALSE)
# # test <- list.files(paste0(dir_data, "res/"), 
# #                    pattern = "res_indicator_knots_knots10_structure_flexible_lambda_4_norm_L1", 
# #                    full.names = FALSE)
# 
# 
# missing_seed1 <- function(j, n){
#   seed_pattern <- paste0(pre_res[j], "_n",n,"_seed(\\d+).RData")
#   seed_found <- grepl(seed_pattern,
#                       all_files)
#   # Extract seed numbers using sub
#   extracted_seeds <- sort(as.numeric(sapply(all_files[seed_found], function(x) {
#     sub("_seed", "", sub(seed_pattern, "\\1", x))
#   })))
#   
#   ref <- 1:1000
#   missing <- ref[! ref %in% extracted_seeds]
#   # return(paste(missing, collapse = ", "))
#   return(missing)
# }
# missing_seed2 <- function(j, n){
#   seed_pattern <- paste0(pre_res[j], "_n",n,"_seed(\\d+).RData")
#   seed_found <- grepl(seed_pattern,
#                       all_files)
#   # Extract seed numbers using sub
#   extracted_seeds <- sort(as.numeric(sapply(all_files[seed_found], function(x) {
#     sub("_seed", "", sub(seed_pattern, "\\1", x))
#   })))
#   
#   ref <- 1:1000
#   missing <- ref[! ref %in% extracted_seeds]
#   return(paste(missing, collapse = ", "))
#   # return(missing)
# }
# 
# # for (ii in 20:22) {
# for (ii in 1:25) {
#   print(pre_res[ii])
#   print(missing_seed2(ii, 100))
#   print(missing_seed2(ii, 300))
#   print(missing_seed2(ii, 500))
#   print(missing_seed2(ii, 1000))
#   # print(c(length(missing_seed1(ii, 100)),
#   #         length(missing_seed1(ii, 300)),
#   #         length(missing_seed1(ii, 500)),
#   #         length(missing_seed1(ii, 1000))))
# }
# 
# check file existing for later added simulations
# tmp <- params3[inds,]
# for (ii in 60:78) {
#   x <- tmp[ii,]
#   z <- file.exists(paste0(dir_data, 
#                           "res/res_",x$G_type,"_knots",x$knots,"_structure_",x$structure,
#                           "_lambda_",x$lambda,"_norm_",x$norm,"_n",x$n_samples,
#                           "_seed", x$sim, ".RData"))
#   # print(x);print(z)
#   print(c(ii,z))
# }
# #-------------------------------------------------------------------------------






#-------------------------------------------------------------------------------
#  sign_As_knotsNA_structure_NA_lambda_NA_norm_NA
#-------------------------------------------------------------------------------
# load all p-value results
lsp_p <- vector(mode = "list", length = length(pre_res))
for (i in seq(pre_res)){
  lsp_n <- list()
  lsp_n[[1]] <- get(load(paste0(dir_data, "pval/", pre_res[i], "_n",100,"_lst_pval.RData")))
  lsp_n[[2]] <- get(load(paste0(dir_data, "pval/", pre_res[i], "_n",300,"_lst_pval.RData")))
  lsp_n[[3]] <- get(load(paste0(dir_data, "pval/", pre_res[i], "_n",500,"_lst_pval.RData")))
  lsp_n[[4]] <- get(load(paste0(dir_data, "pval/", pre_res[i], "_n",1000,"_lst_pval.RData")))
  names(lsp_n) <- paste0("n_", n)
  
  lsp_p[[i]] <- lsp_n
}
names(lsp_p) <- pre_res

saveRDS(lsp_p, file = paste0(dir_save, "list_pval_DGPcc.RDS"))
# saveRDS(lsp_p, file = paste0(dir_save, "list_pval_DGPnull.RDS"))





# all p-values for time points t=10 & 25
typeIerr_summary <- function(lst_p){
  res_pvals <- lapply(seq_along(times), 
                      function(i){
                        sapply(lst_p, function(x){x[[i]]})
                      })
  res_typeIerr <- do.call(rbind, lapply(res_pvals,
                                        function(r){
                                          apply(r, 2, function(x){mean(x<0.05, na.rm = TRUE)})
                                        }))
  rownames(res_typeIerr) <- paste0("t_", times)
  return(res_typeIerr)
}

# apply to 25 scenarios
p_t <- lapply(lsp_p, typeIerr_summary)
# timepoint 1 (t=25)
vec_p_t1 <- round(do.call(c, lapply(p_t, function(x){x[1,]})), digits = 3)






#-------------------------------------
# draw typeI-error plot
#-------------------------------------
dat_plot <- rbind(
  # sign_As
  expand.grid(n_samples = c(100, 300, 500, 1000),   # sample sizes
              G_type = c("sign_As"),
              knots = "Max",
              structure = "monotone",
              lambda = "Monotone",
              norm  = "L1"),
  # sign_knots
  expand.grid(n_samples = c(100, 300, 500, 1000),   # sample sizes
              G_type = c("sign_knots"),
              knots = as.character(c(10, 20, 50)), #c(10, 20, 50),
              structure = "monotone",
              lambda = "Monotone",
              norm  = "L1"),
  # flexible L1
  expand.grid(n_samples = c(100, 300, 500, 1000),   # sample sizes
              G_type = c("indicator_knots"),
              knots = as.character(c(10, 20, 50)),
              lambda = as.character(c(4,6,Inf)),
              structure = c("flexible"),
              norm  = c("L1")),
  # flexible L2
  expand.grid(n_samples = c(100, 300, 500, 1000),   # sample sizes
              G_type = c("indicator_knots"),
              knots = as.character(c(10, 20, 50)),
              lambda = as.character(c(4,6,Inf)),
              structure = c("flexible"),
              norm  = c("L2")),
  # monotone L2
  expand.grid(n_samples = c(100, 300, 500, 1000),   # sample sizes
              G_type = c("indicator_knots"),
              knots = as.character(c(10, 20, 50)),
              lambda = c("Monotone"),
              structure = c("monotone"),
              norm  = c("L2"))
)

# long table with parameters and type I errors (t=10)
dat_plotx <- cbind(dat_plot, 
                   typeIerr = vec_p_t1)
dat_plotx$n_samples <- factor(dat_plotx$n_samples, levels = n)
saveRDS(dat_plotx, file = paste0(dir_save, "dat_plotx_DGPcc.RDS"))
# saveRDS(dat_plotx, file = paste0(dir_save, "dat_plotx_DGPnull.RDS"))

# plot multi-panel plots t=25
p1 <- ggplot(dat_plotx, aes(x = n_samples,
                            y = typeIerr,
                            color = knots,
                            group = knots)) + 
  theme(legend.position = "top") +
  labs(title = "Simulation summary for L1/L2-norm, t=25, DGP=Concave",
       x = "Sample size", 
       y = "Power") +
  geom_line(aes(color=knots), linewidth = 0.8) + 
  geom_point(aes(color=knots, shape=knots), size = 1.2) + 
  geom_hline(yintercept=0.95, color = "black", linetype=2) + 
  facet_grid(cols = vars(factor(lambda, 
                                levels = c("Monotone",4,6,Inf),
                                labels = c("Monotone", 
                                           expression(lambda == 4),
                                           expression(lambda == 6),
                                           expression(lambda == Inf)))),
             rows = vars(factor(norm)),
             scales = "free_y",
             labeller=labeller(treatment = labels))


# save the ggplot for the summary
pdf(file = paste0(dir_save, "summary_by_designs_cc_t25.pdf"), width=15, height=5)
print(p1)
dev.off()




# long table with parameters and type I errors (t=25)
dat_plotx <- cbind(dat_plot, 
                   typeIerr = vec_p_t2)
dat_plotx$n_samples <- factor(dat_plotx$n_samples, levels = n)

# L1
dat_L1 <- dat_plotx[1:(13*length(n)),]
# L2
dat_L2 <- dat_plotx[-c(1:(13*length(n))),]

# plot multi-panel plots t=25
p3 <- ggplot(dat_L1, aes(x = n_samples,
                      y = typeIerr,
                      color = knots,
                      group = knots)) + 
  theme(legend.position = "top") +
  labs(title = "Simulation summary for L1-norm, t=25",
       x = "Sample size", 
       y = "Type I error") +
  geom_line(aes(color=knots), linewidth = 0.8) + 
  geom_point(aes(color=knots, shape=knots), size = 1.2) + 
  geom_hline(yintercept=0.05, color = "black", linetype=2) + 
  facet_grid(cols = vars(factor(lambda, 
                                levels = c(2,4,6,Inf))),
             # rows = vars(factor(structure)),
             scales = "free_y",
             labeller=labeller(treatment = labels))

p4 <- ggplot(dat_L2, aes(x = n_samples,
                      y = typeIerr,
                      color = knots,
                      group = knots)) + 
  theme(legend.position = "top") +
  labs(title = "Simulation summary for L2-norm, t=25",
       x = "Sample size", 
       y = "Type I error") +
  geom_line(aes(color=knots), linewidth = 0.8) + 
  geom_point(aes(color=knots, shape=knots), size = 1.2) + 
  geom_hline(yintercept=0.05, color = "black", linetype=2) + 
  facet_grid(cols = vars(factor(lambda, 
                                levels = c(2,4,6,Inf))),
             rows = vars(factor(structure)),
             scales = "free_y",
             labeller=labeller(treatment = labels))

# save the ggplot for the summary
pdf(file = paste0(dir_save, "summary_by_designs_t25.pdf"), width=10, height=5)
print(p3)
print(p4)
dev.off()



#-------------------------------------
# QQ plot for observed p values 
#-------------------------------------
qqplot <- function (observed_pval, ni, ti) {
  plot(x = -log10(1:length(observed_pval)/length(observed_pval)), 
       y = -log10(sort(observed_pval)), 
       xlim = c(0, 3), ylim = c(0, 4),
       xlab = expression(-log10(rank(p)/1000)),
       ylab = expression(-log10(p)),
       main = paste0("n=", ni, ", t=",ti),
       pch = 20)
  abline(0, 1, col = "red")
}

# alternative way
# qqmath(~-log10(my.pvalues),
#        distribution=function(x){-log10(qunif(1-x))})

# # wald-type p values
# pdf(file = paste0(dir_save, "qqplot_by_n_and_t_waldp.pdf"), width = 4*length(n), height = 4*length(times))
# par(mfrow = c(length(times), length(n)))
# for (i in seq_along(times)) {
#   for (j in seq_along(n)) {
#     qqplot(res_all_by_time[[i]][[j]]$pval_wald, ni=n[j], ti=times[i])
#   }
# }
# dev.off()


# empirical p values
pdf(file = paste0(dir_save, "qqplot_by_n_and_t_sign_As.pdf"), width = 4*length(n), height = 4*length(times))
par(mfrow = c(length(times), length(n)))
for (l in 1:4){
  lst_tmp <- lsp_p[[l]]
  for (i in seq_along(times)) {
    for (j in seq_along(n)) {
      print(pre_res[l])
      qqplot(lst_tmp[[j]][[i]], ni=n[j], ti=times[i])
    }
  }
}
dev.off()









