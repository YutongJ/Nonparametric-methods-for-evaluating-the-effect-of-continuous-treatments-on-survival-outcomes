rm(list = ls())


# --------------
#  library
# --------------
library(ggplot2)
library(patchwork) # for arranging multiple plots



# --------------
#   Path
# --------------
dir <- path_read
dir_fun <- paste0(dir, "R_func/")
dir_save <- path_save


# read plotting data
dat_plotx_null <- readRDS(file = paste0(dir_save, "dat_plotx_DGPnull.RDS"))
dat_plotx_mono <- readRDS(file = paste0(dir_save, "dat_plotx_DGPmono.RDS"))
dat_plotx_cc <- readRDS(file = paste0(dir_save, "dat_plotx_DGPcc.RDS"))


# plot multi-panel plots t=25
p1 <- ggplot(dat_plotx_null, aes(x = n_samples,
                            y = typeIerr,
                            color = knots,
                            group = knots)) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20, angle=0),
        # legend.position = "none",
        legend.position = "bottom",
        legend.text = element_text(size = 20, angle=0),
        legend.title = element_text(size=20),
        title = element_text(size=20)) +
  labs(title = "(A)",
       x = "Sample size", 
       y = "Type I error") +
  geom_line(aes(color=knots), linewidth = 0.8) + 
  geom_point(aes(color=knots, shape=knots), size = 1.2) + 
  geom_hline(yintercept=0.05, color = "black", linetype=2) + 
  facet_grid(cols = vars(factor(lambda, 
                                levels = c("Monotone",4,6,Inf),
                                labels = c("Monotone", 
                                           expression(lambda == 4),
                                           expression(lambda == 6),
                                           expression(lambda == Inf)))),
             rows = vars(factor(norm)),
             scales = "free_y",
             labeller=labeller(treatment = labels))

# plot multi-panel plots t=25
p2 <- ggplot(dat_plotx_mono, aes(x = n_samples,
                            y = typeIerr,
                            color = knots,
                            group = knots)) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20, angle=0),
        # legend.position = "none",
        legend.position = "bottom",
        legend.text = element_text(size = 20, angle=0),
        legend.title = element_text(size=20),
        title = element_text(size=20)) +
  labs(title = "(B)",
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

# plot multi-panel plots t=25
p3 <- ggplot(dat_plotx_cc, aes(x = n_samples,
                            y = typeIerr,
                            color = knots,
                            group = knots)) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20, angle=0),
        legend.position = "bottom",
        legend.text = element_text(size = 20, angle=0),
        legend.title = element_text(size=20),
        title = element_text(size=20)) +
  labs(title = "(C)",
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
pdf(file = paste0(dir_save, "summary_by_designs_all_t25.pdf"), width=20, height=16)
p1 + p2 + p3 + plot_layout(ncol = 1)
dev.off()


