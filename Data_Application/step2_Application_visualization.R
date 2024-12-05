rm(list = ls())

# --------------
#  library
# --------------
# library(hal9001)
library(ggplot2)
library(gridExtra)

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

dat <- readRDS(file = paste0(dir_save, "dat_pooled_onehot.RDS"))


# plot
plotting <- function(dat_plot, k){
  density_x <- density(dat_plot$concentration, 
                       adjust = 0.8,
                       weights = dat$wt / sum(dat$wt))
  par(mar = c(7, 4, 5, 4) + 0.3) 
  # plot(x = dat_plot$concentration,
  #      y = dat_plot$survprob,
  #      type = "n")
  plot(density_x, 
       axes = FALSE, 
       col = "yellow",
       xlab="", ylab = "",main=paste0("Knots=",k, " adjusting for weights & geographic"),
       ylim = c(0, 0.015))
  
  axis(side = 4, at = seq(0, 0.015, 0.003), las =1)       
  mtext("Density", side = 4, line = 3)
  
  polygon(density_x, col = rgb(1, 0.855, 0, alpha = 0.5), border = NA)
  
  
  # set parameter new=True for a new axis 
  par(new = TRUE)          
  plot(x = dat_plot$concentration,
       y = dat_plot$survprob_wt,
       type = "l",
       # pch=16,
       lwd = 2,
       ylim = c(0.8, 1),
       las = 1, 
       xlab = "D61 Concentration", 
       ylab = "Counterfactual survival probability by D595",
       main = "")
  # points(x = dat_plot$concentration,
  #        y = dat_plot$survprob,
  #        pch="*")
  # add p value from our hypothesis test
  pval <- round(dat_plot$p[1], digits=3)
  text(x = 250, y = 0.4, labels = paste0("p = ", pval), cex = 1.6)
  
}


# figure
pdf(file = paste0(dir_save, "Figs_wt_ind_knots_lambda46_adjust_for_weight_geo_survlib_cox_gam_rf.pdf"), width=40, height=20)
par(mfrow = c(2, 4))
for (ind in 1:8) {
  (k <- params$k[ind])
  (w <- params$w[ind])
  (lambda <- params$lambda[ind])
  
  # plotting counterfactual survival probability vs D61 dose
  dat_plot <- readRDS(file = paste0(dir_save, "res/dat_plot_knots",k, "_lambda", lambda,".rds"))
  
  plotting(dat_plot = dat_plot, k = k)
}
dev.off()






##### Including in the draft
# prepare table in the plot
res_vec <- rep(NA, length = nrow(params))
for(ind in 1:8){
  res <- readRDS(file = paste0(dir_save, "res/res_knots",params$k[ind], "_lambda", params$lambda[ind],".rds"))
  res_vec[ind] <- res$res_tests[[1]]$p
}
tab_pval <- cbind(c("Lambda", "N", "p values"),
                  rbind(as.character(params$lambda),
                  as.character(params$k),
                  sprintf("%.2f", res_vec)))


# plot
plotting_draft <- function(dat_plot, 
                           adjust = 0.75,
                           minval_y = 0.8, # minimum value on primary y-axis (for survival probability)
                           maxval_den = 0.015, # maximum value on secondy y-axis (for density)
                           table_data,
                           print_table = FALSE,
                           ...){
  # density for concentration
  density_x <- density(dat_plot$concentration, 
                       adjust = adjust,
                       weights = dat$wt / sum(dat$wt))
  
  # Prepare the density data for plotting
  density_data <- data.frame(
    x = density_x$x,
    y = density_x$y)
  
  # Create combined plot -- background:density, line:survival probability
  p <- ggplot() +
    # Layer 1: Density plot mapped to the secondary y-axis
    geom_polygon(
      data = density_data,
      aes(x = x, y = minval_y+y * ((1-minval_y)/maxval_den)), # Scale density to align with secondary y-axis
      fill = "orange",
      alpha = 0.5,
      color = NA
    ) +
    # Layer 2: Survival probability line plot mapped to the primary y-axis
    geom_line(
      data = dat_plot,
      aes(x = concentration, y = survprob_wt),
      color = "black",
      linewidth = 0.8
    ) +
    # Layer 3: Rug plot for concentration values
    geom_rug(
      data = dat_plot,
      aes(x = concentration),
      sides = "b", # Add rug only on the bottom (x-axis)
      color = "brown3", # Set color for rug marks
      alpha = 0.7
    ) +
    # Primary y-axis for survival probability
    scale_y_continuous(
      name = "Counterfactual survival probability by D595",
      limits = c(minval_y, 1)
    ) +
    # Secondary y-axis for density
    scale_y_continuous(
      name = "Counterfactual survival probability by D595",
      sec.axis = sec_axis(~ (.-minval_y) / ((1-minval_y)/maxval_den), name = "Density"), # Rescale secondary axis
      limits = c(minval_y, 1) # Match primary axis limits after scaling
    ) +
    # X-axis
    scale_x_continuous(
      name = "D61 Concentration"
    ) +
    # Themes and styling
    theme_bw() +
    theme(
      axis.title.y.left = element_text(color = "black"),
      axis.text.y.left = element_text(color = "black"),
      axis.title.y.right = element_text(color = "orange"),
      axis.text.y.right = element_text(color = "orange")
    )
  
  if (print_table){
    # Create plot with caption
    plot_top <- p + ggtitle("")
    
    # Convert the matrix into a tableGrob (table graphic object)
    table_grob <- tableGrob(table_data)
    # table_grob$widths <- unit(rep(1, ncol(table_data)), "null")
    
    # Arrange the plot and the table together
    grid.arrange(
      plot_top,
      table_grob,
      ncol = 1,
      heights = c(2, 1), # Adjust heights to give more space for the plot
      padding = unit(0.25, "lines")
    )
  }else{
    print(p)
  }
}


# plotting counterfactual survival probability vs D61 dose
dat_plot <- readRDS(file = paste0(dir_save, "res/dat_plot_knots",20, "_lambda", 4,".rds"))

# figure
p <- plotting_draft(dat_plot = dat_plot, table_data = tab_pval,
                    adjust = 0.7, minval_y = 0.8, maxval_den = 0.015)

ggsave(paste0(dir_save, "Draft_DA_Fig.svg"), plot = p, width = 8, height = 4)