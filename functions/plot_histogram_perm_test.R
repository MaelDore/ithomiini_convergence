### Function to plot histogram of null distribution during permutation test ####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

### Generate histrogram associated with a permutation test

### Inputs
 # List with obs stats [[1]], obs stats [[1]], mean null distribution [[3]], Q5 null distribution [[4]], Q95 null distribution [[5]], p-value left-sided [[6]], p-value right-sided [[7]]

### Outputs
 # Histogram of the test

### Function

plot_histogram_perm_test <- function(stat_summary, # List of obs stats [[1]],
                                                   # obs stats [[2]],
                                                   # mean null distribution [[3]],
                                                   # Q5 null distribution [[4]],
                                                   # Q95 null distribution [[5]],
                                                   # p-value left-sided [[6]],
                                                   # p-value right-sided [[7]]
                                     metric,       # Metric name
                                     title, cex_title = 1.3,
                                     cex_axis = 1.3, cex_lab = 1.4,
                                     cex_legend = 1.0,
                                     arrow_btm = 10, arrow_top = 60, arrow_adjust = 0, # Adjust arrow position
                                     panel_letter = "", cex_panel_letter = 1.6) # Facet letter
{
  
  # Extract info from test results
  stat_obs <- format(round(stat_summary[[1]], 3), nsmall = 3)
  stat_null <- format(round(stat_summary[[2]], 3), nsmall = 3)
  stat_mean <- format(round(stat_summary[[3]], 3), nsmall = 3)
  Q5 <- format(round(stat_summary[[4]], 3), nsmall = 3)
  Q95 <- format(round(stat_summary[[5]], 3), nsmall = 3)
  p_value_clustering <- stat_summary$p_value_clustering
  p_value_overdispersion <- stat_summary$p_value_overdispersion
  
  # Select the appropriate test
  if(p_value_clustering < p_value_overdispersion)
  {
    p_value <- p_value_clustering
  } else {
    p_value <- p_value_overdispersion
  }
  
  hist(x = c(as.numeric(stat_obs), as.numeric(stat_null)),
       breaks = 30, freq = TRUE, col = "gray", 
       main = title, 
       xlab = metric,
       cex.axis = cex_axis, cex.lab = cex_lab, cex.main = cex_title, lwd = 2)
  arrows(x0 = as.numeric(stat_obs) + arrow_adjust, y0 = arrow_top, x1 = as.numeric(stat_obs) + arrow_adjust, y1 = arrow_btm, length = 0.1, lwd = 3)
  
  # Draw control lines
  abline(v = as.numeric(stat_mean), lty = 2, lwd = 2)
  abline(v = as.numeric(Q5), lty = 2, lwd = 2, col = "red")
  abline(v = as.numeric(Q95), lty = 2, lwd = 2, col = "red")
  
  # Insert quantiles legend
  legend(legend = c(paste0("Mean ", metric, " = ",stat_mean), 
                    paste0("CI 5% = ",Q5), 
                    paste0("CI 95% = ",Q95)),
         x = "topright", inset = c(0.00, 0.00), lty = 2 , lwd = 2, col = c("black", "red"), cex = cex_legend, bty ="n")
  
  # Insert p-value legend
  legend(legend = c(paste0(metric, " obs = ",stat_obs),
                    paste0("p = ", p_value)),
         x = "topleft", inset = c(-0.01, 0.02), xjust = 1, # Use inset to manually adjust position
         cex = cex_legend, bty = "n", bg = "white", box.col = NA)
  
  # Add panel legend
  legend(legend = panel_letter, text.font = 2,
         x = "bottomright", inset = c(0.01, 0.08), xjust = 0.5,
         cex = cex_panel_letter, bty ="n")
  
}

