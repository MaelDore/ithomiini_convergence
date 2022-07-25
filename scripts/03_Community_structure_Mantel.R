##### Script 03: Mantel tests for pairwise Ist ~ Dclim #####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

# Compute pairwise Ist, geographic distances and standardized climatic distances among all communities
# Compute Mantel test for Ist ~ Dclim ; Ist ~ Dgeo ; Ist ~ Dclim + Dgeo
# Plot scatter plot for Ist ~ Dclim using residuals from Ist ~ Dgeo


### Input files

# Summary table of OMUs
# Stack of OMUs probabilities of presence
# Mimicry richness stack for all rings
# Stack of environmental variables

### Output files

# Pairwise Ist, geographic distance, and environmental distance matrices
# Mantel test for Ist ~ Dclim ; Ist ~ Dgeo ; Ist ~ Dclim + Dgeo
# Associated plot of null distri for each
# Scatter plot for Ist ~ Dclim using residuals from Ist ~ Dgeo


# Effacer l'environnement
rm(list = ls())

### 1/ Load directly the mimicry ring richness brick ####
# (each layer = expected local number of species for a mimicry ring)

library(raster)

# Load the mimicry ring richness stack
ring_richness_stack <- readRDS(file = paste0("./input_data/SDM_stacks/All_ring_rich_stack_Jaccard.80.rds"))
ring_richness_brick <- ring_richness_stack*1 # Transform into brick format to get a unique df for data
ring_richness_df <- ring_richness_brick[] # Extract unique df of communities x ring richness

# Load climate stack
env_stack <- readRDS(file = "./input_data/Select_env_15.rds")
climate_stack <- subset(env_stack, 1:4)
names(climate_stack) <- c("Tmean", "Tvar", "Htot", "Hvar")

### 2/ Functions to compute indices

#  Function to compute simpson index = probability to find species from a different ring when two species are drawn at random
simpson <- function (x, na.rm) {
  if (any(is.na(x))) { # Case with some NA values
    y <- NA
  }else{ # Case with all layers with data
    y <- x/sum(x)
    y <- 1-(sum(y*y))
  }
  return(y)
}

### 2/ Loop the analyses N times to show reproducibility ####

set.seed(seed = 54542) # Ensure reproducibility

library(vegan)

# Subsampling for X community to limit spatial autocorrelation, and save computation time
subsampling <- 1000 # Set number of communities to subsample
iterations <- 100 # Set number of iterations of the analyses

Mantel_summary_df <- matrix(nrow = iterations, ncol = 6)
Mantel_summary_df <- data.frame(Mantel_summary_df)
names(Mantel_summary_df) <- c("Dclim_stat", "Dclim_p_value", "Dgeo_stat", "Dgeo_p_value", "Dclim_geo_stat", "Dclim_geo_p_value")

for (k in 1:iterations) 
{
  ### 2.1/ Select subsampled communities ####
  
  sample.index <- sample(x = which(!is.na(ring_richness_stack[[1]]@data@values)), size = subsampling, replace = F)
  
  # Subsample the richness data
  sampled_ring_richness <- ring_richness_df[sample.index,]
  
  ### 2.2/ Compute pairwise_Ist ####
  
  pairwise_Ist <- matrix(ncol = subsampling, nrow = subsampling) # Generate squared matrix to store parwise Ist
  
  for (i in 1:subsampling) # For each community in row
  {
    for (j in 1:subsampling) # For each community in column
    {
      com1 <- sampled_ring_richness[i,] # Get mimicry richnesses of all rings for community 1
      com2 <- sampled_ring_richness[j,] # Get mimicry richnesses of all rings for community 2
      D1 <- simpson(com1) # Compute local simpson index for community 1
      D2 <- simpson(com2) # Compute local simpson index for community 2
      Ds <- mean(c(D1,D2)) # Compute mean simpson index among the pair
      full_com <- apply(X = sampled_ring_richness[c(i,j),], MARGIN = 2, FUN = sum) # Get mimcry richnesses of all rings for the two communities combined
      Dt <- simpson(full_com) # Compute Simpson for both communities merged
      pairwise_Ist[i,j] <- 1-(Ds/Dt) # Compute pairwise Ist
    }
    
    # Show i every 10 iterations
    if (i%%100==0) {cat(paste0("\n",i,"/",subsampling))}
  }
  
  ### 2.3/ Compute pairwise geographical distances ####
  coord <- coordinates(ring_richness_stack[[1]])[sample.index,] # Extract spatial coordinates for subsampled communities
  # plot(coord) # To get an idea of the sampling locations
  pairwise_geodist = geosphere::distm(x = coord)/1000 # Geographic distance on the WGS84 ellipsoid, in km
  
  ### 2.4/ Compute pairwise euclidian climatic distances on standardized climatic variables ####
  sampled_com_env <- climate_stack@data@values[sample.index,] # Extract climatic variables for subsampled communities
  sampled_com_env <- scale(x = sampled_com_env, center = T, scale = T) # Standardize the variables
  pairwise_climdist <- as.matrix(dist(x = sampled_com_env, method = "euclidian")) # Compute euclidian distances
  
  # Save all pairwise distances for this iteration
  save(pairwise_Ist, pairwise_climdist, pairwise_geodist, file = paste0("./outputs/Community_structure/Mantel_tests/Iterations/Pairwise_values_",k,".RData"))
  
  ### 2.5/ Compute Mantel tests ####
  
  # library("vegan")
  
  # load(file = paste0("./outputs/Community_structure/Mantel_tests/Iterations/Pairwise_values_",k,".RData"))
  
  # # Check for normality
  # hist(pairwise_climdist)
  # hist(pairwise_geodist)
  # hist(pairwise_Ist)
  
  # No normality => Spearman rank's test
  
  # 2.5.1/ Ist ~ clim avec rho de Spearman
  Mantel_Spearman_Ist_clim <- mantel(xdis = pairwise_Ist, ydis = pairwise_climdist, method = "spearman", na.rm = T)
  Mantel_Spearman_Ist_clim
  save(Mantel_Spearman_Ist_clim, file = paste0("./outputs/Community_structure/Mantel_tests/Iterations/Resultats_Mantel_Ist_clim_Spearman_",k,".RData"))
  
  cat("Ist ~ clim - Ok")
  
  # 2.5.2/ Ist ~ geo avec rho de Spearman
  Mantel_Spearman_Ist_geo <- mantel(xdis = pairwise_Ist, ydis = pairwise_geodist, method = "spearman", na.rm = T)
  Mantel_Spearman_Ist_geo
  save(Mantel_Spearman_Ist_geo, file = paste0("./outputs/Community_structure/Mantel_tests/Iterations/Resultats_Mantel_Ist_geo_Spearman_",k,".RData"))
  
  cat("Ist ~ geo - Ok")
  
  # 2.5.3/ Mantel Partial Test : Ist ~ Dclim + cov(Dgeo)
  partial_Mantel_Spearman_Ist_clim_geo <- mantel.partial(xdis = pairwise_Ist, ydis = pairwise_climdist, zdis = pairwise_geodist, method = "spearman", na.rm = T)
  partial_Mantel_Spearman_Ist_clim_geo
  save(partial_Mantel_Spearman_Ist_clim_geo, file = paste0("./outputs/Community_structure/Mantel_tests/Iterations/Resultats_Mantel_partiel_Ist_Spearman_",k,".RData"))
  
  cat("Ist ~ clim + cov(geo) - Ok")
  
  ### 2.6/ Save statistics in summary df ####
  
  Mantel_summary_df$Dclim_stat[k] <- round(Mantel_Spearman_Ist_clim$statistic, 3)
  Mantel_summary_df$Dclim_p_value[k] <- round(Mantel_Spearman_Ist_clim$signif, 3)
  Mantel_summary_df$Dgeo_stat[k] <- round(Mantel_Spearman_Ist_geo$statistic, 3)
  Mantel_summary_df$Dgeo_p_value[k] <- round(Mantel_Spearman_Ist_geo$signif, 3)
  Mantel_summary_df$Dclim_geo_stat[k] <- round(partial_Mantel_Spearman_Ist_clim_geo$statistic, 3)
  Mantel_summary_df$Dclim_geo_p_value[k] <- round(partial_Mantel_Spearman_Ist_clim_geo$signif, 3)
  
  save(Mantel_summary_df, file = paste0("./outputs/Community_structure/Mantel_tests/Mantel_summary_df.RData"))
  
  cat(paste0("\n", Sys.time()," ------ Process over for ",k,"/",iterations," iterations ------\n"))
}


##### 3/ Plot Mantel test results for iteration k ####

# Chose the iteration to use as example
k <- 1

# Copy in the main folder files from iteration k to use for the main text plots
file.copy(from = paste0("./outputs/Community_structure/Mantel_tests/Iterations/Resultats_Mantel_Ist_clim_Spearman_",k,".RData"), to = paste0("./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_clim_Spearman.RData"), overwrite = T)
file.copy(from = paste0("./outputs/Community_structure/Mantel_tests/Iterations/Resultats_Mantel_Ist_geo_Spearman_",k,".RData"), to = paste0("./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_geo_Spearman.RData"), overwrite = T)
file.copy(from = paste0("./outputs/Community_structure/Mantel_tests/Iterations/Resultats_Mantel_partiel_Ist_Spearman_",k,".RData"), to = paste0("./outputs/Community_structure/Mantel_tests/Resultats_Mantel_partiel_Ist_Spearman.RData"), overwrite = T)

file.copy(from = paste0("./outputs/Community_structure/Mantel_tests/Iterations/Pairwise_values_",k,".RData"), to = paste0("./outputs/Community_structure/Mantel_tests/Pairwise_values.RData"), overwrite = T)

### 3.1 Plotting functions ####

# Function to plot scatterplot of pairwise distances
plot_pairwise_distances <- function(title, cex_title = 1.3,
                                    y, x, y_lab, x_lab,
                                    y_lim = NULL,
                                    cex_axis = 1.5, cex_lab = 1.6, cex_legend = 1.4,
                                    rho_value, p_value, regression,
                                    panel_letter = "", cex_panel_letter = 2.0)
{
  plot(x = x,y = y, 
       ylim = y_lim,
       main = title, ylab = y_lab, xlab = x_lab,
       cex.axis = cex_axis, cex.lab = cex_lab, cex.main = cex_title, lwd = 1, type = "n", axes = F)
  
  points(x = x, y = y, pch = 16, col = "#00000070")
  
  # Insert blank
  legend(legend = c("              ",
                    "              "),
         x = "topleft", inset = c(0, 0.02), xjust = 1, y.intersp = 1.1, # Use inset to manually adjust position
         cex = cex_legend, bty = "o", bg = "white", box.col = NA)
  
  # Insert rho value legend
  legend(legend = bquote(bold(rho) ~ bold('=') ~ bold(.(rho_value))), text.font = 2,
         x = "topleft", inset = c(-0.01, 0.02), xjust = 1, # Use inset to manually adjust position
         cex = cex_legend, bty = "n", bg = "white", box.col = NA)
  
  # Insert p-value legend
  legend(legend = c("",paste0("p = ",p_value)), text.font = 2,
         x = "topleft", inset = c(-0.01, 0.02), xjust = 1, # Use inset to manually adjust position
         cex = cex_legend, bty = "n", bg = "white", box.col = NA)
  
  axis(side = 1, lwd = 2, cex.axis = cex_axis)
  axis(side = 2, lwd = 2, cex.axis = cex_axis)
  abline(regression, lwd = 3, col = "red")
  
  # Add panel legend
  legend(legend = panel_letter, text.font = 2,
         x = "bottomright", inset = c(-0.03, 0.03), xjust = 0.5,
         cex = cex_panel_letter, bty ="n")
  
}


### Function to plot histogram of null distribution of Spearman's rho during permutation test
plot_histogram_Mantel_test <- function(test, title, cex_title = 1.3,
                                       cex_axis = 1.5, cex_lab = 1.6, cex_legend = 1.4,
                                       arrow_btm, arrow_top, arrow_adjust,
                                       rho_value, p_value,
                                       panel_letter = "", cex_panel_letter = 2.0)
{
  hist(c(test$perm, test$statistic),
       breaks = 30, freq = TRUE, col = "gray", 
       main = title, 
       xlab = expression(paste("Spearman's ", rho)),
       cex.axis = cex_axis, cex.lab = cex_lab, cex.main = cex_title, lwd = 2)
  arrows(x0 = test$statistic + arrow_adjust, y0 = arrow_top, x1 = test$statistic + arrow_adjust, y1 = arrow_btm, length = 0.1, lwd = 3)
  abline(v = mean(c(test$perm, test$statistic)), lty = 2, lwd = 2)
  abline(v = quantile(c(test$perm, test$statistic), 0.95), lty = 2, lwd = 2, col = "red")
  
  # Insert quantiles legend
  legend(legend = c(paste("Mean = 0.000"), 
                    paste0("CI 95% = ", format(round(quantile(c(test$perm, test$statistic), 0.95),3), nsmall = 3))), 
         x = "topright", inset = c(0.02, 0.20), y.intersp = 1.2, lty = 2 , lwd = 2, col = c("black", "red"), cex = cex_legend, bty ="n")
  
  # Insert rho value legend
  legend(legend = bquote(bold(rho) ~ bold('obs =') ~ bold(.(rho_value))), text.font = 2,
         x = "bottomright", inset = c(0.02, 0.41), xjust = 1, # Use inset to manually adjust position
         cex = cex_legend, bty = "n", bg = "white", box.col = NA)
  
  # Insert p-value legend
  legend(legend = c(paste0("p = ",p_value)), text.font = 2,
         x = "bottomright", inset = c(0.02, 0.35), xjust = 1, # Use inset to manually adjust position
         cex = cex_legend, bty = "n", bg = "white", box.col = NA)
  
  # Add panel legend
  legend(legend = panel_letter, text.font = 2,
         x = "topright", inset = c(-0.03, 0.00), xjust = 0.5,
         cex = cex_panel_letter, bty ="n")
  
}


### 3.2/ Load and prepare data ####

load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_clim_Spearman.RData")
load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_geo_Spearman.RData")
load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_partiel_Ist_Spearman.RData")

load(file = "./outputs/Community_structure/Mantel_tests/Pairwise_values.RData")


# Extract vectors of distances, without replicates and NA
pairwise_Ist_vec <- as.dist(pairwise_Ist)[which(!is.na(as.dist(pairwise_Ist)))]
pairwise_climdist_vec <- as.dist(pairwise_climdist)[which(!is.na(as.dist(pairwise_Ist)))]
pairwise_geodist_vec <- as.dist(pairwise_geodist)[which(!is.na(as.dist(pairwise_Ist)))]

# Extract only 1000 pairwise distances for the plot to make it readable
set.seed(seed = 55644) # Ensure reproductibility
sample_indices <- sample(x= seq_along(pairwise_Ist_vec), size = 1000, replace = F)
save(sample_indices, file = "./outputs/Community_structure/Mantel_tests/sample_indices.RData")


### 3.3/ Plot Ist ~ Dclim ####

# Extract stats for legend
rho_value <- format(round(Mantel_Spearman_Ist_clim$statistic, 3), nsmall = 3)
p_value <- format(Mantel_Spearman_Ist_clim$signif, nsmall = 3)

# Compute a linear regression to get coefficients to draw a predict line
reg_clim <- lm(pairwise_Ist_vec ~ pairwise_climdist_vec)
summary(reg_clim)

# Plot
pdf(file = paste0("./graphs/Community_Structure/Ist_Dclim.pdf"), height = 6, width = 8)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

plot_pairwise_distances(title = "Climate structures \n mimicry community composition", cex_title = 1.3, 
                        y = pairwise_Ist_vec[sample_indices], x = pairwise_climdist_vec[sample_indices], 
                        y_lab = bquote('Pairwise' ~I[ST]), x_lab = "Standardized climatic distance", 
                        cex_axis = 1.3, cex_lab = 1.4, 
                        rho_value = rho_value, p_value = p_value, regression = reg_clim)

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()


### 3.4/ Plot Ist ~ Dgeo ####


# Extract stats for legend
rho_value <- format(round(Mantel_Spearman_Ist_geo$statistic, 3), nsmall = 3)
p_value <- format(Mantel_Spearman_Ist_geo$signif, nsmall = 3)

# Compute a linear regression to get coefficients to draw a predict line
reg_geo <- lm(pairwise_Ist_vec ~ pairwise_geodist_vec)
summary(reg_geo)

# Plot
pdf(file = paste0("./graphs/Community_Structure/Ist_Dgeo.pdf"), height = 6, width = 8)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

plot_pairwise_distances(title = "Space structures \n mimicry community composition", cex_title = 1.3, 
                        y = pairwise_Ist_vec[sample_indices], x = pairwise_geodist_vec[sample_indices], 
                        y_lab = bquote('Pairwise' ~I[ST]), x_lab = "Geographic distance [km]", 
                        cex_axis = 1.3, cex_lab = 1.4, 
                        rho_value = rho_value, p_value = p_value, regression = reg_geo)

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()


# 3.5/ Plot Ist ~ Dclim + covar(Dgeo) ####

# Extract stats for legend
rho_value <- format(round(partial_Mantel_Spearman_Ist_clim_geo$statistic, 3), nsmall = 3)
p_value <- format(partial_Mantel_Spearman_Ist_clim_geo$signif, nsmall = 3)

# Compute linear regression of residuals from Ist ~ Dgeo on Dclim
Ist_resid <- residuals(reg_geo, type = "response") # Extract residuals from LM on geographic distances
reg_clim_geo <- lm(Ist_resid ~ pairwise_climdist_vec)
summary(reg_clim_geo)

# Plot
pdf(file = paste0("./graphs/Community_Structure/Ist_Dclim_geo.pdf"), height = 6, width = 8)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

plot_pairwise_distances(title = "Climate structures \n mimicry community composition \n even when space is taken into account", cex_title = 1.3, 
                        y = Ist_resid[sample_indices], x = pairwise_climdist_vec[sample_indices], 
                        y_lab = bquote('Residual Pairwise' ~I[ST]), x_lab = "Standardized climatic distance", 
                        y_lim = c(-0.2, 0.5),
                        cex_axis = 1.3, cex_lab = 1.4, 
                        rho_value = rho_value, p_value = p_value, regression = reg_clim_geo)

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()


### 3.6/ Plot null distri for partial Mantel test Ist ~ Dclim + covar(Dgeo) ####

# Extract stats for legend
rho_value <- format(round(partial_Mantel_Spearman_Ist_clim_geo$statistic, 3), nsmall = 3)
p_value <- format(partial_Mantel_Spearman_Ist_clim_geo$signif, nsmall = 3)

pdf(file = paste0("./graphs/Community_Structure/Ist_Dclim_geo_null_distri.pdf"), height = 6, width = 8)

original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))

plot_histogram_Mantel_test(test = partial_Mantel_Spearman_Ist_clim_geo,
                           title = "Partial Mantel test statistic distribution\n under 999 permutations",
                           cex_title = 1.5, cex_axis = 1.3, cex_lab = 1.4,
                           arrow_btm = 5, arrow_top = 70, arrow_adjust = 0.000,
                           rho_value = rho_value, p_value = p_value)

par(mar = original_int_margins)

dev.off()


##### 4/ Final plot with all 4 sub-plots ####

### Load and prepare data

load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_clim_Spearman.RData")
load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_geo_Spearman.RData")
load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_partiel_Ist_Spearman.RData")

load(file = "./outputs/Community_structure/Mantel_tests/Pairwise_values.RData")


# Extract vectors of distances, without replicates and NA
pairwise_Ist_vec <- as.dist(pairwise_Ist)[which(!is.na(as.dist(pairwise_Ist)))]
pairwise_climdist_vec <- as.dist(pairwise_climdist)[which(!is.na(as.dist(pairwise_Ist)))]
pairwise_geodist_vec <- as.dist(pairwise_geodist)[which(!is.na(as.dist(pairwise_Ist)))]

# Extract only 1000 pairwise distances for the plot to make it readable
set.seed(seed = 55644) # Ensure reproductibility
sample_indices <- sample(x= seq_along(pairwise_Ist_vec), size = 1000, replace = F)
save(sample_indices, file = "./outputs/Community_structure/Mantel_tests/sample_indices.RData")


### Plot

pdf(file = paste0("./graphs/Community_structure/Mantel_all_plots_2.pdf"), height = 10.3, width = 10)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5,5,1.5,1.5), mfrow = c(2,2))

# Panel A = Ist ~ Dclim


# Extract stats for legend
rho_value <- format(round(Mantel_Spearman_Ist_clim$statistic, 3), nsmall = 3)
p_value <- format(Mantel_Spearman_Ist_clim$signif, nsmall = 3)

# Compute a linear regression to get coefficients to draw a predict line
reg_clim <- lm(pairwise_Ist_vec ~ pairwise_climdist_vec)
summary(reg_clim)

plot_pairwise_distances(title = NULL, cex_title = 1.3, 
                        y = pairwise_Ist_vec[sample_indices], x = pairwise_climdist_vec[sample_indices], 
                        y_lab = bquote('Pairwise' ~I[ST]), x_lab = "Standardized climatic distance", 
                        cex_axis = 1.5, cex_lab = 1.6, 
                        rho_value = rho_value, p_value = p_value, regression = reg_clim,
                        panel_letter = "A")

# Panel B = Ist ~ Dgeo

# Extract stats for legend
rho_value <- format(round(Mantel_Spearman_Ist_geo$statistic, 3), nsmall = 3)
p_value <- format(Mantel_Spearman_Ist_geo$signif, nsmall = 3)

# Compute a linear regression to get coefficients to draw a predict line
reg_geo <- lm(pairwise_Ist_vec ~ pairwise_geodist_vec)
summary(reg_geo)

plot_pairwise_distances(title = NULL, cex_title = 1.3, 
                        y = pairwise_Ist_vec[sample_indices], x = pairwise_geodist_vec[sample_indices], 
                        y_lab = bquote('Pairwise' ~I[ST]), x_lab = "Geographic distance [km]", 
                        cex_axis = 1.5, cex_lab = 1.6, 
                        rho_value = rho_value, p_value = p_value, regression = reg_geo,
                        panel_letter = "B")

# Panel C = Ist ~ Dclim + covar(Dgeo)

# Extract stats for legend
rho_value <- format(round(partial_Mantel_Spearman_Ist_clim_geo$statistic, 3), nsmall = 3)
p_value <- format(partial_Mantel_Spearman_Ist_clim_geo$signif, nsmall = 3)

# Compute linear regression of residuals from Ist ~ Dgeo on Dclim
Ist_resid <- residuals(reg_geo, type = "response") # Extract residuals from LM on geographic distances
reg_clim_geo <- lm(Ist_resid ~ pairwise_climdist_vec)
summary(reg_clim_geo)


plot_pairwise_distances(title = NULL, cex_title = 1.3, 
                        y = Ist_resid[sample_indices], x = pairwise_climdist_vec[sample_indices], 
                        y_lab = bquote('Residual Pairwise' ~I[ST]), x_lab = "Standardized climatic distance", 
                        y_lim = c(-0.2, 0.5),
                        cex_axis = 1.5, cex_lab = 1.6, 
                        rho_value = rho_value, p_value = p_value, regression = reg_clim_geo,
                        panel_letter = "C")

# Panel D = NUll distri for Ist ~ Dclim + covar(Dgeo)

# Extract stats for legend
rho_value <- format(round(partial_Mantel_Spearman_Ist_clim_geo$statistic, 3), nsmall = 3)
p_value <- format(partial_Mantel_Spearman_Ist_clim_geo$signif, nsmall = 3)

plot_histogram_Mantel_test(test = partial_Mantel_Spearman_Ist_clim_geo,
                           title = NULL, cex_title = 1.5, 
                           cex_axis = 1.5, cex_lab = 1.6,
                           arrow_btm = 5, arrow_top = 70, arrow_adjust = 0.000,
                           rho_value = rho_value, p_value = p_value,
                           panel_letter = "D")

par(oma = original_ext_margins, mar = original_int_margins, mfrow = c(1,1))        

dev.off()


##### 5/ Replicates analysis #####

### Retrieve results from other replicates and compute summary stats

load(file = paste0("./outputs/Community_structure/Mantel_tests/Mantel_summary_df.RData"))

library(dplyr)

rho_summary <- Mantel_summary_df %>%
  summarize(
    Dclim_mean = round(mean(Dclim_stat), 3),
    Dclim_sd = round(sd(Dclim_stat), 3),
    Dclim_CV = round(sd(Dclim_stat)/mean(Dclim_stat)*100, 1),
    Dclim_2.5 = round(quantile(Dclim_stat, probs = 0.025), 3),
    Dclim_97.5 = round(quantile(Dclim_stat, probs = 0.975), 3),
    Dclim_min = round(min(Dclim_stat), 3),
    Dclim_max = round(max(Dclim_stat), 3),
    Dgeo_mean = round(mean(Dgeo_stat), 3),
    Dgeo_sd = round(sd(Dclim_stat), 3),
    Dgeo_CV = round(sd(Dgeo_stat)/mean(Dgeo_stat)*100, 1),
    Dgeo_2.5 = round(quantile(Dgeo_stat, probs = 0.025), 3),
    Dgeo_97.5 = round(quantile(Dgeo_stat, probs = 0.975), 3),
    Dgeo_min = round(min(Dgeo_stat), 3),
    Dgeo_max = round(max(Dgeo_stat), 3),
    Dclim_geo_mean = round(mean(Dclim_geo_stat), 3),
    Dclim_geo_sd = round(sd(Dclim_geo_stat), 3),
    Dclim_geo_CV = round(sd(Dclim_geo_stat)/mean(Dclim_geo_stat)*100, 1),
    Dclim_geo_2.5 = round(quantile(Dclim_geo_stat, probs = 0.025), 3),
    Dclim_geo_97.5 = round(quantile(Dclim_geo_stat, probs = 0.975), 3),
    Dclim_geo_min = round(min(Dclim_geo_stat), 3),
    Dclim_geo_max = round(max(Dclim_geo_stat), 3))

p_summary <- Mantel_summary_df %>%
  summarize(
    Dclim_mean = round(mean(Dclim_p_value), 3),
    Dclim_sd = round(sd(Dclim_p_value), 3),
    Dclim_CV = round(sd(Dclim_p_value)/mean(Dclim_p_value)*100, 1),
    Dclim_2.5 = round(quantile(Dclim_p_value, probs = 0.025), 3),
    Dclim_97.5 = round(quantile(Dclim_p_value, probs = 0.975), 3),
    Dclim_min = round(min(Dclim_p_value), 3),
    Dclim_max = round(max(Dclim_p_value), 3),
    Dgeo_mean = round(mean(Dgeo_p_value), 3),
    Dgeo_sd = round(sd(Dclim_p_value), 3),
    Dgeo_CV = round(sd(Dgeo_p_value)/mean(Dgeo_p_value)*100, 1),
    Dgeo_2.5 = round(quantile(Dgeo_p_value, probs = 0.025), 3),
    Dgeo_97.5 = round(quantile(Dgeo_p_value, probs = 0.975), 3),
    Dgeo_min = round(min(Dgeo_p_value), 3),
    Dgeo_max = round(max(Dgeo_p_value), 3),
    Dclim_geo_mean = round(mean(Dclim_geo_p_value), 3),
    Dclim_geo_sd = round(sd(Dclim_geo_p_value), 3),
    Dclim_geo_CV = round(sd(Dclim_geo_p_value)/mean(Dclim_geo_p_value)*100, 1),
    Dclim_geo_2.5 = round(quantile(Dclim_geo_p_value, probs = 0.025), 3),
    Dclim_geo_97.5 = round(quantile(Dclim_geo_p_value, probs = 0.975), 3),
    Dclim_geo_min = round(min(Dclim_geo_p_value), 3),
    Dclim_geo_max = round(max(Dclim_geo_p_value), 3))

all_summary <- rbind(rho_summary, p_summary)
row.names(all_summary) <- c("Spearman_rho", "p_values")

write.csv2(x = all_summary, file = "./tables/Mantel_tests_summary.csv")




############# Old plotting script #################

# 
# # Extract vectors of data without NA
# pairwise_Ist_vec <- as.dist(pairwise_Ist)[which(!is.na(as.dist(pairwise_Ist)))]
# pairwise_climdist_vec <- as.dist(pairwise_climdist)[which(!is.na(as.dist(pairwise_Ist)))]
# pairwise_geodist_vec <- as.dist(pairwise_geodist)[which(!is.na(as.dist(pairwise_Ist)))]
# 
# ### Select data to subsample for the plot (499 500 points it's a bit too much)
# # Set seed for reproductibility
# set.seed(23018)
# sample_indices <- sample(x= 1:length(pairwise_Ist_vec), size = 1000, replace = F)
# save(sample_indices, file = "./outputs/Community_structure/Mantel_tests/sample_indices.RData")
# 
# ### 3.1/ Ist ~ Dclim ##
# 
# # Compute a linear regression to get coefficients to draw a predict line
# reg.clim <- lm(pairwise_Ist_vec ~ pairwise_climdist_vec)
# summary(reg.clim)
# 
# round(Mantel_Spearman_Ist_clim$statistic, 3) # Need to write manually the value on the plot # rho = 0.351
# 
# # Plot
# pdf(file = "./graphs/Community_Structure/Ist_Dclim.pdf", height = 6, width = 8)
# 
# original_ext_margins <- par()$oma
# original_int_margins <- par()$mar
# par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))
# 
# plot(pairwise_climdist_vec[sample_indices], pairwise_Ist_vec[sample_indices], 
#      main = "Climate structures \n mimicry community composition", ylab = bquote('Pairwise' ~I[ST]), xlab = "Standardized climatic distance",
#      cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.3, lwd = 1, type = "n", axes = F)
# 
# points(pairwise_climdist_vec[sample_indices], pairwise_Ist_vec[sample_indices], pch = 16, col = "#00000070")
# 
# legend(legend = c("                   ",
#                   "                   "),
#        x = "topleft", inset = c(0, 0.02), xjust = 1, y.intersp = 1.1, # Use inset to manually adjust position
#        cex = 1.2, bty ="o", bg = "white", box.col = NA)
# 
# 
# # legend(legend = c((paste0("\u03C1 = ", # Unicode for small case rho
# #                           round(Mantel_Spearman_Ist_clim$statistic, 3))), 
# #                   paste0("p < 0.001")), 
# #        x = "topleft", inset = c(-0.01, 0.02), xjust = 1, # Use inset to manually adjust position
# #        cex = 1.2, bty ="o", bg = "white", box.col = NA)
# 
# # legend(legend = c(bquote(rho == .(round(Mantel_Spearman_Ist_clim$statistic, 3))),
# #                   paste0("p < 0.001")),
# #        x = "topleft", inset = c(-0.01, 0.02), xjust = 1, # Use inset to manually adjust position
# #        cex = 1.2, bty ="n", bg = "white", box.col = NA)
# 
# legend(legend = c(expression(paste(rho, " = 0.351")),
#                   paste0("p = 0.001")),
#        x = "topleft", inset = c(-0.01, 0.02), xjust = 1, # Use inset to manually adjust position
#        cex = 1.2, bty ="n", bg = "white", box.col = NA)
# 
# axis(side = 1, lwd = 2, cex.axis = 1.3)
# axis(side = 2, lwd = 2, cex.axis = 1.3)
# abline(reg.clim, lwd = 3, col = "red")
# 
# par(oma = original_ext_margins, mar = original_int_margins)
# 
# dev.off()
# 
# ### 3.2/ Ist ~ Dgeo ##
# 
# # Compute a linear regression to get coefficients to draw a predict line
# reg.geo <- lm(pairwise_Ist_vec~pairwise_geodist_vec)
# summary(reg.geo)
# 
# round(Mantel_Spearman_Ist_geo$statistic, 3) # Need to write manually the value on the plot # rho = 0.517
# 
# pdf(file = "./graphs/Community_Structure/Ist_Dgeo.pdf", height = 6, width = 8)
# 
# original_ext_margins <- par()$oma
# original_int_margins <- par()$mar
# par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))
# 
# plot(pairwise_geodist_vec[sample_indices], pairwise_Ist_vec[sample_indices],
#      main = "Space structures \n mimicry community composition", ylab = bquote('Pairwise' ~I[ST]), xlab = "Geographic distance [km]",
#      cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.3, lwd = 1, type = "n", axes = F)
# 
# points(pairwise_geodist_vec[sample_indices], pairwise_Ist_vec[sample_indices], pch = 16, col = "#00000070")
# 
# legend(legend = c("                   ",
#                   "                   "),
#        x = "topleft", inset = c(0, 0.02), xjust = 1, y.intersp = 1.113, # Use inset to manually adjust position
#        cex = 1.2, bty ="o", bg = "white", box.col = NA)
# 
# legend(legend = c(expression(paste(rho, " = 0.517")),
#                   paste0("p = 0.001")),
#        x = "topleft", inset = c(-0.01, 0.02), xjust = 1, # Use inset to manually adjust position
#        cex = 1.2, bty ="n", bg = "white", box.col = NA)
# 
# axis(side = 1, lwd = 2, cex.axis = 1.3)
# axis(side = 2, lwd = 2, cex.axis = 1.3)
# abline(reg.geo, lwd = 3, col = "red")
# 
# par(oma = original_ext_margins, mar = original_int_margins)
# 
# dev.off()
# 
# # 3.3/ Ist ~ Dclim + covar(Dgeo) ##
# 
# Ist_resid <- residuals(reg.geo, type = "response") # Extract residuals from LM on geographic distances
# reg.clim.geo <- lm(Ist_resid ~ pairwise_climdist_vec)
# summary(reg.clim.geo)
# 
# round(partial_Mantel_Spearman_Ist_clim_geo$statistic, 3) # Need to wirte manually the value on the plot # rho = 0.195
# 
# # Plot
# pdf(file = "./graphs/Community_Structure/Ist_Dclim_geo.pdf", height = 6, width = 8)
# 
# original_ext_margins <- par()$oma
# original_int_margins <- par()$mar
# par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))
# 
# plot(pairwise_climdist_vec[sample_indices], Ist_resid[sample_indices], 
#      main = "Climate structures \n mimicry community composition \n even when space is taking into account", ylab = bquote('Residual Pairwise' ~I[ST]), xlab = "Standardized climatic distance",
#      cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.3, lwd = 1, type = "n", axes = F, 
#      ylim = c(-0.2, 0.5))
# 
# points(pairwise_climdist_vec[sample_indices], Ist_resid[sample_indices], pch = 16, col = "#00000070")
# 
# legend(legend = c("                   ",
#                   "                   "),
#        x = "topleft", inset = c(0, 0.02), xjust = 1, y.intersp = 0.9, # Use inset to manually adjust position
#        cex = 1.2, bty ="o", bg = "white", box.col = NA)
# 
# legend(legend = c(expression(paste(rho, " = 0.195")),
#                   paste0("p = 0.001")),
#        x = "topleft", inset = c(-0.01, 0.02), xjust = 1, # Use inset to manually adjust position
#        cex = 1.2, bty ="n", bg = "white", box.col = NA)
# 
# axis(side = 1, lwd = 2, cex.axis = 1.3)
# axis(side = 2, lwd = 2, cex.axis = 1.3)
# abline(reg.clim.geo, lwd = 3, col = "red")
# 
# par(oma = original_ext_margins, mar = original_int_margins)
# 
# dev.off()
# 
# 
# ### 3.4/ Plot null distri for partial Mantel test Ist ~ Dclim + covar(Dgeo) ##
# load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_partiel_Ist_Spearman.RData")
# 
# pdf(file = "./graphs/Community_Structure/Ist_Dclim_geo_null_distri.pdf", height = 6, width = 8)
# 
# original_int_margins <- par()$mar
# par(mar = c(5.1,5,4.1,2.1))
# 
# 
# hist(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic),
#      breaks = 30, freq = TRUE, col = "gray", 
#      main = "Partial Mantel test statistic distribution\n under 999 permutations", 
#      xlab = expression(paste("Spearman's ", rho)),
#      cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5, lwd = 2)
# arrows(partial_Mantel_Spearman_Ist_clim_geo$statistic + 0.000, 80, partial_Mantel_Spearman_Ist_clim_geo$statistic + 0.000, 10, length = 0.1, lwd = 2)
# abline(v = mean(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic)), lty = 2, lwd = 2)
# abline(v = quantile(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic), 0.95), lty = 2, lwd = 2, col = "red")
# 
# legend(legend = c(expression(paste("Mean ",rho, " obs = 0")), 
#                   paste0("CI 95% = ", round(quantile(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic), 0.95),3))), 
#        x = "topright", inset = c(0.02, 0.00), lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")
# 
# legend(legend = c(expression(paste(rho, " obs = 0.195")),
#                   paste0("       p = 0.001")),
#        x = "bottomright", inset = c(0.02, 0.40), xjust = 1, # Use inset to manually adjust position
#        cex = 1.2, bty ="n", bg = "white", box.col = NA)
# 
# # legend(legend = c(paste0("Rho obs = ", round(partial_Mantel_Spearman_Ist_clim_geo$statistic, 3)),
# #                   paste0("                  p = ", 1-(ecdf(x = c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic))(partial_Mantel_Spearman_Ist_clim_geo$statistic)))),
# #        x = "right", cex = 1.2, bty ="n", xjust = 1)
# 
# par(mar = original_int_margins)
# dev.off()
# 
# 
# ##### 4/ Final plot with all 4 sub-plots ##
# 
# load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_clim_Spearman.RData")
# load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_geo_Spearman.RData")
# load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_partiel_Ist_Spearman.RData")
# 
# load(file = "./outputs/Community_structure/Mantel_tests/Pairwise_values.RData")
# 
# # Extract vectors of data without NA
# pairwise_Ist_vec <- as.dist(pairwise_Ist)[which(!is.na(as.dist(pairwise_Ist)))]
# pairwise_climdist_vec <- as.dist(pairwise_climdist)[which(!is.na(as.dist(pairwise_Ist)))]
# pairwise_geodist_vec <- as.dist(pairwise_geodist)[which(!is.na(as.dist(pairwise_Ist)))]
# 
# load(file = "./outputs/Community_structure/Mantel_tests/sample_indices.RData")
# 
# pdf(file = "./graphs/Community_Structure/Mantel_all_plots.pdf", height = 10.3, width = 10)
# 
# original_ext_margins <- par()$oma
# original_int_margins <- par()$mar
# par(oma = c(0,0,0,0), mar = c(5,5,1,1), mfrow = c(2,2))
# 
# # Panel A = Ist ~ Dclim
# 
# # Compute a linear regression to get coefficients to draw a predict line
# reg.clim <- lm(pairwise_Ist_vec ~ pairwise_climdist_vec)
# # summary(reg.clim)
# 
# plot(pairwise_climdist_vec[sample_indices], pairwise_Ist_vec[sample_indices], 
#      # main = "Climate structures \n mimicry community composition",
#      main = NA, 
#      ylab = bquote('Pairwise' ~I[ST]), xlab = "Standardized climatic distance",
#      cex.axis = 1.5, cex.lab = 1.7, cex.main = 1.3, lwd = 1, type = "n", axes = F)
# 
# points(pairwise_climdist_vec[sample_indices], pairwise_Ist_vec[sample_indices], pch = 16, col = "#00000070")
# 
# legend(legend = c("                   ",
#                   "                   "),
#        x = "topleft", inset = c(0, 0.02), xjust = 1, y.intersp = 1.1, # Use inset to manually adjust position
#        cex = 1.2, bty ="o", bg = "white", box.col = NA)
# 
# legend(legend = c(as.expression(bquote(bold(paste(rho, " = 0.351")))),
#                   as.expression(bquote(bold(paste("p = 0.001"))))),
#        x = "topleft", inset = c(-0.03, 0.02), xjust = 1, # Use inset to manually adjust position
#        cex = 1.5, text.font = 2,
#        bty ="n", bg = "white", box.col = NA)
# 
# axis(side = 1, lwd = 2, cex.axis = 1.5)
# axis(side = 2, lwd = 2, cex.axis = 1.5)
# abline(reg.clim, lwd = 3, col = "red")
# 
# legend(legend = as.expression(bquote(bold("A"))), 
#        x = "bottomright", inset = c(0.02, 0.01), xjust = 0.5,
#        cex = 1.6, bty ="n")
# 
# # Panel B = Ist ~ Dgeo
# 
# # Compute a linear regression to get coefficients to draw a predict line
# reg.geo <- lm(pairwise_Ist_vec~pairwise_geodist_vec)
# # summary(reg.geo)
# 
# plot(pairwise_geodist_vec[sample_indices], pairwise_Ist_vec[sample_indices],
#      # main = "Space structures \n mimicry community composition",
#      main = NA, 
#      ylab = bquote('Pairwise' ~I[ST]), xlab = "Geographic distance [km]",
#      cex.axis = 1.5, cex.lab = 1.7, cex.main = 1.3, lwd = 1, type = "n", axes = F)
# 
# points(pairwise_geodist_vec[sample_indices], pairwise_Ist_vec[sample_indices], pch = 16, col = "#00000070")
# 
# legend(legend = c("                   ",
#                   "                   "),
#        x = "topleft", inset = c(0, 0.02), xjust = 1, y.intersp = 1.113, # Use inset to manually adjust position
#        cex = 1.2, bty ="o", bg = "white", box.col = NA)
# 
# legend(legend = c(as.expression(bquote(bold(paste(rho, " = 0.517")))),
#                   as.expression(bquote(bold(paste("p = 0.001"))))),
#        x = "topleft", inset = c(-0.03, 0.02), xjust = 1, # Use inset to manually adjust position
#        cex = 1.5, bty ="n", bg = "white", box.col = NA)
# 
# axis(side = 1, lwd = 2, cex.axis = 1.5)
# axis(side = 2, lwd = 2, cex.axis = 1.5)
# abline(reg.geo, lwd = 3, col = "red")
# 
# legend(legend = as.expression(bquote(bold("B"))), 
#        x = "bottomright", inset = c(0.02, 0.01), xjust = 0.5,
#        cex = 1.6, bty ="n")
# 
# # Panel C = Ist ~ Dclim + covar(Dgeo)
# 
# Ist_resid <- residuals(reg.geo, type = "response") # Extract residuals from LM on geographic distances
# reg.clim.geo <- lm(Ist_resid ~ pairwise_climdist_vec)
# # summary(reg.clim.geo)
# 
# plot(pairwise_climdist_vec[sample_indices], Ist_resid[sample_indices], 
#      # main = "Climate structures \n mimicry community composition \n even when space is taking into account",
#      main = NA, 
#      ylab = bquote('Residual Pairwise' ~I[ST]), xlab = "Standardized climatic distance",
#      cex.axis = 1.5, cex.lab = 1.7, cex.main = 1.3, lwd = 1, type = "n", axes = F, 
#      ylim = c(-0.2, 0.5))
# 
# points(pairwise_climdist_vec[sample_indices], Ist_resid[sample_indices], pch = 16, col = "#00000070")
# 
# legend(legend = c("                   ",
#                   "                   "),
#        x = "topleft", inset = c(0, 0.02), xjust = 1, y.intersp = 0.9, # Use inset to manually adjust position
#        cex = 1.2, bty ="o", bg = "white", box.col = NA)
# 
# legend(legend = c(as.expression(bquote(bold(paste(rho, " = 0.195")))),
#                   as.expression(bquote(bold(paste("p = 0.001"))))),
#        x = "topleft", inset = c(-0.03, 0.02), xjust = 1, # Use inset to manually adjust position
#        cex = 1.5, bty ="n", bg = "white", box.col = NA)
# 
# axis(side = 1, lwd = 2, cex.axis = 1.5)
# axis(side = 2, lwd = 2, cex.axis = 1.5)
# abline(reg.clim.geo, lwd = 3, col = "red")
# 
# text(x = -0.59, y = -0.1, labels = "-0.1", cex = 1.5, srt = 90, xpd = T)
# 
# legend(legend = as.expression(bquote(bold("C"))), 
#        x = "topright", inset = c(0.02, 0.001), xjust = 0.5,
#        cex = 1.6, bty ="n")
# 
# # Panel D = NUll distri for Ist ~ Dclim + covar(Dgeo)
# 
# # Extract stats for legend
# rho_value <- format(round(partial_Mantel_Spearman_Ist_clim_geo$statistic, 3), nsmall = 3)
# p_value <- format(partial_Mantel_Spearman_Ist_clim_geo$signif, nsmall = 3)
# 
# 
# 
# hist(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic),
#      breaks = seq(-0.05, 0.20, 0.05/4),
#      freq = TRUE, col = "gray", 
#      # main = "partial Mantel test statistic distribution\n under 999 permutations", 
#      main = NA, 
#      xlab = expression(paste("Spearman's ", rho)),
#      cex.axis = 1.5, cex.lab = 1.7, cex.main = 1.5, lwd = 2)
# arrows(partial_Mantel_Spearman_Ist_clim_geo$statistic + 0.000, 100, partial_Mantel_Spearman_Ist_clim_geo$statistic + 0.000, 10, length = 0.1, lwd = 2)
# abline(v = mean(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic)), lty = 2, lwd = 2)
# abline(v = quantile(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic), 0.95), lty = 2, lwd = 2, col = "red")
# 
# legend(legend = c(expression(paste("Mean = 0.001")), 
#                   paste0("CI 95% = ", round(quantile(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic), 0.95),3))), 
#        x = "topright", inset = c(0.02, 0.17), lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.5, bty ="n")
# 
# legend(legend = c(as.expression(bquote(bold(paste(rho, " obs = 0.195")))),
#                   as.expression(bquote(bold(paste("       p = 0.001"))))),
#        x = "bottomright", inset = c(0.02, 0.40), xjust = 1, # Use inset to manually adjust position
#        cex = 1.5, bty ="n", bg = "white", box.col = NA)
# 
# legend(legend = as.expression(bquote(bold("D"))), 
#        x = "topright", inset = c(0.02, 0.001), xjust = 0.5,
#        cex = 1.6, bty ="n")
# 
# par(oma = original_ext_margins, mar = original_int_margins, mfrow = c(1,1))        
# 
# dev.off()
# 
# 
