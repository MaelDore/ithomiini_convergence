##### Mantel tests for pairwise Ist ~ Dclim #####

# Compute pairwise Ist, geographic distances and standardized climatic distances among all communities
# Compute Mantel test for Ist ~ Dclim ; Ist ~ Dgeo ; Ist ~ Dclim + Dgeo
# Plot scatter plot for Ist ~ Dclim using residuals from Ist ~ Dgeo


### Input files

# Summary table of OMUs
# Stack of OMUs probabilities of presence
# Mimicry richness stack for all rings
# Stack of environmental variables

### Output files

# Mantel test for Ist ~ Dclim ; Ist ~ Dgeo ; Ist ~ Dclim + Dgeo
# Associated plot of null distri for each
# Scatter plot for Ist ~ Dclim using residuals from Ist ~ Dgeo


# Effacer l'environnement
rm(list = ls())

library(raster)

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

### 1/ Load directly the mimicry ring richness brick ####
# (each layer = expected local number of species for a mimicry ring)

# Load the mimicry ring richness stack
ring_richness_stack <- readRDS(file = paste0("./input_data/SDM_stacks/All_ring_rich_stack_Jaccard.80.rds"))
ring_richness_brick <- ring_richness_stack*1 # Transform into brick format to get a unique df for data
ring_richness_df <- ring_richness_brick[] # Extract unique df of communities x ring richness

# Load climate stack
env_stack <- readRDS(file = "./input_data/Select_env_15.rds")
climate_stack <- subset(env_stack, 1:4)
names(climate_stack) <- c("Tmean", "Tvar", "Htot", "Hvar")

### 2/ Select subsampled communities ####

# Subsampling for N community to limit spatial autocorrelation, and save computation time
set.seed(seed = 54542) # Ensure reproductibility
subsampling <- 1000
sample.index <- sample(x = which(!is.na(ring_richness_stack[[1]]@data@values)), size = subsampling, replace = F)

# Subsample the richness data
sampled_ring_richness <- ring_richness_df[sample.index,]

### 3/ Compute pairwise_Ist ####

pairwise_Ist <- matrix(ncol = subsampling, nrow = subsampling) # Generate squared matrix to store parwise Ist

for (i in 1:subsampling) # For each community in row
{
  for (j in 1:subsampling) # For each community in column
  {
    com1 <- sampled_ring_richness[i,] # Get mimcry richnesses of all rings for community 1
    com2 <- sampled_ring_richness[j,] # Get mimcry richnesses of all rings for community 2
    D1 <- simpson(com1) # Compute local simpson index for community 1
    D2 <- simpson(com2) # Compute local simpson index for community 2
    Ds <- mean(c(D1,D2)) # Compute mean simpson index among the pair
    full_com <- apply(X = sampled_ring_richness[c(i,j),], MARGIN = 2, FUN = sum) # Get mimcry richnesses of all rings for the two communities combined
    Dt <- simpson(full_com) # Compute Simpson for both communities merged
    pairwise_Ist[i,j] <- 1-(Ds/Dt) # Compute pairwise Ist
  }
  
  # Show i every 10 iterations
  if (i%%10==0) {print(i)}
}

### 4/ Compute pairwise geographical distances ####
coord <- coordinates(ring_richness_stack[[1]])[sample.index,] # Extract spatial coordinates for subsampled communities
# plot(coord) # To get an idea of the sampling locations
pairwise_geodist = geosphere::distm(x = coord)/1000 # Geographic distance on the WGS84 ellipsoid, in km

### 5/ Compute pairwise euclidian climatic distances on standardized climatic variables ####
sampled_com_env <- climate_stack@data@values[sample.index,] # Extract climatic variables for subsampled communities
sampled_com_env <- scale(x = sampled_com_env, center = T, scale = T) # Standardize the variables
pairwise_climdist <- as.matrix(dist(x = sampled_com_env, method = "euclidian")) # Compute euclidian distances

# Save all pairwise distances
save(pairwise_Ist, pairwise_climdist, pairwise_geodist, file = "./outputs/Community_structure/Mantel_tests/Pairwise_values.RData")


### 6/ Compute Mantel tests ####

library("vegan")

load(file = "./outputs/Community_structure/Mantel_tests/Pairwise_values.RData")

# Check for normality
hist(pairwise_climdist)
hist(pairwise_geodist)
hist(pairwise_Ist)

# No normality => Spearman rank's test

# 6.1/ Ist ~ clim avec rho de Spearman
Mantel_Spearman_Ist_clim <- mantel(xdis = pairwise_Ist, ydis = pairwise_climdist, method = "spearman", na.rm = T)
Mantel_Spearman_Ist_clim
save(Mantel_Spearman_Ist_clim , file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_clim_Spearman.RData")

# 6.2/ Ist ~ geo avec rho de Spearman
Mantel_Spearman_Ist_geo <- mantel(xdis = pairwise_Ist, ydis = pairwise_geodist, method = "spearman", na.rm = T)
Mantel_Spearman_Ist_geo
save(Mantel_Spearman_Ist_geo, file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_geo_Spearman.RData")

# 6.3/ Mantel Partial Test : Ist ~ Dclim + cov(Dgeo)

partial_Mantel_Spearman_Ist_clim_geo <- mantel.partial(xdis = pairwise_Ist, ydis = pairwise_climdist, zdis = pairwise_geodist, method = "spearman", na.rm = T)
partial_Mantel_Spearman_Ist_clim_geo
save(partial_Mantel_Spearman_Ist_clim_geo, file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_partiel_Ist_Spearman.RData")

##### 7/ Plot Mantel test results ####

### Plot relationship between distances

load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_clim_Spearman.RData")
load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_geo_Spearman.RData")
load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_partiel_Ist_Spearman.RData")

load(file = "./outputs/Community_structure/Mantel_tests/Pairwise_values.RData")

# Extract vectors of data without NA
pairwise_Ist_vec <- as.dist(pairwise_Ist)[which(!is.na(as.dist(pairwise_Ist)))]
pairwise_climdist_vec <- as.dist(pairwise_climdist)[which(!is.na(as.dist(pairwise_Ist)))]
pairwise_geodist_vec <- as.dist(pairwise_geodist)[which(!is.na(as.dist(pairwise_Ist)))]

### Select data to subsample for the plot (499 500 points it's a bit too much)
# Set seed for reproductibility
set.seed(23018)
sample_indices <- sample(x= 1:length(pairwise_Ist_vec), size = 1000, replace = F)
save(sample_indices, file = "./outputs/Community_structure/Mantel_tests/sample_indices.RData")

### 7.1/ Ist ~ Dclim ####

# Compute a linear regression to get coefficients to draw a predict line
reg.clim <- lm(pairwise_Ist_vec ~ pairwise_climdist_vec)
summary(reg.clim)

round(Mantel_Spearman_Ist_clim$statistic, 3) # Need to write manually the value on the plot # rho = 0.351

# Plot
pdf(file = "./graphs/Community_Structure/Ist_Dclim.pdf", height = 6, width = 8)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

plot(pairwise_climdist_vec[sample_indices], pairwise_Ist_vec[sample_indices], 
     main = "Climate structures \n mimicry community composition", ylab = bquote('Pairwise' ~I[ST]), xlab = "Standardized climatic distance",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.3, lwd = 1, type = "n", axes = F)

points(pairwise_climdist_vec[sample_indices], pairwise_Ist_vec[sample_indices], pch = 16, col = "#00000070")

legend(legend = c("                   ",
                  "                   "),
       x = "topleft", inset = c(0, 0.02), xjust = 1, y.intersp = 1.1, # Use inset to manually adjust position
       cex = 1.2, bty ="o", bg = "white", box.col = NA)


# legend(legend = c((paste0("\u03C1 = ", # Unicode for small case rho
#                           round(Mantel_Spearman_Ist_clim$statistic, 3))), 
#                   paste0("p < 0.001")), 
#        x = "topleft", inset = c(-0.01, 0.02), xjust = 1, # Use inset to manually adjust position
#        cex = 1.2, bty ="o", bg = "white", box.col = NA)

# legend(legend = c(bquote(rho == .(round(Mantel_Spearman_Ist_clim$statistic, 3))),
#                   paste0("p < 0.001")),
#        x = "topleft", inset = c(-0.01, 0.02), xjust = 1, # Use inset to manually adjust position
#        cex = 1.2, bty ="n", bg = "white", box.col = NA)

legend(legend = c(expression(paste(rho, " = 0.351")),
                  paste0("p = 0.001")),
       x = "topleft", inset = c(-0.01, 0.02), xjust = 1, # Use inset to manually adjust position
       cex = 1.2, bty ="n", bg = "white", box.col = NA)

axis(side = 1, lwd = 2, cex.axis = 1.3)
axis(side = 2, lwd = 2, cex.axis = 1.3)
abline(reg.clim, lwd = 3, col = "red")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()

### 7.2/ Ist ~ Dgeo ####

# Compute a linear regression to get coefficients to draw a predict line
reg.geo <- lm(pairwise_Ist_vec~pairwise_geodist_vec)
summary(reg.geo)

round(Mantel_Spearman_Ist_geo$statistic, 3) # Need to write manually the value on the plot # rho = 0.517

pdf(file = "./graphs/Community_Structure/Ist_Dgeo.pdf", height = 6, width = 8)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

plot(pairwise_geodist_vec[sample_indices], pairwise_Ist_vec[sample_indices],
     main = "Space structures \n mimicry community composition", ylab = bquote('Pairwise' ~I[ST]), xlab = "Geographic distance [km]",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.3, lwd = 1, type = "n", axes = F)

points(pairwise_geodist_vec[sample_indices], pairwise_Ist_vec[sample_indices], pch = 16, col = "#00000070")

legend(legend = c("                   ",
                  "                   "),
       x = "topleft", inset = c(0, 0.02), xjust = 1, y.intersp = 1.113, # Use inset to manually adjust position
       cex = 1.2, bty ="o", bg = "white", box.col = NA)

legend(legend = c(expression(paste(rho, " = 0.517")),
                  paste0("p = 0.001")),
       x = "topleft", inset = c(-0.01, 0.02), xjust = 1, # Use inset to manually adjust position
       cex = 1.2, bty ="n", bg = "white", box.col = NA)

axis(side = 1, lwd = 2, cex.axis = 1.3)
axis(side = 2, lwd = 2, cex.axis = 1.3)
abline(reg.geo, lwd = 3, col = "red")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()

# 7.3/ Ist ~ Dclim + covar(Dgeo) ####

Ist_resid <- residuals(reg.geo, type = "response") # Extract residuals from LM on geographic distances
reg.clim.geo <- lm(Ist_resid ~ pairwise_climdist_vec)
summary(reg.clim.geo)

round(partial_Mantel_Spearman_Ist_clim_geo$statistic, 3) # Need to wirte manually the value on the plot # rho = 0.195

# Plot
pdf(file = "./graphs/Community_Structure/Ist_Dclim_geo.pdf", height = 6, width = 8)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

plot(pairwise_climdist_vec[sample_indices], Ist_resid[sample_indices], 
     main = "Climate structures \n mimicry community composition \n even when space is taking into account", ylab = bquote('Residual Pairwise' ~I[ST]), xlab = "Standardized climatic distance",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.3, lwd = 1, type = "n", axes = F, 
     ylim = c(-0.2, 0.5))

points(pairwise_climdist_vec[sample_indices], Ist_resid[sample_indices], pch = 16, col = "#00000070")

legend(legend = c("                   ",
                  "                   "),
       x = "topleft", inset = c(0, 0.02), xjust = 1, y.intersp = 0.9, # Use inset to manually adjust position
       cex = 1.2, bty ="o", bg = "white", box.col = NA)

legend(legend = c(expression(paste(rho, " = 0.195")),
                  paste0("p = 0.001")),
       x = "topleft", inset = c(-0.01, 0.02), xjust = 1, # Use inset to manually adjust position
       cex = 1.2, bty ="n", bg = "white", box.col = NA)

axis(side = 1, lwd = 2, cex.axis = 1.3)
axis(side = 2, lwd = 2, cex.axis = 1.3)
abline(reg.clim.geo, lwd = 3, col = "red")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()


### 7.4/ Plot null distri for partial Mantel test Ist ~ Dclim + covar(Dgeo) ####
load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_partiel_Ist_Spearman.RData")

pdf(file = "./graphs/Community_Structure/Ist_Dclim_geo_null_distri.pdf", height = 6, width = 8)

original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))


hist(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic),
     breaks = 30, freq = TRUE, col = "gray", 
     main = "partial Mantel test statistic distribution\n under 999 permutations", 
     xlab = expression(paste("Spearman's ", rho)),
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5, lwd = 2)
arrows(partial_Mantel_Spearman_Ist_clim_geo$statistic + 0.000, 80, partial_Mantel_Spearman_Ist_clim_geo$statistic + 0.000, 10, length = 0.1, lwd = 2)
abline(v = mean(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic)), lty = 2, lwd = 2)
abline(v = quantile(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic), 0.95), lty = 2, lwd = 2, col = "red")

legend(legend = c(expression(paste("Mean ",rho, " obs = 0")), 
                  paste0("CI 95% = ", round(quantile(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic), 0.95),3))), 
       x = "topright", inset = c(0.02, 0.00), lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")

legend(legend = c(expression(paste(rho, " obs = 0.195")),
                  paste0("       p = 0.001")),
       x = "bottomright", inset = c(0.02, 0.40), xjust = 1, # Use inset to manually adjust position
       cex = 1.2, bty ="n", bg = "white", box.col = NA)

# legend(legend = c(paste0("Rho obs = ", round(partial_Mantel_Spearman_Ist_clim_geo$statistic, 3)),
#                   paste0("                  p = ", 1-(ecdf(x = c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic))(partial_Mantel_Spearman_Ist_clim_geo$statistic)))),
#        x = "right", cex = 1.2, bty ="n", xjust = 1)

par(mar = original_int_margins)
dev.off()


##### 8/ Final plot with all 4 sub-plots ####

load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_clim_Spearman.RData")
load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_Ist_geo_Spearman.RData")
load(file = "./outputs/Community_structure/Mantel_tests/Resultats_Mantel_partiel_Ist_Spearman.RData")

load(file = "./outputs/Community_structure/Mantel_tests/Pairwise_values.RData")

# Extract vectors of data without NA
pairwise_Ist_vec <- as.dist(pairwise_Ist)[which(!is.na(as.dist(pairwise_Ist)))]
pairwise_climdist_vec <- as.dist(pairwise_climdist)[which(!is.na(as.dist(pairwise_Ist)))]
pairwise_geodist_vec <- as.dist(pairwise_geodist)[which(!is.na(as.dist(pairwise_Ist)))]

load(file = "./outputs/Community_structure/Mantel_tests/sample_indices.RData")

pdf(file = "./graphs/Community_Structure/Mantel_all_plots.pdf", height = 10.3, width = 10)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5,5,1,1), mfrow = c(2,2))

# Panel A = Ist ~ Dclim

# Compute a linear regression to get coefficients to draw a predict line
reg.clim <- lm(pairwise_Ist_vec ~ pairwise_climdist_vec)
# summary(reg.clim)

plot(pairwise_climdist_vec[sample_indices], pairwise_Ist_vec[sample_indices], 
     # main = "Climate structures \n mimicry community composition",
     main = NA, 
     ylab = bquote('Pairwise' ~I[ST]), xlab = "Standardized climatic distance",
     cex.axis = 1.5, cex.lab = 1.7, cex.main = 1.3, lwd = 1, type = "n", axes = F)

points(pairwise_climdist_vec[sample_indices], pairwise_Ist_vec[sample_indices], pch = 16, col = "#00000070")

legend(legend = c("                   ",
                  "                   "),
       x = "topleft", inset = c(0, 0.02), xjust = 1, y.intersp = 1.1, # Use inset to manually adjust position
       cex = 1.2, bty ="o", bg = "white", box.col = NA)

legend(legend = c(as.expression(bquote(bold(paste(rho, " = 0.351")))),
                  as.expression(bquote(bold(paste("p = 0.001"))))),
       x = "topleft", inset = c(-0.03, 0.02), xjust = 1, # Use inset to manually adjust position
       cex = 1.5, text.font = 2,
       bty ="n", bg = "white", box.col = NA)

axis(side = 1, lwd = 2, cex.axis = 1.5)
axis(side = 2, lwd = 2, cex.axis = 1.5)
abline(reg.clim, lwd = 3, col = "red")

legend(legend = as.expression(bquote(bold("A"))), 
       x = "bottomright", inset = c(0.02, 0.01), xjust = 0.5,
       cex = 1.6, bty ="n")

# Panel B = Ist ~ Dgeo

# Compute a linear regression to get coefficients to draw a predict line
reg.geo <- lm(pairwise_Ist_vec~pairwise_geodist_vec)
# summary(reg.geo)

plot(pairwise_geodist_vec[sample_indices], pairwise_Ist_vec[sample_indices],
     # main = "Space structures \n mimicry community composition",
     main = NA, 
     ylab = bquote('Pairwise' ~I[ST]), xlab = "Geographic distance [km]",
     cex.axis = 1.5, cex.lab = 1.7, cex.main = 1.3, lwd = 1, type = "n", axes = F)

points(pairwise_geodist_vec[sample_indices], pairwise_Ist_vec[sample_indices], pch = 16, col = "#00000070")

legend(legend = c("                   ",
                  "                   "),
       x = "topleft", inset = c(0, 0.02), xjust = 1, y.intersp = 1.113, # Use inset to manually adjust position
       cex = 1.2, bty ="o", bg = "white", box.col = NA)

legend(legend = c(as.expression(bquote(bold(paste(rho, " = 0.517")))),
                  as.expression(bquote(bold(paste("p = 0.001"))))),
       x = "topleft", inset = c(-0.03, 0.02), xjust = 1, # Use inset to manually adjust position
       cex = 1.5, bty ="n", bg = "white", box.col = NA)

axis(side = 1, lwd = 2, cex.axis = 1.5)
axis(side = 2, lwd = 2, cex.axis = 1.5)
abline(reg.geo, lwd = 3, col = "red")

legend(legend = as.expression(bquote(bold("B"))), 
       x = "bottomright", inset = c(0.02, 0.01), xjust = 0.5,
       cex = 1.6, bty ="n")

# Panel C = Ist ~ Dclim + covar(Dgeo)

Ist_resid <- residuals(reg.geo, type = "response") # Extract residuals from LM on geographic distances
reg.clim.geo <- lm(Ist_resid ~ pairwise_climdist_vec)
# summary(reg.clim.geo)

plot(pairwise_climdist_vec[sample_indices], Ist_resid[sample_indices], 
     # main = "Climate structures \n mimicry community composition \n even when space is taking into account",
     main = NA, 
     ylab = bquote('Residual Pairwise' ~I[ST]), xlab = "Standardized climatic distance",
     cex.axis = 1.5, cex.lab = 1.7, cex.main = 1.3, lwd = 1, type = "n", axes = F, 
     ylim = c(-0.2, 0.5))

points(pairwise_climdist_vec[sample_indices], Ist_resid[sample_indices], pch = 16, col = "#00000070")

legend(legend = c("                   ",
                  "                   "),
       x = "topleft", inset = c(0, 0.02), xjust = 1, y.intersp = 0.9, # Use inset to manually adjust position
       cex = 1.2, bty ="o", bg = "white", box.col = NA)

legend(legend = c(as.expression(bquote(bold(paste(rho, " = 0.195")))),
                  as.expression(bquote(bold(paste("p = 0.001"))))),
       x = "topleft", inset = c(-0.03, 0.02), xjust = 1, # Use inset to manually adjust position
       cex = 1.5, bty ="n", bg = "white", box.col = NA)

axis(side = 1, lwd = 2, cex.axis = 1.5)
axis(side = 2, lwd = 2, cex.axis = 1.5)
abline(reg.clim.geo, lwd = 3, col = "red")

text(x = -0.59, y = -0.1, labels = "-0.1", cex = 1.5, srt = 90, xpd = T)

legend(legend = as.expression(bquote(bold("C"))), 
       x = "topright", inset = c(0.02, 0.001), xjust = 0.5,
       cex = 1.6, bty ="n")

# Panel D = NUll distri for Ist ~ Dclim + covar(Dgeo)

hist(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic),
     breaks = seq(-0.05, 0.20, 0.05/4),
     freq = TRUE, col = "gray", 
     # main = "partial Mantel test statistic distribution\n under 999 permutations", 
     main = NA, 
     xlab = expression(paste("Spearman's ", rho)),
     cex.axis = 1.5, cex.lab = 1.7, cex.main = 1.5, lwd = 2)
arrows(partial_Mantel_Spearman_Ist_clim_geo$statistic + 0.000, 100, partial_Mantel_Spearman_Ist_clim_geo$statistic + 0.000, 10, length = 0.1, lwd = 2)
abline(v = mean(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic)), lty = 2, lwd = 2)
abline(v = quantile(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic), 0.95), lty = 2, lwd = 2, col = "red")

legend(legend = c(expression(paste("Mean = 0.001")), 
                  paste0("CI 95% = ", round(quantile(c(partial_Mantel_Spearman_Ist_clim_geo$perm, partial_Mantel_Spearman_Ist_clim_geo$statistic), 0.95),3))), 
       x = "topright", inset = c(0.02, 0.17), lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.5, bty ="n")

legend(legend = c(as.expression(bquote(bold(paste(rho, " obs = 0.195")))),
                  as.expression(bquote(bold(paste("       p = 0.001"))))),
       x = "bottomright", inset = c(0.02, 0.40), xjust = 1, # Use inset to manually adjust position
       cex = 1.5, bty ="n", bg = "white", box.col = NA)

legend(legend = as.expression(bquote(bold("D"))), 
       x = "topright", inset = c(0.02, 0.001), xjust = 0.5,
       cex = 1.6, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins, mfrow = c(1,1))        

dev.off()
