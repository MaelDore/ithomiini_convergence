##### Community composition analysis #####

### Preparation ###

# Effacer l'environnement
rm(list = ls())

library(raster)
library(rangeBuilder)

### 1/ Load functions and files ####

# Function to aggregate probabilities among several OMUs of the same species
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x)
  return(y) 
}

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

# From Dore et al., 2020 - Global diversity patterns in Ithomiine butterflies in the face of human impacts ###

# Load summary table for species list
load(file = paste0("./input_data/Summary_tables/list.sp.RData"))
# Load summary table for OMU list
load(file = paste0("./input_data/Summary_tables/list.unit.RData"))

# Load the OMU probability stack
OMU_proba_stack <- readRDS(file = paste0("./input_data/SDM_stacks/All_OMU_stack_Jaccard.80.rds"))
# # Load the species probability stack
# sp_proba_stack <- readRDS(file = paste0("./input_data/SDM_stacks/All_sp_proba_stack_Jaccard.80.rds"))
# Load the mimicry ring richness stack
ring_richness_stack <- readRDS(file = paste0("./input_data/SDM_stacks/All_ring_rich_stack_Jaccard.80.rds"))

# Color palet for plot
pal_bl_red_Mannion <- readRDS(file = "./input_data/Map_stuff/pal_bl_red_Mannion.rds")

# Load mask for continent borders
continent_mask <- readRDS(file = "./input_data/Map_stuff/continent_mask_15.rds")
crop_mask_shp <- readRDS(file = "./input_data/Map_stuff/crop_mask_shp_15.rds")
# bg_mask_pixel <- readRDS(file = "./input_data/Map_stuff/bg_mask_pixel.rds")
bg_mask <- readRDS(file = "./input_data/Map_stuff/bg_mask.rds")

##### 2/ Compute observed Ist value for whole communities #####

### 2.1/ Shift to raster brick to extract data under a single df for all layers ####
ring_richness_brick <- ring_richness_stack*1

### 2.2/ Compute Community diversity index/map = Dk ####
Dk <- calc(ring_richness_stack, fun = simpson, na.rm = T)*1
saveRDS(Dk, file = "./outputs/Community_Structure/Dk_map.rds")

### Plot Simpson's Mimicry Diversity Map
Dk <- readRDS(file = "./outputs/Community_Structure/Dk_map.rds")
Dk_map <- merge(Dk, continent_mask)
# Add 0 values in terrestrial pixel with no Ithomiini

pdf(file = paste0("./maps/Community_structure/Dk.pdf"), height = 6.3, width = 6.5)

image(Dk_map, col = pal_bl_red_Mannion, main = "Simpson's Mimicry Diversity Map", 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
addRasterLegend(Dk_map, locs = seq(0, 0.9, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.1, label = "Simpson's\n diversity")

dev.off()


### 2.3/ Compute Mean community diversity = Ds = 0.762 ####
Global_Ds <- mean(Dk@data@values, na.rm = T)

### 2.4/ Compute Total diversity index = Dt = for all communities merged in one = 0.911 ####
mimicry.abundances <- NA # To store each mimicry global expected richness (i.e., nb of species of this ring among all communities, duplicates included)
for (i in 1:nlayers(ring_richness_stack)) {
  mimicry.abundances[i] <- sum(ring_richness_stack[[i]]@data@values, na.rm = T)
}
# hist(mimicry.abundances)
Global_Dt <- simpson(mimicry.abundances)


### 2.5/ Compute Global Ist = 0.164
Global_Ist <- 1-(Global_Ds/Global_Dt)

# Save stuff
save(Global_Ds, Global_Dt, Global_Ist, file = "./outputs/Community_Structure/Global.Ist.RData")


##### 3/ Generate new virtual community matrices under null hypothesis of no effect of mimicry ring on species presence ####
# analogous to DeVries et al.'s [1999] and Hill's [2010] test of mimicry structure across microhabitats

mimicry.list <- as.character(unique(list.unit$Mimicry.model)) # 44 Mimicry rings

## Start the loop to shuffle mimicry patterns among OMUs
Ist_null <- NA
for (k in 1:1000) # 1000 permutations
{ 
  # k <- 1
  
  ## Shuffle mimicry ring among units
  shuffle.list.unit <- list.unit
  shuffle.list.unit$Mimicry.model <- sample(as.character(shuffle.list.unit$Mimicry.model))
  
  ## Generate the new Mimicry rings Richness Stack with random attribution of species to mimicry ring
  
  simul.ring_richness_stack <- Dk # 1st layer of the stack to remove later
  for (i in 1:length(mimicry.list)) # Per mimicry ring
  { 
    # i <- 1
    
    ring <- as.character(mimicry.list[i]) # Load mimicry ring name
    
    temp_ring_stack <- subset(x = OMU_proba_stack, subset = which(shuffle.list.unit$Mimicry.model == ring), drop = F)
    
    simul_ring_richness <- calc(temp_ring_stack, fun = sum)*1
    
    simul.ring_richness_stack <- addLayer(simul.ring_richness_stack, simul_ring_richness)
    
    # print(i)
    
  }
  simul.ring_richness_stack <- dropLayer(simul.ring_richness_stack, i = 1)*1
  names(simul.ring_richness_stack) <- mimicry.list
  
  save(simul.ring_richness_stack, file = paste0("./outputs/Community_Structure/Permutations/Ist/Ring_Richness_Stack_simul_",k,".RData"))
  
  # plot(ring_richness_stack)        # Original  
  # plot(simul.ring_richness_stack)  # Simulated
  
  
  ## Compute Community diversity index/map = Dk
  simul.Dk <- calc(simul.ring_richness_stack, fun = simpson, na.rm = T)*1
  # Plot Simulated Simpson's Mimicry Diversity Map
  # plot(simul.Dk, col = pal_bl_red_Mannion, main = paste0("Simulated Simpson's Mimicry Diversity Map n°", k))
  
  ## Compute index for all communities
  
  # Mean community diversity = Ds
  Simul_global_Ds <- mean(simul.Dk@data@values, na.rm = T)
  
  # Total diversity index = Dt
  mimicry.abundances <- NA # To store each mimicry global expected richness
  for (i in 1:nlayers(simul.ring_richness_stack)) 
  {
    mimicry.abundances[i] <- sum(simul.ring_richness_stack[[i]]@data@values, na.rm = T)
  }
  # hist(mimicry.abundances)
  Simul_global_Dt <- simpson(mimicry.abundances)
  
  # Global Ist
  Simul_global_Ist <- 1-(Simul_global_Ds/Simul_global_Dt)
  
  ## Save stuff
  Ist_null[k] <- Simul_global_Ist
  
  save(Ist_null, file = "./outputs/Community_Structure/Ist_null.RData")
  saveRDS(Ist_null, file = "./outputs/Community_Structure/Ist_null.rds")
  
  cat(paste0(Sys.time(), " - Simul n°", k,"\n"))
}

summary(Ist_null)

##### 4/ Plot distri of global null Ist ####

load(file = "./outputs/Community_Structure/Ist_null.RData")
load(file = "./outputs//Community_Structure/Global.Ist.RData")

# Final plot

pdf(file = "./graphs/Community_Structure/Global_Ist_null.pdf", height = 5.3, width = 6.5)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

hist(c(Ist_null, Global_Ist), 50, freq = TRUE, col = "gray", xlab = bquote(~I[ST]), 
     main = bquote('Distribution of' ~I[ST]~ '\nunder null Hypothesis'),
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5, lwd = 2)
arrows(Global_Ist - 0.00055, 75, Global_Ist- 0.00055, 10, length = 0.1, lwd = 2) # Draw arrow above Ist obs
abline(v = mean(c(Ist_null, Global_Ist)), lty = 2, lwd = 2) # Add vertical line for mean value
abline(v = quantile(c(Ist_null, Global_Ist), 0.95), lty = 2, lwd = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("   Mean = ", format(round(mean(c(Ist_null, Global_Ist)),3), nsmall = 3)), 
                  paste0("CI 95% = ", round(quantile(c(Ist_null, Global_Ist), 0.95),3))), 
       x = "topright", inset = c(0, 0.17), xjust = 1, # Use inset to manually adjust position
       cex = 1.2, # text.width = c(0.9, 1),
       lty = 2, lwd = 2, col = c("black", "red"), bty ="n")
# legend(legend = c(paste0("          Ist obs = ", round(Global_Ist, 3)),
#                   paste0("                   p = 0.001")),
#        x = "right", cex = 1.2, bty ="n", xjust = 1)
legend(legend = bquote("           " ~ I[ST] ~ 'obs =' ~ .(round(Global_Ist, 3))),
       x = "right", cex = 1.2, bty ="n", xjust = 1)
legend(legend = c("",
                  paste0("                   p = 0.001")),
       x = "right", cex = 1.2, bty ="n", xjust = 1, y.intersp = 1.8)

legend(legend = as.expression(bquote(bold("A"))), 
       x = "topright", inset = c(0.001, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()



########################################## Approach with Bray-Curtis indices on comimetic pairs of species #############################

# Effacer l'environnement
rm(list = ls())

library(raster)

### 1/ Load directly the unit probability stack ####
OMU_proba_stack <- readRDS(file = paste0("./input_data/SDM_stacks/All_OMU_stack_Jaccard.80.rds"))

### 2/ Generate matrix communities * OMU for probabilities of presence ####
OMU_proba_brick <- OMU_proba_stack*1
com_unit_mat_783 <- NA
for (i in 1:nlayers(OMU_proba_brick)) 
{
  if (i==1){
    com_unit_mat_783 <- OMU_proba_brick[[i]]@data@values
  }else{
    unit.values <- OMU_proba_brick[[i]]@data@values
    com_unit_mat_783 <- cbind(com_unit_mat_783, unit.values)
  }
  if (i %% 10 == 0) {print(i)}
}
names(com_unit_mat_783) <- names(OMU_proba_brick)
rm(OMU_proba_brick)

saveRDS(com_unit_mat_783, file = "./input_data/SDM_stacks/com_unit_mat_783.rds")

# Filter out communities with NA
real_com_index <- NA
for (i in 1:nrow(com_unit_mat_783)) 
{
  com.row <- com_unit_mat_783[i,]
  real_com_index[i] <- (!any(is.na(com.row)))&(!any(is.nan(com.row)))
}
Filtered_com_unit_mat_783 <- com_unit_mat_783[real_com_index,]
nrow(Filtered_com_unit_mat_783) # 24986 community with no NA

saveRDS(Filtered_com_unit_mat_783, file = "./input_data/SDM_stacks/Filtered_com_unit_mat_783.rds")

# Load directly the matrix of community * units probability presence
com_unit_mat_783 <- readRDS(file = "./input_data/SDM_stacks/Filtered_com_unit_mat_783.rds")

### 3/ Compute Bray-Curtis indices ####

### 3.1/ Compute Bray-Curtis index for all pairs of units ####
library(vegan)

dim(com_unit_mat_783)
unit_BC_dist <- vegdist(x = t(com_unit_mat_783), method = "bray") # Compute dissimilarities between rows

saveRDS(unit_BC_dist, file = "./outputs/Community_Structure/unit_BC_dist.rds")
save(unit_BC_dist, file = "./outputs/Community_Structure/unit_BC_dist.RData")

# Load directly the vector of BD distances among all units
unit_BC_dist <- readRDS(file = "./outputs/Community_Structure/unit_BC_dist.rds")
OMU_proba_stack <- readRDS(file = "./input_data/SDM_stacks/All_OMU_stack_Jaccard.80.rds")

# hist(unit_BC_dist)

unit_BC_dist_mat <- as.matrix(unit_BC_dist)

mimicry.list <- as.character(unique(list.unit$Mimicry.model))

### 3.2/ Compute mean BC per mimicry rings ####
mean_BC <- NA
for (i in 1:length(mimicry.list)) # Per mimicry rings
{ 
  ring <- mimicry.list[i]
  
  # Get names and indices of all units/OMUs for this ring
  tags <- as.character(list.unit$Tag.model[list.unit$Mimicry.model == ring]) # Get names of all OMUs in the ring
  index <- which(names(OMU_proba_stack) %in% tags) # Get index of layer for these OMUs
  
  if (length(tags) == 1) # Case with only one OMU in the mimicry ring. Impossible to compute pairwise distances
  {
    mat.index <- data.frame(matrix(nrow = 2, ncol = 0)) # Empty df of pairs indices (no pairs of OMUs)
    save(mat.index, file = paste0("./outputs/Community_Structure/BC_mimicry_ring_mat/mat.index.",ring,".RData"))
    BC <- c() # Empty vector of Bray-Curtis values
    save(BC, file = paste0("./outputs/Community_Structure/BC_mimicry_ring_mat/BC.",ring,".RData"))
    mean_BC[i] <- NA
  } else { # Case with at least 2 OMUs.
    mat.index <- combn(x = index, m = 2, FUN = c) # Get all combinations possibles for pairs of OMUs indices
    save(mat.index, file = paste0("./outputs/Community_Structure/BC_mimicry_ring_mat/mat.index.",ring,".RData"))
    BC <- NA # Initiate the vector used to store all BC values
    for (j in 1:ncol(mat.index)) # For each pair of OMUs
    {
      BC[j] <- unit_BC_dist_mat[mat.index[1,j], mat.index[2,j]] # Extract Bray-Curtis index from the complete matrix of BC indices
    }
    save(BC, file = paste0("./outputs/Community_Structure/BC_mimicry_ring_mat/BC.",ring,".RData")) # Save the vector of BC for all pairs of OMUs in this ring
    mean_BC[i]  <- mean(BC) # Compute mean BC for all pairs of this mimicry ring and store it in final vector
  }
  
  cat(paste0(Sys.time(), " - ", ring," - n°",i, " on ",length(mimicry.list),"\n"))
  
}
names(mean_BC) <- mimicry.list # Associate mean BC value with ring name
save(mean_BC, file = paste0("./outputs/Community_Structure/mean_BC.RData")) # Save final vector with mean BC values per ring
saveRDS(mean_BC, file = paste0("./outputs/Community_Structure/mean_BC.rds")) # Save final vector with mean BC values per ring

# Load mean BC of mimicry rings
load(file = paste0("./outputs/Community_Structure/mean_BC.RData"))
mean_BC

### 3.3/ Compute mean BC for all comimics, all rings taken into account ####

# Retrieve all co-mimic coordinates

all_mimic_mat_index <- data.frame(matrix(ncol = 0, nrow=2)) # Generate empty df to store indices of pairs of comimics
for (i in 1:length(mimicry.list))  # Per mimetic ring
{ 
  ring <- mimicry.list[i]
  load(file = paste0("./outputs/Community_Structure/BC_mimicry_ring_mat/mat.index.",ring,".RData")) # Load the matrix of indices of pairsof comimics for this ring
  all_mimic_mat_index <- cbind(all_mimic_mat_index,mat.index) # Merge them all in one df
}
dim(all_mimic_mat_index) # 14 659 pairs of co-mimic OMUs

save(all_mimic_mat_index, file = paste0("./outputs/Community_Structure/BC_mimicry_ring_mat/all_mimic_mat_index.RData"))
saveRDS(all_mimic_mat_index, file = paste0("./outputs/Community_Structure/BC_mimicry_ring_mat/all_mimic_mat_index.rds"))

# Retrieve all co-mimic BC values and non-co-mimic BC values at the same time
BC_mimic <- NA # Initiate final vector to store BC values for comimics
unit_BC_dist_mat_no_mimic <- unit_BC_dist_mat # Copy matrix of all BC indices for all pairs of OMUs
for (j in 1:ncol(all_mimic_mat_index)) # For all pairs of comimics
{
  BC_mimic[j] <- unit_BC_dist_mat[all_mimic_mat_index[1,j],all_mimic_mat_index[2,j]] # Extract BC index for this pair of comimics
  unit_BC_dist_mat_no_mimic[all_mimic_mat_index[1,j],all_mimic_mat_index[2,j]] <- NA # Remove BC value of comimics from matrix of non-comimics pairs
  if (j %% 1000 == 0) {print(j)}
}
length(BC_mimic) # 14 659 pairs of co-mimic OMUs

save(BC_mimic, file = paste0("./outputs/Community_Structure/BC_mimic.RData")) # Save final vector with all BC values for comimics
saveRDS(BC_mimic, file = paste0("./outputs/Community_Structure/BC_mimic.rds")) # Save final vector with all BC values for comimics

# Retrieve all non-mimic BC values
BC_no_mimic <- unit_BC_dist_mat_no_mimic[upper.tri(unit_BC_dist_mat_no_mimic)] # Extract only one side of the triangle to avoid duplicate values of pairs
BC_no_mimic <- na.omit(BC_no_mimic) # Remove NA (the comimics pairs)
length(BC_no_mimic) # 291 494 pairs of non-mimic units

save(BC_no_mimic, file = paste0("./outputs/Community_Structure/BC_no_mimic.RData")) # Save final vector with all BC values for non-comimics
saveRDS(BC_no_mimic, file = paste0("./outputs/Community_Structure/BC_no_mimic.rds")) # Save final vector with all BC values for non-comimics


# Compute mean obs BC for each group
load(file = "./outputs/Community_Structure/unit_BC_dist.RData") # All BC index (distance format)
load(file = "./outputs/Community_Structure/BC_mimic.RData") # Only mimic pairs (vector format)
load(file = "./outputs/Community_Structure/BC_no_mimic.RData") # Only non mimic pairs (vector format)

Global_mean_BC <- mean(unit_BC_dist) ; Global_mean_BC # 0.949
Global_mimic_mean_BC <- mean(BC_mimic) ; Global_mimic_mean_BC # 0.896
Global_no.mimic_mean_BC <- mean(BC_no_mimic) ; Global_no.mimic_mean_BC # 0.952

# Save all mean BC values for all pairs, only co-mimics, only non-comimics
save(Global_mean_BC, Global_mimic_mean_BC, Global_no.mimic_mean_BC, file = "./outputs/Community_Structure/All_Global_BC.RData")

# Quick ugly boxplots
All_BC <- c(BC_mimic, BC_no_mimic)
Status <- as.factor(c(rep("Mimic", length(BC_mimic)), rep("Non-Mimic", length(BC_no_mimic))))
boxplot(All_BC ~ Status)

library(ggplot2)

boxplot_df <- data.frame(All_BC, Status, stringsAsFactors = T)
gg_boxplot <- ggplot(data = boxplot_df, aes(x = Status, y = All_BC)) +
  geom_violin() +
  coord_trans(y = "exp")
print(gg_boxplot)

# Histogram
# Plot distri of BC for mimic pairs
pdf(file = paste0("./graphs/Community_structure/hist_BC_mimic.pdf"), height = 6.3, width = 6.5)
hist(BC_mimic, xlab = "Bray-Curtis index", main = "Bray-Curtis indices of comimic pairs")
abline(v = mean(BC_mimic), col = "red", lty = 2, lwd = 2)
legend(legend = c(paste0("Mean = ", round(mean(BC_mimic, na.rm = T),3)), 
                  paste0("CI 5% = ", round(quantile(BC_mimic, 0.05),3)),
                  paste0("CI 95% = ", round(quantile(BC_mimic, 0.95),3))),
       x = "topleft", cex = 1, bty ="n") 
dev.off()


##### 4/ Generate new virtual community matrices under null hypothesis of no effect of mimicry ring on species presence ####
# analogous to DeVries et al.'s [1999] and Hill's [2010] test of mimicry structure across microhabitats

# Load directly the vector of BC distances among all units
unit_BC_dist <- readRDS(file = "./outputs/Community_Structure/unit_BC_dist.rds")
OMU_proba_stack <- readRDS(file = "./input_data/SDM_stacks/All_OMU_stack_Jaccard.80.rds")

unit_BC_dist_mat <- as.matrix(unit_BC_dist)

mimicry.list <- as.character(unique(list.unit$Mimicry.model))


## Start the loop for permutations
BC_mimic_null <- BC_no_mimic_null <- NA # Create vectors to store simulated mean BC values for co-mimics and non-comimics
mean_BC_null <- matrix(ncol = length(mimicry.list), nrow = 0) # Create matrix to store mean simulated BC values for each mimicry ring for each simulation
for (k in 1:1000) # 1000 simulations/permutations
{ 
  # k <- 1
  
  ## Suffle mimicry ring among units
  shuffle.list.unit <- list.unit
  shuffle.list.unit$Mimicry.model <- sample(as.character(shuffle.list.unit$Mimicry.model))
  
  # # Check if number of unit per ring is preserved
  # table(list.unit$Mimicry.model)
  # table(shuffle.list.unit$Mimicry.model)
  
  ## Generate the new Mimicry rings Richness Stack with random attribution of OMUs to mimicry ring
  
  # Mimicry list
  mimicry.list <- as.character(unique(list.unit$Mimicry.model)) # 44 Mimicry rings
  
  mean_BC <- NA
  for (i in 1:length(mimicry.list)) # Per mimicry rings
  { 
    # i <-  1
    ring <- mimicry.list[i]
    
    # Get names and indices of all units/OMUs for this ring
    tags <- as.character(list.unit$Tag.model[shuffle.list.unit$Mimicry.model == ring])
    index <- which(names(OMU_proba_stack) %in% tags)
    
    if (length(tags)==1) # Case with only one OMU in the mimicry ring. Impossible to compute pairwise distances
    { 
      mat.index <- data.frame(matrix(nrow = 2, ncol = 0)) # Empty df of pairs coordinates
      save(mat.index, file = paste0("./outputs/Community_Structure/Permutations/BC_mimicry_ring_mat/simul_mat.index.",ring,".RData")) # Save temp mat indices df to be used later to merge for all rings
      mean_BC[i] <- NA
    } else { # Case with at least 2 OMUs
      mat.index <- combn(x = index, m = 2, FUN = c) # Get all combinations possibles for pairs of OMUs indices
      save(mat.index, file = paste0("./outputs/Community_Structure/Permutations/BC_mimicry_ring_mat/simul_mat.index.",ring,".RData")) # Save temp mat indices df to be used later to merge for all rings
      BC <- NA # Initiate the vector used to store all BC values
      for(j in 1:ncol(mat.index)) # For each pair of OMUs
      {
        BC[j] <- unit_BC_dist_mat[mat.index[1,j],mat.index[2,j]] # Extract Bray-Curtis index from the complete matrix of BC indices
      }
      mean_BC[i]  <- mean(BC) # Compute mean BC for all pairs of this mimicry ring and store it in final vector
    }
    
    # print(i)
  }
  names(mean_BC) <- mimicry.list
  mean_BC_null <- rbind(mean_BC_null,mean_BC) # Store the mean BC vector into the matrix of simulations, as a row
  
  # Retrieve all mimic coordinates
  all_mimic_mat_index <- data.frame(matrix(ncol = 0, nrow = 2)) # Generate empty df to store indices of pairs of comimics
  for (i in 1:length(mimicry.list)) { # Per mimicry rings
    ring <- mimicry.list[i]
    load(file =  paste0("./outputs/Community_Structure/Permutations/BC_mimicry_ring_mat/simul_mat.index.",ring,".RData")) # Load the matrix of indices of pairs of comimics for this ring
    all_mimic_mat_index <- cbind(all_mimic_mat_index,mat.index) # Merge them all in one df
  }
  # dim(all_mimic_mat_index) # 14 659 pairs of mimic units
  
  # Retrieve all mimic values and non-co-mimic BC values at the same time
  BC_mimic <- NA # Initiate final vector to store BC values for comimics
  unit_BC_dist_mat_no_mimic <- unit_BC_dist_mat # Copy matrix of all BC indices for all pairs of OMUs
  for(j in 1:ncol(all_mimic_mat_index)) # For all pairs of comimics
  {
    BC_mimic[j] <- unit_BC_dist_mat[all_mimic_mat_index[1,j],all_mimic_mat_index[2,j]] # Extract BC index for this pair of comimics
    unit_BC_dist_mat_no_mimic[all_mimic_mat_index[1,j],all_mimic_mat_index[2,j]] <- NA # Remove BC value of comimics from matrix of non-comimics pairs
  }
  # length(BC_mimic) # 14 659 pairs of mimic units
  
  # Retrieve all non-mimic values
  BC_no_mimic <- unit_BC_dist_mat_no_mimic[upper.tri(unit_BC_dist_mat_no_mimic)] # Extract only one side of the triangle to avoid duplicate values of pairs
  BC_no_mimic <- na.omit(BC_no_mimic) # Remove NA (the comimics pairs)
  # length(BC_no_mimic) # 291 494 pairs of non-mimic units
  
  # Save global mean computation into final vectors
  BC_mimic_null[k] <- mean(BC_mimic)
  BC_no_mimic_null[k] <- mean(BC_no_mimic)
  
  cat(paste0(Sys.time(), " - Simul n°", k," out of 1000\n"))
  save(mean_BC_null, BC_mimic_null, BC_no_mimic_null, file = paste0("./outputs/Community_Structure/Permutations/All_simul_mean_BC.RData"))
}

summary(mean_BC_null)     # Per mimicry ring
summary(BC_mimic_null)    # For all comimic pairs
summary(BC_no_mimic_null) # For all non-comimic pairs

##### 5/ Plot distri of null mean BC ####

### 5.1/ Plot for all comimics ####

load(file = "./outputs/Community_Structure/All_Global_BC.RData") # Load mean BC obs values 
load(file = "./outputs/Community_Structure/Permutations/All_simul_mean_BC.RData") # Load null mean BC values

pdf(file = "./graphs/Community_Structure/BC_mimic_null.pdf", height = 5.3, width = 6.5)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

hist(x = c(BC_mimic_null, Global_mimic_mean_BC), breaks = 40, freq = TRUE, col = "gray",
     xlim = c(0.89, 0.96),
     ylim = c(0, 250),
     main = "Distribution of mean pairwise BC indices\nof co-mimetic species\nunder null Hypothesis",
     xlab = "Mean Bray-Curtis Index",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)

arrows(Global_mimic_mean_BC + 0.0001, 60, Global_mimic_mean_BC+ 0.0001, 10, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = mean(c(BC_mimic_null, Global_mimic_mean_BC)), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(c(BC_mimic_null, Global_mimic_mean_BC), 0.05), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("Mean = ", round(mean(c(BC_mimic_null, Global_mimic_mean_BC)),3)), 
                  paste0("CI 5% = ", round(quantile(c(BC_mimic_null, Global_mimic_mean_BC), 0.95),3))), 
       x = "topleft", inset = c(0, 0.17), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty = "n")
legend(legend = c(paste0("Mean BC obs = ", round(Global_mimic_mean_BC, 3)),
                  paste0("p = 0.001")),
       x = "left", inset = c(-0.05, 0),
       cex = 1.2, bty ="n", xjust = 1)

legend(legend = as.expression(bquote(bold("B"))), 
       x = "topleft", inset = c(-0.03, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()



### 5.2/ Plot distri of null BC per mimicry ring co-mimics ####

load(file = paste0("./outputs/Community_Structure/mean_BC.RData")) # Load mean BC obs value for each ring
load(file = "./outputs/Community_Structure/Permutations/All_simul_mean_BC.RData") # Load null mean BC values

# Build summary table for mimicry ring BC results at the same time
BC_ring_summary_table <- as.data.frame(matrix(ncol = 8, nrow = 44, data = NA))
names(BC_ring_summary_table) <- c("ring", "N_units", "N_pairs", "BC_obs", "mean_BC", "BC_2.5", "BC_97.5", "p_value")

for (i in 1:length(mimicry.list))  # Per mimicry ring
{
  # i <- 1
  
  BC_ring_summary_table$ring[i] <- ring <- mimicry.list[i] # Get ring name
  BC_ring_summary_table$N_units[i] <- N_units <- sum(list.unit$Mimicry.model == ring) # Get number of OMUs
  BC_ring_summary_table$N_pairs[i] <- N_pairs <- N_units*(N_units-1)/2 # Get number of pairs
  
  if (is.na(mean_BC[i]))# Case for ring with no pair of OMUs (only one OMU/species)
  {
    # Plot BC null distri for this ring (i.e. NA message)
    pdf(file = paste0("./graphs/Community_Structure/Per_mimicry_ring/BC_mimic_null_",ring,".pdf"), height = 5.3, width = 6.5)
    
    plot(1:100,1:100, type = "n", main = paste0("Distribution of mean pairwise Bray-Curtis indices\nof ", ring, " OMUs\nunder null Hypothesis"))
    text(x = 50, y = 50, labels = "Only one OMU for this mimicry ring \n No pair available for index computation")
    
    dev.off()
    
  } else { # Case for ring with pair(s) of OMUs 
    
    BC_ring_summary_table$BC_obs[i] <- round(mean_BC[i],3) # Get mean BC obs
    BC_ring_summary_table$mean_BC[i] <- round(mean(c(mean_BC_null[,i], mean_BC[i]), na.rm = T),3) # Get mean of mean BC null from simulations
    BC_ring_summary_table$BC_2.5[i] <- round(quantile(c(mean_BC_null[,i], mean_BC[i]), 0.025),3) # Get 2.5% quantile of mean BC null from simulations
    BC_ring_summary_table$BC_97.5[i] <- round(quantile(c(mean_BC_null[,i], mean_BC[i]), 0.975),3) # Get 97.5% quantile of mean BC null from simulations
    BC_ring_summary_table$p_value[i] <- round(ecdf(x = c(mean_BC_null[,i], mean_BC[i]))(mean_BC[i]),3) # Get p-value from simulations
    
    save(BC_ring_summary_table, file = "./outputs/Community_Structure/BC_ring_summary_table.Rdata")
    
    # Plot BC null distri for this ring
    pdf(file = paste0("./graphs/Community_Structure/Per_mimicry_ring/BC_mimic_null_",ring,".pdf"), height = 5.3, width = 6.5)
    
    hist(c(mean_BC_null[,i], mean_BC[i]), 
         breaks = seq(from = floor(floor(min(min(mean_BC_null[,i]), mean_BC[i])*1000)/2)/500, to = ceiling(ceiling(max(max(mean_BC_null[,i]), mean_BC[i])*1000)/2)/500, by = 0.002), 
         xlab = "Mean pairwise Bray-Curtis index", 
         main = paste0("Distribution of mean pairwise Bray-Curtis indices\nof ", ring, " OMUs\nunder null Hypothesis"))
    
    # arrows(mean_BC[i] + 0.0001, 60, mean_BC[i] + 0.0001, 10, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
    abline(v = mean_BC[i], lwd = 2, lty = 2, col = "red") # Add vertical line for mean BC obs value
    
    abline(v = mean(c(mean_BC_null[,i], mean_BC[i]), na.rm = T), lwd = 2, lty = 2) # Add vertical line for mean value
    abline(v = quantile(c(mean_BC_null[,i], mean_BC[i]), 0.025, na.rm = T), lwd = 2, lty = 2, col = "grey30") # Add vertical line for 2.5% value
    abline(v = quantile(c(mean_BC_null[,i], mean_BC[i]), 0.975, na.rm = T), lwd = 2, lty = 2, col = "grey30") # Add vertical line for 97.5% value
    
    legend(legend = c(paste0("N units = ", N_units), 
                      paste0("N pairs = ", N_pairs)),
           x = "topright", xjust = 1, cex = 1, bty ="n")
    
    legend(legend = c(paste0("Mean = ", round(mean(c(mean_BC_null[,i], mean_BC[i]), na.rm = T),3)), 
                      paste0("CI 2.5% = ", round(quantile(c(mean_BC_null[,i], mean_BC[i]), 0.025),3)),
                      paste0("CI 97.5% = ", round(quantile(c(mean_BC_null[,i], mean_BC[i]), 0.975),3))),
           x = "topleft", cex = 1, bty ="n") 
    
    legend(legend = paste0("Mean BC obs = ", round(mean_BC[i], 3), "\n p = ", round(ecdf(x = mean_BC_null[,i])(mean_BC[i]),3)), 
           inset = c(0, 0.2), text.col = "red", x = "bottomleft", cex = 1, bty ="n")
    
    dev.off()
    
    
  }
  
  cat(paste0(i, " out of 44 - ",ring, " - Done \n"))
}

View(BC_ring_summary_table)


##### 6/ Final plot with both Ist and BC #####

# Remove titles

load(file = "./outputs/Community_Structure/Ist_null.RData")
load(file = "./outputs//Community_Structure/Global.Ist.RData")

# Final plot

pdf(file = "./graphs/Community_Structure/Ist_&_BC.pdf", height = 5.3, width = 14)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4), mfrow = (c(1,2)))

hist(c(Ist_null, Global_Ist), 50, freq = TRUE, col = "gray", xlab = bquote(~I[ST]), 
     main = "",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5, lwd = 2)
arrows(Global_Ist - 0.00055, 75, Global_Ist- 0.00055, 10, length = 0.1, lwd = 2) # Draw arrow above Ist obs
abline(v = mean(c(Ist_null, Global_Ist)), lty = 2, lwd = 2) # Add vertical line for mean value
abline(v = quantile(c(Ist_null, Global_Ist), 0.95), lty = 2, lwd = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("   Mean = ", format(round(mean(c(Ist_null, Global_Ist)),3), nsmall = 3)), 
                  paste0("CI 95% = ", round(quantile(c(Ist_null, Global_Ist), 0.95),3))), 
       x = "topright", inset = c(0, 0.17), xjust = 1, # Use inset to manually adjust position
       cex = 1.2, # text.width = c(0.9, 1),
       lty = 2, lwd = 2, col = c("black", "red"), bty ="n")
# legend(legend = c(paste0("          Ist obs = ", round(Global_Ist, 3)),
#                   paste0("                   p = 0.001")),
#        x = "right", cex = 1.2, bty ="n", xjust = 1)
legend(legend = bquote("           " ~ I[ST] ~ 'obs =' ~ .(round(Global_Ist, 3))),
       x = "right", cex = 1.2, bty ="n", xjust = 1)
legend(legend = c("",
                  paste0("                   p = 0.001")),
       x = "right", cex = 1.2, bty ="n", xjust = 1, y.intersp = 1.8)

legend(legend = as.expression(bquote(bold("A"))), 
       x = "topright", inset = c(0.001, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

hist(x = c(BC_mimic_null, Global_mimic_mean_BC), breaks = 40, freq = TRUE, col = "gray",
     xlim = c(0.89, 0.96),
     ylim = c(0, 250),
     main = "Distribution of mean pairwise BC indices\nof co-mimetic species\nunder null Hypothesis",
     xlab = "Mean Bray-Curtis Index",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)

arrows(Global_mimic_mean_BC + 0.0001, 60, Global_mimic_mean_BC+ 0.0001, 10, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = mean(c(BC_mimic_null, Global_mimic_mean_BC)), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(c(BC_mimic_null, Global_mimic_mean_BC), 0.05), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("Mean = ", round(mean(c(BC_mimic_null, Global_mimic_mean_BC)),3)), 
                  paste0("CI 5% = ", round(quantile(c(BC_mimic_null, Global_mimic_mean_BC), 0.95),3))), 
       x = "topleft", inset = c(0, 0.17), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty = "n")
legend(legend = c(paste0("Mean BC obs = ", round(Global_mimic_mean_BC, 3)),
                  paste0("p = 0.001")),
       x = "left", inset = c(-0.05, 0),
       cex = 1.2, bty ="n", xjust = 1)

legend(legend = as.expression(bquote(bold("B"))), 
       x = "topleft", inset = c(-0.03, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins, mfrow = (c(1,1)))

dev.off()






##################  Mantel distances for pairwise Ist ~ Dclim  #########################

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

# # Load species richness map
# sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))

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
Pairwise_geodist = geosphere::distm(x = coord)/1000 # Geographic distance on the WGS84 ellipsoid, in km

### 5/ Compute pairwise euclidian climatic distances on standardized climatic variables ####
sampled_com_env <- climate_stack@data@values[sample.index,] # Extract climatic variables for subsampled communities
sampled_com_env <- scale(x = sampled_com_env, center = T, scale = T) # Standardize the variables
Pairwise_climdist <- as.matrix(dist(x = sampled_com_env, method = "euclidian")) # Compute euclidian distances

# Save all pairwise distances
save(pairwise_Ist, Pairwise_climdist, Pairwise_geodist, file = "./outputs/Niche_evolution/Mantel_tests/Pairwise_values.RData")

load(file = paste0(internal.wd, "/Community_Structure/Pairwise_values.RData"))

### 6/ Compute Mantel tests ####

library("vegan")

# Check for normality
hist(Pairwise_climdist)
hist(Pairwise_geodist)
hist(pairwise_Ist)

# No normality => Spearman rank's test

# 6.1/ Ist ~ clim avec rho de Spearman
Mantel_Spearman_Ist_clim <- mantel(xdis = pairwise_Ist, ydis = Pairwise_climdist, method = "spearman", na.rm = T)
Mantel_Spearman_Ist_clim
save(Mantel_Spearman_Ist_clim , file = "./outputs/Niche_evolution/Mantel_tests/Resultats_Mantel_Ist_clim_Spearman.RData")

# 6.2/ Ist ~ geo avec rho de Spearman
Mantel_Spearman_Ist_geo <- mantel(xdis = pairwise_Ist, ydis = Pairwise_geodist, method = "spearman", na.rm = T)
Mantel_Spearman_Ist_geo
save(Mantel_Spearman_Ist_geo, file = "./outputs/Niche_evolution/Mantel_tests/Resultats_Mantel_Ist_geo_Spearman.RData")

# 6.3/ Mantel Partial Test : Ist ~ Dclim + cov(Dgeo)

partial_Mantel_Spearman_Ist_clim_geo <- mantel.partial(xdis = pairwise_Ist, ydis = Pairwise_climdist, zdis = Pairwise_geodist, method = "spearman", na.rm = T)
partial_Mantel_Spearman_Ist_clim_geo
save(partial_Mantel_Spearman_Ist_clim_geo, file = "./outputs/Niche_evolution/Mantel_tests/Resultats_Mantel_partiel_Ist_Spearman.RData")


# Plot de la distri de la stat
load(file = paste0(internal.wd, "/Community_Structure/Resultats_Mantel_partiel_Ist_Spearman.RData"))
tiff(filename = paste0(internal.wd,"/Community_Structure/Plots/Distri_Mantel_Spearman.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
hist(c(Mantel_Spearman$perm, Mantel_Spearman$statistic), 30, freq = TRUE, col = "gray", 
     main = "Mantel test statistic distribution\n under 999 permutations", xlab = "Spearman's rho",
     cex.axis = 1.7, cex.main = 1, cex.lab = 1.7)
arrows(Mantel_Spearman$statistic + 0.0051, 80, Mantel_Spearman$statistic+ 0.0051, 10, length = 0.1, lwd = 2)
abline(v = mean(c(Mantel_Spearman$perm, Mantel_Spearman$statistic)), lty = 2, lwd = 2)
abline(v = quantile(c(Mantel_Spearman$perm, Mantel_Spearman$statistic), 0.95), lty = 2, lwd = 2, col = "red")
legend(legend = c(paste0("Mean = ", round(mean(c(Mantel_Spearman$perm, Mantel_Spearman$statistic)),3)), 
                  paste0("CI 95% = ", round(quantile(c(Mantel_Spearman$perm, Mantel_Spearman$statistic), 0.95),3))), 
       x = "topright", lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")
legend(legend = c(paste0("Rho obs = ", round(Mantel_Spearman$statistic, 3)),
                  paste0("                  p = ", 1-(ecdf(x = c(Mantel_Spearman$perm, Mantel_Spearman$statistic))(Mantel_Spearman$statistic)))),
       x = "right", cex = 1.2, bty ="n", xjust = 1)
par(mar = original_int_margins)
dev.off()


## Plot de la relation avec 1000 points tir?s au hasard parmi les n(n-1)/2 valeurs de distances

load(file = "./outputs/Niche_evolution/Mantel_tests/Resultats_Mantel_Ist_clim_Spearman.RData")
load(file = "./outputs/Niche_evolution/Mantel_tests/Resultats_Mantel_Ist_geo_Spearman.RData")
load(file = "./outputs/Niche_evolution/Mantel_tests/Resultats_Mantel_partiel_Ist_Spearman.RData")

# Extraction des vecteurs
pairwise_Ist_vec <- as.dist(pairwise_Ist)[which(!is.na(as.dist(pairwise_Ist)))]
pairwise_climdist_vec <- as.dist(Pairwise_climdist)[which(!is.na(as.dist(pairwise_Ist)))]
pairwise_geodist_vec <- as.dist(Pairwise_geodist)[which(!is.na(as.dist(pairwise_Ist)))]

# Ist ~ Dclim
reg.clim <- lm(pairwise_Ist_vec~pairwise_climdist_vec)
summary(reg.clim)
random.sample <- sample(x= 1: length(pairwise_Ist_vec), size = 1000, replace = F)
tiff(filename = paste0(internal.wd,"/Community_Structure/Plots/Ist_Dclim.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
plot(pairwise_climdist_vec[random.sample], pairwise_Ist_vec[random.sample], 
     main = "Climate structures \n mimicry community composition", ylab = "Pairwise Ist", xlab = "Standardized climatic distance",
     cex.axis = 1.7, cex.main = 1, axis.args=list(cex.axis=1.7), cex.lab = 1.7)
abline(reg.clim, lwd = 2, col = "red")
# text(x = 7, y = 0.4, labels = paste0("Rho = ", round(Mantel_Spearman_Ist_clim$statistic,3), "\n p < 0.001 "), col = "red")
par(mar = original_int_margins)
dev.off()

# Ist ~ Dgeo
reg.geo <- lm(pairwise_Ist_vec~pairwise_geodist_vec)
summary(reg.geo)
random.sample <- sample(x= 1: length(pairwise_Ist_vec), size = 1000, replace = F)
tiff(filename = paste0(internal.wd,"/Community_Structure/Plots/Ist_Dgeo.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
plot(pairwise_geodist_vec[random.sample], pairwise_Ist_vec[random.sample],
     main = "Space structures \n mimicry community composition", ylab = "Pairwise Ist", xlab = "Geographic distance",
     cex.axis = 1.7, cex.main = 1, axis.args=list(cex.axis=1.7), cex.lab = 1.7)
abline(reg.geo, lwd = 2, col = "red")
# text(x = 2000, y = 0.65, labels = paste0("Rho = ", round(Mantel_Spearman$statistic,3), "\n p < 0.001 "), col = "red")
par(mar = original_int_margins)
dev.off()

# Ist ~ Dclim + covar(Dgeo) (r?gression sur les r?sidus)
Ist_resid <- residuals(reg.geo, type = "response")
reg.clim.geo <- lm(Ist_resid~pairwise_climdist_vec)
summary(reg.clim.geo)
random.sample <- sample(x= 1: length(pairwise_Ist_vec), size = 1000, replace = F)
tiff(filename = paste0(internal.wd,"/Community_Structure/Plots/Ist_Dclim_geo.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
plot(pairwise_climdist_vec[random.sample], Ist_resid[random.sample], 
     main = "Climate structures \n mimicry community composition \n even when space is taking into account", ylab = "Residual Pairwise Ist", xlab = "Standardized climatic distance",
     cex.axis = 1.7, cex.main = 1, axis.args=list(cex.axis=1.7), cex.lab = 1.7)
abline(reg.clim.geo, lwd = 2, col = "red")
# text(x = 8, y = 0.3, labels = paste0("Rho = ", round(Mantel_Spearman_Ist_geo$statistic,3), "\n p < 0.001 "), col = "red")
par(mar = original_int_margins)
dev.off()
