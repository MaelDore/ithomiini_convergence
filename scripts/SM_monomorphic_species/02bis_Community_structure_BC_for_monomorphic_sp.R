##### Script 02bis: Community structure with Bray-Curtis indices on comimetic pairs of OMUs #####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

### Focus only on monomorphic/species => 1 OMU = 1 species

# Compute BC indices for all OMUs
# Compute mean BC for mimics only and non-mimics only
# Compute mean BC per mimicry rings
# Generate null distribution with randomization of mimicry patterns
# Plot test results

### Input files

# Summary table of OMUs
# Stack of OMUs probabilities of presence

### Output files

# Matrix of probabilities of presence for OMUs x communities
# BC indices for all pairs of OMUs, only comimics, only non-comimics, and per rings, and associated means
# Null distribution of mean values
# Plot of null distribution and test for all rings, and per rings
# Summary table for each ring of BC values


### Preparation ###

# Effacer l'environnement
rm(list = ls())

library(raster)
library(tidyverse)

### 1/ Load stuff ####

# Load summary table for OMU list
load(file = paste0("./input_data/Summary_tables/list.unit.RData"))
# Load OMU/unit probability stack 
OMU_proba_stack <- readRDS(file = paste0("./input_data/SDM_stacks/All_OMU_stack_Jaccard.80.rds"))

### 2/ Keep only OMU from monomorphic species

# Identify monomorphic the 178 species and associated OMUs
mono_sp_names <- names(table(as.factor(list.unit$Sp_full)))[table(as.factor(list.unit$Sp_full)) == 1]
mono_indices <- which(list.unit$Sp_full %in% mono_sp_names)
saveRDS(object = mono_indices, file = "./outputs/Community_structure/Mono_178/mono_indices.rds")

# Extract list of OMU
list.unit_mono_178 <- list.unit[mono_indices,]
saveRDS(object = list.unit_mono_178, file = "./input_data/Summary_tables/list.unit_mono_178.rds")

# Generate mimicry ring list
mimicry.list_mono_178 <- as.character(unique(list.unit_mono_178$Mimicry.model))
saveRDS(object = mimicry.list_mono_178, file = "./outputs/Community_structure/Mono_178/mimicry.list_mono_178.rds")

# Extract raster layers
OMU_proba_stack_mono_178 <- raster::subset(OMU_proba_stack, mono_indices)
saveRDS(object = OMU_proba_stack_mono_178, file = "./outputs/Community_structure/Mono_178/OMU_proba_stack_mono_178.rds")


### 2/ Extract matrix communities * OMU for probabilities of presence ####

Filtered_com_unit_mat_783 <- readRDS(file = "./input_data/SDM_stacks/Filtered_com_unit_mat_783.rds")
Filtered_com_unit_mat_178 <- Filtered_com_unit_mat_783[, mono_indices]

saveRDS(Filtered_com_unit_mat_178, file = "./input_data/SDM_stacks/Filtered_com_unit_mat_178.rds")

# Load directly the matrix of community * units probability presence
com_unit_mat_178 <- readRDS(file = "./input_data/SDM_stacks/Filtered_com_unit_mat_178.rds")

### 3/ Extract Bray-Curtis observed values ####

### 3.1/ Extract Bray-Curtis index for all pairs of units of monomorphic species ####

library(vegan)

# Load vector of BD distances among all units
unit_BC_dist_783 <- readRDS(file = "./outputs/Community_Structure/unit_BC_dist.rds")

# Convert to matrix and extract data for monomorphic species
unit_BC_dist_mat_783 <- as.matrix(unit_BC_dist_783)
unit_BC_dist_mat_178 <- unit_BC_dist_mat_783[mono_indices, mono_indices]

# Convert back to distances
unit_BC_dist_mono_178 <- as.dist(unit_BC_dist_mat_178)

# Save
saveRDS(unit_BC_dist_mono_178, file = "./outputs/Community_Structure/Mono_178/unit_BC_dist_mono_178.rds")

# Load vector of BD distances among all units of monomorphic species
unit_BC_dist_mono_178 <- readRDS(file = "./outputs/Community_Structure/Mono_178/unit_BC_dist_mono_178.rds")

# hist(unit_BC_dist_mono_178)


### 3.2/ Compute mean BC per mimicry rings ####

mean_BC_mono_178 <- NA
for (i in 1:length(mimicry.list_mono_178)) # Per mimicry rings
{ 
  ring <- mimicry.list_mono_178[i]
  
  # Get names and indices of all units/OMUs for this ring
  tags <- as.character(list.unit_mono_178$Tag.model[list.unit_mono_178$Mimicry.model == ring]) # Get names of all OMUs in the ring
  index <- which(names(OMU_proba_stack_mono_178) %in% tags) # Get index of layer for these OMUs
  
  if (length(tags) == 1) # Case with only one OMU in the mimicry ring. Impossible to compute pairwise distances
  {
    mat.index <- data.frame(matrix(nrow = 2, ncol = 0)) # Empty df of pairs indices (no pairs of OMUs)
    save(mat.index, file = paste0("./outputs/Community_Structure/Mono_178/BC_mimicry_ring_mat/mat.index.",ring,".RData"))
    BC <- c() # Empty vector of Bray-Curtis values
    save(BC, file = paste0("./outputs/Community_Structure/Mono_178/BC_mimicry_ring_mat/BC.",ring,".RData"))
    mean_BC_mono_178[i] <- NA
  } else { # Case with at least 2 OMUs.
    mat.index <- combn(x = index, m = 2, FUN = c) # Get all combinations possibles for pairs of OMUs indices
    save(mat.index, file = paste0("./outputs/Community_Structure/Mono_178/BC_mimicry_ring_mat/mat.index.",ring,".RData"))
    BC <- NA # Initiate the vector used to store all BC values
    for (j in 1:ncol(mat.index)) # For each pair of OMUs
    {
      BC[j] <- unit_BC_dist_mat_178[mat.index[1,j], mat.index[2,j]] # Extract Bray-Curtis index from the complete matrix of BC indices
    }
    # save(BC, file = paste0("./outputs/Community_Structure/BC_mimicry_ring_mat/BC.",ring,".RData")) # Save the vector of BC for all pairs of OMUs in this ring
    mean_BC_mono_178[i]  <- mean(BC) # Compute mean BC for all pairs of this mimicry ring and store it in final vector
  }
  
  cat(paste0(Sys.time(), " - ", ring," - n°",i, " on ",length(mimicry.list_mono_178),"\n"))
  
}
mean_BC_mono_178
names(mean_BC_mono_178) <- mimicry.list_mono_178 # Associate mean BC value with ring name
saveRDS(mean_BC_mono_178, file = paste0("./outputs/Community_Structure/Mono_178/mean_BC_mono_178.rds")) # Save final vector with mean BC values per ring

# Load mean BC of mimicry rings
mean_BC_mono_178 <- readRDS(file = paste0("./outputs/Community_Structure/Mono_178/mean_BC_mono_178.rds"))
mean_BC_mono_178

### 3.3/ Compute mean BC for all comimics, all rings taken into account ####

# Retrieve all co-mimic coordinates

all_mimic_mat_index_mono_178 <- data.frame(matrix(ncol = 0, nrow = 2)) # Generate empty df to store indices of pairs of comimics
for (i in 1:length(mimicry.list_mono_178))  # Per mimetic ring
{ 
  ring <- mimicry.list_mono_178[i]
  load(file = paste0("./outputs/Community_Structure/Mono_178/BC_mimicry_ring_mat/mat.index.",ring,".RData")) # Load the matrix of indices of pairs of comimics for this ring
  all_mimic_mat_index_mono_178 <- cbind(all_mimic_mat_index_mono_178,mat.index) # Merge them all in one df
}
dim(all_mimic_mat_index_mono_178) # 1,281 pairs of co-mimic OMUs in monomorphic species

saveRDS(all_mimic_mat_index_mono_178, file = paste0("./outputs/Community_Structure/Mono_178/BC_mimicry_ring_mat/all_mimic_mat_index_mono_178.rds"))

# Retrieve all co-mimic BC values and non-co-mimic BC values at the same time
BC_mimic_mono_178 <- NA # Initiate final vector to store BC values for comimics
unit_BC_dist_mat_178_no_mimic <- unit_BC_dist_mat_178 # Copy matrix of all BC indices for pairs of OMUs in monomorphic species
for (j in 1:ncol(all_mimic_mat_index_mono_178)) # For all pairs of comimics
{
  BC_mimic_mono_178[j] <- unit_BC_dist_mat_178[all_mimic_mat_index_mono_178[1,j],all_mimic_mat_index_mono_178[2,j]] # Extract BC index for this pair of comimics
  unit_BC_dist_mat_178_no_mimic[all_mimic_mat_index_mono_178[1,j],all_mimic_mat_index_mono_178[2,j]] <- NA # Remove BC value of comimics from matrix of non-comimics pairs
  if (j %% 100 == 0) {print(j)}
}
length(BC_mimic_mono_178) # 1,281 pairs of co-mimic OMUs

saveRDS(BC_mimic_mono_178, file = paste0("./outputs/Community_Structure/Mono_178/BC_mimic_mono_178.rds")) # Save final vector with all BC values for comimics

# Retrieve all non-mimic BC values
BC_no_mimic_mono_178 <- unit_BC_dist_mat_178_no_mimic[upper.tri(unit_BC_dist_mat_178_no_mimic)] # Extract only one side of the triangle to avoid duplicate values of pairs
BC_no_mimic_mono_178 <- na.omit(BC_no_mimic_mono_178) # Remove NA (the comimics pairs)
length(BC_no_mimic_mono_178) # 14,472 pairs of non-mimic units

saveRDS(BC_no_mimic_mono_178, file = paste0("./outputs/Community_Structure/Mono_178/BC_no_mimic_mono_178.rds")) # Save final vector with all BC values for non-comimics


# Compute mean obs BC for each group
load(file = "./outputs/Community_Structure/Mono_178/unit_BC_dist_mono_178.RData") # All BC index (distance format)
load(file = "./outputs/Community_Structure/Mono_178/BC_mimic_mono_178.RData") # Only mimic pairs (vector format)
load(file = "./outputs/Community_Structure/Mono_178/BC_no_mimic_mono_178.RData") # Only non mimic pairs (vector format)

Global_mean_BC_mono_178 <- mean(unit_BC_dist_mono_178) ; Global_mean_BC_mono_178 # 0.953
Global_mimic_mean_BC_mono_178 <- mean(BC_mimic_mono_178) ; Global_mimic_mean_BC_mono_178 # 0.916
Global_no.mimic_mean_BC_mono_178 <- mean(BC_no_mimic_mono_178) ; Global_no.mimic_mean_BC_mono_178 # 0.957

# Save all mean BC values for all pairs, only co-mimics, only non-comimics
save(Global_mean_BC_mono_178, Global_mimic_mean_BC_mono_178, Global_no.mimic_mean_BC_mono_178, file = "./outputs/Community_Structure/Mono_178/All_Global_BC_mono_178.RData")

# Quick ugly boxplots
All_BC <- c(BC_mimic_mono_178, BC_no_mimic_mono_178)
Status <- as.factor(c(rep("Mimic", length(BC_mimic_mono_178)), rep("Non-Mimic", length(BC_no_mimic_mono_178))))
boxplot(All_BC ~ Status)

library(ggplot2)

boxplot_df <- data.frame(All_BC, Status, stringsAsFactors = T)
gg_boxplot <- ggplot(data = boxplot_df, aes(x = Status, y = All_BC)) +
  geom_violin() +
  coord_trans(y = "exp")
print(gg_boxplot)

# Histogram
# Plot distri of BC for mimic pairs
pdf(file = paste0("./graphs/Community_structure/Mono_178/hist_BC_mimic_mono_178.pdf"), height = 6.3, width = 6.5)
hist(BC_mimic_mono_178, xlab = "Bray-Curtis index", main = "Bray-Curtis indices of comimic pairs")
abline(v = mean(BC_mimic_mono_178), col = "red", lty = 2, lwd = 2)
legend(legend = c(paste0("Mean = ", round(mean(BC_mimic_mono_178, na.rm = T),3)), 
                  paste0("CI 5% = ", round(quantile(BC_mimic_mono_178, 0.05),3)),
                  paste0("CI 95% = ", round(quantile(BC_mimic_mono_178, 0.95),3))),
       x = "topleft", cex = 1, bty ="n") 
dev.off()


##### 4/ Generate new virtual community matrices under null hypothesis of no effect of mimicry ring on OMUs presence ####
# analogous to DeVries et al.'s [1999] and Hill's [2010] test of mimicry structure across microhabitats

# Load directly the vector of BC distances among all units
unit_BC_dist_mono_178 <- readRDS(file = "./outputs/Community_Structure/Mono_178/unit_BC_dist_mono_178.rds")
OMU_proba_stack_mono_178 <- readRDS(file = "./outputs/Community_structure/Mono_178/OMU_proba_stack_mono_178.rds")

unit_BC_dist_mat_178 <- as.matrix(unit_BC_dist_mono_178)

mimicry.list_mono_178 <- as.character(unique(list.unit_mono_178$Mimicry.model))


## Start the loop for permutations
BC_mimic_mono_178_null <- BC_no_mimic_mono_178_null <- NA # Create vectors to store simulated mean BC values for co-mimics and non-comimics
mean_BC_null_mono_178 <- matrix(ncol = length(mimicry.list_mono_178), nrow = 0) # Create matrix to store mean simulated BC values for each mimicry ring for each simulation
for (k in 1:1000) # 1000 simulations/permutations
{ 
  # k <- 1
  
  ## Suffle mimicry ring among units
  shuffle.list.unit <- list.unit_mono_178
  shuffle.list.unit$Mimicry.model <- sample(as.character(shuffle.list.unit$Mimicry.model))
  
  # # Check if number of unit per ring is preserved
  # table(list.unit_mono_178$Mimicry.model)
  # table(shuffle.list.unit$Mimicry.model)
  
  ## Generate the new Mimicry rings Richness Stack with random attribution of OMUs to mimicry ring
  
  # Mimicry list
  mimicry.list_mono_178 <- as.character(unique(list.unit_mono_178$Mimicry.model)) # 30 Mimicry rings
  
  mean_BC <- NA
  for (i in 1:length(mimicry.list_mono_178)) # Per mimicry rings
  { 
    # i <-  1
    ring <- mimicry.list_mono_178[i]
    
    # Get names and indices of all units/OMUs for this ring
    tags <- as.character(list.unit_mono_178$Tag.model[shuffle.list.unit$Mimicry.model == ring])
    index <- which(names(OMU_proba_stack_mono_178) %in% tags)
    
    if (length(tags)==1) # Case with only one OMU in the mimicry ring. Impossible to compute pairwise distances
    { 
      mat.index <- data.frame(matrix(nrow = 2, ncol = 0)) # Empty df of pairs coordinates
      save(mat.index, file = paste0("./outputs/Community_Structure/Mono_178/Permutations/BC_mimicry_ring_mat/simul_mat.index.",ring,".RData")) # Save temp mat indices df to be used later to merge for all rings
      mean_BC[i] <- NA
    } else { # Case with at least 2 OMUs
      mat.index <- combn(x = index, m = 2, FUN = c) # Get all combinations possibles for pairs of OMUs indices
      save(mat.index, file = paste0("./outputs/Community_Structure/Mono_178/Permutations/BC_mimicry_ring_mat/simul_mat.index.",ring,".RData")) # Save temp mat indices df to be used later to merge for all rings
      BC <- NA # Initiate the vector used to store all BC values
      for(j in 1:ncol(mat.index)) # For each pair of OMUs
      {
        BC[j] <- unit_BC_dist_mat_178[mat.index[1,j],mat.index[2,j]] # Extract Bray-Curtis index from the complete matrix of BC indices
      }
      mean_BC[i]  <- mean(BC) # Compute mean BC for all pairs of this mimicry ring and store it in final vector
    }
    
    # print(i)
  }
  names(mean_BC) <- mimicry.list_mono_178
  mean_BC_null_mono_178 <- rbind(mean_BC_null_mono_178, mean_BC) # Store the mean BC vector into the matrix of simulations, as a row
  
  # Retrieve all mimic coordinates
  all_mimic_mat_index_mono_178 <- data.frame(matrix(ncol = 0, nrow = 2)) # Generate empty df to store indices of pairs of comimics
  for (i in 1:length(mimicry.list_mono_178)) { # Per mimicry rings
    ring <- mimicry.list_mono_178[i]
    load(file =  paste0("./outputs/Community_Structure/Mono_178/Permutations/BC_mimicry_ring_mat/simul_mat.index.",ring,".RData")) # Load the matrix of indices of pairs of comimics for this ring
    all_mimic_mat_index_mono_178 <- cbind(all_mimic_mat_index_mono_178, mat.index) # Merge them all in one df
  }
  # dim(all_mimic_mat_index_mono_178) # 14 659 pairs of mimic units
  
  # Retrieve all mimic values and non-co-mimic BC values at the same time
  BC_mimic_mono_178 <- NA # Initiate final vector to store BC values for comimics
  unit_BC_dist_mat_178_no_mimic <- unit_BC_dist_mat_178 # Copy matrix of all BC indices for all pairs of OMUs
  for(j in 1:ncol(all_mimic_mat_index_mono_178)) # For all pairs of comimics
  {
    BC_mimic_mono_178[j] <- unit_BC_dist_mat_178[all_mimic_mat_index_mono_178[1,j], all_mimic_mat_index_mono_178[2,j]] # Extract BC index for this pair of comimics
    unit_BC_dist_mat_178_no_mimic[all_mimic_mat_index_mono_178[1,j], all_mimic_mat_index_mono_178[2,j]] <- NA # Remove BC value of comimics from matrix of non-comimics pairs
  }
  # length(BC_mimic_mono_178) # 14 659 pairs of mimic units
  
  # Retrieve all non-mimic values
  BC_no_mimic_mono_178 <- unit_BC_dist_mat_178_no_mimic[upper.tri(unit_BC_dist_mat_178_no_mimic)] # Extract only one side of the triangle to avoid duplicate values of pairs
  BC_no_mimic_mono_178 <- na.omit(BC_no_mimic_mono_178) # Remove NA (the comimics pairs)
  # length(BC_no_mimic_mono_178) # 291 494 pairs of non-mimic units
  
  # Save global mean computation into final vectors
  BC_mimic_mono_178_null[k] <- mean(BC_mimic_mono_178)
  BC_no_mimic_mono_178_null[k] <- mean(BC_no_mimic_mono_178)
  
  cat(paste0(Sys.time(), " - Simul n°", k," out of 1000\n"))
  save(mean_BC_null_mono_178, BC_mimic_mono_178_null, BC_no_mimic_mono_178_null, file = paste0("./outputs/Community_Structure/Mono_178/Permutations/All_simul_mean_BC_mono_178.RData"))
}

summary(mean_BC_null_mono_178)     # Per mimicry ring
summary(BC_mimic_mono_178_null)    # For all comimic pairs
summary(BC_no_mimic_mono_178_null) # For all non-comimic pairs

##### 5/ Plot distri of null mean BC ####

### 5.1/ Plot for all comimics ####

load(file = "./outputs/Community_Structure/Mono_178/All_Global_BC_mono_178.RData") # Load mean BC obs values 
load(file = "./outputs/Community_Structure/Mono_178/Permutations/All_simul_mean_BC_mono_178.RData") # Load null mean BC values

pdf(file = "./graphs/Community_Structure/Mono_178/BC_mimic_mono_178_null.pdf", height = 5.3, width = 6.5)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

hist(x = c(BC_mimic_mono_178_null, Global_mimic_mean_BC_mono_178), breaks = 40, freq = TRUE, col = "gray",
     xlim = c(0.91, 0.97),
     ylim = c(0, 100),
     main = "Distribution of mean pairwise BC indices\nof co-mimetic monomorphic species\nunder null Hypothesis",
     xlab = "Mean Bray-Curtis Index",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)

arrows(Global_mimic_mean_BC_mono_178 + 0.0005, 40, Global_mimic_mean_BC_mono_178 + 0.0005, 5, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = mean(c(BC_mimic_mono_178_null, Global_mimic_mean_BC_mono_178)), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(c(BC_mimic_mono_178_null, Global_mimic_mean_BC_mono_178), 0.05), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("Mean = ", format(round(mean(c(BC_mimic_mono_178_null, Global_mimic_mean_BC_mono_178)),3), nsmall = 3)), 
                  paste0("CI 5% = ", round(quantile(c(BC_mimic_mono_178_null, Global_mimic_mean_BC_mono_178), 0.05),3))), 
       x = "topleft", inset = c(0, 0.17), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty = "n")
legend(legend = c(paste0("Mean BC obs = ", round(Global_mimic_mean_BC_mono_178, 3)),
                  paste0("p = 0.001")),
       x = "left", inset = c(0.00, -0.05),
       cex = 1.2, bty ="n", xjust = 1)

legend(legend = as.expression(bquote(bold("B"))), 
       x = "topleft", inset = c(-0.03, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()



### 5.2/ Plot distri of null BC per mimicry ring co-mimics ####

mean_BC_mono_178 <- readRDS(file = paste0("./outputs/Community_Structure/Mono_178/mean_BC_mono_178.rds")) # Load mean BC obs value for each ring
load(file = "./outputs/Community_Structure/Mono_178/Permutations/All_simul_mean_BC_mono_178.RData") # Load null mean BC values

# Build summary table for mimicry ring BC results at the same time
BC_ring_summary_table_mono_178 <- as.data.frame(matrix(ncol = 8, nrow = length(mimicry.list_mono_178), data = NA))
names(BC_ring_summary_table_mono_178) <- c("ring", "N_units", "N_pairs", "BC_obs", "mean_BC", "BC_2.5", "BC_97.5", "p_value")

for (i in 1:length(mimicry.list_mono_178))  # Per mimicry ring
{
  # i <- 1
  
  BC_ring_summary_table_mono_178$ring[i] <- ring <- mimicry.list_mono_178[i] # Get ring name
  BC_ring_summary_table_mono_178$N_units[i] <- N_units <- sum(list.unit_mono_178$Mimicry.model == ring) # Get number of OMUs
  BC_ring_summary_table_mono_178$N_pairs[i] <- N_pairs <- N_units*(N_units-1)/2 # Get number of pairs
  
  if (is.na(mean_BC_mono_178[i]))# Case for ring with no pair of OMUs (only one OMU/species)
  {
    # Plot BC null distri for this ring (i.e. NA message)
    pdf(file = paste0("./graphs/Community_Structure/Mono_178/Per_mimicry_ring/BC_mimic_mono_178_null_",ring,".pdf"), height = 5.3, width = 6.5)
    
    plot(1:100,1:100, type = "n", main = paste0("Distribution of mean pairwise Bray-Curtis indices\nof ", ring, " OMUs\nunder null Hypothesis"))
    text(x = 50, y = 50, labels = "Only one OMU for this mimicry ring \n No pair available for index computation")
    
    dev.off()
    
  } else { # Case for ring with pair(s) of OMUs 
    
    BC_ring_summary_table_mono_178$BC_obs[i] <- round(mean_BC_mono_178[i],3) # Get mean BC obs
    BC_ring_summary_table_mono_178$mean_BC[i] <- round(mean(c(mean_BC_null_mono_178[,i], mean_BC_mono_178[i]), na.rm = T),3) # Get mean of mean BC null from simulations
    BC_ring_summary_table_mono_178$BC_2.5[i] <- round(quantile(c(mean_BC_null_mono_178[,i], mean_BC_mono_178[i]), 0.025),3) # Get 2.5% quantile of mean BC null from simulations
    BC_ring_summary_table_mono_178$BC_97.5[i] <- round(quantile(c(mean_BC_null_mono_178[,i], mean_BC_mono_178[i]), 0.975),3) # Get 97.5% quantile of mean BC null from simulations
    BC_ring_summary_table_mono_178$p_value[i] <- round(ecdf(x = c(mean_BC_null_mono_178[,i], mean_BC_mono_178[i]))(mean_BC_mono_178[i]),3) # Get p-value from simulations
    
    save(BC_ring_summary_table_mono_178, file = "./outputs/Community_Structure/BC_ring_summary_table_mono_178.Rdata")
    
    # Plot BC null distri for this ring
    pdf(file = paste0("./graphs/Community_Structure/Mono_178/Per_mimicry_ring/BC_mimic_mono_178_null_",ring,".pdf"), height = 5.3, width = 6.5)
    
    hist(c(mean_BC_null_mono_178[,i], mean_BC_mono_178[i]), 
         breaks = seq(from = floor(floor(min(min(mean_BC_null_mono_178[,i]), mean_BC_mono_178[i])*1000)/2)/500, to = ceiling(ceiling(max(max(mean_BC_null_mono_178[,i]), mean_BC_mono_178[i])*1000)/2)/500, by = 0.002), 
         xlab = "Mean pairwise Bray-Curtis index", 
         main = paste0("Distribution of mean pairwise Bray-Curtis indices\nof ", ring, " OMUs\nunder null Hypothesis"))
    
    # arrows(mean_BC_mono_178[i] + 0.0001, 60, mean_BC_mono_178[i] + 0.0001, 10, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
    abline(v = mean_BC_mono_178[i], lwd = 2, lty = 2, col = "red") # Add vertical line for mean BC obs value
    
    abline(v = mean(c(mean_BC_null_mono_178[,i], mean_BC_mono_178[i]), na.rm = T), lwd = 2, lty = 2) # Add vertical line for mean value
    abline(v = quantile(c(mean_BC_null_mono_178[,i], mean_BC_mono_178[i]), 0.025, na.rm = T), lwd = 2, lty = 2, col = "grey30") # Add vertical line for 2.5% value
    abline(v = quantile(c(mean_BC_null_mono_178[,i], mean_BC_mono_178[i]), 0.975, na.rm = T), lwd = 2, lty = 2, col = "grey30") # Add vertical line for 97.5% value
    
    legend(legend = c(paste0("N units = ", N_units), 
                      paste0("N pairs = ", N_pairs)),
           x = "topright", xjust = 1, cex = 1, bty ="n")
    
    legend(legend = c(paste0("Mean = ", round(mean(c(mean_BC_null_mono_178[,i], mean_BC_mono_178[i]), na.rm = T),3)), 
                      paste0("CI 2.5% = ", round(quantile(c(mean_BC_null_mono_178[,i], mean_BC_mono_178[i]), 0.025),3)),
                      paste0("CI 97.5% = ", round(quantile(c(mean_BC_null_mono_178[,i], mean_BC_mono_178[i]), 0.975),3))),
           x = "topleft", cex = 1, bty ="n") 
    
    legend(legend = paste0("Mean BC obs = ", round(mean_BC_mono_178[i], 3), "\n p = ", round(ecdf(x = mean_BC_null_mono_178[,i])(mean_BC_mono_178[i]),3)), 
           inset = c(0, 0.2), text.col = "red", x = "bottomleft", cex = 1, bty ="n")
    
    dev.off()
    
    
  }
  
  cat(paste0(i, " out of ", length(mimicry.list_mono_178), " - ",ring, " - Done \n"))
}

### 6/ Export summary tables ####

### 6.1/ Export table for monomorphic species

# Sort in alphabetic order
BC_ring_summary_table_mono_178 <- arrange(BC_ring_summary_table_mono_178, ring)

View(BC_ring_summary_table_mono_178)

# Save
saveRDS(BC_ring_summary_table_mono_178, file = "./tables/BC_ring_summary_table_mono_178.rds")
write.csv2(BC_ring_summary_table_mono_178, file = "./tables/BC_ring_summary_table_mono_178.csv")

### 6.2/ Add entries for mimicry rings absent from monomorphic species

BC_ring_summary_table_all_OMU <- read.csv2(file = "./tables/BC_ring_summary_table.csv", row.names = 1)
BC_ring_summary_table_all_OMU <- arrange(BC_ring_summary_table_all_OMU, ring)

all_rings <-  BC_ring_summary_table_all_OMU[, "ring", F]

BC_ring_summary_table_mono_178_with_null <- left_join(all_rings, BC_ring_summary_table_mono_178)
BC_ring_summary_table_mono_178_with_null$N_units[is.na(BC_ring_summary_table_mono_178_with_null$N_units)] <- 0
BC_ring_summary_table_mono_178_with_null$N_pairs[is.na(BC_ring_summary_table_mono_178_with_null$N_pairs)] <- 0

# Save
saveRDS(BC_ring_summary_table_mono_178, file = "./tables/BC_ring_summary_table_mono_178.rds")
write.csv2(BC_ring_summary_table_mono_178, file = "./tables/BC_ring_summary_table_mono_178.csv")

### 6.3/ Fuse with the complete table for all OMUs

BC_ring_summary_table_both <- BC_ring_summary_table_mono_178
names(BC_ring_summary_table_both) <- c("ring", "N_units_mono", "N_pairs_mono", "BC_obs_mono", "mean_BC_mono", "BC_2.5_mono", "BC_97.5_mono", "p_value_mono")

BC_ring_summary_table_both <- left_join(BC_ring_summary_table_all_OMU, BC_ring_summary_table_both)

# Save
saveRDS(BC_ring_summary_table_both, file = "./tables/BC_ring_summary_table_both.rds")
write.csv2(BC_ring_summary_table_both, file = "./tables/BC_ring_summary_table_both.csv")
