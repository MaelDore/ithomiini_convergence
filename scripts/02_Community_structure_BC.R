##### Script 02: Community structure with Bray-Curtis indices on comimetic pairs of OMUs #####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

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

# Final plot with null distri for Ist and mean BC for comimics


### Preparation ###

# Effacer l'environnement
rm(list = ls())

library(raster)

### 1/ Load stuff ####

# Load summary table for OMU list
load(file = paste0("./input_data/Summary_tables/list.unit.RData"))
# Load OMU/unit probability stack 
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


##### 4/ Generate new virtual community matrices under null hypothesis of no effect of mimicry ring on OMUs presence ####
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
     ylim = c(0, 400),
     main = "Distribution of mean pairwise BC indices\nof co-mimetic OMUs\nunder null Hypothesis",
     xlab = "Mean Bray-Curtis Index",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)

arrows(Global_mimic_mean_BC + 0.0007, 150, Global_mimic_mean_BC+ 0.0007, 10, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = mean(c(BC_mimic_null, Global_mimic_mean_BC)), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(c(BC_mimic_null, Global_mimic_mean_BC), 0.05), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("Mean = ", format(round(mean(c(BC_mimic_null, Global_mimic_mean_BC)),3), nsmall = 3)), 
                  paste0("CI 5% = ", round(quantile(c(BC_mimic_null, Global_mimic_mean_BC), 0.05),3))), 
       x = "topleft", inset = c(0, 0.17), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty = "n")
legend(legend = c(paste0("Mean BC obs = ", round(Global_mimic_mean_BC, 3)),
                  paste0("p = 0.001")),
       x = "left", inset = c(0.00, -0.05),
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
write.csv2(BC_ring_summary_table, file = "./tables/BC_ring_summary_table.csv")

##### 6/ Final plot with both Ist and BC #####

# Remove titles

load(file = "./outputs/Community_Structure/Ist_null.RData")
load(file = "./outputs//Community_Structure/Global.Ist.RData")

load(file = "./outputs/Community_Structure/All_Global_BC.RData") # Load mean BC obs values 
load(file = "./outputs/Community_Structure/Permutations/All_simul_mean_BC.RData") # Load null mean BC values

round(Global_mimic_mean_BC, 3) # Mean BC obs = 0.896
round(Global_Ist, 3) # Ist obs = 0.164

# Final plot

pdf(file = "./graphs/Community_Structure/Ist_&_BC.pdf", height = 5.5, width = 10)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,5,1,1), mfrow = c(1,2))

# Panel A = Ist

hist(c(Ist_null, Global_Ist), 30, freq = TRUE, col = "gray", xlab = bquote(~I[ST]), 
     main = "",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5, lwd = 2)
arrows(Global_Ist - 0.001, 175, Global_Ist- 0.001, 10, length = 0.1, lwd = 2) # Draw arrow above Ist obs
abline(v = mean(c(Ist_null, Global_Ist)), lty = 2, lwd = 2) # Add vertical line for mean value
abline(v = quantile(c(Ist_null, Global_Ist), 0.95), lty = 2, lwd = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("   Mean = ", format(round(mean(c(Ist_null, Global_Ist)),3), nsmall = 3)), 
                  paste0("CI 95% = ", round(quantile(c(Ist_null, Global_Ist), 0.95),3))), 
       x = "topright", inset = c(0, 0.17), xjust = 1, # Use inset to manually adjust position
       cex = 1.2, # text.width = c(0.9, 1),
       lty = 2, lwd = 2, col = c("black", "red"), bty ="n")

# legend(legend = bquote("           " ~ I[ST] ~ 'obs =' ~ .(round(Global_Ist, 3))),
#        x = "right", cex = 1.2, bty ="n", xjust = 1)
# legend(legend = c("",
#                   paste0("                   p = 0.001")),
#        x = "right", cex = 1.2, bty ="n", xjust = 1, y.intersp = 1.8)

# legend(legend = c(as.expression(bquote(bold(paste(~ I[ST] ~ " obs = 0.164")))),
#                   as.expression(bquote(bold(paste("           p = 0.001"))))),
#        x = "bottomright", inset = c(0.02, 0.40), xjust = 1, # Use inset to manually adjust position
#        cex = 1.2, bty ="n", bg = "white", box.col = NA)

legend(legend = as.expression(bquote(bold(paste(~ I[ST] ~ " obs = 0.164")))),
       x = "bottomright", inset = c(0.02, 0.45), xjust = 1, # Use inset to manually adjust position
       cex = 1.2, bty ="n", bg = "white", box.col = NA)
legend(legend = as.expression(bquote(bold(paste("p = 0.001")))),
       x = "bottomright", inset = c(0.02, 0.40), xjust = 1, # Use inset to manually adjust position
       cex = 1.2, bty ="n", bg = "white", box.col = NA)

legend(legend = as.expression(bquote(bold("A"))), 
       x = "topright", inset = c(0.05, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

# Panel B = BC Indices

hist(x = c(BC_mimic_null, Global_mimic_mean_BC), 
     breaks = seq(0.89, 0.96, 1/300),
     freq = TRUE, col = "gray",
     xlim = c(0.89, 0.96),
     ylim = c(0, 500),
     main = "",
     xlab = "Mean Bray-Curtis Index",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)

arrows(Global_mimic_mean_BC - 0.0012, 195, Global_mimic_mean_BC - 0.0012, 10, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = mean(c(BC_mimic_null, Global_mimic_mean_BC)), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(c(BC_mimic_null, Global_mimic_mean_BC), 0.05), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("Mean = ", format(round(mean(c(BC_mimic_null, Global_mimic_mean_BC)),3), nsmall = 3)),
                  paste0("CI 5% = ", round(quantile(c(BC_mimic_null, Global_mimic_mean_BC), 0.05),3))), 
       x = "topleft", inset = c(0.05, 0.17), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty = "n")

# legend(legend = c(paste0("Mean BC obs = ", round(Global_mimic_mean_BC, 3)),
#                   paste0("p = 0.001")),
#        x = "bottomleft", inset = c(-0.02, 0.40),
#        cex = 1.2, bty ="n", xjust = 1)

legend(legend = c(as.expression(bquote(bold(paste("Mean BC obs = 0.896")))),
                  as.expression(bquote(bold(paste("p = 0.001"))))),
       x = "bottomleft", inset = c(-0.02, 0.40),
       cex = 1.2, bty ="n", xjust = 1)

legend(legend = as.expression(bquote(bold("B"))), 
       x = "topleft", inset = c(-0.03, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins, mfrow = c(1,1))

dev.off()


##### 7/ Illustration of BC overlap between OMU #####

# Select the OMU to illustrate

load(file = paste0("./input_data/Summary_tables/list.unit.RData"))

unit_BC_dist <- readRDS(file = "./outputs/Community_Structure/unit_BC_dist.rds")
unit_BC_dist_mat <- as.matrix(unit_BC_dist)

table(unit_BC_dist_mat[which(list.unit$Tag.model == "Dircenna.dero.DILUCIDA"), ])
table(unit_BC_dist_mat[which(list.unit$Tag.model == "Dircenna.jemina.DILUCIDA"), ])

list.unit$Tag.model[which((unit_BC_dist_mat[60, ] > 0.93) & (unit_BC_dist_mat[60, ] < 0.97) & (list.unit$initial_model_type == "complete"))]

which(list.unit$Tag.model == "Dircenna.jemina.DILUCIDA")
which(list.unit$Tag.model == "Dircenna.dero.DILUCIDA")
which(list.unit$Tag.model == "Oleria.amalda.LERIDA")
which(list.unit$Tag.model == "Mechanitis.mazaeus.MAELUS")

unit_BC_dist_mat[60, 58]
unit_BC_dist_mat[60, 499]
unit_BC_dist_mat[60, 372]

# Load OMU/unit probability stack 
OMU_proba_stack <- readRDS(file = paste0("./input_data/SDM_stacks/All_OMU_stack_Jaccard.80.rds"))

# Load packages
library(raster)
library(rgdal)

# Load stuff for maps
pal_bl_red_Mannion <- readRDS(file = "./input_data/Map_stuff/pal_bl_red_Mannion.rds")

grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds")
country_borders_Mollweide <- readRDS(file = "./input_data/Map_stuff/country_borders_Mollweide.rds")

# Function to project raster into Mollweide
Mollweide_projection <- function(x) # Raster to project
{
  new_map <- projectRaster(from = x, 
                           method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                           crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                           alignOnly = F)
  return(new_map)
}

# Function to retrieve type of model and parse name
extract_OMU_name_and_pattern <- function(unit)
{
  load(file = "./input_data/list.models.RData")
  load(file = "./input_data/list.rings.RData")
  
  ID <- which(list.models$Tag.model == unit)
  species_name <- paste(list.models$Genus_new[ID], list.models$Species_new[ID])
  
  # Update pattern name
  pattern <- as.character(list.models$Mimicry.model[ID])
  pattern <- list.rings$Mimicry_new[which(list.rings$Mimicry == pattern)]
  pattern <- paste0("(", pattern, ")")

  return(c(species_name, pattern))
}


# Function to plot a unit map
{
  map_unit_Mollweide <- function(x,                                       # Raster to map
                                 unit_pts,                                # Shp with data points of ocurrences
                                 color_palette = pal_bl_red_Mannion,   # Color palette
                                 main_title,                           # Main title
                                 main_title_cex = 1.5,                 # Main title size
                                 
                                 xlim = c(-4600, 4600),    # Limit of plot on x-axis (Longitude)
                                 ylim = c(-4450, 3400),    # Limit of plot on y-axis (Latitude)
                                 axis_cex = 1.4,           # Axes size
                                 
                                 xlab = "",                # X-axis label
                                 ylab = "",                # Y-axis label
                                 x_axis_breaks = c(-3930, -2170, -420, 1500, 3050),            # X-axis tick breaks
                                 y_axis_breaks = c(-3650, -2450, -1220, 0, 1230, 2445),        # Y-axis tick breaks
                                 x_axis_labels = c("120°W", "100°W", "80°W", "60°W", "40°W"),      # X-axis tick labels
                                 y_axis_labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"),  # Y-axis tick labels
                                 
                                 occ_cex = 0.5,         # Size of occurrences
                                 occ_col = "#00000080",  # Color of occurrences
                                 
                                 legend_title,             # Legend title
                                 legend_title_cex = 1.4,   # Legend title size
                                 legend_title_x = -3550,   # Legend title x position
                                 legend_title_y = 430,     # Legend title y position
                                 legend_cex = 1.4,         # Legend size
                                 legend_breaks,            # Legend tick positions
                                 legend_location = c(-4100, -3800, -3950, 0),  # Legend position
                                 
                                 scale_bar_position = c(-2600, -4000),  # Scale bar position
                                 
                                 arrow_scale = 0.45,           # North arrow size
                                 arrow_padin = c(0.15, 0.15),  # North arrow position adjustement
                                 
                                 facet_letter = "",                  # Small case letter for facet
                                 facet_letter_col = "black",         # Color of case letter for facet
                                 facet_letter_cex = 2.2,             # Size of small case letter for facet
                                 facet_letter_inset = c(0, -0.008))  # Position adjustment of small case letter for facet
  
  {
    # Plot raster background without axis
    image(x, col = color_palette,
          xlim = xlim, ylim = ylim, axes = F,
          xlab = xlab, ylab = ylab)
    title(main = main_title, cex.main = main_title_cex, line = 1.3)
    
    # Generate axes with manual positioning of ticks
    axis(side = 1, at = x_axis_breaks, labels = x_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1, gap.axis = 0, padj = 0.5)
    axis(side = 2, at = y_axis_breaks, labels = y_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1, gap.axis = 0)
    
    # Add background, borders and graticules
    plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
    plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
    plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
    plot(country_borders_Mollweide, lwd = 1, border = "#00000030", col = NA, add = T)
    
    # Add occurrence data
    plot(unit_pts, add = TRUE, col = occ_col, pch = 16, cex = occ_cex)
    
    # Add occurrence legend
    legend(legend = NA, pch = 16, cex = 1, x = -2700, y = -2600, bty = "n")
    text(labels = "Occurrences", cex = 1.15, x = -1200, y = -2945, bty = "n", font = 2)
    
    
    # Add scale bar in legend
    scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = scale_bar_position, label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1.2)
    prettymapr::addnortharrow(pos = "bottomright", scale = arrow_scale, padin = arrow_padin, text.col = "#00000000")
    rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
    rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
    graphics::text(x = legend_title_x, y = legend_title_y, font = 2, cex = legend_title_cex, label = legend_title)
    
    # Add facet letter
    legend(legend = facet_letter, x = "bottomright", bty = "n", text.col = facet_letter_col,
           text.font = 2, cex = facet_letter_cex, inset = facet_letter_inset)
    
  }
}

pdf(file = paste0("./maps/Community_Structure/BC_illustration.pdf"), height = 10, width = 10)

par(mfrow = c(2,2))
original_int_margins <- par()$mar
par(mar = c(3, 3, 3, 1)) # b l t r

OMU_names <- c("Dircenna.jemina.DILUCIDA", "Dircenna.dero.DILUCIDA", "Oleria.amalda.LERIDA", "Mechanitis.mazaeus.MAELUS")
BC_list <- c("", "0.35", "0.75", "0.95")

for (i in seq_along(OMU_names))
{
  # i <- 1
  
  unit <- OMU_names[i]
  
  # Load unit occurrence points
  load(paste0("./input_data/Species_data/occurrences_", unit,".RData")) # Load Spatial object with occurrences
  Mollweide_shp_projection(unit.points)
  
  # Retrieve type of model and parse name
  unit_spec <- extract_OMU_name_and_pattern(unit)
  unit_title <- bquote(bolditalic(.(unit_spec[1]))~bold(.(unit_spec[2])))
  
  # Adjust map projection and range
  map <- Mollweide_projection(OMU_proba_stack[[unit]])
  map@data@max <- 1
  
  # Map the OMU habitat suitability
  map_unit_Mollweide(x = map,
                     unit_pts = unit.points_Mollweide,
                     color_palette = c("#EDEDED", tmaptools::get_brewer_pal("RdYlGn", n = 199, plot = F)),
                     main_title = unit_title,
                     scale_bar_position = c(-2400, -4000),
                     legend_breaks = seq(0, 1, 0.2),
                     legend_title = "Habitat\n   suitability",
                     legend_title_x = -3400,
                     legend_title_y = 850)
  
  BC_value <- BC_list[i]
  
  if (i > 1)
  {
    legend(x = 1750, y = 3130, legend = "              ", cex = 1.4, bg = "white", box.lwd = 2, box.col = "grey40")
    legend(x = 1125, y = 3200, legend = bquote(bold('BC ='~.(BC_value))), cex = 1.6, bty = "n")
    
  }
  
}

par(mar = original_int_margins)
par(mfrow = c(1,1))

dev.off()


# Panel A: "Dircenna.jemina.DILUCIDA"
OMU_proba_stack[["Dircenna.jemina.DILUCIDA"]]


# Panel B: "Dircenna.dero.DILUCIDA" ; BC = 0.35

# Panel C: "Oleria.amalda.LERIDA" ; BC = 0.75

# Panel D: "Mechanitis.mazaeus.MAELUS" ; BC = 0.95