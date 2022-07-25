##### Script 06bis: Test fit of Evolutionary models, test for mean Climatic distance and phyloMANOVA #####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

### Focus only on monomorphic/species => 1 OMU = 1 species

### MCD analyses limited to the 139 monomorphic species in the phylogeny
### phyloMANOVA limited to mimicry ring with N >= 5 => 106 monomorphic species/OMUs

# Test for convergence in the climatic niche via mean climatic distance of comimics
# Run Phylogenetic MANOVA to test for divergence among rings for the climatic niche, accounting for the phylogeny 


### Input files

# Summary tables of 139/106 monomorphic species included in the phylogeny
# Phylogeny of the Ithomiini (Chazot et al., 2019)
# pPCA climatic values of the 139/106 OMUS

### Output files

# Simulated values for the 139/106 OMUs in the pPCA-space following the best evolutionary model
# Comimicry matrix for the 139/106 OMUs
# Global MCD for comimics, and per mimicry rings + null distribution
# Plot null distribution for test of MCD, global and per mimicry ring
# Plot null distribution of Wilks' lambda and pseudo-F of the phyloMANOVA


# Clean environment
rm(list = ls())

##### 1/ Load stuff  #####

# Load summary table for the 719 OMUs in the phylogeny
list.unit_phyl_order_719 <- readRDS(file = paste0("./outputs/Niche_evolution/list.unit_phyl_order_719.rds"))

# Load summary table for the 619 OMUs in the phylogeny and in big rings
reduced.list.unit_phyl_order_619 <- readRDS(file = paste0("./outputs/Niche_evolution/reduced.list.unit_phyl_order_619.rds"))

### Load phylogeny, unit list, and pPCA env data for evolutionary models and MCD analyses (include all species in the phylogeny)

# 139 monomorphic species/OMUs in the phylogeny for MCD analyses
phylo.Ithomiini_mono_139 <- readRDS(file = paste0("./input_data/Phylogenies/phylo.Ithomiini_mono_139.rds"))
list.unit_phyl_order_mono_139 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/list.unit_phyl_order_mono_139.rds"))

# Extract observed pPCA env data for the 139 monomorphic species
load(file = paste0("./outputs/Niche_evolution/Evol_simul/pPC.env_units_719.RData"))
pPC.env_units_mono_139 <- pPC.env_units_719[match(phylo.Ithomiini_mono_139$tip.label, row.names(pPC.env_units_719)), ]
saveRDS(pPC.env_units_mono_139, file = paste0("./outputs/Niche_evolution/Mono_sp/pPC.env_units_mono_139.rds"))

identical(phylo.Ithomiini_mono_139$tip.label, row.names(list.unit_phyl_order_mono_139))
identical(phylo.Ithomiini_mono_139$tip.label, row.names(pPC.env_units_mono_139))

### Load phylogeny and unit list for the phyloMANOVA analyses (include only rings with N >= 5)

# 106 monomorphic species/OMUs in the phylogeny, in big rings (N >= 5) for phylo MANOVA
phylo.Ithomiini_mono_106 <- readRDS(file = paste0("./input_data/Phylogenies/phylo.Ithomiini_mono_106.rds"))
reduced.list.unit_phyl_order_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/reduced.list.unit_phyl_order_mono_106.rds"))
pPC.env_units_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/pPC.env_units_mono_106.rds"))

identical(phylo.Ithomiini_mono_106$tip.label, row.names(reduced.list.unit_phyl_order_mono_106))
identical(phylo.Ithomiini_mono_106$tip.label, row.names(pPC.env_units_mono_106))


##### 2/ Extract simulated data and comimicry matrices #####

### 2.1/ Load simulated data for the 719 OMUs in the phylogeny ###

intraspecific_choice <- "_from_obs" # Data-informed. Costumed within-species variance based on observations, for each species
load(paste0("./outputs/Niche_evolution/Evol_simul/Sim_clim_2pPCA_OMUs_719",intraspecific_choice,".RData"))

### 2.2/ Extract the 139 OMUs/monomorphic species included in the MCD analyses ####

OMUs_139_in_719_indices <- match(list.unit_phyl_order_mono_139$Tag.model, list.unit_phyl_order_719$Tag.model)

Sim_clim_2pPCA_OMUs_139 <- list() # Create final list to store results
for (i in 1:length(Sim_clim_2pPCA_OMUs_719)) # Per simulation
{
  # i <- 1
  
  # Get the data for the 719 OMUs
  pPCA_simul_temp <- Sim_clim_2pPCA_OMUs_719[[i]]
  
  # Extract only the 139 OMUs and store them in the final list
  Sim_clim_2pPCA_OMUs_139 [[i]] <- pPCA_simul_temp[OMUs_139_in_719_indices, ]
  
  if(i %% 100 == 0)
  {
    cat(paste0("Simulation n° ",i, "/999\n"))
  }
  
}

saveRDS(Sim_clim_2pPCA_OMUs_139, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Sim_clim_2pPCA_OMUs_139",intraspecific_choice,".rds"))


### 2.3/ Extract the 106 OMUs included in the perMANOVA and phlyoMANOVA analyses ####

OMUs_106_in_719_indices <- match(reduced.list.unit_phyl_order_mono_106$Tag.model, list.unit_phyl_order_719$Tag.model)

Sim_clim_2pPCA_OMUs_106 <- list() # Create final list to store results
for (i in 1:length(Sim_clim_2pPCA_OMUs_719)) # Per simulation
{
  # i <- 1
  
  # Get the data for the 719 OMUs
  pPCA_simul_temp <- Sim_clim_2pPCA_OMUs_719[[i]]
  
  # Extract only the 106 OMUs and store them in the final list
  Sim_clim_2pPCA_OMUs_106 [[i]] <- pPCA_simul_temp[OMUs_106_in_719_indices, ]
  
  if(i %% 100 == 0)
  {
    cat(paste0("Simulation n° ",i, "/999\n"))
  }
  
}

saveRDS(Sim_clim_2pPCA_OMUs_106, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Sim_clim_2pPCA_OMUs_106",intraspecific_choice,".rds"))


### 2.4/ Extract co-mimicry matrix for the 139 OMUs in the BC and MCD analyses ####

## Incidence matrix: 1 = comic, 0 = non-comimic. Used to quickly compute mean for comimimetic OMUs.

# Load data for the 719 OMUs in big rings
list.unit_phyl_order_719 <- readRDS(file = paste0("./outputs/Niche_evolution/list.unit_phyl_order_719.rds"))
load(file = paste0("./outputs/Niche_evolution/comimicry_matrix_units.RData"))

identical(list.unit_phyl_order_719$Tag.model, colnames(comimicry_matrix_units))

# Find matches for the 139 OMUs
OMUs_139_in_719_indices <- match(list.unit_phyl_order_mono_139$Tag.model, list.unit_phyl_order_719$Tag.model)

# Extract comimicry matrix for the 139 OMUs
comimicry_matrix_units_mono_139 <- comimicry_matrix_units[OMUs_139_in_719_indices, OMUs_139_in_719_indices]

identical(colnames(comimicry_matrix_units_mono_139), list.unit_phyl_order_mono_139$Tag.model)

saveRDS(comimicry_matrix_units_mono_139, file = paste0("./outputs/Niche_evolution/Mono_sp/comimicry_matrix_units_mono_139.rds"))


##### 3/ Test for convergence in the climatic niche via mean climatic distance of comimics #####

### Analyses on the 139 monomorphic species ###

### 3.1/ Compute observed MCD ####

# Load comimicry matrix for the 139 monomorphic species
comimicry_matrix_units_mono_139 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/comimicry_matrix_units_mono_139.rds"))

# Load observed pPCA env data
pPC.env_units_mono_139 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/pPC.env_units_mono_139.rds"))

# Compute the pairwise climatic distance matrix between the 139 OMUs based on pPCA variables
pairwise_climdist_pPCA_mono_139 <- dist(x = pPC.env_units_mono_139, method = "euclidian")
pairwise_climdist_pPCA_mono_139_mat <- as.matrix(pairwise_climdist_pPCA_mono_139)
save(pairwise_climdist_pPCA_mono_139, pairwise_climdist_pPCA_mono_139_mat, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/pairwise_climdist_pPCA_mono_139.RData"))

length(pairwise_climdist_pPCA_mono_139)
length(as.dist(comimicry_matrix_units_mono_139))

# Compute mean climatic distances
Global_MCD_obs <- mean(pairwise_climdist_pPCA_mono_139)
Global_MCD_obs # Global MCD for all pairs of monomorphic species = 4.69

Comimic_MCD_obs <- weighted.mean(x = pairwise_climdist_pPCA_mono_139, w = as.dist(comimicry_matrix_units_mono_139))
Comimic_MCD_obs # Global MCD only for pairs of comimics = 4.03

# Need to standardized the comimic MCD because some simulations may create a more dilated or reduced climatic space which make comparison of mean distances biased
Comimic_MCD_obs_std <- Comimic_MCD_obs/Global_MCD_obs
Comimic_MCD_obs_std # Standardized MCD obs = 0.860

save(Global_MCD_obs, Comimic_MCD_obs, Comimic_MCD_obs_std, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_obs_stats_mono_139.RData"))

load(paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_obs_stats_mono_139.RData"))


### 3.2/ Compute MCD null distri ####

# Choose the within-species variation model
# intraspecific_choice <- ""
# intraspecific_choice <- "_null"
intraspecific_choice <- "_from_obs"

# Load simulated env data
Sim_clim_2pPCA_OMUs_139 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Sim_clim_2pPCA_OMUs_139",intraspecific_choice,".rds"))

# Loop per simulations
Global_MCD_null <- Comimic_MCD_null <- Comimic_MCD_std_null <-  NA # Initiate vectors to store results
for (i in 1:length(Sim_clim_2pPCA_OMUs_139)) # Per simulation
{
  # i <- 1
  
  # Compute the pairwise climatic distance matrix between the 139 OMUs based on simulated pPCA variables
  pairwise_climdist_pPCA_139_simul <- dist(x = Sim_clim_2pPCA_OMUs_139[[i]], method = "euclidian")
  
  # Compute mean climatic distances
  Global_MCD_simul <- mean(pairwise_climdist_pPCA_139_simul) # For all pairs of OMUs
  Comimic_MCD_simul <- weighted.mean(x = pairwise_climdist_pPCA_139_simul, w = as.dist(comimicry_matrix_units_mono_139)) # For comimics only
  # Need to standardized the comimic MCD because some simulations may create a more dilated or reduced climatic space which make comparison of mean distances biased
  Comimic_MCD_std_simul <- Comimic_MCD_simul/Global_MCD_simul # Standardized MCD
                                     
  # Store results                             
  Global_MCD_null[i] <- Global_MCD_simul
  Comimic_MCD_null[i] <- Comimic_MCD_simul
  Comimic_MCD_std_null[i] <- Comimic_MCD_std_simul
  
  save(Global_MCD_null, Comimic_MCD_null, Comimic_MCD_std_null, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_null_stats_139",intraspecific_choice,".RData"))
  
  if(i %% 10 == 0) {cat(paste0("Simulation n° ",i, "/999\n"))}
}
save(Global_MCD_null, Comimic_MCD_null, Comimic_MCD_std_null, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_null_stats_139",intraspecific_choice,".RData"))

summary(Comimic_MCD_std_null)


### 3.3/ Plot the distri of the stats for all comimics ####

# Choose the within-species variation model
# intraspecific_choice <- ""
# intraspecific_choice <- "_null"
intraspecific_choice <- "_from_obs"

# Load data
load(paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_obs_stats_mono_139.RData"))
load(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_null_stats_139",intraspecific_choice,".RData"))

# Plot
pdf(file = paste0("./graphs/Niche_evolution/Mono_sp/Evol_simul/Comimics_MCD_null_139",intraspecific_choice,".pdf"), height = 6, width = 7)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

hist(x = c(Comimic_MCD_std_null, Comimic_MCD_obs_std), breaks = 40, freq = TRUE, col = "gray",
     xlim = c(0.80, 1.15),
     ylim = c(0, 60),
     main = "Distribution of the Mean Climatic Distance\nof co-mimetic monomorphic species\nunder neutral evolution",
     # main = "",
     xlab = "Standardized Mean pairwise Climatic Distance",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)

arrows(Comimic_MCD_obs_std + 0.002, 23, Comimic_MCD_obs_std + 0.002, 3, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = mean(c(Comimic_MCD_std_null, Comimic_MCD_obs_std)), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(c(Comimic_MCD_std_null, Comimic_MCD_obs_std), 0.05), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("Mean = ", format(round(mean(c(Comimic_MCD_std_null, Comimic_MCD_obs_std)),3), nsmall = 3)), 
                  paste0("CI 5% = ", round(quantile(c(Comimic_MCD_std_null, Comimic_MCD_obs_std), 0.05),3))), 
       x = "topleft", inset = c(0.02, 0.05), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty = "n")
legend(legend = c(paste0("MCD obs = ", round(Comimic_MCD_obs_std, 3)),
                  paste0("p = 0.001")),
       x = "left", inset = c(-0.02, -0.05),
       cex = 1.2, bty ="n", xjust = 1)

# legend(legend = as.expression(bquote(bold("A"))), 
#        x = "topright", inset = c(0.05, 0.001), xjust = 0.5,
#        cex = 1.3, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()


### 3.4/ Compute MCD obs and null distri per mimicry ring ####

load(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_null_stats_139",intraspecific_choice,".RData"))
load(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/pairwise_climdist_pPCA_mono_139.RData"))

# Check if the order is the same
identical(colnames(pairwise_climdist_pPCA_mono_139_mat), list.unit_phyl_order_mono_139$Tag.model)

# Generate list of mimicry rings 
ring_list <- as.character(unique(list.unit_phyl_order_mono_139$Mimicry.model))

### Compute MCD obs per mimicry ring

MCD_per_ring_obs <- MCD_std_per_ring_obs <- rep(NA, length(ring_list)) # Generate empty vector to store MPD obs and its standardized version

for (i in 1:length(ring_list)) # Per mimicry ring
{
  # i <- 10
  
  ring <- ring_list[i]
  
  # Get indices of rows/columns associated of OMUs for this ring
  ring_indices <- which(list.unit_phyl_order_mono_139$Mimicry.model == ring)
  
  if(length(ring_indices) != 1) # Need computation only if there are more than 1 OMU, otherwise, no MCD possible
  {
    # Extract only the pairwise distances for this ring
    pairwise_dist_ring <- as.dist(pairwise_climdist_pPCA_mono_139_mat[ring_indices, ring_indices])
    
    MCD_per_ring_obs[i] <- mean(pairwise_dist_ring) # Compute MCD for this ring
    MCD_std_per_ring_obs[i] <- MCD_per_ring_obs[i]/Global_MCD_obs # Compute standardized MCD for this ring
      
  }
  print(i)
}
names(MCD_per_ring_obs) <- names(MCD_std_per_ring_obs) <- ring_list
MCD_std_per_ring_obs
save(MCD_per_ring_obs, MCD_std_per_ring_obs, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_per_ring_obs_mono_139.RData"))


### Compute MCD per ring for all the simulations

# Load MCD null stats to extract global MCD for each simulation. Used for standardization.
load(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_null_stats_139",intraspecific_choice,".RData"))

MCD_per_ring_simul <- NA # Initiate vectors to store results
MCD_per_ring_null <- MCD_std_per_ring_null <- matrix(data = NA, nrow = length(Sim_clim_2pPCA_OMUs_139), ncol = length(ring_list)) # Initiate matrix to store results
for (i in 1:length(Sim_clim_2pPCA_OMUs_139)) # Per simulation
{
  # i <- 1
  
  # Compute the pairwise climatic distance matrix between the 139 OMUs based on simulated pPCA variables
  pairwise_climdist_pPCA_139_simul <- as.matrix(dist(x = Sim_clim_2pPCA_OMUs_139[[i]], method = "euclidian"))
  
  for (j in 1:length(ring_list)) # Per mimicry ring
  {
    # j <- 1
    
    ring <- ring_list[j]
    
    # Get indices of rows/columns associated of OMUs for this ring
    ring_indices <- which(list.unit_phyl_order_mono_139$Mimicry.model == ring)
    
    if(length(ring_indices) != 1) # Need computation only if there are more than 1 OMU, otherwise, no MCD possible
    {
      # Extract only the pairwise distances for this ring
      pairwise_dist_ring <- as.dist(pairwise_climdist_pPCA_139_simul[ring_indices, ring_indices])
      
      # Compute MCD for this ring
      MCD_per_ring_simul[j] <- mean(pairwise_dist_ring) 

    } else { # Provide NA if not enough OMUs to form pairs
      MCD_per_ring_simul[j] <- NA
    }
  }
  
  # Store results for each simulation                        
  MCD_per_ring_null[i, ] <- MCD_per_ring_simul
  MCD_std_per_ring_null[i, ] <- MCD_per_ring_simul/Global_MCD_null[i] # Compute standardized MCD for this simulation

  save(MCD_per_ring_null, MCD_std_per_ring_null, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_per_ring_null_stats_139",intraspecific_choice,".RData"))
  
  if(i %% 10 == 0) {cat(paste0("Simulation n° ",i, "/999\n"))}
}
colnames(MCD_per_ring_null) <- colnames(MCD_std_per_ring_null) <- ring_list
save(MCD_per_ring_null, MCD_std_per_ring_null, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_per_ring_null_stats_139",intraspecific_choice,".RData"))

summary(MCD_std_per_ring_null)


### 3.5/ Plot MCD null distri per mimicry ring and generate summary table ####

load(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_per_ring_obs_mono_139.RData"))
load(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_per_ring_null_stats_139",intraspecific_choice,".RData")) 

MCD_std_per_ring_obs
MCD_std_per_ring_null

# Reorder in alphabetic order

alphabetic_order <- order(names(MCD_std_per_ring_obs))
MCD_std_per_ring_obs <- MCD_std_per_ring_obs[alphabetic_order]
MCD_std_per_ring_null <- MCD_std_per_ring_null[ , alphabetic_order]

ring_list <- names(MCD_std_per_ring_obs)


MCD_ring_summary_table_mono_139 <- as.data.frame(matrix(ncol = 9, nrow = length(ring_list), data = NA))
names(MCD_ring_summary_table_mono_139) <- c("ring", "N_units", "N_pairs", "MCD_obs", "mean_MCD", "MCD_2.5", "MCD_97.5", "p_value", "pattern")


for (i in 1:length(ring_list)) # Per mimicry ring
{
  # i <- 2
  
  MCD_ring_summary_table_mono_139$ring[i] <- ring <- ring_list[i]
  MCD_ring_summary_table_mono_139$N_units[i] <- N_units <- sum(list.unit_phyl_order_mono_139$Mimicry.model == ring)
  MCD_ring_summary_table_mono_139$N_pairs[i] <- N_pairs <- N_units*(N_units-1)/2
  
  if (is.na(MCD_std_per_ring_obs[i])) # Case for ring with only one OMUs. No pairs. No MCD.
  {
    pdf(file = paste0("./graphs/Niche_evolution/Mono_sp/Evol_simul/Per_ring/MCD_null_",ring,".pdf"), height = 6, width = 7)
    
    plot(1:100,1:100, type = "n", xlab = "Standardized Mean pairwise Climatic Distance",
         main = paste0("Distribution of Mean Climatic Distance \n of ", ring, " monomorphic species \n under neutral evolution"))
    text(x = 50, y = 50, labels = "Only one monomorphic species for this mimicry ring \n No pair available for index computation")
    
    dev.off()
    
  } else 
  {
    mean_val <- round(mean(c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i])),3)
    
    if (MCD_std_per_ring_obs[i] < mean_val)  # Case for significant signal for convergence
    {
      pattern <- "convergence"
      p_value <- round(ecdf(x = c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i]))(MCD_std_per_ring_obs[i]),3)
    } else  # Case for significant signal for divergence
    {   
      pattern <- "divergence"
      p_value <- round(1 - ecdf(x = c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i]))(MCD_std_per_ring_obs[i]),3)
    }
    
    MCD_ring_summary_table_mono_139$MCD_obs[i] <- round(MCD_std_per_ring_obs[i],3)
    MCD_ring_summary_table_mono_139$mean_MCD[i]  <- mean_val
    MCD_ring_summary_table_mono_139$MCD_2.5[i] <- round(quantile(c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i]), 0.025),3)
    MCD_ring_summary_table_mono_139$MCD_97.5[i] <- round(quantile(c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i]), 0.975),3)
    MCD_ring_summary_table_mono_139$p_value[i] <- p_value
    MCD_ring_summary_table_mono_139$pattern[i] <- pattern    
    
    
    histo.save <- hist(MCD_std_per_ring_null[,i],
                       breaks = 30,
                       plot = F)
    
    pdf(file = paste0("./graphs/Niche_evolution/Mono_sp/Evol_simul/Per_ring/MCD_null_",ring,".pdf"), height = 6, width = 7)
    
    hist(c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i]), 
         breaks = 30,
         col = "gray", xlab = "Standardized Mean pairwise Climatic Distance", 
         main = paste0("Distribution of the Mean Climatic Distance\nof co-mimetic monomoprhic species\nunder neutral evolution"),
         cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)
    arrows(MCD_std_per_ring_obs[i], max(histo.save$counts)/3, MCD_std_per_ring_obs[i], max(histo.save$counts)/30, length = 0.1, lwd = 2)
    abline(v = mean(c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i])), lty = 2, lwd = 2)
    abline(v = quantile(c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i]), 0.025, na.rm = T), lty = 2, lwd = 2, col = "red")
    abline(v = quantile(c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i]), 0.975, na.rm = T), lty = 2, lwd = 2, col = "red")
    legend(legend = c(paste0("N units = ", N_units), 
                      paste0("N pairs = ", N_pairs)),
           x = "topright", cex = 1, bty ="n")
    
    legend(legend = c(paste0("Mean = ", mean_val), 
                      paste0("CI 2.5% = ", round(quantile(c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i]), 0.025, na.rm = T),3)),
                      paste0("CI 97.5% = ", round(quantile(c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i]), 0.975, na.rm = T),3))), 
           x = "topleft", cex = 1, bty ="n")
    
    legend(legend = c(paste0("MCD obs = ", round(MCD_std_per_ring_obs[i], 3)),
                      paste0("p = ", p_value)),
           x = "left", cex = 1, bty ="n")
    
    dev.off()
  }
  
  save(MCD_ring_summary_table_mono_139, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_ring_summary_table_mono_139.Rdata"))
  
  cat(paste0("N° ",i, "/",length(ring_list)," - ",ring, " - Done \n"))
}

View(MCD_ring_summary_table_mono_139)


### 3.6/ Export summary tables for MCD analyses ####

## 3.6.1/ Export table for monomorphic species

# Sort in alphabetic order
MCD_ring_summary_table_mono_139 <- arrange(MCD_ring_summary_table_mono_139, ring)

View(MCD_ring_summary_table_mono_139)

# Save
saveRDS(MCD_ring_summary_table_mono_139, file = "./tables/MCD_ring_summary_table_mono_139.rds")
write.csv2(MCD_ring_summary_table_mono_139, file = "./tables/MCD_ring_summary_table_mono_139.csv")

## 3.6.2/ Add entries for mimicry rings absent from monomorphic species

MCD_ring_summary_table_all_OMU <- read.csv2(file = "./tables/MCD_ring_summary_table_719.csv", row.names = 1)
MCD_ring_summary_table_all_OMU <- arrange(MCD_ring_summary_table_all_OMU, ring)

all_rings <- MCD_ring_summary_table_all_OMU[, "ring", F]

MCD_ring_summary_table_mono_139_with_null <- left_join(all_rings, MCD_ring_summary_table_mono_139)
MCD_ring_summary_table_mono_139_with_null$N_units[is.na(MCD_ring_summary_table_mono_139_with_null$N_units)] <- 0
MCD_ring_summary_table_mono_139_with_null$N_pairs[is.na(MCD_ring_summary_table_mono_139_with_null$N_pairs)] <- 0

View(MCD_ring_summary_table_mono_139_with_null)

# Save
saveRDS(MCD_ring_summary_table_mono_139, file = "./tables/MCD_ring_summary_table_mono_139.rds")
write.csv2(MCD_ring_summary_table_mono_139, file = "./tables/MCD_ring_summary_table_mono_139.csv")

### 6.3/ Fuse with the complete table for all OMUs

MCD_ring_summary_table_both <- MCD_ring_summary_table_mono_139
names(MCD_ring_summary_table_both) <- c("ring", "N_units_mono", "N_pairs_mono", "BC_obs_mono", "mean_BC_mono", "BC_2.5_mono", "BC_97.5_mono", "p_value_mono", "pattern_mono")

MCD_ring_summary_table_both <- left_join(MCD_ring_summary_table_all_OMU, MCD_ring_summary_table_both)

View(MCD_ring_summary_table_both)

# Save
saveRDS(MCD_ring_summary_table_both, file = "./tables/MCD_ring_summary_table_both.rds")
write.csv2(MCD_ring_summary_table_both, file = "./tables/MCD_ring_summary_table_both.csv")



### Could run the same MCD analyses for the 106 monomorphic species included in the perMANOVA and phlyoMANOVA analyses 


##### 5/ phlyoMANOVA to test for divergence of climatic niche among rings #####

# Choose the within-species variation model
# intraspecific_choice <- ""
# intraspecific_choice <- "_null"
intraspecific_choice <- "_from_obs"

# Need to use only the 106 OMUs included in rings with N >= 10

# Load OMU summary table
reduced.list.unit_phyl_order_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/reduced.list.unit_phyl_order_mono_106.rds"))
# Load simulated env data
Sim_clim_2pPCA_OMUs_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Sim_clim_2pPCA_OMUs_106",intraspecific_choice,".rds"))
# Load observed env data
pPC.env_units_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/pPC.env_units_mono_106.rds"))

View(reduced.list.unit_phyl_order_mono_106)
str(Sim_clim_2pPCA_OMUs_106)
View(pPC.env_units_mono_106)

identical(rownames(pPC.env_units_mono_106), reduced.list.unit_phyl_order_mono_106$Tag.model)

### 5.1/ Compute MANOVA on observed data and extract the Wilks' summary stat and the pseudo-F ####

# Not use geiger::aov.phlyo because it does not allow to use data that have already been simulated
# library(geiger)
# ?aov.phylo

MANOVA_obs_pPCA_mono_106 <- manova(pPC.env_units_mono_106[,1:2] ~ reduced.list.unit_phyl_order_mono_106$Mimicry.model)
summary(MANOVA_obs_pPCA_mono_106, test="Wilks") # Global test
summary.aov(MANOVA_obs_pPCA_mono_106) # Test per response variable (the two pPCA-axis)
Wilks_lambda_obs_pPCA_mono_106 <- summary(MANOVA_obs_pPCA_mono_106, test="Wilks")$stats[,2][1] # Save the wilk's lambda
Pseudo_F_obs_pPCA_mono_106 <- summary(MANOVA_obs_pPCA_mono_106, test="Wilks")$stats[,3][1] # Save the associated pseudo-F used for the test (not the phylo test, the regular MANOVA test)
saveRDS(Wilks_lambda_obs_pPCA_mono_106, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Wilks_lambda_obs_pPCA_mono_106.rds"))
saveRDS(Pseudo_F_obs_pPCA_mono_106, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Pseudo_F_obs_pPCA_mono_106.rds"))


### 5.2/ Compute MANOVA on each simulation and extract the Wilks' summary stat and the pseudo-F ####

Wilks_lambda_null_pPCA_mono_106 <- Pseudo_F_null_pPCA_mono_106 <- NA
for (i in 1:length(Sim_clim_2pPCA_OMUs_106)) # Per simulation
{
  MANOVA_simul_pPCA_mono_106 <- manova(as.matrix(Sim_clim_2pPCA_OMUs_106[[i]]) ~ reduced.list.unit_phyl_order_mono_106$Mimicry.model)
  Wilks_lambda_null_pPCA_mono_106[i] <- summary(MANOVA_simul_pPCA_mono_106, test="Wilks")$stats[,2][1] # Save the wilk's lambda
  Pseudo_F_null_pPCA_mono_106[i] <- summary(MANOVA_simul_pPCA_mono_106, test="Wilks")$stats[,3][1] # Save the associated pseudo-F used for the test (not the phylo test, the regular MANOVA test)
  
  if(i %% 100 == 0) {cat(paste0("Simulation n° ",i, "/",length(Sim_clim_2pPCA_OMUs_106),"\n"))}
}

summary(Wilks_lambda_null_pPCA_mono_106)
summary(Pseudo_F_null_pPCA_mono_106)

saveRDS(Wilks_lambda_null_pPCA_mono_106, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Wilks_lambda_null_pPCA_mono_106",intraspecific_choice,".rds"))
saveRDS(Pseudo_F_null_pPCA_mono_106, file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Pseudo_F_null_pPCA_mono_106",intraspecific_choice,".rds"))

### 5.3/ Plot the distri of the Wilks lambda ####

# Load observed and simulated Wilks lambda
Wilks_lambda_obs_pPCA_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Wilks_lambda_obs_pPCA_mono_106.rds"))
Wilks_lambda_null_pPCA_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Wilks_lambda_null_pPCA_mono_106",intraspecific_choice,".rds"))


pdf(file = paste0("./graphs/Niche_evolution/Mono_sp/Evol_simul/MANOVA_phylo_pPCA_Wilks",intraspecific_choice,".pdf"), height = 6, width = 7)

hist(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106), 30, freq = TRUE, col = "gray", 
     main = bquote("Phylogenetic MANOVA \n Distribution of Wilks' lambda\nunder neutral evolution"),
     xlab = bquote("Wilks'" ~lambda),
     xlim = c(0.3, 1.0),
     ylim = c(0, 150),
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5, lwd = 2)

arrows(Wilks_lambda_obs_pPCA_mono_106 - 0.002, 50, Wilks_lambda_obs_pPCA_mono_106 - 0.002, 5, length = 0.1, lwd = 2)
abline(v = mean(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106)), lty = 2, lwd = 2)
abline(v = quantile(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106), 0.05), lty = 2, lwd = 2, col = "red")

legend(legend = c(paste0("Mean = ", round(mean(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106)),3)), 
                  paste0("CI 5% = ", round(quantile(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106), 0.05),3))), 
       x = "topleft", inset = c(0.05, 0.00), lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")

legend(legend = bquote(lambda ~ "obs =" ~ .(format(round(Wilks_lambda_obs_pPCA_mono_106, 3), nsmall = 3))),
       x = "bottomleft", cex = 1.2, bty ="n",
       inset = c(0.02, 0.43))
legend(legend = paste0("p = ", round(ecdf(x = c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106))(Wilks_lambda_obs_pPCA_mono_106),3)),
       x = "bottomleft", cex = 1.2, bty ="n", y.intersp = 2.1,
       inset = c(0.02, 0.35))

# legend(legend = as.expression(bquote(bold("B"))), 
#        x = "topright", inset = c(0.05, 0.001), xjust = 0.5,
#        cex = 1.3, bty ="n")

dev.off()

### 5.4/ Plot the distri of the pseudo-F ####

# Load observed and simulated pseudo-F
Pseudo_F_obs_pPCA_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Pseudo_F_obs_pPCA_mono_106.rds"))
Pseudo_F_null_pPCA_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Pseudo_F_null_pPCA_mono_106",intraspecific_choice,".rds"))


pdf(file = paste0("./graphs/Niche_evolution/Mono_sp/Evol_simul/MANOVA_phylo_pPCA_pseudo_F",intraspecific_choice,".pdf"), height = 6, width = 7)

hist(c(Pseudo_F_null_pPCA_mono_106, Pseudo_F_obs_pPCA_mono_106), 50, freq = TRUE, col = "gray", 
     main = "Phylogenetic MANOVA\n Distribution of pseudo-F\n under neutral evolution",
     xlab = "pseudo-F",
     xlim = c(0, 9),
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5, lwd = 2)

arrows(Pseudo_F_obs_pPCA_mono_106 + 0.08, 60, Pseudo_F_obs_pPCA_mono_106 + 0.08, 5, length = 0.1, lwd = 2)
abline(v = mean(c(Pseudo_F_null_pPCA_mono_106, Pseudo_F_obs_pPCA_mono_106)), lty = 2, lwd = 2)
abline(v = quantile(c(Pseudo_F_null_pPCA_mono_106, Pseudo_F_obs_pPCA_mono_106), 0.95), lty = 2, lwd = 2, col = "red")

legend(legend = c(paste0("Mean = ", round(mean(c(Pseudo_F_null_pPCA_mono_106, Pseudo_F_obs_pPCA_mono_106)),2)), 
                  paste0("CI 95% = ", round(quantile(c(Pseudo_F_null_pPCA_mono_106, Pseudo_F_obs_pPCA_mono_106), 0.95),2))), 
       x = "topright", inset = c(0.02, 0.15), lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")

legend(legend = paste0("pseudo-F obs = ", round(Pseudo_F_obs_pPCA_mono_106, 2)),
       x = "bottomright", cex = 1.2, bty ="n",
       inset = c(0.00, 0.44))
legend(legend = paste0("p = 0.001"),
       x = "bottomright", cex = 1.2, bty ="n", y.intersp = 2.1,
       inset = c(0.00, 0.36))

# legend(legend = as.expression(bquote(bold("B"))), 
#        x = "topright", inset = c(0.05, 0.001), xjust = 0.5,
#        cex = 1.3, bty ="n")

dev.off()


##### 6/ Plot both MCD and Wilks' lambda #####

# Choose the within-species variation model
# intraspecific_choice <- ""
# intraspecific_choice <- "_null"
intraspecific_choice <- "_from_obs"

# Load data for MCD plot
load(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_per_ring_obs_mono_139.RData"))
load(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/MCD_per_ring_null_stats_139",intraspecific_choice,".RData")) 

# Load data for Wilk's lambda plot
Wilks_lambda_obs_pPCA_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Wilks_lambda_obs_pPCA_mono_106.rds"))
Wilks_lambda_null_pPCA_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/Evol_simul/Wilks_lambda_null_pPCA_mono_106",intraspecific_choice,".rds"))

round(Comimic_MCD_obs_std, 3) # MCD obs = 0.86
format(round(Wilks_lambda_obs_pPCA_mono_106, 3), nsmall = 3) # Lambda obs = 0.352


### 6.1/ Plot Wilk's lambda (A) and MCD (B) ####

pdf(file = paste0("./graphs/Niche_evolution/Mono_sp/Evol_simul/Wilks_lambda_&_MCD_both_plots",intraspecific_choice,".pdf"), height = 6.4, width = 14)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,5,5.5,1.5), mfrow = c(1,2))

# Panel A: Wilk's lambda from phylogenetic MANOVA

hist(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106), 30, freq = TRUE, col = "gray", 
     main = bquote("Phylogenetic MANOVA \n Distribution of Wilks' lambda\nunder neutral evolution"),
     xlab = bquote("Wilks'" ~lambda),
     xlim = c(0.3, 1.0),
     ylim = c(0, 150),
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.3, lwd = 2)

arrows(Wilks_lambda_obs_pPCA_mono_106 - 0.002, 50, Wilks_lambda_obs_pPCA_mono_106 - 0.002, 5, length = 0.1, lwd = 2)
abline(v = mean(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106)), lty = 2, lwd = 2)
abline(v = quantile(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106), 0.05), lty = 2, lwd = 2, col = "red")

legend(legend = c(paste0("Mean = ", round(mean(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106)),3)), 
                  paste0("CI 5% = ", round(quantile(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106), 0.05),3))), 
       x = "topleft", inset = c(0.02, 0.20), lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")

legend(legend = bquote(lambda ~ "obs =" ~ .(format(round(Wilks_lambda_obs_pPCA_mono_106, 3), nsmall = 3))),
       x = "bottomleft", cex = 1.2, bty ="n",
       inset = c(0.01, 0.43))
legend(legend = paste0("p = ", round(ecdf(x = c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106))(Wilks_lambda_obs_pPCA_mono_106),3)),
       x = "bottomleft", cex = 1.2, bty ="n", y.intersp = 2.1,
       inset = c(0.01, 0.35))

legend(legend = as.expression(bquote(bold("A"))), 
       x = "topleft", inset = c(-0.03, -0.03), xjust = 0.5,
       cex = 1.8, bty ="n")


# Panel B: Mean Climatic Distances (MCD)

hist(x = c(Comimic_MCD_std_null, Comimic_MCD_obs_std), breaks = 30, freq = TRUE, col = "gray",
     xlim = c(0.85, 1.15),
     ylim = c(0, 120),
     main = "Distribution of the Mean Climatic Distance\nof co-mimetic monomorphic species\nunder neutral evolution",
     # main = "",
     xlab = "Standardized Mean pairwise Climatic Distance",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.3, lwd = 2)

arrows(Comimic_MCD_obs_std + 0.004, 45, Comimic_MCD_obs_std + 0.004, 3, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = mean(c(Comimic_MCD_std_null, Comimic_MCD_obs_std)), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(c(Comimic_MCD_std_null, Comimic_MCD_obs_std), 0.05), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("Mean = ", format(round(mean(c(Comimic_MCD_std_null, Comimic_MCD_obs_std)),3), nsmall = 3)), 
                  paste0("CI 5% = ", round(quantile(c(Comimic_MCD_std_null, Comimic_MCD_obs_std), 0.05),3))), 
       x = "topleft", inset = c(0.02, 0.20), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty = "n")
legend(legend = c(paste0("MCD obs = ", round(Comimic_MCD_obs_std, 3)),
                  paste0("p = 0.001")),
       x = "left", inset = c(-0.02, -0.05),
       cex = 1.2, bty ="n", xjust = 1)

legend(legend = as.expression(bquote(bold("B"))), 
       x = "topleft", inset = c(-0.03, -0.03), xjust = 0.5,
       cex = 1.8, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins, mfrow = c(1,1))

dev.off()


### 6.2/ Same but with MCD as panel A and Wilk's lambda as panel B ####

pdf(file = paste0("./graphs/Niche_evolution/Mono_sp/Evol_simul/MCD_&_Wilks_lambda_both_plots",intraspecific_choice,".pdf"), height = 6.4, width = 14)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,5,5.5,1.5), mfrow = c(1,2))

hist(x = c(Comimic_MCD_std_null, Comimic_MCD_obs_std), breaks = 30, freq = TRUE, col = "gray",
     xlim = c(0.85, 1.15),
     ylim = c(0, 120),
     main = "Distribution of the Mean Climatic Distance\nof co-mimetic monomorphic species\nunder neutral evolution",
     # main = "",
     xlab = "Standardized Mean pairwise Climatic Distance",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.3, lwd = 2)

arrows(Comimic_MCD_obs_std + 0.004, 45, Comimic_MCD_obs_std + 0.004, 3, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = mean(c(Comimic_MCD_std_null, Comimic_MCD_obs_std)), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(c(Comimic_MCD_std_null, Comimic_MCD_obs_std), 0.05), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("Mean = ", format(round(mean(c(Comimic_MCD_std_null, Comimic_MCD_obs_std)),3), nsmall = 3)), 
                  paste0("CI 5% = ", round(quantile(c(Comimic_MCD_std_null, Comimic_MCD_obs_std), 0.05),3))), 
       x = "topleft", inset = c(0.02, 0.20), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty = "n")
legend(legend = c(paste0("MCD obs = ", round(Comimic_MCD_obs_std, 3)),
                  paste0("p = 0.001")),
       x = "left", inset = c(-0.02, -0.05),
       cex = 1.2, bty ="n", xjust = 1)

legend(legend = as.expression(bquote(bold("A"))),
       x = "topleft", inset = c(-0.03, -0.03), xjust = 0.5,
       cex = 1.8, bty ="n")

# Panel B: Wilk's lambda from phylogenetic MANOVA

hist(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106), 30, freq = TRUE, col = "gray", 
     main = bquote("Phylogenetic MANOVA \n Distribution of Wilks' lambda\nunder neutral evolution"),
     xlab = bquote("Wilks'" ~lambda),
     xlim = c(0.3, 1.0),
     ylim = c(0, 150),
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.3, lwd = 2)

arrows(Wilks_lambda_obs_pPCA_mono_106 - 0.002, 50, Wilks_lambda_obs_pPCA_mono_106 - 0.002, 5, length = 0.1, lwd = 2)
abline(v = mean(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106)), lty = 2, lwd = 2)
abline(v = quantile(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106), 0.05), lty = 2, lwd = 2, col = "red")

legend(legend = c(paste0("Mean = ", round(mean(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106)),3)), 
                  paste0("CI 5% = ", round(quantile(c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106), 0.05),3))), 
       x = "topleft", inset = c(0.02, 0.20), lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")

legend(legend = bquote(lambda ~ "obs =" ~ .(format(round(Wilks_lambda_obs_pPCA_mono_106, 3), nsmall = 3))),
       x = "bottomleft", cex = 1.2, bty ="n",
       inset = c(0.01, 0.43))
legend(legend = paste0("p = ", round(ecdf(x = c(Wilks_lambda_null_pPCA_mono_106, Wilks_lambda_obs_pPCA_mono_106))(Wilks_lambda_obs_pPCA_mono_106),3)),
       x = "bottomleft", cex = 1.2, bty ="n", y.intersp = 2.1,
       inset = c(0.01, 0.35))

legend(legend = as.expression(bquote(bold("B"))),
       x = "topleft", inset = c(-0.03, -0.03), xjust = 0.5,
       cex = 1.8, bty ="n")


par(oma = original_ext_margins, mar = original_int_margins, mfrow = c(1,1))

dev.off()
