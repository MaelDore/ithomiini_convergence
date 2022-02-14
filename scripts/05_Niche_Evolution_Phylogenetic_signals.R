##### Script 05: Phylogenetic signals #####

# Test for phylogenetic signal in the evolution of mimicry patterns AND climatic niche

### Input files

# Summary tables of OMUs and species
# Phylogeny of the Ithomiini (Chazot et al., 2019), with and without OMUs
# 


### Output files

# Summary table for climatic variable with only for the 339 species in the phylogeny
# Co-mimicry matrix for OMUs (who is comimic of who?)
# Global MPD for comimics, and per mimicry rings + null distribution
# Plot null distri for test of MPD, global and per mimcry ring
# Plot Kmult null distri for test for phylogenetic signal in climatic niche
# Plot with both MPD and Kmult


# Effacer l'environnement
rm(list = ls())

##### 1/ Prepare phylogeny and env data ####

# Load phylogeny, sp list and unit list
load(file = paste0("./input_data/Phylogenies/Final_phylogeny.RData"))
load(file = paste0("./input_data/Phylogenies/Final_units_phylogeny.RData"))
load(file = paste0("./input_data/Summary_tables/list.sp.RData"))
load(file = paste0("./input_data/Summary_tables/list.unit.RData"))

# Extract only the 339 species included in the phylogeny from list.sp and list.unit
list.sp <- list.sp[list.sp$Sp_full %in% phylo.Ithomiini$tip.label,] # 339 species
list.unit <- list.unit[list.unit$Sp_full %in% phylo.Ithomiini$tip.label,] # 719 units

# Reorder list.unit following Phylogeny to match rows in the phenetic distance matrix with env data in the summary table
order.index <- NA
for (i in 1:length(phylo.Ithomiini.units$tip.label)) {
  unit <- as.character(phylo.Ithomiini.units$tip.label[i])
  order.index[i] <- which(list.unit$Tag.model==unit)
}

list.unit_phyl_order <- list.unit[order.index,]
row.names(list.unit_phyl_order) <- as.character(list.unit_phyl_order$Tag.model)

save(list.unit_phyl_order, file = paste0("./outputs/Niche_evolution/list.unit_phyl_order.RData"))

# Get environmental data for units/OMUs
unit.env.table_719 <- list.unit_phyl_order[,c("bio1","bio4","bio12","bio15")]
names(unit.env.table_719) <- c("Tmean", "Tvar", "Hmean", "Hvar")
save(unit.env.table_719, file = paste0("./outputs/Niche_evolution/unit.env.table_719.RData"))

# Check the order is the same
identical(row.names(unit.env.table_719), phylo.Ithomiini.units$tip.label)

# Get environmental data for species
sp_env_table_339 <- list.unit_phyl_order[,c("Sp_full","bio1","bio4","bio12","bio15")]
sp_env_table_339 <- sp_env_table_339 %>% 
  group_by(Sp_full) %>%
  summarise(Tmean = mean(bio1),
            Tvar = mean(bio4),
            Hmean = mean(bio12),
            Hvar = mean(bio15))

# Reorder list.unit following Phylogeny to match rows in the phenetic distance matrix with env data in the summary table
order.index <- NA
for (i in 1:length(phylo.Ithomiini$tip.label)) {
  sp <- as.character(phylo.Ithomiini$tip.label[i])
  order.index[i] <- which(sp_env_table_339$Sp_full == sp)
}

sp_env_table_339 <- sp_env_table_339[order.index,]
row.names(sp_env_table_339) <- as.character(sp_env_table_339$Sp_full)

# Check the order is the same
identical(sp_env_table_339$Sp_full, phylo.Ithomiini$tip.label)

save(sp_env_table_339, file = paste0("./outputs/Niche_evolution/sp_env_table_339.RData"))

##### 2/ Test for phylogenetic signal in mimicry patterns with MPD : global and per ring #####

######## Version on OMUs ########

### 2.1/ Compute Phylogenetic/patristic distance matrix ####
library(phytools)

phylo.dist.mat_units <- cophenetic.phylo(x = phylo.Ithomiini.units)

### 2.2/ Compute co-mimicry matrix for OMUs ####
comimicry_matrix_units <- matrix(nrow = nrow(list.unit_phyl_order), ncol = nrow(list.unit_phyl_order), data = 0)
for (i in 1:nrow(list.unit_phyl_order))
{
  ring_1 <- as.character(list.unit_phyl_order$Mimicry.model)[i]
  
  for (j in 1:nrow(list.unit_phyl_order)) 
  {
    ring_2 <- as.character(list.unit_phyl_order$Mimicry.model)[j]
    
    if (ring_1 == ring_2) 
    {
      comimicry_matrix_units[i,j] <- 1
    }
  }
  if (i %% 10 == 0) {print(i)}
}
save(comimicry_matrix_units, file = paste0("./outputs/Niche_evolution/comimicry_matrix_units.RData"))

load(file = paste0("./outputs/Niche_evolution/comimicry_matrix_units.RData"))

hist(as.dist(phylo.dist.mat_units))

### 2.3/ Compute the global observed MPD weighted by mimicry similarity ####
mean(as.dist(phylo.dist.mat_units)) # Global MPD for all pairs of OMUs = 37.17 Ma
Global_MPD_obs <- weighted.mean(x = as.dist(phylo.dist.mat_units), w = as.dist(comimicry_matrix_units)) # Global MPD only for pairs of comimics = 34.43 Ma

save(Global_MPD_obs, file = paste0("./outputs/Niche_evolution/Phylo_signal/Global_MPD_obs.RData"))

### 2.4/ Compute the MPD per mimicry ring ####

mimic.list <- as.character(unique(list.unit_phyl_order$Mimicry.model))

MPD_per_mimic_obs <- rep(NA, length(mimic.list)) # Generate empty vector to store MPD obs

for (i in 1:length(mimic.list)) # Per mimic ring
{
  # i <-  3
  
  ring <- mimic.list[i]
  
  # Get indices of rows/columns associated of OMUs for this ring
  ring_indices <- which(list.unit_phyl_order$Mimicry.model == ring)
  
  if(length(ring_indices) != 1) # Need computation only if there are more than 1 OMU, otherwise, no MPD possible
  {
    # Extract only the pairwise distances for this ring
    pairwise_dist_ring <- as.dist(phylo.dist.mat_units[ring_indices, ring_indices])
    
    MPD_per_mimic_obs[i] <- mean(pairwise_dist_ring)
  }
  print(i)
}
names(MPD_per_mimic_obs) <- mimic.list
MPD_per_mimic_obs
save(MPD_per_mimic_obs, file = paste0("./outputs/Niche_evolution/Phylo_signal/MPD_per_mimic_obs.RData"))


### 2.5/ Simulate MPD under null hypothesis by permutation ####

## Start the loop
Global_MPD_null <-  NA
MPD_per_mimic_null <- data.frame(matrix(nrow = 999, ncol = length(mimic.list))) # Generate df to store results
for (l in 1:999) # 999 perumations/simulations
{
  # l <- 1
  
  ## Suffle mimicry ring among units in the comimicry matrix
  shuffle.indices <- sample(x = 1:nrow(list.unit_phyl_order), size = nrow(list.unit_phyl_order), replace = F)
  comimicry_matrix_units_shuffle <- comimicry_matrix_units[shuffle.indices, shuffle.indices]
  
  # Compute the global MPD
  Global_MPD_simul <- weighted.mean(x = as.dist(phylo.dist.mat_units), w = as.dist(comimicry_matrix_units_shuffle)) # Global MPD only for pairs of comimics = 34.43 Ma
  Global_MPD_null[l] <- Global_MPD_simul
  
  save(Global_MPD_null, file = paste0("./outputs/Niche_evolution/Phylo_signal/Global_MPD_null.RData"))
  
  
  # Compute the MPD per ring
  MPD_per_mimic_simul <- rep(NA, length(mimic.list)) # Generate empty vector to store MPD obs
  
  for (i in 1:length(mimic.list)) # Per mimic ring
  {
    # i <-  3
    
    ring <- mimic.list[i]
    
    # Get indices of rows/columns associated of OMUs for this ring in the real dataset
    ring_indices <- which(list.unit_phyl_order$Mimicry.model == ring)
    
    if (length(ring_indices) != 1) # Need computation only if there are more than 1 OMU, otherwise, no MPD possible
    {
      shuffle_ring_indices <- shuffle.indices[ring_indices] # Get the new indices of the OMU in the shuffled comimicry matrix
      
      # Extract only the pairwise distances for this ring
      pairwise_dist_ring <- as.dist(phylo.dist.mat_units[shuffle_ring_indices, shuffle_ring_indices])
      
      MPD_per_mimic_simul[i] <- mean(pairwise_dist_ring)
    }
  }
  
  MPD_per_mimic_null[l, ] <- MPD_per_mimic_simul
  
  save(MPD_per_mimic_null, file = paste0("./outputs/Niche_evolution/Phylo_signal/MPD_per_mimic_null.RData")) 
  cat(paste0(Sys.time(), " - Simul n°", l, " out of 999\n"))
}
names(MPD_per_mimic_null) <- mimic.list
save(MPD_per_mimic_null, file = paste0("./outputs/Niche_evolution/Phylo_signal/MPD_per_mimic_null.RData")) 


### 2.6/ Plot null distri for global comimics MPD ####

load(file = paste0("./outputs/Niche_evolution/Phylo_signal/Global_MPD_obs.RData"))
load(file = paste0("./outputs/Niche_evolution/Phylo_signal/Global_MPD_null.RData"))

pdf(file = "./graphs/Niche_evolution/Phylo_signal/Global_MPD_null.pdf", height = 6, width = 7)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

hist(x = c(Global_MPD_null, Global_MPD_obs), breaks = 40, freq = TRUE, col = "gray",
     # xlim = c(34.5, 38.0),
     # ylim = c(0, 200),
     main = "Distribution of the Mean pairwise Phylogenetic Distance \nof co-mimetic OMUs \nunder the null hypothesis",
     # main = "",
     xlab = "Mean pairwise Phylogenetic Distance [My]",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)

arrows(Global_MPD_obs + 0.015, 60, Global_MPD_obs+ 0.015, 5, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = mean(c(Global_MPD_null, Global_MPD_obs)), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(c(Global_MPD_null, Global_MPD_obs), 0.05), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("Mean = ", round(mean(c(Global_MPD_null, Global_MPD_obs)),2)), 
                  paste0("CI 5% = ", round(quantile(c(Global_MPD_null, Global_MPD_obs), 0.05),2))), 
       x = "topleft", inset = c(0.02, 0.05), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty = "n")
legend(legend = c(paste0("MPD obs = ", round(Global_MPD_obs, 2)),
                  paste0("p = 0.001")),
       x = "left", inset = c(-0.04, -0.05),
       cex = 1.2, bty ="n", xjust = 1)

legend(legend = as.expression(bquote(bold("A"))), 
       x = "topright", inset = c(0.05, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()



### 2.7/ Plot null distri MPD per mimic ring ####

load(file = paste0("./outputs/Niche_evolution/Phylo_signal/MPD_per_mimic_obs.RData"))
load(file = paste0("./outputs/Niche_evolution/Phylo_signal/MPD_per_mimic_null.RData")) 

MPD_per_mimic_obs
MPD_per_mimic_null

# Reorder in alphabetic order

alphabetic_order <- order(names(MPD_per_mimic_obs))
MPD_per_mimic_obs <- MPD_per_mimic_obs[alphabetic_order]
MPD_per_mimic_null <- MPD_per_mimic_null[ , alphabetic_order]

mimicry.list <- names(MPD_per_mimic_obs)

# Create summary table for phylogenetic signal of mimicry rings
MPD_ring_summary_table <- as.data.frame(matrix(ncol = 9, nrow = 44, data = NA))
names(MPD_ring_summary_table) <- c("ring", "N_units", "N_pairs", "MPD_obs", "mean_MPD", "MPD_2.5", "MPD_97.5", "p_value", "pattern")

for (i in 1:length(mimicry.list)) 
{
  # i <- 2
  
  MPD_ring_summary_table$ring[i] <- ring <- mimicry.list[i]
  MPD_ring_summary_table$N_units[i] <- N_units <- sum(list.unit$Mimicry.model == ring)
  MPD_ring_summary_table$N_pairs[i] <- N_pairs <- N_units*(N_units-1)/2
  
  if (is.na(MPD_per_mimic_obs[i]))  # Case for ring with only one OMUs. No pairs. No MPD.
  {
    pdf(file = paste0("./graphs/Niche_evolution/Phylo_signal/Per_ring/MPD_null_",ring,".pdf"), height = 6, width = 7)
    plot(1:100,1:100, type = "n", xlab = "Mean pairwise Phylogenetic Distance [My]", ylab = "",
         main = paste0("Distribution of Mean Phylogenetic Distance \n of ", ring, " OMUs \n under null Hypothesis"))
    text(x = 50, y = 50, labels = "Only one OMU for this mimicry ring \n No pair available for index computation")
    dev.off()
    
  } else 
  {
    mean_val <- round(mean(MPD_per_mimic_null[,i]),3)
    if (MPD_per_mimic_obs[i] < mean_val)  # Case for significant signal for clustering
    {
      pattern <- "clustering"
      p_value <- round(ecdf(x = MPD_per_mimic_null[,i])(MPD_per_mimic_obs[i]),3)
    } else {   # Case for significant signal for overdispersal
      pattern <- "overdispersal"
      p_value <- round(1 - ecdf(x = MPD_per_mimic_null[,i])(MPD_per_mimic_obs[i]),3)
    }
    
    MPD_ring_summary_table$MPD_obs[i] <- round(MPD_per_mimic_obs[i],3)
    MPD_ring_summary_table$mean_MPD[i] <- mean_val
    MPD_ring_summary_table$MPD_2.5[i] <- round(quantile(MPD_per_mimic_null[,i], 0.025, na.rm = T),3)
    MPD_ring_summary_table$MPD_97.5[i] <- round(quantile(MPD_per_mimic_null[,i], 0.975, na.rm = T),3)
    MPD_ring_summary_table$p_value[i] <- p_value
    MPD_ring_summary_table$pattern[i] <- pattern
    
    histo.save <- hist(MPD_per_mimic_null[,i],
                       breaks = seq(from = floor(min(min(MPD_per_mimic_null[,i], na.rm = T), MPD_per_mimic_obs[i])), to = ceiling(max(max(MPD_per_mimic_null[,i], na.rm = T), MPD_per_mimic_obs[i])), by = 0.5),
                       plot = F)
    
    pdf(file = paste0("./graphs/Niche_evolution/Phylo_signal/Per_ring/MPD_null_",ring,".pdf"), height = 6, width = 7)
    
    hist(MPD_per_mimic_null[,i], 
         breaks = seq(from = floor(min(min(MPD_per_mimic_null[,i], na.rm = T), MPD_per_mimic_obs[i])), to = ceiling(max(max(MPD_per_mimic_null[,i], na.rm = T), MPD_per_mimic_obs[i])), by = 0.5),
         col = "gray", xlab = "Mean Phylogenetic Distance [My]", 
         main = paste0("Distribution of the Mean Phylogenetic distances \n of ",ring," OMUs \n under the null Hypothesis"),
         cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)
    arrows(MPD_per_mimic_obs[i], max(histo.save$counts)/3, MPD_per_mimic_obs[i], max(histo.save$counts)/30, length = 0.1, lwd = 2)
    abline(v = mean(MPD_per_mimic_null[,i]), lty = 2, lwd = 2)
    abline(v = quantile(MPD_per_mimic_null[,i], 0.025, na.rm = T), lty = 2, lwd = 2, col = "red")
    abline(v = quantile(MPD_per_mimic_null[,i], 0.975, na.rm = T), lty = 2, lwd = 2, col = "red")
    legend(legend = c(paste0("N units = ", N_units), 
                      paste0("N pairs = ", N_pairs)),
           x = "topright", cex = 1, bty ="n")
    
    legend(legend = c(paste0("Mean = ", mean_val), 
                      paste0("CI 2.5% = ", round(quantile(MPD_per_mimic_null[,i], 0.025, na.rm = T),3)),
                      paste0("CI 97.5% = ", round(quantile(MPD_per_mimic_null[,i], 0.975, na.rm = T),3))), 
           x = "topleft", cex = 1, bty ="n")
    
    legend(legend = c(paste0("MPD obs = ", round(MPD_per_mimic_obs[i], 3)),
                      paste0("p = ", p_value)),
           x = "left", cex = 1, bty ="n")
    
    dev.off()
  }
  save(MPD_ring_summary_table, file = paste0("./outputs/Niche_evolution/MPD_ring_summary_table.RData"))
  
  cat(paste0("N° ",i, "/",length(mimicry.list)," - ",ring, " - Done \n"))
}

View(MPD_ring_summary_table)
write.csv2(MPD_ring_summary_table, file = "./tables/MPD_ring_summary_table.csv")





# ######## Version on species ########
# 
# ## Compute Phylogenetic/patristic distance matrix
# phylo.dist.mat <- cophenetic.phylo(x = phylo.Ithomiini)
# 
# ## Compute matrix of pairwise weights as Jaccard index (i.e., proportion of mimicry ring in commun out of all different mimicry rings of the two species)
# # Generate incidence table of species and mimicry rings
# mimic.incidence.table <- matrix(nrow = 339, ncol = 44)
# mimic.list <- as.character(unique(list.models$Mimicry.model))
# for (i in 1:nrow(list.sp)) { # 339 sp = rows
#   sp <- as.character(list.sp$Sp_full[i])
#   for (j in 1:length(mimic.list)) { # 44 mimicry rings = columns
#     mimic <- mimic.list[j]
#     temp <- as.character(list.models$Mimicry.model[list.models$Sp_full==sp]) # Extract mimicry ring of all units for the i species
#     mimic.incidence.table[i,j] <- as.numeric(any(temp==mimic)) # 0 if mimicry ring not present ; 1 if found
#   }
# }
# save(mimic.incidence.table, file = paste0(internal.wd,"/Niche_evolution/Mimicry_Incidence_Table.RData"))
# 
# # Compute the Jaccard index between pairs of species
# library(ade4)
# jaccard.index <- dist.binary(df = mimic.incidence.table, method = 1)
# # Computed as distance, so need to be back-transformed to the similarity index
# jaccard.index <- as.matrix(jaccard.index)
# jaccard.index <- jaccard.index*jaccard.index
# jaccard.index <- 1-jaccard.index
# mimicry.similarity.weights <- as.dist(jaccard.index)
# save(mimicry.similarity.weights, file = paste0(internal.wd,"/Niche_evolution/Jaccard_index.RData"))
# 
# # Compute the global observed MPD weighted by mimicry similarity
# mean(as.dist(phylo.dist.mat)) # MPD globale toutes paires (sans pond?ration par la similarit? de pattern mim?tiques) = 37.01 Ma
# Global_MPD_obs <- weighted.mean(x = as.dist(phylo.dist.mat), w = mimicry.similarity.weights) # 33.68 Ma
# 
# save(Global_MPD_obs, file = paste0(internal.wd,"/Niche_evolution/Global_MPD_obs.RData"))
# 
# ## Compute the MPD per mimicry ring
# MPD_per_mimic_obs <- rep(NA, length(mimic.list))
# for (i in 1:length(mimic.list)) {
#   # Compute weights focusing only on the mimicry ring targeted => weight = proportion of the targeted mimicry ring among all the units of the pair of species
#   mimic <- mimic.list[i]
#   mimic.weights <- matrix(nrow = nrow(mimic.incidence.table), ncol = nrow(mimic.incidence.table))
#   for (j in 1:nrow(mimic.incidence.table)) {
#     for (k in 1:nrow(mimic.incidence.table)) {
#       if (mimic.incidence.table[j,i]+mimic.incidence.table[k,i]==2) { # Case when the two species share the targeted mimicry ring
#         sp1 <- which(mimic.incidence.table[j,]==1) # Index of mimicry ring of sp1
#         sp2 <- which(mimic.incidence.table[k,]==1) # Index of mimicry ring of sp2
#         tot_rings <- length(union(sp1,sp2)) # Number of different mimicry ring in total among the two species
#         mimic.weights[j,k] <- 1/tot_rings # Weight is the proportion of the targeted mimicry ring among all the units of the two species
#       }else{ # Case when the two species don't share the targeted mimicry ring
#         mimic.weights[j,k] <- 0 # Weight is null
#       }
#     }
#   }
#   save(mimic.weights, file = paste0(internal.wd,"/Niche_evolution/mimic.weights_,",mimic,".RData"))
#   
#   # Compute the global observed MPD weighted by mimicry similarity for one mimicry ring
#   MPD_obs <- weighted.mean(x = as.dist(phylo.dist.mat), w = as.dist(mimic.weights)) # 33.68 Ma
#   if (!is.nan(MPD_obs)) { # When mimicry ring is found only in one species, the computation of the weighted mean will lead to Inf (NaN)
#     MPD_per_mimic_obs[i] <- MPD_obs # Store value only if it is not NaN
#   }
#   print(i)
# }
# names(MPD_per_mimic_obs) <- mimic.list
# MPD_per_mimic_obs
# save(MPD_per_mimic_obs, file = paste0(internal.wd,"/Niche_evolution/MPD_per_mimic_obs.RData"))
# 
# 
# ### Simulate MPD under null hypothesis by permutation
# 
# ## Start the loop
# Global_MPD_null <-  NA
# MPD_per_mimic_null <- data.frame(matrix(nrow = 1000, ncol = length(mimic.list)))
# for (l in 1:1000) {
#   
#   ## Suffle mimicry ring among units
#   shuffle.list.unit <- list.models
#   shuffle.list.unit$Mimicry.model <- sample(as.character(shuffle.list.unit$Mimicry.model))
#   
#   ## Compute matrix of pairwise weights as Jaccard index (i.e., proportion of mimicry ring in commun out of all different mimicry rings of the two species)
#   # Generate incidence table of species and mimicry rings
#   simul.mimic.incidence.table <- matrix(nrow = 339, ncol = 44)
#   mimic.list <- as.character(unique(list.models$Mimicry.model))
#   for (i in 1:nrow(list.sp)) { # 339 sp = rows
#     sp <- as.character(list.sp$Sp_full[i])
#     for (j in 1:length(mimic.list)) { # 44 mimicry rings = columns
#       mimic <- mimic.list[j]
#       temp <- as.character(shuffle.list.unit$Mimicry.model[shuffle.list.unit$Sp_full==sp]) # Extract mimicry ring of all units for the i species
#       simul.mimic.incidence.table[i,j] <- as.numeric(any(temp==mimic)) # 0 if mimicry ring not present ; 1 if found
#     }
#   }
#   # Compute the Jaccard index between pairs of species
#   jaccard.index <- dist.binary(df = simul.mimic.incidence.table, method = 1)
#   # Computed as distance, so need to be back-transformed to the similarity index
#   jaccard.index <- as.matrix(jaccard.index)
#   jaccard.index <- jaccard.index*jaccard.index
#   jaccard.index <- 1-jaccard.index
#   simul.mimicry.similarity.weights <- as.dist(jaccard.index)
#   
#   # Compute the global observed MPD weighted by mimicry similarity
#   mean(as.dist(phylo.dist.mat)) # MPD globale toutes paires (sans pond?ration par la similarit? de pattern mim?tiques) = 37.01 Ma
#   Global_MPD_simul <- weighted.mean(x = as.dist(phylo.dist.mat), w = simul.mimicry.similarity.weights) # 33.68 Ma
#   Global_MPD_null[l] <- Global_MPD_simul
#   
#   save(Global_MPD_null, file = paste0(internal.wd,"/Niche_evolution/Global_MPD_null.RData"))
#   
#   ## Compute the MPD per mimicry ring
#   MPD_per_mimic_simul <- rep(NA, length(mimic.list))
#   for (i in 1:length(mimic.list)) {
#     # Compute weights focusing only on the mimicry ring targeted => weight = proportion of the targeted mimicry ring among all the units of the pair of species
#     mimic <- mimic.list[i]
#     simul.mimic.weights <- matrix(nrow = nrow(simul.mimic.incidence.table), ncol = nrow(simul.mimic.incidence.table))
#     for (j in 1:nrow(simul.mimic.incidence.table)) {
#       for (k in 1:nrow(simul.mimic.incidence.table)) {
#         if (simul.mimic.incidence.table[j,i]+simul.mimic.incidence.table[k,i]==2) { # Case when the two species share the targeted mimicry ring
#           sp1 <- which(simul.mimic.incidence.table[j,]==1) # Index of mimicry ring of sp1
#           sp2 <- which(simul.mimic.incidence.table[k,]==1) # Index of mimicry ring of sp2
#           tot_rings <- length(union(sp1,sp2)) # Number of different mimicry ring in total among the two species
#           simul.mimic.weights[j,k] <- 1/tot_rings # Weight is the proportion of the targeted mimicry ring among all the units of the two species
#         }else{ # Case when the two species don't share the targeted mimicry ring
#           simul.mimic.weights[j,k] <- 0 # Weight is null
#         }
#       }
#     }
#     
#     # Compute the global observed MPD weighted by mimicry similarity for one mimicry ring
#     MPD_simul <- weighted.mean(x = as.dist(phylo.dist.mat), w = as.dist(simul.mimic.weights)) # 33.68 Ma
#     if (!is.nan(MPD_simul)) { # When mimicry ring is found only in one species, the computation of the weighted mean will lead to Inf (NaN)
#       MPD_per_mimic_simul[i] <- MPD_simul # Store value only if it is not NaN
#     }
#   }
#   MPD_per_mimic_null[l,] <- MPD_per_mimic_simul
#   save(MPD_per_mimic_null, file = paste0(internal.wd,"/Niche_evolution/MPD_per_mimic_null.RData")) 
#   cat(paste0(Sys.time(), " - Simul n?", l, "\n"))
# }
# names(MPD_per_mimic_null) <- mimic.list
# 
# 
# # Plot de l'histo global
# load(file = paste0(internal.wd,"/Niche_evolution/Global_MPD_obs.RData"))
# load(file = paste0(internal.wd,"/Niche_evolution/Global_MPD_null.RData"))
# 
# tiff(filename = paste0(internal.wd,"/Niche_evolution/Plots/Global_MPD_null.tiff"))
# original_int_margins <- par()$mar
# par(mar = c(5.1,5,4.1,2.1))
# hist(c(Global_MPD_null, Global_MPD_obs), 30, freq = TRUE, col = "gray", 
#      main = "Distribution of the Mean pairwise Phylogenetic Distance \n under the null hypothesis", 
#      xlab = "Mean pairwise Phylogenetic Distance (My)",
#      cex.axis = 1.7, cex.main = 1, cex.lab = 1.7)
# arrows(Global_MPD_obs + 0.0001, 80, Global_MPD_obs+ 0.0001, 10, length = 0.1, lwd = 2)
# abline(v = mean(c(Global_MPD_null, Global_MPD_obs)), lty = 2)
# abline(v = quantile(c(Global_MPD_null, Global_MPD_obs), 0.05), lty = 2, col = "red")
# legend(legend = c(paste0("Mean = ", round(mean(c(Global_MPD_null, Global_MPD_obs)),3)), 
#                   paste0("CI 5% = ", round(quantile(c(Global_MPD_null, Global_MPD_obs), 0.05),3))), 
#        x = "topleft", lty = 2, lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")
# legend(legend = c(paste0("MPD obs = ", round(Global_MPD_obs, 3)),
#                   paste0("         p = 0.001")),
#        x = "left", cex = 1, bty ="n", xjust = 1)
# par(mar = original_int_margins)
# dev.off()
# 
# 
# # Plot des histo per mimic ring
# load(file = paste0(internal.wd,"/Niche_evolution/MPD_per_mimic_obs.RData"))
# load(file = paste0(internal.wd,"/Niche_evolution/MPD_per_mimic_null.RData"))
# 
# mimicry.list <- as.character(unique(list.models$Mimicry.model))
# 
# MPD_ring_summary_table <- as.data.frame(matrix(ncol = 8, nrow = 44, data = NA))
# names(MPD_ring_summary_table) <- c("ring", "N_units", "N_pairs", "MPD_obs", "mean_MPD", "MPD_2.5", "MPD_97.5", "p_value")
# 
# for (i in 1:length(mimicry.list)) {
#   MPD_ring_summary_table$ring[i] <- mimic <- mimicry.list[i]
#   MPD_ring_summary_table$N_units[i] <- N_units <- sum(list.models$Mimicry.model==mimic)
#   MPD_ring_summary_table$N_pairs[i] <- N_pairs <- N_units*(N_units-1)/2
#   if (is.na(MPD_per_mimic_obs[i])) {
#     jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/MPD_null_",mimic,".jpeg"), quality =100)
#     plot(1:100,1:100, type = "n", main = paste0("Distribution of Mean Phylogenetic Distance \n of ", mimic, " species \n under null Hypothesis"))
#     text(x = 50, y = 50, labels = "Only one unit for this mimicry ring \n No pair available for index computation")
#     dev.off()
#   }else{
#     MPD_ring_summary_table$MPD_obs[i] <- round(MPD_per_mimic_obs[i],3)
#     MPD_ring_summary_table$mean_MPD[i] <- round(mean(MPD_per_mimic_null[,i], na.rm = T),3)
#     MPD_ring_summary_table$MPD_2.5[i] <- round(quantile(MPD_per_mimic_null[,i], 0.025, na.rm = T),3)
#     MPD_ring_summary_table$MPD_97.5[i] <- round(quantile(MPD_per_mimic_null[,i], 0.975, na.rm = T),3)
#     MPD_ring_summary_table$p_value[i] <- round(ecdf(x = c(MPD_per_mimic_null[,i], MPD_per_mimic_obs[i]))(MPD_per_mimic_obs[i]),3)
#     jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/MPD_null_",mimic,".jpeg"), quality =100)
#     histo.save <- hist(MPD_per_mimic_null[,i], breaks = seq(from = floor(min(min(MPD_per_mimic_null[,i], na.rm = T), MPD_per_mimic_obs[i])), to = ceiling(max(max(MPD_per_mimic_null[,i], na.rm = T), MPD_per_mimic_obs[i])), by = 1), col = "gray", xlab = "Mean Phylogenetic Distance", main = paste0("Distribution of the Mean Phylogenetic distances \n of ",mimic," species \n under the null Hypothesis"))
#     hist(MPD_per_mimic_null[,i], breaks = seq(from = floor(min(min(MPD_per_mimic_null[,i], na.rm = T), MPD_per_mimic_obs[i])), to = ceiling(max(max(MPD_per_mimic_null[,i], na.rm = T), MPD_per_mimic_obs[i])), by = 1), col = "gray", xlab = "Mean Phylogenetic Distance", main = paste0("Distribution of the Mean Phylogenetic distances \n of ",mimic," species \n under the null Hypothesis"))
#     arrows(MPD_per_mimic_obs[i], max(histo.save$counts)/3, MPD_per_mimic_obs[i], max(histo.save$counts)/30, length = 0.1, lwd = 2)
#     abline(v = mean(MPD_per_mimic_null[,i]), lty = 2, lwd = 1)
#     abline(v = quantile(MPD_per_mimic_null[,i], 0.025, na.rm = T), lty = 2, col = "red")
#     abline(v = quantile(MPD_per_mimic_null[,i], 0.975, na.rm = T), lty = 2, col = "red")
#     legend(legend = c(paste0("N units = ", N_units), 
#                       paste0("N pairs = ", N_pairs)),
#            x = "topright", cex = 1, bty ="n")
#     legend(legend = c(paste0("Mean = ", round(mean(MPD_per_mimic_null[,i]),3)), 
#                       paste0("CI 2.5% = ", round(quantile(MPD_per_mimic_null[,i], 0.025, na.rm = T),3)),
#                       paste0("CI 97.5% = ", round(quantile(MPD_per_mimic_null[,i], 0.975, na.rm = T),3))), 
#            x = "topleft", cex = 1, bty ="n")
#     legend(legend = c(paste0("MPD obs = ", round(MPD_per_mimic_obs[i], 3)),
#                       paste0("p = ", round(ecdf(x = MPD_per_mimic_null[,i])(MPD_per_mimic_obs[i]),3))),
#            x = "left", cex = 1, bty ="n")
#     dev.off()
#   }
#   save(MPD_ring_summary_table, file = paste0(internal.wd,"/Niche_evolution/MPD_ring_summary_table.RData"))
#   cat(paste0(i, " - ",mimic, " - Done \n"))
# }
# 
# View(MPD_ring_summary_table)
# 



##### 3/ Test for phylogenetic signal in the evolution of climatic niche (Kmult) #####

# Must be done at species level because the matrix of phylogenetic distances need to be inverted

source(paste0("./scripts/Test.Kmult.R")) # From Adams, 2014

load(file = paste0("./input_data/Phylogenies/Final_phylogeny.RData"))
load(file = paste0("./outputs/Niche_evolution/sp_env_table_339.RData"))

# Check the order is the same
identical(row.names(sp_env_table_339), phylo.Ithomiini$tip.label)
sp_env_mat_339 <- as.matrix(sp_env_table_339[,-1]) ; row.names(sp_env_mat_339) <- sp_env_table_339$Sp_full

# Compute the Kmult test
test.Kmult <- Test.Kmult(x = sp_env_mat_339, phy = phylo.Ithomiini, iter = 999)
save(test.Kmult, file = paste0("./outputs/Niche_evolution/Phylo_signal/test.Kmult.Rdata"))
test.Kmult$phy.signal ; test.Kmult$pvalue

# Plot null distri
pdf(file = "./graphs/Niche_evolution/Phylo_signal/Kmult_test.pdf", height = 6, width = 7)

hist(test.Kmult$K.simul.val, 30, freq = TRUE, col = "gray", 
     main = bquote('Distribution of' ~K[mult]~ 'of the climatic niche\nunder null Hypothesis'),
     xlab = bquote(~K[mult]~ 'index'),
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5, lwd = 2)

arrows(test.Kmult$phy.signal, 50, test.Kmult$phy.signal, 5, length = 0.1, lwd = 2)
abline(v = mean(test.Kmult$K.simul.val), lty = 2, lwd = 2)
abline(v = quantile(test.Kmult$K.simul.val, 0.95), lty = 2, lwd = 2, col = "red")

legend(legend = c(paste0("Mean = ", round(mean(test.Kmult$K.simul.val),3)), 
                  paste0("CI 95% = ", round(quantile(test.Kmult$K.simul.val, 0.95),3))), 
       x = "topleft", inset = c(0.02, 0.00), lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")

legend(legend = bquote("  " ~K[mult]~ "obs =" ~ .(format(round(test.Kmult$phy.signal, 3), nsmall = 3))),
       x = "right", cex = 1.2, bty ="n",
       inset = c(0.03, 0.00))
legend(legend = c("",
                  paste0("       p = ", test.Kmult$pvalue)),
       x = "right", cex = 1.2, bty ="n", y.intersp = 2.1,
       inset = c(0.11, 0.00))

legend(legend = as.expression(bquote(bold("B"))), 
       x = "topright", inset = c(0.05, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

dev.off()


##### 4/ Plot both MPD and Kmult #####

load(file = paste0("./outputs/Niche_evolution/Phylo_signal/Global_MPD_obs.RData"))
load(file = paste0("./outputs/Niche_evolution/Phylo_signal/Global_MPD_null.RData"))

load(file = paste0("./outputs/Niche_evolution/Phylo_signal/test.Kmult.Rdata"))

round(Global_MPD_obs, 2) # MPD obs = 34.43
format(round(test.Kmult$phy.signal, 3), nsmall = 3) # Kmult obs = 0.120

pdf(file = "./graphs/Niche_evolution/Phylo_signal/MPD_&_Kmult_both_plots.pdf", height = 6.4, width = 14)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,5,1.5,1), mfrow = c(1,2))

hist(x = c(Global_MPD_null, Global_MPD_obs), 
     breaks = 20, 
     freq = TRUE, col = "gray",
     # xlim = c(34.5, 38.0),
     # ylim = c(0, 200),
     # main = "Distribution of the Mean pairwise Phylogenetic Distance \nof co-mimetic OMUs \nunder the null hypothesis",
     main = "",
     xlab = "Mean pairwise Phylogenetic Distance [My]",
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.2, lwd = 2)

arrows(Global_MPD_obs + 0.055, 120, Global_MPD_obs+ 0.055, 5, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = mean(c(Global_MPD_null, Global_MPD_obs)), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(c(Global_MPD_null, Global_MPD_obs), 0.05), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("Mean = ", round(mean(c(Global_MPD_null, Global_MPD_obs)),2)), 
                  paste0("CI 5% = ", round(quantile(c(Global_MPD_null, Global_MPD_obs), 0.05),2))), 
       x = "topleft", inset = c(0.00, 0.15), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.6, bty = "n")

legend(legend = c(as.expression(bquote(bold("MPD obs = 34.43"))),
                  as.expression(bquote(bold("p = 0.001")))),
       x = "bottomleft", inset = c(-0.035, 0.40),
       cex = 1.6, bty ="n", xjust = 1)

legend(legend = as.expression(bquote(bold("A"))), 
       x = "topleft", inset = c(-0.03, -0.03), xjust = 0.5,
       cex = 1.8, bty ="n")

hist(c(test.Kmult$K.simul.val, test.Kmult$phy.signal),
     breaks = seq(0, 0.22, 0.01),
     freq = TRUE, col = "gray", 
     # main = bquote('Distribution of' ~K[mult]~ 'of the climatic niche\nunder null Hypothesis'),
     main = "",
     xlim = c(0.00, 0.22),
     xlab = bquote(~K[mult]~ 'index'),
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.2, lwd = 2)

arrows(test.Kmult$phy.signal, 95, test.Kmult$phy.signal, 5, length = 0.1, lwd = 2)
abline(v = mean(test.Kmult$K.simul.val), lty = 2, lwd = 2)
abline(v = quantile(test.Kmult$K.simul.val, 0.95), lty = 2, lwd = 2, col = "red")

legend(legend = c(paste0("Mean = ", round(mean(test.Kmult$K.simul.val),3)), 
                  paste0("CI 95% = ", round(quantile(test.Kmult$K.simul.val, 0.95),3))), 
       x = "topright", inset = c(0.02, 0.15), lty = 2 , lwd = 2,
       col = c("black", "red"), cex = 1.6, bty ="n")

# legend(legend = bquote("  " ~K[mult]~ "obs =" ~ .(format(round(test.Kmult$phy.signal, 3), nsmall = 3))),
#        x = "bottomright", cex = 1.3, bty ="n",
#        inset = c(0.018, 0.435))
# legend(legend = c("",
#                   paste0("p = ", test.Kmult$pvalue)),
#        x = "bottomright", cex = 1.3, bty ="n", y.intersp = 2.1,
#        inset = c(0.142, 0.35))

legend(legend = as.expression(bquote(bold(K[mult]~ "obs = 0.120"))),
       x = "bottomright", cex = 1.6, bty ="n",
       inset = c(0.10, 0.435))
legend(legend = c("",
                  as.expression(bquote(bold(paste("p = 0.013"))))),
       x = "bottomright", cex = 1.6, bty ="n", y.intersp = 2.1,
       inset = c(0.26, 0.345))

legend(legend = as.expression(bquote(bold("B"))), 
       x = "topright", inset = c(0.05, -0.03), xjust = 0.5,
       cex = 1.9, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins, mfrow = c(1,1))

dev.off()

