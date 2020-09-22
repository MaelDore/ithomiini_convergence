##### Convergence de l'?volution de la niche climatique chez les co-mimes #####

### Preparation ###

# Effacer l'environnement
rm(list = ls())

library(raster)

# Set Wd
internal.wd <- "F:/Documents/Etudes/Fac/Cours/TROPIMUNDO/Stage_S4/Projet_Papillons/R_codes"

# Load phylogeny, sp list and unit list
load(file = paste0(internal.wd, "/Phylogenies/Final_phylogeny.RData"))
load(file = paste0(internal.wd,"/list.sp.RData"))
load(file = paste0(internal.wd,"/list.models.RData"))

# Modeling resolution
res <-  "5"

# Extract only the 339 species included in the phylogeny from list.sp, list.models and the sp.stack
list.sp <- list.sp[list.sp$Sp_full %in% phylo.Ithomiini$tip.label,] # 339 species
list.models <- list.models[list.models$Sp_full %in% phylo.Ithomiini$tip.label,] # 719 units


# Extraction des donn?es climatiques moyennes des units
Unit_mean_env_data <- data.frame() # Create df to store final data
for ( i in 1:nrow(list.models)){  # Per unit
  unit <- as.character(list.models$Tag.model[i])
  load(file = paste0(internal.wd,"/Species_data/",res,"/Env_Stacks/Env_stack_",unit,".RData")) # load env. stack of the unit
  load(paste0(internal.wd,"/Species_data/",res,"/Spatial_Points_Objects/occurrences_", unit,".RData")) # load occurrences spatial points of the species
  all.occ.env.data <- extract(x = unit.env, y = unit.points) # Extract env. data for each occurrence
  mean.env.data <- apply(X = all.occ.env.data, MARGIN = 2, FUN = mean, na.rm = T) # Compute mean for each variable
  Unit_mean_env_data <- rbind(Unit_mean_env_data, mean.env.data) # Store results in the final df
  print(i)
}
names(Unit_mean_env_data) <- names(mean.env.data) # Retrieve names of env. variables

list.models <- cbind(list.models, Unit_mean_env_data) # Bind to original list of units available in the phylogeny
save(list.models, file = paste0(internal.wd,"/unit.env.table.RData"))

# Extraction des donn?es climatiques moyennes des esp?ces
Sp_mean_env_data <- data.frame() # Create final df to store data
for (j in 1:nrow(list.sp)) { # Per species
  sp <- as.character(list.sp$Sp_full[j])
  list.unit.per.sp <- list.models[list.models$Sp_full==sp,] # Extract only units of this sp
  all.occ.env.data <- data.frame() # Generate clean df to store data for the new species
  for ( i in 1:nrow(list.unit.per.sp)){ # Per units
    unit <- as.character(list.unit.per.sp$Tag.model[i])
    load(file = paste0(internal.wd,"/Species_data/",res,"/Env_Stacks/Env_stack_",unit,".RData")) # load env. stack
    load(paste0(internal.wd,"/Species_data/",res,"/Spatial_Points_Objects/occurrences_", unit,".RData")) # load occurrences spatial points of the species
    patch <- extract(x = unit.env, y = unit.points)
    all.occ.env.data <- rbind(all.occ.env.data, patch) # Store env data of occurrences of all units of the species
  }
  mean.env.data <- apply(X = all.occ.env.data, MARGIN = 2, FUN = mean, na.rm = T) # Compute mean for each variable
  Sp_mean_env_data <- rbind(Sp_mean_env_data, mean.env.data) # Store final data
  print(j)
}
names(Sp_mean_env_data) <- names(mean.env.data) # Retrieve env. data names

list.sp <- cbind(list.sp, Sp_mean_env_data) # Bind with sp list df

# Reorder list.sp following Phylogeny
order.index <- NA
for (i in 1:length(phylo.Ithomiini$tip.label)) {
  sp <- as.character(phylo.Ithomiini$tip.label[i]) 
  order.index[i] <- which(list.sp$Sp_full==sp)
}

list.sp <- list.sp[order.index,]
row.names(list.sp) <- as.character(list.sp$Sp_full)

save(list.sp, file = paste0(internal.wd,"/sp.env.table.RData"))

# Load the env. variables df
load(file = paste0(internal.wd,"/sp.env.table.RData"))


###### Approche multivari?e ? 5 dimensions ######

### MPD entre co-mimic = signal phylog?n?tique des cercles mim?tiques

## Compute Phylogenetic/patristic distance matrix
phylo.dist.mat <- cophenetic.phylo(x = phylo.Ithomiini)

## Compute matrix of pairwise weights as Jaccard index (i.e., proportion of mimicry ring in commun out of all different mimicry rings of the two species)
# Generate incidence table of species and mimicry rings
mimic.incidence.table <- matrix(nrow = 339, ncol = 44)
mimic.list <- as.character(unique(list.models$Mimicry.model))
for (i in 1:nrow(list.sp)) { # 339 sp = rows
  sp <- as.character(list.sp$Sp_full[i])
  for (j in 1:length(mimic.list)) { # 44 mimicry rings = columns
    mimic <- mimic.list[j]
    temp <- as.character(list.models$Mimicry.model[list.models$Sp_full==sp]) # Extract mimicry ring of all units for the i species
    mimic.incidence.table[i,j] <- as.numeric(any(temp==mimic)) # 0 if mimicry ring not present ; 1 if found
  }
}
save(mimic.incidence.table, file = paste0(internal.wd,"/Niche_evolution/Mimicry_Incidence_Table.RData"))

# Compute the Jaccard index between pairs of species
library(ade4)
jaccard.index <- dist.binary(df = mimic.incidence.table, method = 1)
# Computed as distance, so need to be back-transformed to the similarity index
jaccard.index <- as.matrix(jaccard.index)
jaccard.index <- jaccard.index*jaccard.index
jaccard.index <- 1-jaccard.index
mimicry.similarity.weights <- as.dist(jaccard.index)
save(mimicry.similarity.weights, file = paste0(internal.wd,"/Niche_evolution/Jaccard_index.RData"))

# Compute the global observed MPD weighted by mimicry similarity
mean(as.dist(phylo.dist.mat)) # MPD globale toutes paires (sans pond?ration par la similarit? de pattern mim?tiques) = 37.01 Ma
Global_MPD_obs <- weighted.mean(x = as.dist(phylo.dist.mat), w = mimicry.similarity.weights) # 33.68 Ma

save(Global_MPD_obs, file = paste0(internal.wd,"/Niche_evolution/Global_MPD_obs.RData"))

## Compute the MPD per mimicry ring
MPD_per_mimic_obs <- rep(NA, length(mimic.list))
for (i in 1:length(mimic.list)) {
  # Compute weights focusing only on the mimicry ring targeted => weight = proportion of the targeted mimicry ring among all the units of the pair of species
  mimic <- mimic.list[i]
  mimic.weights <- matrix(nrow = nrow(mimic.incidence.table), ncol = nrow(mimic.incidence.table))
  for (j in 1:nrow(mimic.incidence.table)) {
    for (k in 1:nrow(mimic.incidence.table)) {
      if (mimic.incidence.table[j,i]+mimic.incidence.table[k,i]==2) { # Case when the two species share the targeted mimicry ring
        sp1 <- which(mimic.incidence.table[j,]==1) # Index of mimicry ring of sp1
        sp2 <- which(mimic.incidence.table[k,]==1) # Index of mimicry ring of sp2
        tot_rings <- length(union(sp1,sp2)) # Number of different mimicry ring in total among the two species
        mimic.weights[j,k] <- 1/tot_rings # Weight is the proportion of the targeted mimicry ring among all the units of the two species
      }else{ # Case when the two species don't share the targeted mimicry ring
        mimic.weights[j,k] <- 0 # Weight is null
      }
    }
  }
  save(mimic.weights, file = paste0(internal.wd,"/Niche_evolution/mimic.weights_,",mimic,".RData"))
  
  # Compute the global observed MPD weighted by mimicry similarity for one mimicry ring
  MPD_obs <- weighted.mean(x = as.dist(phylo.dist.mat), w = as.dist(mimic.weights)) # 33.68 Ma
  if (!is.nan(MPD_obs)) { # When mimicry ring is found only in one species, the computation of the weighted mean will lead to Inf (NaN)
    MPD_per_mimic_obs[i] <- MPD_obs # Store value only if it is not NaN
  }
  print(i)
}
names(MPD_per_mimic_obs) <- mimic.list
MPD_per_mimic_obs
save(MPD_per_mimic_obs, file = paste0(internal.wd,"/Niche_evolution/MPD_per_mimic_obs.RData"))


### Simulate MPD under null hypothesis by permutation

## Start the loop
Global_MPD_null <-  NA
MPD_per_mimic_null <- data.frame(matrix(nrow = 1000, ncol = length(mimic.list)))
for (l in 1:1000) {
  
  ## Suffle mimicry ring among units
  shuffle.list.unit <- list.models
  shuffle.list.unit$Mimicry.model <- sample(as.character(shuffle.list.unit$Mimicry.model))

  ## Compute matrix of pairwise weights as Jaccard index (i.e., proportion of mimicry ring in commun out of all different mimicry rings of the two species)
  # Generate incidence table of species and mimicry rings
  simul.mimic.incidence.table <- matrix(nrow = 339, ncol = 44)
  mimic.list <- as.character(unique(list.models$Mimicry.model))
  for (i in 1:nrow(list.sp)) { # 339 sp = rows
    sp <- as.character(list.sp$Sp_full[i])
    for (j in 1:length(mimic.list)) { # 44 mimicry rings = columns
      mimic <- mimic.list[j]
      temp <- as.character(shuffle.list.unit$Mimicry.model[shuffle.list.unit$Sp_full==sp]) # Extract mimicry ring of all units for the i species
      simul.mimic.incidence.table[i,j] <- as.numeric(any(temp==mimic)) # 0 if mimicry ring not present ; 1 if found
    }
  }
  # Compute the Jaccard index between pairs of species
  jaccard.index <- dist.binary(df = simul.mimic.incidence.table, method = 1)
  # Computed as distance, so need to be back-transformed to the similarity index
  jaccard.index <- as.matrix(jaccard.index)
  jaccard.index <- jaccard.index*jaccard.index
  jaccard.index <- 1-jaccard.index
  simul.mimicry.similarity.weights <- as.dist(jaccard.index)
  
  # Compute the global observed MPD weighted by mimicry similarity
  mean(as.dist(phylo.dist.mat)) # MPD globale toutes paires (sans pond?ration par la similarit? de pattern mim?tiques) = 37.01 Ma
  Global_MPD_simul <- weighted.mean(x = as.dist(phylo.dist.mat), w = simul.mimicry.similarity.weights) # 33.68 Ma
  Global_MPD_null[l] <- Global_MPD_simul
  
  save(Global_MPD_null, file = paste0(internal.wd,"/Niche_evolution/Global_MPD_null.RData"))
  
  ## Compute the MPD per mimicry ring
  MPD_per_mimic_simul <- rep(NA, length(mimic.list))
  for (i in 1:length(mimic.list)) {
    # Compute weights focusing only on the mimicry ring targeted => weight = proportion of the targeted mimicry ring among all the units of the pair of species
    mimic <- mimic.list[i]
    simul.mimic.weights <- matrix(nrow = nrow(simul.mimic.incidence.table), ncol = nrow(simul.mimic.incidence.table))
    for (j in 1:nrow(simul.mimic.incidence.table)) {
      for (k in 1:nrow(simul.mimic.incidence.table)) {
        if (simul.mimic.incidence.table[j,i]+simul.mimic.incidence.table[k,i]==2) { # Case when the two species share the targeted mimicry ring
          sp1 <- which(simul.mimic.incidence.table[j,]==1) # Index of mimicry ring of sp1
          sp2 <- which(simul.mimic.incidence.table[k,]==1) # Index of mimicry ring of sp2
          tot_rings <- length(union(sp1,sp2)) # Number of different mimicry ring in total among the two species
          simul.mimic.weights[j,k] <- 1/tot_rings # Weight is the proportion of the targeted mimicry ring among all the units of the two species
        }else{ # Case when the two species don't share the targeted mimicry ring
          simul.mimic.weights[j,k] <- 0 # Weight is null
        }
      }
    }
    
    # Compute the global observed MPD weighted by mimicry similarity for one mimicry ring
    MPD_simul <- weighted.mean(x = as.dist(phylo.dist.mat), w = as.dist(simul.mimic.weights)) # 33.68 Ma
    if (!is.nan(MPD_simul)) { # When mimicry ring is found only in one species, the computation of the weighted mean will lead to Inf (NaN)
      MPD_per_mimic_simul[i] <- MPD_simul # Store value only if it is not NaN
    }
  }
  MPD_per_mimic_null[l,] <- MPD_per_mimic_simul
  save(MPD_per_mimic_null, file = paste0(internal.wd,"/Niche_evolution/MPD_per_mimic_null.RData")) 
  cat(paste0(Sys.time(), " - Simul n?", l, "\n"))
}
names(MPD_per_mimic_null) <- mimic.list


# Plot de l'histo global
load(file = paste0(internal.wd,"/Niche_evolution/Global_MPD_obs.RData"))
load(file = paste0(internal.wd,"/Niche_evolution/Global_MPD_null.RData"))

tiff(filename = paste0(internal.wd,"/Niche_evolution/Plots/Global_MPD_null.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
hist(c(Global_MPD_null, Global_MPD_obs), 30, freq = TRUE, col = "gray", 
     main = "Distribution of the Mean pairwise Phylogenetic Distance \n under the null hypothesis", 
     xlab = "Mean pairwise Phylogenetic Distance (My)",
     cex.axis = 1.7, cex.main = 1, cex.lab = 1.7)
arrows(Global_MPD_obs + 0.0001, 80, Global_MPD_obs+ 0.0001, 10, length = 0.1, lwd = 2)
abline(v = mean(c(Global_MPD_null, Global_MPD_obs)), lty = 2)
abline(v = quantile(c(Global_MPD_null, Global_MPD_obs), 0.05), lty = 2, col = "red")
legend(legend = c(paste0("Mean = ", round(mean(c(Global_MPD_null, Global_MPD_obs)),3)), 
                  paste0("CI 5% = ", round(quantile(c(Global_MPD_null, Global_MPD_obs), 0.05),3))), 
       x = "topleft", lty = 2, lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")
legend(legend = c(paste0("MPD obs = ", round(Global_MPD_obs, 3)),
                  paste0("         p = 0.001")),
       x = "left", cex = 1, bty ="n", xjust = 1)
par(mar = original_int_margins)
dev.off()


# Plot des histo per mimic ring
load(file = paste0(internal.wd,"/Niche_evolution/MPD_per_mimic_obs.RData"))
load(file = paste0(internal.wd,"/Niche_evolution/MPD_per_mimic_null.RData"))

mimicry.list <- as.character(unique(list.models$Mimicry.model))

MPD_ring_summary_table <- as.data.frame(matrix(ncol = 8, nrow = 44, data = NA))
names(MPD_ring_summary_table) <- c("ring", "N_units", "N_pairs", "MPD_obs", "mean_MPD", "MPD_2.5", "MPD_97.5", "p_value")

for (i in 1:length(mimicry.list)) {
  MPD_ring_summary_table$ring[i] <- mimic <- mimicry.list[i]
  MPD_ring_summary_table$N_units[i] <- N_units <- sum(list.models$Mimicry.model==mimic)
  MPD_ring_summary_table$N_pairs[i] <- N_pairs <- N_units*(N_units-1)/2
  if (is.na(MPD_per_mimic_obs[i])) {
    jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/MPD_null_",mimic,".jpeg"), quality =100)
    plot(1:100,1:100, type = "n", main = paste0("Distribution of Mean Phylogenetic Distance \n of ", mimic, " species \n under null Hypothesis"))
    text(x = 50, y = 50, labels = "Only one unit for this mimicry ring \n No pair available for index computation")
    dev.off()
  }else{
    MPD_ring_summary_table$MPD_obs[i] <- round(MPD_per_mimic_obs[i],3)
    MPD_ring_summary_table$mean_MPD[i] <- round(mean(MPD_per_mimic_null[,i], na.rm = T),3)
    MPD_ring_summary_table$MPD_2.5[i] <- round(quantile(MPD_per_mimic_null[,i], 0.025, na.rm = T),3)
    MPD_ring_summary_table$MPD_97.5[i] <- round(quantile(MPD_per_mimic_null[,i], 0.975, na.rm = T),3)
    MPD_ring_summary_table$p_value[i] <- round(ecdf(x = c(MPD_per_mimic_null[,i], MPD_per_mimic_obs[i]))(MPD_per_mimic_obs[i]),3)
    jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/MPD_null_",mimic,".jpeg"), quality =100)
    histo.save <- hist(MPD_per_mimic_null[,i], breaks = seq(from = floor(min(min(MPD_per_mimic_null[,i], na.rm = T), MPD_per_mimic_obs[i])), to = ceiling(max(max(MPD_per_mimic_null[,i], na.rm = T), MPD_per_mimic_obs[i])), by = 1), col = "gray", xlab = "Mean Phylogenetic Distance", main = paste0("Distribution of the Mean Phylogenetic distances \n of ",mimic," species \n under the null Hypothesis"))
    hist(MPD_per_mimic_null[,i], breaks = seq(from = floor(min(min(MPD_per_mimic_null[,i], na.rm = T), MPD_per_mimic_obs[i])), to = ceiling(max(max(MPD_per_mimic_null[,i], na.rm = T), MPD_per_mimic_obs[i])), by = 1), col = "gray", xlab = "Mean Phylogenetic Distance", main = paste0("Distribution of the Mean Phylogenetic distances \n of ",mimic," species \n under the null Hypothesis"))
    arrows(MPD_per_mimic_obs[i], max(histo.save$counts)/3, MPD_per_mimic_obs[i], max(histo.save$counts)/30, length = 0.1, lwd = 2)
    abline(v = mean(MPD_per_mimic_null[,i]), lty = 2, lwd = 1)
    abline(v = quantile(MPD_per_mimic_null[,i], 0.025, na.rm = T), lty = 2, col = "red")
    abline(v = quantile(MPD_per_mimic_null[,i], 0.975, na.rm = T), lty = 2, col = "red")
    legend(legend = c(paste0("N units = ", N_units), 
                      paste0("N pairs = ", N_pairs)),
           x = "topright", cex = 1, bty ="n")
    legend(legend = c(paste0("Mean = ", round(mean(MPD_per_mimic_null[,i]),3)), 
                      paste0("CI 2.5% = ", round(quantile(MPD_per_mimic_null[,i], 0.025, na.rm = T),3)),
                      paste0("CI 97.5% = ", round(quantile(MPD_per_mimic_null[,i], 0.975, na.rm = T),3))), 
           x = "topleft", cex = 1, bty ="n")
    legend(legend = c(paste0("MPD obs = ", round(MPD_per_mimic_obs[i], 3)),
                      paste0("p = ", round(ecdf(x = MPD_per_mimic_null[,i])(MPD_per_mimic_obs[i]),3))),
           x = "left", cex = 1, bty ="n")
    dev.off()
  }
  save(MPD_ring_summary_table, file = paste0(internal.wd,"/Niche_evolution/MPD_ring_summary_table.RData"))
  cat(paste0(i, " - ",mimic, " - Done \n"))
}

View(MPD_ring_summary_table)


### Signal phylog?n?tique dans la niche climatique (Kmult)

source(paste0(internal.wd, "./Test.Kmult.R")) # From Adams, 2014
load(file = paste0(internal.wd,"/sp.env.table.RData"))

test.Kmult <- Test.Kmult(x = list.sp[,c("bio1", "bio3", "bio4", "bio12", "bio15")], phy =  phylo.Ithomiini, iter = 999)
save(test.Kmult, file = paste0(internal.wd,"/Niche_evolution/test.Kmult.Rdata"))
test.Kmult$phy.signal ; test.Kmult$pvalue


jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/Kmult_test.jpeg"), quality =100)
hist(test.Kmult$K.simul.val, 30, freq = TRUE, col = "gray", main = "Distribution of Kmult of the climatic niche \n under the null hypothesis", xlab = "Kmult index")
arrows(test.Kmult$phy.signal, 50, test.Kmult$phy.signal, 5, length = 0.1, lwd = 2)
abline(v = mean(test.Kmult$K.simul.val), lty = 2)
abline(v = quantile(test.Kmult$K.simul.val, 0.95), lty = 2, col = "red")
legend(legend = c(paste0("    Mean = ", round(mean(test.Kmult$K.simul.val),3)), 
                  paste0("  CI 95% = ", round(quantile(test.Kmult$K.simul.val, 0.95),3))), 
       x = "topright", cex = 1, bty ="n")
legend(legend = c(paste0("Kmult obs = ", round(test.Kmult$phy.signal, 3)),
                  paste0("       p = ", test.Kmult$pvalue)),
       x = "right", cex = 1, bty ="n")
dev.off()

# Plot de l'histo global
load(file = paste0(internal.wd,"/Niche_evolution/test.Kmult.Rdata"))

tiff(filename = paste0(internal.wd,"/Niche_evolution/Plots/Kmult_test.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
hist(c(test.Kmult$K.simul.val, test.Kmult$phy.signal), 30, freq = TRUE, col = "gray", 
     main = "Distribution of Kmult of the climatic niche \n under the null hypothesis", 
     xlab = "Kmult index",
     cex.axis = 1.7, cex.main = 1, cex.lab = 1.7)
arrows(test.Kmult$phy.signal + 0.0021, 50, test.Kmult$phy.signal+ 0.0021, 10, length = 0.1, lwd = 2)
abline(v = mean(c(test.Kmult$K.simul.val, test.Kmult$phy.signal)), lty = 2)
abline(v = quantile(c(test.Kmult$K.simul.val, test.Kmult$phy.signal), 0.95), lty = 2, col = "red")
legend(legend = c(paste0("Mean = ", round(mean(c(test.Kmult$K.simul.val, test.Kmult$phy.signal)),3)), 
                  paste0("CI 95% = ", round(quantile(c(test.Kmult$K.simul.val, test.Kmult$phy.signal), 0.95),3))), 
       x = "topleft", lty = 2, lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")
legend(legend = c(paste0("Kmult obs = ", round(test.Kmult$phy.signal, 3)),
                  paste0("         p = ", test.Kmult$pvalue)),
       x = "right", cex = 1, bty ="n", xjust = 1)
par(mar = original_int_margins)
dev.off()


##### Simulation of character evolution following evolution models

library(ape)

??corMatrix
?rTraitMult

library(geiger)

?fitContinuous

library(adephylo)

?ppca

## Models fitting and comparison

library(motmot.2.0)

?transformPhylo.ML

load(file = paste0(internal.wd,"/sp.env.table.RData"))

sp.env.table <- list.sp[,c("bio1","bio3","bio4","bio12","bio15")]
row.names(sp.env.table) <- list.sp$Sp_full

# V?rifier que les esp?ces sont bien ordonn?es de la m?me mani?re
list.sp$Sp_full
phylo.Ithomiini$tip.label

hist(sp.env.table[,1])
hist(sp.env.table[,2])
hist(sp.env.table[,3])
hist(sp.env.table[,4])
hist(sp.env.table[,5])

?transformPhylo.ML

par(mfcol=(c(1,1)))

# Create summary table
Models.df <- data.frame(Model = character(), Likelihood = numeric(), AICc = numeric(), param1 = numeric(), param2 = numeric())
Models.df$Model <- as.character(Models.df$Model)
# Mod?le Brownien
bm.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(sp.env.table), model="bm")
bm.ml
Models.df[1,] <-  c("BM", round(bm.ml$logLikelihood, 3), round(bm.ml$AICc, 3), NA, NA)

# Mod?le avec Pagel's lambda
lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(sp.env.table), model="lambda", profilePlot=T)
lambda.ml
Models.df[2,] <-  c("Lambda", round(lambda.ml$MaximumLikelihood, 3), round(lambda.ml$AICc, 3), round(lambda.ml$Lambda[1,1],3), NA)
p.value <- 1 - pchisq(lambda.ml$MaximumLikelihood - bm.ml$logLikelihood, 1) ; p.value
bm.ml$AICc - lambda.ml$AICc

# Pagel's lambda univari? pour voir
var <- 1:2
test <- as.matrix(sp.env.table[,var]) ; row.names(test) <- list.sp$Sp_full 
lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=test, model="lambda", profilePlot=T)
lambda.ml

# Test du signal phylog?n?tique en multivari?e avec lambda de Pagel
tree0 <- transformPhylo(phy = phylo.Ithomiini, model = "lambda", y = as.matrix(sp.env.table), lambda = 0) # Transformation de l'arbre sous Pagel's lambda = 0
lambda0.ml <- transformPhylo.ML(phy=tree0, y=as.matrix(sp.env.table), model="bm") # Estimation du meilleur mod?le en BM sur l'arbre transform?
p.value <- 1 - pchisq(lambda.ml$MaximumLikelihood - lambda0.ml$logLikelihood, 1) ; p.value # LRT via mod?les emboit?s

# Mod?le avec Pagel's Kappa
kappa.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(sp.env.table), model="kappa", profilePlot=T)
p.value <- 1 - pchisq(kappa.ml$MaximumLikelihood - bm.ml$logLikelihood, 1) ; p.value
bm.ml$AICc - kappa.ml$AICc
Models.df[3,] <-  c("Kappa", round(kappa.ml$MaximumLikelihood, 3), round(kappa.ml$AICc, 3), round(kappa.ml$Kappa[1,1],3), NA)


# Mod?le avec Pagel's lambda et Pagel's Kappa
kappa.lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(sp.env.table), model="kappa", lambdaEst = T, profilePlot=T)
kappa.lambda.ml
Models.df[4,] <-  c("Kappa - Lambda", round(kappa.lambda.ml$MaximumLikelihood, 3), round(kappa.lambda.ml$AICc, 3), round(kappa.lambda.ml$Kappa[1,1],3), round(kappa.lambda.ml$lambda,3))
p.value <- 1 - pchisq(kappa.lambda.ml$MaximumLikelihood - bm.ml$logLikelihood, 2) ; p.value
bm.ml$AICc - kappa.lambda.ml$AICc
p.value <- 1 - pchisq(kappa.lambda.ml$MaximumLikelihood - kappa.ml$MaximumLikelihood, 1) ; p.value
kappa.ml$AICc - kappa.lambda.ml$AICc
p.value <- 1 - pchisq(kappa.lambda.ml$MaximumLikelihood - lambda.ml$MaximumLikelihood, 1) ; p.value
lambda.ml$AICc - kappa.lambda.ml$AICc

save(Models.df, file = paste0(internal.wd,"/Niche_evolution/Models.df.RData"))

# Mod?les bonus
delta.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(sp.env.table), model="delta", profilePlot=T)
ou.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(sp.env.table), model="OU", profilePlot=T)
acdc.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(sp.env.table), model="ACDC", profilePlot=T)


## Traits simulation under null hypothesis

# Meilleur mod?le : Lambda = 0 soit pas de signal. Peu cr?dible et peu utile d'utiliser la phylog?nie pour smuler des valeurs dans ce cas.



##### Approche bivari?e apr?s ACP classique #####

load(file = paste0(internal.wd,"/sp.env.table.RData"))

sp.env.table <- list.sp[,c("bio1","bio3","bio4","bio12","bio15")]
row.names(sp.env.table) <- list.sp$Sp_full

# V?rifier que les esp?ces sont bien ordonn?es de la m?me mani?re
list.sp$Sp_full
phylo.Ithomiini$tip.label

# Compute PCA to extract new axis

library(ade4) 

ACP <-  dudi.pca(df = sp.env.table, center = T, scale = T, scannf = F, nf = 2) # G?n?re l'ACP et le graph des ?boulis directement depuis les donn?es brutes
round(ACP$eig/sum(ACP$eig)*100,3) # % Variance expliqu?e
cumsum(round(ACP$eig/sum(ACP$eig)*100,3)) # % Variance expliqu?e cumulative
# 93% de la variance expliqu?e avec les 2 premiers axes, on ne garde que ces deux axes !

s.corcircle(ACP$co,xax=1,yax=2) # Pour tracer le cercle des correlations entre variables d'origine et CP. On peut choisir quelles CP sont repr?sent?es sur les axes x et y

PC.env <- ACP$li # Extraction des nouvelles variables

hist(PC.env[,1])
hist(PC.env[,2])

?transformPhylo.ML

par(mfcol=(c(1,1)))

# Create summary table
Models.PC.df <- data.frame(Model = character(), Likelihood = numeric(), AICc = numeric(), param1 = numeric(), param2 = numeric())
Models.PC.df$Model <- as.character(Models.PC.df$Model)
# Mod?le Brownien
bm.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.env), model="bm")
bm.ml
Models.PC.df[1,] <-  c("BM", round(bm.ml$logLikelihood, 3), round(bm.ml$AICc, 3), NA, NA)

# Mod?le avec Pagel's lambda
lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.env), model="lambda", profilePlot=T)
lambda.ml
Models.PC.df[2,] <-  c("Lambda", round(lambda.ml$MaximumLikelihood, 3), round(lambda.ml$AICc, 3), round(lambda.ml$Lambda[1,1],3), NA)
p.value <- 1 - pchisq(lambda.ml$MaximumLikelihood - bm.ml$logLikelihood, 1) ; p.value
bm.ml$AICc - lambda.ml$AICc

# Pagel's lambda univari? pour voir
var <- 1
test <- as.matrix(PC.env[,var]) ; row.names(test) <- list.sp$Sp_full 
test.lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=test, model="lambda", profilePlot=T)
test.lambda.ml

# Test du signal phylog?n?tique en multivari?e avec lambda de Pagel
tree0 <- transformPhylo(phy = phylo.Ithomiini, model = "lambda", y = as.matrix(PC.env), lambda = 0) # Transformation de l'arbre sous Pagel's lambda = 0
lambda0.ml <- transformPhylo.ML(phy=tree0, y=as.matrix(PC.env), model="bm") # Estimation du meilleur mod?le en BM sur l'arbre transform?
p.value <- 1 - pchisq(lambda.ml$MaximumLikelihood - lambda0.ml$logLikelihood, 1) ; p.value # LRT via mod?les emboit?s

# Mod?le avec Pagel's Kappa
kappa.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.env), model="kappa", profilePlot=T)
p.value <- 1 - pchisq(kappa.ml$MaximumLikelihood - bm.ml$logLikelihood, 1) ; p.value
bm.ml$AICc - kappa.ml$AICc
Models.PC.df[3,] <-  c("Kappa", round(kappa.ml$MaximumLikelihood, 3), round(kappa.ml$AICc, 3), round(kappa.ml$Kappa[1,1],3), NA)

# Mod?le avec Pagel's lambda et Pagel's Kappa
kappa.lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.env), model="kappa", lambdaEst = T, profilePlot=T)
kappa.lambda.ml
Models.PC.df[4,] <-  c("Kappa - Lambda", round(kappa.lambda.ml$MaximumLikelihood, 3), round(kappa.lambda.ml$AICc, 3), round(kappa.lambda.ml$Kappa[1,1],3), round(kappa.lambda.ml$lambda,3))
p.value <- 1 - pchisq(kappa.lambda.ml$MaximumLikelihood - bm.ml$logLikelihood, 2) ; p.value
bm.ml$AICc - kappa.lambda.ml$AICc
p.value <- 1 - pchisq(kappa.lambda.ml$MaximumLikelihood - kappa.ml$MaximumLikelihood, 1) ; p.value
kappa.ml$AICc - kappa.lambda.ml$AICc
p.value <- 1 - pchisq(kappa.lambda.ml$MaximumLikelihood - lambda.ml$MaximumLikelihood, 1) ; p.value
lambda.ml$AICc - kappa.lambda.ml$AICc
# Diff?rence non-significative, on garde le mod?le le plus simple soit le mod?le avec simplement un lambda.

save(Models.PC.df, file = paste0(internal.wd,"/Niche_evolution/Models.PC.df.RData"))

# Mod?les bonus
delta.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.env), model="delta", profilePlot=T)
ou.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.env), model="OU", profilePlot=T)
acdc.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.env), model="ACDC", profilePlot=T)


## Traits simulation under null hypothesis

library(phylocurve)

best.model <- lambda.ml

?sim.traits # Besoin d'une matrice de variance-covariance entre les traits (obtenue depuis motmot)

Sim.res <- sim.traits(tree = phylo.Ithomiini, v = best.model$brownianVariance, anc = best.model$root.state, model = "lambda", parameters = list(lambda = best.model$Lambda[1,1]), nsim = 999, return.type = "matrix")

str(Sim.res)
plot(Sim.res$tree) # Arbre original
plot(Sim.res$sim_tree) # Arbre apr?s tranformation via lambda = 0.47

Sim.PC <- Sim.res$trait_data # Liste des simulations

## Stat computation and distribution plot

# Compute the climatic distance matrix
Pairwise_climdist <- dist(x = PC.env, method = "euclidian")

# Load the Jaccard index similarities to use as weights
load(file = paste0(internal.wd,"/Niche_evolution/Jaccard_index.RData"))

# Compute the observed mean climatic distance (MCD) weighted by mimicry similarity
MCD_obs <- weighted.mean(x = Pairwise_climdist, w = mimicry.similarity.weights) # 2.175
# Standardize the MCD on full species mean pairwise climatic distance
MCD_fullsp <- mean(Pairwise_climdist) # 2.702
MCD_obs <- MCD_obs/MCD_fullsp # 0.805

# Compute the stat for all the simulation
MCD_null <-  NA
for (i in 1:length(Sim.PC)) {
  Pairwise_climdist_sim <- dist(x = Sim.PC[[i]], method = "euclidian")
  MCD_sim <- weighted.mean(x = Pairwise_climdist_sim, w = mimicry.similarity.weights)/mean(Pairwise_climdist_sim)
  MCD_null[i] <- MCD_sim
}
save(MCD_obs, MCD_null, file = paste0(internal.wd,"/Niche_evolution/MCD_values.RData"))

summary(MCD_null)

# Plot the distri of the stats

jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/MCD.classic.PCA_test.jpeg"), quality =100)
hist(c(MCD_null, MCD_obs), 30, freq = TRUE, col = "gray", main = "Distribution of the Mean Climatic Distance \n of co-mimetic species \n under the null hypothesis", xlab = "Standardized Mean Climatic Distance")
arrows(MCD_obs, 50, MCD_obs, 5, length = 0.1, lwd = 2)
abline(v = mean(c(MCD_null, MCD_obs)), lty = 2)
abline(v = quantile(c(MCD_null, MCD_obs), 0.05), lty = 2, col = "red")
legend(legend = c(paste0("Mean = ", round(mean(c(MCD_null, MCD_obs)),3)), 
                  paste0("CI 5% = ", round(quantile(c(MCD_null, MCD_obs), 0.05),3))), 
       x = "topleft", cex = 1, bty ="n")
legend(legend = c(paste0("MCD obs = ", round(MCD_obs, 3)),
                  paste0("p = ", ecdf(x = c(MCD_null, MCD_obs))(MCD_obs))),
       x = "left", cex = 1, bty ="n")
dev.off()





##### Approche bivari?e apr?s pPCA (car Revell dit que l'on ne peut pas faire de PCA sur des donn?es phylo. sans prendre en compte la phylog?nie) #####

load(file = paste0(internal.wd,"/sp.env.table.RData"))

sp.env.table <- list.sp[,c("bio1","bio3","bio4","bio12","bio15")]
row.names(sp.env.table) <- list.sp$Sp_full

# V?rifier que les esp?ces sont bien ordonn?es de la m?me mani?re
list.sp$Sp_full
phylo.Ithomiini$tip.label

# Compute PCA to extract new axis

source(paste0(internal.wd,"/Niche_evolution/Revell_phyl_pca.R"))
library(ape)

C <- vcv(phylo.Ithomiini)
X <- as.matrix(sp.env.table)

PCA.Revell <- Revell_phyl_pca(C, X, mode = "corr")
diag(PCA.Revell$Eval) # eigenvalues
PCA.Revell$Evec # eigenvectors
PCA.Revell$S # scores
PCA.Revell$L # PC loadings

round(diag(PCA.Revell$Eval)/sum(PCA.Revell$Eval),3) # % Variance expliqu?e
round(cumsum(diag(PCA.Revell$Eval))/sum(PCA.Revell$Eval),3) # % Variance expliqu?e cumulative
# 91% de la variance expliqu?e avec les 2 premiers axes, on ne garde que ces deux axes !

library(ade4)
s.corcircle(PCA.Revell$L, xax=1, yax=2, label = names(sp.env.table)) # Pour tracer le cercle des correlations entre variables d'origine et CP. On peut choisir quelles CP sont repr?sent?es sur les axes x et y

PC.Revell.env <- PCA.Revell$S[,1:2]
save(PC.Revell.env, file = paste0(internal.wd, "/Niche_evolution/PC.Revell.env.RData"))

load(file = paste0(internal.wd, "/Niche_evolution/PC.Revell.env.RData"))

hist(PC.Revell.env[,1])
hist(PC.Revell.env[,2])

?transformPhylo.ML

par(mfcol=(c(1,1)))

# Create summary table
Models.PC.Revell.df <- data.frame(Model = character(), Likelihood = numeric(), AICc = numeric(), param1 = numeric(), param2 = numeric())
Models.PC.Revell.df$Model <- as.character(Models.PC.Revell.df$Model)
# Mod?le Brownien
Revell.bm.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="bm")
Revell.bm.ml
Models.PC.Revell.df[1,] <-  c("BM", round(Revell.bm.ml$logLikelihood, 3), round(Revell.bm.ml$AICc, 3), NA, NA)

# Mod?le avec Pagel's lambda
Revell.lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="lambda", profilePlot=T)
Revell.lambda.ml
Models.PC.Revell.df[2,] <-  c("Lambda", round(Revell.lambda.ml$MaximumLikelihood, 3), round(Revell.lambda.ml$AICc, 3), round(Revell.lambda.ml$Lambda[1,1],3), NA)
p.value <- 1 - pchisq(Revell.lambda.ml$MaximumLikelihood - Revell.bm.ml$logLikelihood, 1) ; p.value
Revell.bm.ml$AICc - Revell.lambda.ml$AICc

# Pagel's lambda univari? pour voir
var <- 2
test <- as.matrix(PC.Revell.env[,var]) ; row.names(test) <- list.sp$Sp_full 
test.lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=test, model="lambda", profilePlot=T)
test.lambda.ml

# Test du signal phylog?n?tique en multivari?e avec lambda de Pagel
tree0 <- transformPhylo(phy = phylo.Ithomiini, model = "lambda", y = as.matrix(PC.Revell.env), lambda = 0) # Transformation de l'arbre sous Pagel's lambda = 0
Revell.lambda0.ml <- transformPhylo.ML(phy=tree0, y=as.matrix(PC.Revell.env), model="bm") # Estimation du meilleur mod?le en BM sur l'arbre transform?
p.value <- 1 - pchisq(Revell.lambda.ml$MaximumLikelihood - Revell.lambda0.ml$logLikelihood, 1) ; p.value # LRT via mod?les emboit?s

# Mod?le avec Pagel's Kappa
Revell.kappa.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="kappa", profilePlot=T)
p.value <- 1 - pchisq(Revell.kappa.ml$MaximumLikelihood - Revell.bm.ml$logLikelihood, 1) ; p.value
Revell.bm.ml$AICc - Revell.kappa.ml$AICc
Models.PC.Revell.df[3,] <-  c("Kappa", round(Revell.kappa.ml$MaximumLikelihood, 3), round(Revell.kappa.ml$AICc, 3), round(Revell.kappa.ml$Kappa[1,1],3), NA)

# Mod?le avec Pagel's lambda et Pagel's Kappa
Revell.kappa.lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="kappa", lambdaEst = T, profilePlot=T)
Revell.kappa.lambda.ml
Models.PC.Revell.df[4,] <-  c("Kappa - Lambda", round(Revell.kappa.lambda.ml$MaximumLikelihood, 3), round(Revell.kappa.lambda.ml$AICc, 3), round(Revell.kappa.lambda.ml$Kappa[1,1],3), round(Revell.kappa.lambda.ml$lambda,3))
p.value <- 1 - pchisq(Revell.kappa.lambda.ml$MaximumLikelihood - Revell.bm.ml$logLikelihood, 2) ; p.value
Revell.bm.ml$AICc - Revell.kappa.lambda.ml$AICc
p.value <- 1 - pchisq(Revell.kappa.lambda.ml$MaximumLikelihood - Revell.kappa.ml$MaximumLikelihood, 1) ; p.value
Revell.kappa.ml$AICc - Revell.kappa.lambda.ml$AICc
p.value <- 1 - pchisq(Revell.kappa.lambda.ml$MaximumLikelihood - Revell.lambda.ml$MaximumLikelihood, 1) ; p.value
Revell.lambda.ml$AICc - Revell.kappa.lambda.ml$AICc
# Diff?rence non-significative, on garde le mod?le le plus simple soit le mod?le avec simplement un lambda.

save(Models.PC.Revell.df, file = paste0(internal.wd,"/Niche_evolution/Models.PC.Revell.df.RData"))
load(file = paste0(internal.wd,"/Niche_evolution/Models.PC.Revell.df.RData"))

# Mod?les bonus
delta.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="delta", profilePlot=T)
ou.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="OU", profilePlot=T)
acdc.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="ACDC", profilePlot=T)


## Traits simulation under null hypothesis

library(phylocurve)

best.model <- Revell.lambda.ml
save(best.model, file = paste0(internal.wd, "/Niche_evolution/best.model.for.pPC_evol"))

load(file = paste0(internal.wd, "/Niche_evolution/best.model.for.pPC_evol"))

?sim.traits # Besoin d'une matrice de variance-covariance entre les traits (obtenue depuis motmot)

Sim.Revell.res <- sim.traits(tree = phylo.Ithomiini, v = best.model$brownianVariance, anc = best.model$root.state, model = "lambda", parameters = list(lambda = best.model$Lambda[1,1]), nsim = 999, return.type = "matrix")

str(Sim.Revell.res)
plot(Sim.Revell.res$tree) # Arbre original
plot(Sim.Revell.res$sim_tree) # Arbre apr?s tranformation via lambda = 0.47

Sim.Revell.PC <- Sim.Revell.res$trait_data # Liste des simulations

## Stat computation and distribution plot

# Compute the climatic distance matrix
Revell_Pairwise_climdist <- dist(x = PC.Revell.env, method = "euclidian")

# Load the Jaccard index similarities to use as weights
load(file = paste0(internal.wd,"/Niche_evolution/Jaccard_index.RData"))

# Compute the observed mean climatic distance (MCD) weighted by mimicry similarity
MCD.Revell_obs <- weighted.mean(x = Revell_Pairwise_climdist, w = mimicry.similarity.weights) # 3.544
# Standardize the MCD on full species mean pairwise climatic distance
MCD.Revell_fullsp <- mean(Revell_Pairwise_climdist) # 4.374
MCD.Revell_obs <- MCD.Revell_obs/MCD.Revell_fullsp # 0.810

# Compute the stat for all the simulation
MCD.Revell_null <-  NA
for (i in 1:length(Sim.Revell.PC)) {
  Revell_Pairwise_climdist_sim <- dist(x = Sim.Revell.PC[[i]], method = "euclidian")
  MCD.Revell_sim <- weighted.mean(x = Revell_Pairwise_climdist_sim, w = mimicry.similarity.weights)/mean(Revell_Pairwise_climdist_sim)
  MCD.Revell_null[i] <- MCD.Revell_sim
}
save(MCD.Revell_obs, MCD.Revell_null, file = paste0(internal.wd,"/Niche_evolution/MCD.Revell_values.RData"))

summary(MCD.Revell_null)

# Plot de la distri de la stat
load(file = paste0(internal.wd,"/Niche_evolution/MCD.Revell_values.RData"))

tiff(filename = paste0(internal.wd,"/Niche_evolution/Plots/MCD.pPCA_test.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
hist(c(MCD.Revell_null, MCD.Revell_obs), 30, freq = TRUE, col = "gray", 
     main = "Distribution of the Mean Climatic Distance \n of co-mimetic species \n under the null hypothesis", 
     xlab = "Standardized Mean Climatic Distance",
     cex.axis = 1.7, cex.main = 1, cex.lab = 1.7)
arrows(MCD.Revell_obs + 0.0051, 50, MCD.Revell_obs+ 0.0051, 10, length = 0.1, lwd = 2)
abline(v = mean(c(MCD.Revell_null, MCD.Revell_obs)), lty = 2, lwd = 2)
abline(v = quantile(c(MCD.Revell_null, MCD.Revell_obs), 0.05), lty = 2, lwd = 2, col = "red")
legend(legend = c(paste0("Mean = ", round(mean(c(MCD.Revell_null, MCD.Revell_obs)),3)), 
                  paste0("CI 5% = ", round(quantile(c(MCD.Revell_null, MCD.Revell_obs), 0.05),3))), 
       x = "topleft", lty = 2, lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")
legend(legend = c(paste0("\n\n\nMCD obs = ", round(MCD.Revell_obs, 3)),
                  paste0("p = ", test.Kmult$pvalue)),
       x = "left", cex = 1.2, bty ="n", xjust = 1)
par(mar = original_int_margins)
dev.off()


### Stats computation and plot per mimicry ring !

# Generate list of mimicry rings 
mimic.list <- as.character(unique(list.models$Mimicry.model))

# Generate incidence table of species and mimicry rings to be sure the species are in the good order in the incidence table
mimic.incidence.table <- matrix(nrow = 339, ncol = 44)
mimic.list <- as.character(unique(list.models$Mimicry.model))
for (i in 1:nrow(list.sp)) { # 339 sp = rows
  sp <- as.character(list.sp$Sp_full[i])
  for (j in 1:length(mimic.list)) { # 44 mimicry rings = columns
    mimic <- mimic.list[j]
    temp <- as.character(list.models$Mimicry.model[list.models$Sp_full==sp]) # Extract mimicry ring of all units for the i species
    mimic.incidence.table[i,j] <- as.numeric(any(temp==mimic)) # 0 if mimicry ring not present ; 1 if found
  }
}

# Loop 
MCD_per_mimic_obs <- rep(NA, length(mimic.list))
MCD_per_mimic_null <- NA
for (i in 1:length(mimic.list)) {
  mimic <- mimic.list[i]
  # Compute the Jaccard index similarities to use as weights, specific to the mimicry ring
  # Weight = proportion of the targeted mimicry ring among all the units of the pair of species
  mimic.weights <- matrix(nrow = nrow(mimic.incidence.table), ncol = nrow(mimic.incidence.table))
  for (j in 1:nrow(mimic.incidence.table)) {
    for (k in 1:nrow(mimic.incidence.table)) {
      if (mimic.incidence.table[j,i]+mimic.incidence.table[k,i]==2) { # Case when the two species share the targeted mimicry ring
        sp1 <- which(mimic.incidence.table[j,]==1) # Index of mimicry ring of sp1
        sp2 <- which(mimic.incidence.table[k,]==1) # Index of mimicry ring of sp2
        tot_rings <- length(union(sp1,sp2)) # Number of different mimicry ring in total among the two species
        mimic.weights[j,k] <- 1/tot_rings # Weight is the proportion of the targeted mimicry ring among all the units of the two species
      }else{ # Case when the two species don't share the targeted mimicry ring
        mimic.weights[j,k] <- 0 # Weight is null
      }
    }
  }
  
 
  # Compute the climatic distance matrix on observed data
  Revell_Pairwise_climdist <- dist(x = PC.Revell.env, method = "euclidian")
  # Compute the observed mean climatic distance (MCD) weighted by mimicry similarity
  mimic.MCD.Revell_obs <- weighted.mean(x = Revell_Pairwise_climdist, w = as.dist(mimic.weights))
  # Standardize the MCD on full species mean pairwise climatic distance
  MCD.Revell_fullsp <- mean(Revell_Pairwise_climdist) 
  mimic.MCD.Revell_obs <- mimic.MCD.Revell_obs/MCD.Revell_fullsp
  MCD_per_mimic_obs[i] <- mimic.MCD.Revell_obs # Store MCD_obs for this mimicry ring
 
  ### Compute the null distribution using the simulated traits
  # Compute the stat for all the simulation
  mimic.MCD.Revell_null <-  NA
  for (j in 1:length(Sim.Revell.PC)) {
    Revell_Pairwise_climdist_sim <- dist(x = Sim.Revell.PC[[j]], method = "euclidian")
    mimic.MCD.Revell_sim <- weighted.mean(x = Revell_Pairwise_climdist_sim, w = as.dist(mimic.weights))/mean(Revell_Pairwise_climdist_sim)
    mimic.MCD.Revell_null[j] <- mimic.MCD.Revell_sim
  }
  
  if (i == 1) { # Store the 1st mimic MCD null distri in a new table
    MCD_per_mimic_null <- mimic.MCD.Revell_null
  }else{ # Add the next mimic MCD null distri to the table
    MCD_per_mimic_null <- cbind(MCD_per_mimic_null,mimic.MCD.Revell_null)
  }
  
  save(MCD_per_mimic_obs, MCD_per_mimic_null, file = paste0(internal.wd,"/Niche_evolution/MCD.Revell_per_mimic_values.RData"))
  print(i)
}
names(MCD_per_mimic_obs) <- names(MCD_per_mimic_null) <- mimic.list
save(MCD_per_mimic_obs, MCD_per_mimic_null, file = paste0(internal.wd,"/Niche_evolution/MCD.Revell_per_mimic_values.RData"))


MCD_ring_summary_table <- as.data.frame(matrix(ncol = 8, nrow = 44, data = NA))
names(MCD_ring_summary_table) <- c("ring", "N_units", "N_pairs", "MCD_obs", "mean_MCD", "MCD_2.5", "MCD_97.5", "p_value")

### Plot

load(file = paste0(internal.wd,"/Niche_evolution/MCD.Revell_per_mimic_values.RData"))

for (i in 1:length(mimic.list)) {
  MCD_ring_summary_table$ring[i] <- mimic <- mimic.list[i]
  MCD_ring_summary_table$N_units[i] <- N_units <- sum(list.models$Mimicry.model==mimic)
  MCD_ring_summary_table$N_pairs[i] <- N_pairs <- N_units*(N_units-1)/2
  
  if (is.nan(MCD_per_mimic_obs[i])) {
    jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/MCD_per_mimic/MCD.pPCA_test_",mimic,".jpeg"), quality = 100)
    plot(1:100,1:100, type = "n", main = paste0("Distribution of Mean Climatic Distance \n of ", mimic, " species \n under null Hypothesis"))
    text(x = 50, y = 50, labels = "Only one unit for this mimicry ring \n No pair available for index computation")
    dev.off()
  }else{
    MCD_ring_summary_table$MCD_obs[i] <- round(MCD_per_mimic_obs[i],3)
    MCD_ring_summary_table$mean_MCD[i] <- round(mean(MCD_per_mimic_null[,i], na.rm = T),3)
    MCD_ring_summary_table$MCD_2.5[i] <- round(quantile(MCD_per_mimic_null[,i], 0.025),3)
    MCD_ring_summary_table$MCD_97.5[i] <- round(quantile(MCD_per_mimic_null[,i], 0.975),3)
    MCD_ring_summary_table$p_value[i] <- round(ecdf(x = MCD_per_mimic_null[,i])(MCD_per_mimic_obs[i]),3)
    jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/MCD_per_mimic/MCD.pPCA_test_,",mimic,".jpeg"), quality =100)
    hist(c(MCD_per_mimic_null[,i], MCD_per_mimic_obs[i]), 30, freq = TRUE, col = "gray", main = paste0("Distribution of the Mean Climatic Distance \n of ",mimic, " co-mimetic species \n under the null hypothesis"), xlab = "Standardized Mean Climatic Distance")
    arrows(MCD_per_mimic_obs[i]+0.004, 50, MCD_per_mimic_obs[i]+0.004, 5, length = 0.1, lwd = 2)
    abline(v = mean(c(MCD_per_mimic_null[,i], MCD_per_mimic_obs[i])), lty = 2)
    abline(v = quantile(c(MCD_per_mimic_null[,i], MCD_per_mimic_obs[i]), 0.025), lty = 2, col = "red")
    abline(v = quantile(c(MCD_per_mimic_null[,i], MCD_per_mimic_obs[i]), 0.975), lty = 2, col = "red")
    legend(legend = c(paste0("N units = ", N_units), 
                      paste0("N pairs = ", N_pairs)),
           x = "topright", cex = 1, bty ="n")
    legend(legend = c(paste0("Mean = ", round(mean(c(MCD_per_mimic_null[,i], MCD_per_mimic_obs[i])),3)),
                      paste0("CI 2.5% = ", round(quantile(c(MCD_per_mimic_null[,i], MCD_per_mimic_obs[i]), 0.025),3)),
                      paste0("CI 97.5% = ", round(quantile(c(MCD_per_mimic_null[,i], MCD_per_mimic_obs[i]), 0.975),3))), 
           x = "topleft", cex = 1, bty ="n")
    legend(legend = c(paste0("MCD obs = ", round(MCD_per_mimic_obs[i], 3)),
                      paste0("p = ", ecdf(x = c(MCD_per_mimic_null[,i], MCD_per_mimic_obs[i]))(MCD_per_mimic_obs[i]))),
           x = "left", cex = 1, bty ="n")
    dev.off()
    save(MCD_ring_summary_table, file = paste0(internal.wd,"/Niche_evolution/MCD_ring_summary_table.Rdata"))
  }
  cat(paste0(i, " - ",mimic, " - Done \n"))
}

View(MCD_ring_summary_table)






##### Approche univari? apr?s Revell's pPCA #####

load(file = paste0(internal.wd,"/sp.env.table.RData"))

sp.env.table <- list.sp[,c("bio1","bio3","bio4","bio12","bio15")]
row.names(sp.env.table) <- list.sp$Sp_full

# V?rifier que les esp?ces sont bien ordonn?es de la m?me mani?re
list.sp$Sp_full
phylo.Ithomiini$tip.label

# Compute PCA to extract new axis

source(paste0(internal.wd,"/Niche_evolution/Revell_phyl_pca.R"))
library(ape)

C <- vcv(phylo.Ithomiini)
X <- as.matrix(sp.env.table)

PCA.Revell <- Revell_phyl_pca(C, X, mode = "corr")
diag(PCA.Revell$Eval) # eigenvalues
PCA.Revell$Evec # eigenvectors
PCA.Revell$S # scores
PCA.Revell$L # PC loadings

round(diag(PCA.Revell$Eval)/sum(PCA.Revell$Eval),3) # % Variance expliqu?e
round(cumsum(diag(PCA.Revell$Eval))/sum(PCA.Revell$Eval),3) # % Variance expliqu?e cumulative
# 91% de la variance expliqu?e avec les 2 premiers axes, on ne garde que ces deux axes !

library(ade4)
s.corcircle(PCA.Revell$L, xax=1, yax=2, label = names(sp.env.table)) # Pour tracer le cercle des correlations entre variables d'origine et CP. On peut choisir quelles CP sont repr?sent?es sur les axes x et y

PC.univar.env <- PCA.Revell$S[,1:2]
PC1.univar.env <- PCA.Revell$S[,1]
PC2.univar.env <- PCA.Revell$S[,2]

hist(PC1.univar.env)
hist(PC2.univar.env)

?transformPhylo.ML

par(mfcol=(c(1,1)))

# Create summary table
Models.PC.univar.df <- data.frame(PC = character(), Model = character(), Likelihood = numeric(), AICc = numeric(), param1 = numeric(), param2 = numeric())
Models.PC.univar.df$PC <- as.character(Models.PC.univar.df$PC)
Models.PC.univar.df$Model <- as.character(Models.PC.univar.df$Model)

## Axe 1

# Mod?le Brownien
PC1.univar.bm.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC1.univar.env), model="bm")
PC1.univar.bm.ml
Models.PC.univar.df[1,] <-  c("PC1", "BM", round(PC1.univar.bm.ml$logLikelihood, 3), round(PC1.univar.bm.ml$AICc, 3), NA, NA)

# Mod?le avec Pagel's lambda
PC1.univar.lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC1.univar.env), model="lambda", profilePlot=T)
PC1.univar.lambda.ml
Models.PC.univar.df[2,] <-  c("PC1", "Lambda", round(PC1.univar.lambda.ml$MaximumLikelihood, 3), round(PC1.univar.lambda.ml$AICc, 3), round(PC1.univar.lambda.ml$Lambda[1,1],3), NA)
p.value <- 1 - pchisq(PC1.univar.lambda.ml$MaximumLikelihood - PC1.univar.bm.ml$logLikelihood, 1) ; p.value
PC1.univar.bm.ml$AICc - PC1.univar.lambda.ml$AICc

# Test du signal phylog?n?tique avec lambda de Pagel
tree0 <- transformPhylo(phy = phylo.Ithomiini, model = "lambda", y = as.matrix(PC1.univar.env), lambda = 0) # Transformation de l'arbre sous Pagel's lambda = 0
lambda0.ml <- transformPhylo.ML(phy=tree0, y=as.matrix(PC1.univar.env), model="bm") # Estimation du meilleur mod?le en BM sur l'arbre transform?
p.value <- 1 - pchisq(PC1.univar.lambda.ml$MaximumLikelihood - lambda0.ml$logLikelihood, 1) ; p.value # LRT via mod?les emboit?s
library(phytools)
lambda.test <- phylosig(tree = phylo.Ithomiini, x = PC1.univar.env, method = "lambda", test = T)
lambda.test

# Mod?le avec Pagel's Kappa
PC1.univar.kappa.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC1.univar.env), model="kappa", profilePlot=T)
p.value <- 1 - pchisq(PC1.univar.kappa.ml$MaximumLikelihood - PC1.univar.bm.ml$logLikelihood, 1) ; p.value
PC1.univar.bm.ml$AICc - PC1.univar.kappa.ml$AICc
Models.PC.univar.df[3,] <-  c("PC1", "Kappa", round(PC1.univar.kappa.ml$MaximumLikelihood, 3), round(PC1.univar.kappa.ml$AICc, 3), round(PC1.univar.kappa.ml$Kappa[1,1],3), NA)

# Mod?le avec Pagel's lambda et Pagel's Kappa
PC1.univar.kappa.lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC1.univar.env), model="kappa", lambdaEst = T, profilePlot=T)
PC1.univar.kappa.lambda.ml
Models.PC.univar.df[4,] <-  c("PC1", "Kappa - Lambda", round(PC1.univar.kappa.lambda.ml$MaximumLikelihood, 3), round(PC1.univar.kappa.lambda.ml$AICc, 3), round(PC1.univar.kappa.lambda.ml$Kappa[1,1],3), round(PC1.univar.kappa.lambda.ml$lambda,3))
p.value <- 1 - pchisq(PC1.univar.kappa.lambda.ml$MaximumLikelihood - PC1.univar.bm.ml$logLikelihood, 2) ; p.value
PC1.univar.bm.ml$AICc - PC1.univar.kappa.lambda.ml$AICc
p.value <- 1 - pchisq(PC1.univar.kappa.lambda.ml$MaximumLikelihood - PC1.univar.kappa.ml$MaximumLikelihood, 1) ; p.value
PC1.univar.kappa.ml$AICc - PC1.univar.kappa.lambda.ml$AICc
p.value <- 1 - pchisq(PC1.univar.kappa.lambda.ml$MaximumLikelihood - PC1.univar.lambda.ml$MaximumLikelihood, 1) ; p.value
PC1.univar.lambda.ml$AICc - PC1.univar.kappa.lambda.ml$AICc
# Diff?rence non-significative, on garde le mod?le le plus simple soit le mod?le avec simplement un lambda.

save(Models.PC.univar.df, file = paste0(internal.wd,"/Niche_evolution/Models.PC.univar.df.RData"))

# Mod?les bonus
delta.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC1.univar.env), model="delta", profilePlot=T)
ou.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC1.univar.env), model="OU", profilePlot=T)
acdc.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC1.univar.env), model="ACDC", profilePlot=T)

## Axe 2

# Mod?le Brownien
PC2.univar.bm.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC2.univar.env), model="bm")
PC2.univar.bm.ml
Models.PC.univar.df[5,] <-  c("PC2", "BM", round(PC2.univar.bm.ml$logLikelihood, 3), round(PC2.univar.bm.ml$AICc, 3), NA, NA)

# Mod?le avec Pagel's lambda
PC2.univar.lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC2.univar.env), model="lambda", profilePlot=T)
PC2.univar.lambda.ml
Models.PC.univar.df[6,] <-  c("PC2", "Lambda", round(PC2.univar.lambda.ml$MaximumLikelihood, 3), round(PC2.univar.lambda.ml$AICc, 3), round(PC2.univar.lambda.ml$Lambda[1,1],3), NA)
p.value <- 1 - pchisq(PC2.univar.lambda.ml$MaximumLikelihood - PC2.univar.bm.ml$logLikelihood, 1) ; p.value
PC2.univar.bm.ml$AICc - PC2.univar.lambda.ml$AICc

# Test du signal phylog?n?tique avec lambda de Pagel
tree0 <- transformPhylo(phy = phylo.Ithomiini, model = "lambda", y = as.matrix(PC2.univar.env), lambda = 0) # Transformation de l'arbre sous Pagel's lambda = 0
lambda0.ml <- transformPhylo.ML(phy=tree0, y=as.matrix(PC2.univar.env), model="bm") # Estimation du meilleur mod?le en BM sur l'arbre transform?
p.value <- 1 - pchisq(PC2.univar.lambda.ml$MaximumLikelihood - lambda0.ml$logLikelihood, 1) ; p.value # LRT via mod?les emboit?s
library(phytools)
lambda.test <- phylosig(tree = phylo.Ithomiini, x = PC2.univar.env, method = "lambda", test = T)
lambda.test

# Mod?le avec Pagel's Kappa
PC2.univar.kappa.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC2.univar.env), model="kappa", profilePlot=T)
PC2.univar.kappa.ml
p.value <- 1 - pchisq(PC2.univar.kappa.ml$MaximumLikelihood - PC2.univar.bm.ml$logLikelihood, 1) ; p.value
PC2.univar.bm.ml$AICc - PC2.univar.kappa.ml$AICc
Models.PC.univar.df[7,] <-  c("PC2", "Kappa", round(PC2.univar.kappa.ml$MaximumLikelihood, 3), round(PC2.univar.kappa.ml$AICc, 3), round(PC2.univar.kappa.ml$Kappa[1,1],3), NA)

# Mod?le avec Pagel's lambda et Pagel's Kappa
PC2.univar.kappa.lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC2.univar.env), model="kappa", lambdaEst = T, profilePlot=T)
PC2.univar.kappa.lambda.ml
Models.PC.univar.df[8,] <-  c("PC2", "Kappa - Lambda", round(PC2.univar.kappa.lambda.ml$MaximumLikelihood, 3), round(PC2.univar.kappa.lambda.ml$AICc, 3), round(PC2.univar.kappa.lambda.ml$Kappa[1,1],3), round(PC2.univar.kappa.lambda.ml$lambda,3))
p.value <- 1 - pchisq(PC2.univar.kappa.lambda.ml$MaximumLikelihood - PC2.univar.bm.ml$logLikelihood, 2) ; p.value
PC2.univar.bm.ml$AICc - PC2.univar.kappa.lambda.ml$AICc
p.value <- 1 - pchisq(PC2.univar.kappa.lambda.ml$MaximumLikelihood - PC2.univar.kappa.ml$MaximumLikelihood, 1) ; p.value
PC2.univar.kappa.ml$AICc - PC2.univar.kappa.lambda.ml$AICc
p.value <- 1 - pchisq(PC2.univar.kappa.lambda.ml$MaximumLikelihood - PC2.univar.lambda.ml$MaximumLikelihood, 1) ; p.value
PC2.univar.lambda.ml$AICc - PC2.univar.kappa.lambda.ml$AICc
# Diff?rence significative, on garde le mod?le le plus complexe soit un kappa-lambda... mais c'est un peu pourri pour un lambda = 0

save(Models.PC.univar.df, file = paste0(internal.wd,"/Niche_evolution/Models.PC.univar.df.RData"))

# Mod?les bonus
delta.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC2.univar.env), model="delta", profilePlot=T)
ou.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC2.univar.env), model="OU", profilePlot=T)
acdc.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC2.univar.env), model="ACDC", profilePlot=T)


## Traits simulation under null hypothesis

library(phylocurve)

# Axe 1

best.model.PC1 <- PC1.univar.lambda.ml

Sim.PC1.univar.res <- sim.traits(tree = phylo.Ithomiini, v = best.model.PC1$brownianVariance, anc = best.model.PC1$root.state, model = "lambda", parameters = list(lambda = best.model.PC1$Lambda[1,1]), nsim = 999, return.type = "matrix")

str(Sim.univar.res)
plot(Sim.PC1.univar.res$tree) # Arbre original
plot(Sim.PC1.univar.res$sim_tree) # Arbre apr?s tranformation via lambda = 0.21

Sim.PC1.univar.data <- Sim.PC1.univar.res$trait_data # Liste des simulations

# Axe 2

best.model.PC2 <- PC2.univar.kappa.lambda.ml

# La fonction ne permet pas d'utiliser un mod?le kap^pa-lambda donc on transforme l'arbre avec un lambda avant d'appliquer le Kappa
tree0 <- transformPhylo(phy = phylo.Ithomiini, model = "lambda", y = as.matrix(PC2.univar.env), lambda = 0)
Sim.PC2.univar.res <- sim.traits(tree = tree0, v = best.model.PC2$brownianVariance, anc = best.model.PC2$root.state, model = "kappa", parameters = list(kappa = best.model.PC2$Kappa[1,1]), nsim = 999, return.type = "matrix")

str(Sim.univar.res)
plot(phylo.Ithomiini) # Arbre original
plot(Sim.PC2.univar.res$tree)
plot(Sim.PC2.univar.res$sim_tree) # Arbre apr?s transformation via lambda = 0 et kappa = 0.072

Sim.PC2.univar.data <- Sim.PC2.univar.res$trait_data # Liste des simulations

# Merging simulated PC values
Sim.univar.PC <- list(rep(NA, length(Sim.PC1.univar.data)))
for (i in 1:length(Sim.PC1.univar.data)) {
  Sim.univar.PC[[i]] <- cbind(Sim.PC1.univar.data[[i]],Sim.PC2.univar.data[[i]])
}

## Stat computation and distribution plot

# Compute the climatic distance matrix
univar_Pairwise_climdist <- dist(x = PC.univar.env, method = "euclidian")

# Load the Jaccard index similarities to use as weights
load(file = paste0(internal.wd,"/Niche_evolution/Jaccard_index.RData"))

# Compute the observed mean climatic distance (MCD) weighted by mimicry similarity
MCD.univar_obs <- weighted.mean(x = univar_Pairwise_climdist, w = mimicry.similarity.weights) # 3.544
# Standardize the MCD on full species mean pairwise climatic distance
MCD.univar_fullsp <- mean(univar_Pairwise_climdist) # 4.374
MCD.univar_obs <- MCD.univar_obs/MCD.univar_fullsp # 0.810

# Compute the stat for all the simulation
MCD.univar_null <-  NA
for (i in 1:length(Sim.univar.PC)) {
  univar_Pairwise_climdist_sim <- dist(x = Sim.univar.PC[[i]], method = "euclidian")
  MCD.univar_sim <- weighted.mean(x = univar_Pairwise_climdist_sim, w = mimicry.similarity.weights)/mean(univar_Pairwise_climdist_sim)
  MCD.univar_null[i] <- MCD.univar_sim
}
save(MCD.univar_obs, MCD.univar_null, file = paste0(internal.wd,"/Niche_evolution/MCD.univar_values.RData"))

summary(MCD.univar_null)

# Plot the distri of the stats

jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/MCD.univar_test.jpeg"), quality =100)
hist(c(MCD.univar_null, MCD.univar_obs), 30, freq = TRUE, col = "gray", main = "Distribution of the Mean Climatic Distance \n of co-mimetic species \n under the null hypothesis", xlab = "Standardized Mean Climatic Distance")
arrows(MCD.univar_obs+0.004, 50, MCD.univar_obs+0.004, 5, length = 0.1, lwd = 2)
abline(v = mean(c(MCD.univar_null, MCD.univar_obs)), lty = 2)
abline(v = quantile(c(MCD.univar_null, MCD.univar_obs), 0.05), lty = 2, col = "red")
legend(legend = c(paste0("Mean = ", round(mean(c(MCD.univar_null, MCD.univar_obs)),3)), 
                  paste0("CI 5% = ", round(quantile(c(MCD.univar_null, MCD.univar_obs), 0.05),3))), 
       x = "topleft", cex = 1, bty ="n")
legend(legend = c(paste0("MCD obs = ", round(MCD.univar_obs, 3)),
                  paste0("p = ", ecdf(x = c(MCD.univar_null, MCD.univar_obs))(MCD.univar_obs))),
       x = "left", cex = 1, bty ="n")
dev.off()

##### Approche OU multi-optima avec PhylogeneticEM #####

library(PhylogeneticEM)

### Version ? 5 variables

# Load data
sp.env.table <- list.sp[,c("bio1","bio3","bio4","bio12","bio15")]
row.names(sp.env.table) <- list.sp$Sp_full

# V?rifier que les esp?ces sont bien ordonn?es de la m?me mani?re
list.sp$Sp_full
phylo.Ithomiini$tip.label

### Inf?rence des param?tres du mod?le OU ? partir des donn?es (simul?es) et de l'arbre

## Run algorithm
?PhyloEM
PhyloEM.5 <- PhyloEM(phylo = phylo.Ithomiini,
               Y_data = t(sp.env.table),
               process = "scOU",                   ## scalar OU model = OU model with a single selection strength parameter
               random.root = TRUE,                 ## Root is stationary (true model)
               stationary.root = TRUE,
               alpha = NULL,                       ## Default = a grid of 10 alpha values
               K_max = 80,                         ## Maximal default number of shifts is 18 for a phylogeny of this size (see ?PhyloEM)
               parallel_alpha = T,                 ## This can be set to TRUE for
               Ncores = 6)                         ## parallel computations

PhyloEM.5
save(PhyloEM.5, file = paste0(internal.wd, "/Niche_evolution/PhyloEM.5.RData"))
str(PhyloEM.5, max.level = 2)

# Plot estimated shifts on the phylogeny
plot(PhyloEM.5)
# Same but with another method
plot(PhyloEM.5, params = params_process(PhyloEM.5, method.selection = "DDSE"))
# Same but with fixed parameters we can chose
plot(PhyloEM.5, params = params_process(PhyloEM.5, K = 6, alpha = 0.3329135))

# Plot selection criterion for the number of shifts (K). Best = minimum value
plot_criterion(PhyloEM.5)
# Le meilleur fit est obtenu sans aucun shift...

# To rerun the analysis with some fixed parameters and extract the results in general
params_process(PhyloEM.5, K = 6)

# When there are several equivalent solutions 
params_8 <- params_process(PhyloEM.5, K = 8)
plot(equivalent_shifts(tree, params_8))


### Version ? 2 PC apr?s pPCA

# Load data
load(file = paste0(internal.wd, "/Niche_evolution/PC.Revell.env.RData"))

### Inf?rence des param?tres du mod?le OU ? partir des donn?es (simul?es) et de l'arbre

## Run algorithm
?PhyloEM
PhyloEM.pPC <- PhyloEM(phylo = phylo.Ithomiini,
                     Y_data = t(PC.Revell.env),
                     process = "scOU",                   ## scalar OU model = OU model with a single selection strength parameter
                     random.root = TRUE,                 ## Root is stationary (true model)
                     stationary.root = TRUE,
                     alpha = NULL,                       ## Default = a grid of 10 alpha values
                     K_max = 80,                         ## Maximal default number of shifts is 18 for a phylogeny of this size (see ?PhyloEM)
                     parallel_alpha = T,                 ## This can be set to TRUE for
                     Ncores = 6)                         ## parallel computations

PhyloEM.pPC
save(PhyloEM.pPC, file = paste0(internal.wd, "/Niche_evolution/PhyloEM.pPC.RData"))
str(PhyloEM.pPC, max.level = 2)

load(file = paste0(internal.wd, "/Niche_evolution/PhyloEM.pPC.RData"))

# Plot estimated shifts on the phylogeny
plot(PhyloEM.pPC)
# Same but with another method
plot(PhyloEM.pPC, params = params_process(PhyloEM.pPC, method.selection = "DDSE"))
# Same but with fixed parameters we can chose
plot(PhyloEM.pPC, params = params_process(PhyloEM.pPC, K = 6, alpha = 0.3329135))

# Plot selection criterion for the number of shifts (K). Best = minimum value
plot_criterion(PhyloEM.pPC)
# Le meilleur fit est obtenu sans aucun shift...

# To rerun the analysis with some fixed parameters and extract the results in general
params_process(PhyloEM.pPC, K = 6)

# When there are several equivalent solutions 
params_8 <- params_process(PhyloEM.pPC, K = 8)
plot(equivalent_shifts(tree, params_8))



##### Approche OU multi-optima avec SURFACE #####

### Version 5 variables

# Load data
load(file = paste0(internal.wd, "/Phylogenies/Final_phylogeny.RData"))
load(file = paste0(internal.wd,"/sp.env.table.RData"))
sp.env.table <- list.sp[,c("bio1","bio3","bio4","bio12","bio15")]
row.names(sp.env.table) <- list.sp$Sp_full

# V?rifier que les esp?ces sont bien ordonn?es de la m?me mani?re
list.sp$Sp_full
phylo.Ithomiini$tip.label

library(surface)

### Prepare data
phylo.Ithomiini <- nameNodes(phylo.Ithomiini) # label nodes if they are not

?convertTreeData
olist <- convertTreeData(tree = phylo.Ithomiini, dat = sp.env.table) # Conversion au format de Surface
otree <- olist[[1]]; odata <- olist[[2]] # Conversion au format de Surface

### Lauch forward and backward process at the same time
# No need for transformed input here, put directly the original tree with labeled nods and original df !!!
Surface_5_res <- runSurface(tree = phylo.Ithomiini, dat = sp.env.table, exclude = 0, aic_threshold = 0, max_steps = NULL,
                          verbose = T, plotaic = T, error_skip = T, sample_shifts = F, only_best = TRUE)



### Lauch the forward step algorithm which repeatedly call the function addRegime,
### to add selective regime (i.e., shifts), iteratively => a step = 1 addition of a selective regime
?surfaceForward
Surface_5_fwd <- surfaceForward(otree, odata,
                      starting_list = NULL, # To provide a start for the analysis (either from a previous analysis or a pre-build model with some shift already included)
                      aic_threshold = 0, # Delta AICc needed to accept a new model as better than the previous best model from the previous step
                      exclude = 0, # Portion of worst models to discard at each step
                      max_steps = NULL, # To set a maximum number of regimes/shifts to add even if AIC is still improving
                      sample_shifts = F, # To use a sample of previous nearly best models rather than only the best model to explore next step
                      sample_threshold = 2, # When sample_shifts = T ; AICc threshold to define the "best models" sample use to explore the next step
                      error_skip = T, # To skip models in error that could crash the run
                      verbose = T, # To write progress in the console
                      plotaic = T, # To plot interactively the changes in AIC of the best model at each step
                      save_steps = T, filename = paste0(internal.wd,"/Niche_evolution/Surface_Backups/Run1.RData")) # To save in a back-up file the state at each new step in case of crash

# Explore results of forward algorithm
Surface_5_fwd[[i]] # Best model of the step i of the forward algorithm
k <- length(Surface_5_fwd) # Number of shifts found (i.e, step at which the algorithm stopped)
fsum <- surfaceSummary(Surface_5_fwd) ; fsum
fsum$n_steps # Number of iterations = number of regimes/shifts added
fsum$lnls # LnL of each traits at each step
fsum$aics # Evolution of AICc at each step
fsum$shifts # Localisation of shifts on the tree (nod labels) for the best model
fsum$n_regimes_seq # Summary of regime strucutre evolution at each step
fsum$n_regimes # Summary of regime strucutre at the final step
# k = number of regime shifts (counting the initial one)
# kprime = number of optimums (i.e., regimes) => some can be reach by independant shifts
# deltak = k-kprime = measure of the levels of convergence
# c = number of shifts to a convergent regime
# kprime_conv = number of optimum reached by multiple shifts
# kprime_nonconv = number of optimum reached by a unique shift
fsum$alpha # Value of the parameter representing the strength of the selection, for each trait
fsum$sigma_squared # Value of the parameter representing the evolutionnary rate (i.e., level of divergence per dt) => only parameter of a BM
fsum$theta # Matrix of the optimums per trait

### Lauch the backward step algorithm which repeatedly calls the function collapseRegimes,
### to identify cases where the same (or very similar) regimes are found independently on different branches of the tree,
### and where the model simplification obtained by collapsing them into single regimes results in a further improvement in the AICc
?surfaceBackward
Surface_5_bwd <- surfaceBackward(otree, odata, 
                       starting_model = Surface_5_fwd[[k]], # The final model of the forward phase, to attempt regime collapses on.
                       aic_threshold = 0, # Change in AICc needed to accept the collapse
                       max_steps = NULL, # To set a maximum number of regimes/shifts to collapse even if AIC is still improving
                       only_best = TRUE, # To define if only one pair of regime is collapsed at each step. Otherwise, several compatible pairs can be collapse at once
                       sample_shifts = F, # To use a sample of previous nearly best models rather than only the best model to explore next step
                       sample_threshold = 2, # When sample_shifts = T ; AICc threshold to define the "best models" sample use to explore the next step
                       verbose = T, # To write progress in the console
                       plotaic = T, # To plot interactively the changes in AIC of the best model at each step
                       save_steps = F, filename = "") # To save in a back-up file the state at each new step in case of crash

# Explore results of backward algorithm
Surface_5_bwd[[i]] # Best model of the step i of the backward algorithm
Surface_5_bwd[[i]]$fit# Pour extraire le model et l'utiliser dans des simulations, comparaisons, ...
k_final <- length(Surface_5_bwd) # Number of steps needed to reach the best model (i.e, step at which the algorithm stopped)
bsum <- surfaceSummary(Surface_5_bwd) ; bsum
bsum$n_steps # Number of iterations = number of regime shifts collapsed to a common regime
bsum$lnls # LnL of each traits at each step
bsum$aics # Evolution of AICc at each step
bsum$shifts # Localisation of shifts on the tree (nod labels) for the best model
bsum$n_regimes_seq # Summary of regime strucutre evolution at each step
bsum$n_regimes # Summary of regime strucutre at the final step
# k = number of regime shifts (counting the initial one)
# kprime = number of optimums (i.e., regimes) => some can be reach by independant shifts
# deltak = k-kprime = measure of the levels of convergence
# c = number of shifts to a convergent regime
# kprime_conv = number of optimum reached by multiple shifts
# kprime_nonconv = number of optimum reached by a unique shift
bsum$alpha # Value of the parameter representing the strength of the selection, for each trait
bsum$sigma_squared # Value of the parameter representing the evolutionnary rate (i.e., level of divergence per dt) => only parameter of a BM
bsum$theta # Matrix of the optimums per trait

# k = number of regime shifts (counting the initial one)
# kprime = number of optimums (some can be reach by independant shifts)
# deltak = k-kprime = measure of the levels of convergence
# c = number of shifts to a convergent regime
# kprime_conv = number of optimum reached by multiple shifts
# kprime_nonconv = number of optimum reached by a unique shift


# Plot of regimes on the tree. Convergent regimes are plotted under the same color. Non-convergent regimes remain black
surfaceTreePlot(tree, Surface_5_bwd[[k_final]], labelshifts = T)

# Plot of species in the trait space using color of regimes
?surfaceTraitPlot
par(mfrow=c(1,2), mai=c(0.8,0.8,0.2,0.2))
# convcol = T => species colored by regimes/optimums
# convcol = F => species colored by regime shifts
surfaceTraitPlot(dat, Surface_5_bwd[[k_final]], whattraits = c(1,2))
surfaceTraitPlot(dat, Surface_5_bwd[[k_final]], whattraits = c(3,2))
par(mfrow=c(1,1))

# Generate simpler models to test for significance via LRT ou AICc
bm <- startingModel(otree, odata, brownian=TRUE) # BM
ou1 <- startingModel(otree,odata) # OU avec regime unique (starting point de l'algorithm forward)
H12 <- startingModel(otree,odata,shifts=c("26"="H1","13"="H1","5"="H2","19"="H2")) # Model OU avec 2 r?gimes dont on connait les positions des shifts (ex : corr?l? avec l'habitat)
?getAIC

# Plot the evolution of the AICc during the forward and backward run
surfaceAICPlot(Surface_5_fwd, Surface_5_bwd)
# Add AICc of fixed models on the plot to compare
abline(h=bm[[1]]$aic,lty="longdash")
abline(h=H12[[1]]$aic,lty="longdash")
text(c(6,6),c(bm[[1]]$aic, ou1[[1]]$aic, H12[[1]]$aic)-2,c("BM","OU1","H12"),cex=0.5)

### Version ? 2 pPC

# Load data
load(file = paste0(internal.wd, "/Niche_evolution/PC.Revell.env.RData"))

library(surface)

### Prepare data
phylo.Ithomiini <- nameNodes(phylo.Ithomiini) # label nodes if they are not

olist <- convertTreeData(tree = phylo.Ithomiini, dat = data.frame(PC.Revell.env)) # Conversion au format de Surface
otree <- olist[[1]]; odata <- olist[[2]] # Conversion au format de Surface

plot(otree)

### Lauch forward and backward process at the same time
# No need for transformed input here, put directly the original tree with labeled nods and original df !!!
?runSurface
Surface_pPC_res <- runSurface(tree = phylo.Ithomiini, dat = data.frame(PC.Revell.env), exclude = 0, aic_threshold = 0, max_steps = NULL,
                          verbose = T, plotaic = T, error_skip = T, sample_shifts = F, only_best = TRUE)

save(Surface_pPC_res, file = paste0(internal.wd, "/Niche_evolution/Surface_pPC_res.RData"))
# load(file = paste0(internal.wd, "/Niche_evolution/Surface_pPC_res.RData"))

Surface_pPC_fwd <- Surface_pPC_res$fwd
Surface_pPC_bwd <- Surface_pPC_res$bwd

### Lauch the forward step algorithm which repeatedly call the function addRegime,
### to add selective regime (i.e., shifts), iteratively => a step = 1 addition of a selective regime
?surfaceForward
Surface_pPC_fwd <- surfaceForward(otree, odata,
                      starting_list = NULL, # To provide a start for the analysis (either from a previous analysis or a pre-build model with some shift already included)
                      aic_threshold = 0, # Delta AICc needed to accept a new model as better than the previous best model from the previous step
                      exclude = 0, # Portion of worst models to discard at each step
                      max_steps = NULL, # To set a maximum number of regimes/shifts to add even if AIC is still improving
                      sample_shifts = F, # To use a sample of previous nearly best models rather than only the best model to explore next step
                      sample_threshold = 2, # When sample_shifts = T ; AICc threshold to define the "best models" sample use to explore the next step
                      error_skip = T, # To skip models in error that could crash the run
                      verbose = T, # To write progress in the console
                      plotaic = T, # To plot interactively the changes in AIC of the best model at each step
                      save_steps = T, filename = paste0(internal.wd,"/Niche_evolution/Surface_Backups/Run1.RData")) # To save in a back-up file the state at each new step in case of crash

# Explore results of forward algorithm
Surface_pPC_fwd[[i]] # Best model of the step i of the forward algorithm
k <- length(Surface_pPC_fwd) # Number of shifts found (i.e, step at which the algorithm stopped)
fsum <- surfaceSummary(Surface_pPC_fwd) ; fsum
fsum$n_steps # Number of iterations = number of regimes/shifts added
fsum$lnls # LnL of each traits at each step
fsum$aics # Evolution of AICc at each step
fsum$shifts # Localisation of shifts on the tree (nod labels) for the best model
fsum$n_regimes_seq # Summary of regime strucutre evolution at each step
fsum$n_regimes # Summary of regime strucutre at the final step
# k = number of regime shifts (counting the initial one)
# kprime = number of optimums (i.e., regimes) => some can be reach by independant shifts
# deltak = k-kprime = measure of the levels of convergence
# c = number of shifts to a convergent regime
# kprime_conv = number of optimum reached by multiple shifts
# kprime_nonconv = number of optimum reached by a unique shift
fsum$alpha # Value of the parameter representing the strength of the selection, for each trait
fsum$sigma_squared # Value of the parameter representing the evolutionnary rate (i.e., level of divergence per dt) => only parameter of a BM
fsum$theta # Matrix of the optimums per trait

### Lauch the backward step algorithm which repeatedly calls the function collapseRegimes,
### to identify cases where the same (or very similar) regimes are found independently on different branches of the tree,
### and where the model simplification obtained by collapsing them into single regimes results in a further improvement in the AICc
?surfaceBackward
Surface_pPC_bwd <- surfaceBackward(otree, odata, 
                       starting_model = Surface_pPC_fwd[[k]], # The final model of the forward phase, to attempt regime collapses on.
                       aic_threshold = 0, # Change in AICc needed to accept the collapse
                       max_steps = NULL, # To set a maximum number of regimes/shifts to collapse even if AIC is still improving
                       only_best = TRUE, # To define if only one pair of regime is collapsed at each step. Otherwise, several compatible pairs can be collapse at once
                       sample_shifts = F, # To use a sample of previous nearly best models rather than only the best model to explore next step
                       sample_threshold = 2, # When sample_shifts = T ; AICc threshold to define the "best models" sample use to explore the next step
                       verbose = T, # To write progress in the console
                       plotaic = T, # To plot interactively the changes in AIC of the best model at each step
                       save_steps = F, filename = "") # To save in a back-up file the state at each new step in case of crash

# Explore results of backward algorithm
Surface_pPC_bwd[[i]] # Best model of the step i of the backward algorithm
Surface_pPC_bwd[[i]]$fit# Pour extraire le model et l'utiliser dans des simulations, comparaisons, ...
k_final <- length(Surface_pPC_bwd) # Number of steps needed to reach the best model (i.e, step at which the algorithm stopped)
bsum <- surfaceSummary(Surface_pPC_bwd) ; bsum
bsum$n_steps # Number of iterations = number of regime shifts collapsed to a common regime
bsum$lnls # LnL of each traits at each step
bsum$aics # Evolution of AICc at each step
bsum$shifts # Localisation of shifts on the tree (nod labels) for the best model
bsum$n_regimes_seq # Summary of regime strucutre evolution at each step
bsum$n_regimes # Summary of regime strucutre at the final step
# k = number of regime shifts (counting the initial one)
# kprime = number of optimums (i.e., regimes) => some can be reach by independant shifts
# deltak = k-kprime = measure of the levels of convergence
# c = number of shifts to a convergent regime
# kprime_conv = number of optimum reached by multiple shifts
# kprime_nonconv = number of optimum reached by a unique shift
bsum$alpha # Value of the parameter representing the strength of the selection, for each trait
bsum$sigma_squared # Value of the parameter representing the evolutionnary rate (i.e., level of divergence per dt) => only parameter of a BM
bsum$theta # Matrix of the optimums per trait

# k = number of regime shifts (counting the initial one)
# kprime = number of optimums (some can be reach by independant shifts)
# deltak = k-kprime = measure of the levels of convergence
# c = number of shifts to a convergent regime
# kprime_conv = number of optimum reached by multiple shifts
# kprime_nonconv = number of optimum reached by a unique shift


# Plot of regimes on the tree. Convergent regimes are plotted under the same color. Non-convergent regimes remain black
par(mfrow=c(1,1))
surfaceTreePlot(tree = phylo.Ithomiini, Surface_pPC_bwd[[k_final]], labelshifts = T)

# Plot of species in the trait space using color of regimes
?surfaceTraitPlot
# convcol = T => species colored by regimes/optimums
# convcol = F => species colored by regime shifts
surfaceTraitPlot(dat = data.frame(PC.Revell.env), Surface_pPC_bwd[[k_final]], whattraits = c(1,2))

# Generate simpler models to test for significance via LRT ou AICc
bm <- startingModel(otree, odata, brownian=TRUE) # BM
ou1 <- startingModel(otree,odata) # OU avec regime unique (starting point de l'algorithm forward)
H12 <- startingModel(otree,odata,shifts=c("26"="H1","13"="H1","5"="H2","19"="H2")) # Model OU avec 2 r?gimes dont on connait les positions des shifts (ex : corr?l? avec l'habitat)
?getAIC

# Plot the evolution of the AICc during the forward and backward run
surfaceAICPlot(Surface_pPC_fwd, Surface_pPC_bwd)
# Add AICc of fixed models on the plot to compare
abline(h=bm[[1]]$aic,lty="longdash")
abline(h=H12[[1]]$aic,lty="longdash")
text(c(6,6),c(bm[[1]]$aic, ou1[[1]]$aic, H12[[1]]$aic)-2,c("BM","OU1","H12"),cex=0.5)


