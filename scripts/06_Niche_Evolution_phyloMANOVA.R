##### Script 06: Test fit of Evolutionary models, test for mean Climatic distance and phyloMANOVA #####

# Fit the best neutral evolutionary models
# Test for convergence in the climatic niche via mean climatic distance of comimics
# Run Phylogenetic MANOVA to test for divergence among rings for the climatic niche, accounting for the phylogeny 

### Input files

# Summary tables of species and OMUs, included in the phylogeny or not
# Phylogeny of the Ithomiini (Chazot et al., 2019), with and without OMUs
# pPCA climatic values of the species and OMUS


### Output files

# Summary table for evolutionary models fit to the climatic data, for the 4 bioclimatic variables, and in the pPCA-space (2 axis)
# Best evolutionary model and its parameters
# Heatmap of kappa-lambda model likelihood
# Simulated values for the 719/619 OMUs in the pPCA-space following the best evolutionnary model
# Comimcry matrix for the 619 OMUs
# Global MCD for comimics, and per mimicry rings + null distribution
# Plot null distribution for test of MCD, global and per mimicry ring
# Plot null distribution of Wilks' lambda and pseudo-F of the phyloMANOVA



# Clean environment
rm(list = ls())

##### 1/ Load stuff  #####

### Load phylogenies, sp list and unit list for evolutionary models and MCD analyses (include all species in the phylogeny)

# 339 species
load(file = paste0("./input_data/Phylogenies/Final_phylogeny.RData"))
load(file = paste0("./outputs/Niche_evolution/sp_env_table_339.RData"))

identical(phylo.Ithomiini$tip.label, row.names(sp_env_table_339))

# 719 OMUs
load(file = paste0("./input_data/Phylogenies/Final_units_phylogeny.RData"))
load(file = paste0("./outputs/Niche_evolution/list.unit_phyl_order_719.RData"))
load(file = paste0("./outputs/Niche_evolution/unit.env.table_719.RData"))

identical(phylo.Ithomiini.units$tip.label, list.unit_phyl_order_719$Tag.model)
identical(phylo.Ithomiini.units$tip.label, row.names(unit.env.table_719))

### Load phylogenies, sp list and unit list for the phyloMANOVA analyses (include only rings with N >= 10)

# 322 species
load(file = paste0("input_data/Phylogenies/phylo.Ithomiini_permMANOVA.RData"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/sp_env_table.RData"))

identical(phylo.Ithomiini_permMANOVA$tip.label, row.names(sp_env_table))

# 619 OMUs
load(file = paste0("./input_data/Phylogenies/reduced.phylo.Ithomiini.units.RData"))
load(file = paste0("./outputs/Niche_evolution/reduced.list.unit_phyl_order.RData"))
load(file = paste0("./outputs/Niche_evolution/unit.env.table_619.RData"))

identical(reduced.phylo.Ithomiini.units$tip.label, row.names(unit.env.table_619))
identical(reduced.phylo.Ithomiini.units$tip.label, reduced.list.unit_phyl_order$Tag.model)


##### 2/ Test fit of evolutionary models #####

### 2.1/ On the 4 climatic variables ####

# install.packages(pkgs = "./packages/motmot.2.0_1.1.2.tar.gz", type = "source")

library(motmot.2.0)

?transformPhylo.ML

load(file = paste0("./outputs/Niche_evolution/sp_env_table_339.RData"))
sp_env_mat_339 <- as.matrix(sp_env_table_339[,-1]) ; row.names(sp_env_mat_339) <- sp_env_table_339$Sp_full

# Create summary table for evolutionary models comparison
Evol_mods_summary_df_4var <- data.frame(Model = character(), Likelihood = numeric(), AICc = numeric(), lambda = numeric(), kappa = numeric())
Evol_mods_summary_df_4var$Model <- as.character(Evol_mods_summary_df_4var$Model)

# Brownian Motion model (neutral evolution with phylogenetic signal)
bm.ml <- transformPhylo.ML(phy = phylo.Ithomiini, y = sp_env_mat_339, model = "bm")
bm.ml
Evol_mods_summary_df_4var[1,] <-  c("BM", round(bm.ml$logLikelihood, 3), round(bm.ml$AICc, 3), NA, NA)

# Pagel's lambda model (modulate phylogenetic signal)
lambda.ml <- transformPhylo.ML(phy = phylo.Ithomiini, y = sp_env_mat_339, model = "lambda", profilePlot=T)
lambda.ml
Evol_mods_summary_df_4var[2,] <-  c("Lambda", round(lambda.ml$MaximumLikelihood, 3), round(lambda.ml$AICc, 3), round(lambda.ml$Lambda[1,1],3), NA)
p.value <- 1 - pchisq(2*(lambda.ml$MaximumLikelihood - bm.ml$logLikelihood), 1) ; p.value
bm.ml$AICc - lambda.ml$AICc

# Univariate Pagel's lambda for information
var <- "Tmean"
var <- "Tvar"
var <- "Hmean"
var <- "Hvar"

univariate_test <- as.matrix(sp_env_table_339[,var]) ; row.names(univariate_test) <-  sp_env_table_339$Sp_full 
lambda.ml_univ <- transformPhylo.ML(phy = phylo.Ithomiini, y = univariate_test, model = "lambda", profilePlot=T)
lambda.ml_univ 

# Phylogenetic signal in most variables

# Test for phylogenetic signal based on Pagel's lambda
tree0 <- transformPhylo(phy = phylo.Ithomiini, model = "lambda", y = sp_env_mat_339, lambda = 0) # Get the transformed tree with Pagel's lambda = 0
lambda0.ml <- transformPhylo.ML(phy = tree0, y = sp_env_mat_339, model = "bm") # Estimate best model on BM on transformed tree
p.value <- 1 - pchisq(2*(lambda.ml$MaximumLikelihood - lambda0.ml$logLikelihood), 1) ; p.value # LRT

# Pagel's kappa model (Punctuated equilibrium)
kappa.ml <- transformPhylo.ML(phy = phylo.Ithomiini, y = sp_env_mat_339, model = "kappa", profilePlot=T)
p.value <- 1 - pchisq(2*(kappa.ml$MaximumLikelihood - bm.ml$logLikelihood), 1) ; p.value
bm.ml$AICc - kappa.ml$AICc
Evol_mods_summary_df_4var[3,] <-  c("Kappa", round(kappa.ml$MaximumLikelihood, 3), round(kappa.ml$AICc, 3), NA, round(kappa.ml$Kappa[1,1],3))


# Pagel's lambda-kappa
kappa.lambda.ml <- transformPhylo.ML(phy = phylo.Ithomiini, y = sp_env_mat_339, model="kappa", lambdaEst = T, profilePlot=T)
kappa.lambda.ml
Evol_mods_summary_df_4var[4,] <-  c("Kappa - Lambda", round(kappa.lambda.ml$MaximumLikelihood, 3), round(kappa.lambda.ml$AICc, 3), round(kappa.lambda.ml$lambda,3), round(kappa.lambda.ml$Kappa[1,1],3))
p.value <- 1 - pchisq(2*(kappa.lambda.ml$MaximumLikelihood - bm.ml$logLikelihood), 2) ; p.value
bm.ml$AICc - kappa.lambda.ml$AICc
p.value <- 1 - pchisq(2*(kappa.lambda.ml$MaximumLikelihood - kappa.ml$MaximumLikelihood), 1) ; p.value
kappa.ml$AICc - kappa.lambda.ml$AICc
p.value <- 1 - pchisq(2*(kappa.lambda.ml$MaximumLikelihood - lambda.ml$MaximumLikelihood), 1) ; p.value
lambda.ml$AICc - kappa.lambda.ml$AICc

save(Evol_mods_summary_df_4var, file = paste0("./outputs/Niche_evolution/Phylo_signal/Evol_mods_summary_df_4var.RData"))

# No phylogenetic signal, despite signal in most variables individually. Could be due to model misspecification due to the number of dimensions (Adams & Collyer, 2018)
# Try on a reduced number of dimensions after applying Revell's pPCA.

### 2.2/ On the pPC-axis ####

### Run Revell's PCA on all species in the phylogeny
source(paste0("./scripts/Revell_phyl_pca.R"))

# Create input matrices
C <- vcv(phylo.Ithomiini) # Variance-covariance matrix based on species phylogeny
X <- sp_env_mat_339 # Environmental data matrix

# Check the rows in the two matrices match each other
identical(row.names(C), row.names(X))

PCA.Revell_339 <- Revell_phyl_pca(C, X, mode = "corr")
diag(PCA.Revell_339$Eval) # eigenvalues
PCA.Revell_339$Evec # eigenvectors used to project new objects in the new pPCA space
PCA.Revell_339$S # scores = coordinates of objects in the new space
PCA.Revell_339$L # PC loadings
PCA.Revell_339$ancestral # Standardized ancestral states of traits
PCA.Revell_339$evolVCV # Evolutionary variance/covariance matrix
PCA.Revell_339$evol_corr # Evolutionary correlation matrix
PCA.Revell_339$X # Standardized trait matrix

PCvar <- round(diag(PCA.Revell_339$Eval)/sum(PCA.Revell_339$Eval),3) ; PCvar # % Explained variance per axis
cumVar <- round(cumsum(diag(PCA.Revell_339$Eval))/sum(PCA.Revell_339$Eval),3) ; cumVar # % Cumulative explained variance
# 97.2% of variance explained with the 2 first axes!

save(PCA.Revell_339, PCvar, file = paste0("./outputs/Niche_evolution/Evol_simul/PCA.Revell_339.RData"))

# Extract coordinates for species in the new pPCA climatic space
pPC.env_sp_339 <- PCA.Revell_339$S[,1:2]
save(pPC.env_sp_339, file = paste0("./outputs/Niche_evolution/Evol_simul/pPC.env_sp_339.RData"))

### Fit models on pPC climatic data
library(motmot.2.0)

?transformPhylo.ML

load(file = paste0("./outputs/Niche_evolution/Evol_simul/pPC.env_sp_339.RData"))
identical(phylo.Ithomiini$tip.label, row.names(pPC.env_sp_339))

# Create summary table for evolutionary models comparison
Evol_mods_summary_df_2pPC <- data.frame(Model = character(), Likelihood = numeric(), AICc = numeric(), lambda = numeric(), kappa = numeric())
Evol_mods_summary_df_2pPC$Model <- as.character(Evol_mods_summary_df_2pPC$Model)

# Brownian Motion model (neutral evolution with phylogenetic signal)
bm.ml <- transformPhylo.ML(phy = phylo.Ithomiini, y = pPC.env_sp_339, model = "bm")
bm.ml
Evol_mods_summary_df_2pPC[1,] <-  c("BM", round(bm.ml$logLikelihood, 3), round(bm.ml$AICc, 3), NA, NA)

# Pagel's lambda model (modulate phylogenetic signal)
lambda.ml <- transformPhylo.ML(phy = phylo.Ithomiini, y = pPC.env_sp_339, model = "lambda", profilePlot=T)
lambda.ml
Evol_mods_summary_df_2pPC[2,] <-  c("Lambda", round(lambda.ml$MaximumLikelihood, 3), round(lambda.ml$AICc, 3), round(lambda.ml$Lambda[1,1],3), NA)
p.value <- 1 - pchisq(2*(lambda.ml$MaximumLikelihood - bm.ml$logLikelihood), 1) ; p.value
bm.ml$AICc - lambda.ml$AICc

# Test for phylogenetic signal based on Pagel's lambda
tree0 <- transformPhylo(phy = phylo.Ithomiini, model = "lambda", y = pPC.env_sp_339, lambda = 0) # Get the transformed tree with Pagel's lambda = 0
lambda0.ml <- transformPhylo.ML(phy = tree0, y = pPC.env_sp_339, model = "bm") # Estimate best model on BM on transformed tree
p.value <- 1 - pchisq(2*(lambda.ml$MaximumLikelihood - lambda0.ml$logLikelihood), 1) ; p.value # LRT: Khi = 20.57 ; df = 1 ; p-value < 0.001.


# Pagel's kappa model (Punctuated equilibrium)
kappa.ml <- transformPhylo.ML(phy = phylo.Ithomiini, y = pPC.env_sp_339, model = "kappa", profilePlot=T)
p.value <- 1 - pchisq(2*(kappa.ml$MaximumLikelihood - bm.ml$logLikelihood), 1) ; p.value
bm.ml$AICc - kappa.ml$AICc
Evol_mods_summary_df_2pPC[3,] <-  c("Kappa", round(kappa.ml$MaximumLikelihood, 3), round(kappa.ml$AICc, 3), NA, round(kappa.ml$Kappa[1,1],3))


# Pagel's lambda-kappa
kappa.lambda.ml <- transformPhylo.ML(phy = phylo.Ithomiini, y = pPC.env_sp_339, model="kappa", lambdaEst = T, profilePlot=T)
kappa.lambda.ml
Evol_mods_summary_df_2pPC[4,] <-  c("Kappa - Lambda", round(kappa.lambda.ml$MaximumLikelihood, 3), round(kappa.lambda.ml$AICc, 3), round(kappa.lambda.ml$lambda,3), round(kappa.lambda.ml$Kappa[1,1],3))
p.value <- 1 - pchisq(2*(kappa.lambda.ml$MaximumLikelihood - bm.ml$logLikelihood), 2) ; p.value
bm.ml$AICc - kappa.lambda.ml$AICc
p.value <- 1 - pchisq(2*(kappa.lambda.ml$MaximumLikelihood - kappa.ml$MaximumLikelihood), 1) ; p.value
kappa.ml$AICc - kappa.lambda.ml$AICc
p.value <- 1 - pchisq(2*(kappa.lambda.ml$MaximumLikelihood - lambda.ml$MaximumLikelihood), 1) ; p.value
lambda.ml$AICc - kappa.lambda.ml$AICc

save(Evol_mods_summary_df_2pPC, file = paste0("./outputs/Niche_evolution/Phylo_signal/Evol_mods_summary_df_2pPC.RData"))
write.csv2(Evol_mods_summary_df_2pPC, file = "./tables/Evol_mods_summary_df_2pPC.csv")

### 2.3/ Explore Kappa-Lambda likelihood landscape ####

# Generate the vectors of parameters to explore
lambda <- seq(0, 1, 0.01) ; kappa <- rev(seq(0, 1, 0.01))
# Generate the matrix of results
lk_landscape <- matrix(data = NA, nrow = length(kappa), ncol = length(lambda))

# Loop to compute likelihood for every combination of lambda and kappa
for (i in 1:nrow(lk_landscape)) {
  for (j in 1:ncol(lk_landscape)) {
    tree_transfo <- transformPhylo(phy = phylo.Ithomiini, model = "lambda", lambda = lambda[j]) # Apply lambda first
    tree_transfo <- transformPhylo(phy = tree_transfo, model = "kappa", kappa = kappa[i]) # Then kappa
    kappa.lambda.ml_test <- transformPhylo.ML(phy = tree_transfo, y = pPC.env_sp_339, model="bm") # Estimate best model in BM on transformed tree
    lk_landscape[i,j] <- kappa.lambda.ml_test$logLikelihood 
  }
  print(i)
}
save(lk_landscape, file = paste0("./outputs/Niche_evolution/Phylo_signal/lk_landscape.RData"))

# Find position of the optimum
pos <- arrayInd(which.max(lk_landscape), dim(lk_landscape))
# Extract optimum combination of parameters
lambda.max <- lambda[pos[2]]
kappa.max <- kappa[pos[1]]

hist(lk_landscape[])

# Plot the heatmap properly
pdf(file = "./graphs/Niche_evolution/Phylo_signal/Heatmap_Kappa_Lambda.pdf", height = 6, width = 6)

original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,4))

image(z = t(lk_landscape[nrow(lk_landscape):1,]), xlab = "Lambda", ylab = "Kappa", 
      zlim = c(-1550,-1500), 
      col = hcl.colors(25, "YlOrRd", rev = T),
      main = "Heatmap of Likelihood of Lambda-Kappa models",
      cex.axis = 1.7, cex.lab = 1.8, cex.main = 1)
points(y = kappa.max, x = lambda.max, cex = 2, col = "black", pch = 16)
points(y = kappa.max, x = lambda.max, cex = 1.2, col = "red", pch = 16)

par(mar = original_int_margins)

dev.off()

# Plot the associated scale using raster

library(raster)

pdf(file = "./graphs/Niche_evolution/Phylo_signal/Scale_Heatmap_Kappa_Lambda.pdf", height = 6, width = 7)

original_ext_margins <- par()$oma
par(oma = c(0,0,0,5))

plot(raster(lk_landscape), zlim = c(-1550,-1500), col = hcl.colors(25, "YlOrRd", rev = T), axis.args=list(cex.axis=1.4),
     legend.args = list(text = "              Likelihood", side = 3, 
                        font = 2, line = 1.5, cex = 1.4))

par(oma = original_ext_margins)
dev.off()


# Compute manually best Kappa-Lambda model
tree_maxML <- transformPhylo(phy = phylo.Ithomiini, model = "lambda", lambda = lambda.max) # Apply lambda first
tree_maxML <- transformPhylo(phy = tree_maxML, model = "kappa", kappa = kappa.max) # Then kappa
kappa.lambda.ml_max <- transformPhylo.ML(phy = tree_maxML, y = pPC.env_sp_339, model="bm") # Estimate best model in BM on transformed tree
kappa.lambda.ml_max$logLikelihood ; kappa.lambda.ml$MaximumLikelihood  # Compare to automatic process

kappa.lambda.ml_max$AICc <- kappa.lambda.ml_max$AICc + 4 # Add the penalization associated with the use of 2 parameters that should be included in the AICc

# Actualize summary table with the results from "manual" exploration of Kappa-Lambda landscape
Evol_mods_summary_df_2pPC[4,] <-  c("Kappa - Lambda", round(kappa.lambda.ml_max$logLikelihood, 3), round(kappa.lambda.ml_max$AICc, 3), lambda.max, kappa.max)
save(Evol_mods_summary_df_2pPC, file = paste0("./outputs/Niche_evolution/Phylo_signal/Evol_mods_summary_df_2pPC.RData"))
write.csv2(Evol_mods_summary_df_2pPC, file = "./tables/Evol_mods_summary_df_2pPC.csv")

# Test for significance with LRT
p.value <- 1 - pchisq(2*(kappa.lambda.ml_max$logLikelihood - bm.ml$logLikelihood), 2) ; p.value
bm.ml$AICc - kappa.lambda.ml_max$AICc
p.value <- 1 - pchisq(2*(kappa.lambda.ml_max$logLikelihood - kappa.ml$MaximumLikelihood), 1) ; p.value
kappa.ml$AICc - kappa.lambda.ml_max$AICc

# Is the addition of kappa significant?
p.value <- 1 - pchisq(2*(kappa.lambda.ml_max$logLikelihood - lambda.ml$MaximumLikelihood), 1) ; p.value # LRT: Khi = 4.18 ; df = 1 ; p-value = 0.041.
lambda.ml$AICc - kappa.lambda.ml_max$AICc # Delta AICc = 2.21

# In the end, we chose the lambda model as the best model we can be confident with its parameter estimate.

best_model_2pPCA <- lambda.ml
lamdba_best_model_2pPCA <- best_model_2pPCA$Lambda[1,1]

save(best_model_2pPCA, lamdba_best_model_2pPCA, file = paste0("./outputs/Niche_evolution/Phylo_signal/best_model_2pPCA.RData"))


##### 3/ Traits simulation and comimicry matrices #####

### 3.1/ Traits simulation under null hypothesis ####

load(file = paste0("./outputs/Niche_evolution/Phylo_signal/best_model_2pPCA.RData"))

# install.packages(pkgs = "./packages/phylocurve_2.1.1.tar.gz", type = "source")

library(phylocurve)

?sim.traits # Need to provide the model type, the parameters, the evolutionary covariance matrix and ancestral states estimated from motmot.2.0

# Simulation done at species level since the evolutionary model could not be evaluated for OMUs (vcv matrix not invertible due to duplicate in the phylogenetic tree)
# Since OMUs of the same species occupy the same place in the tree, they will get a value simulated for "within-species observations"
# Thus simulation are conducted with the number of repetitions = max nb of OMUs per species (nreps = 8), to be able to randomly attribute a simulated value to each OMU

OMUs_counts_table <- table(list.unit_phyl_order_719$Sp_full)[order(table(list.unit_phyl_order_719$Sp_full), decreasing = T)] ; OMUs_counts_table
max_nb_OMUs <- OMUs_counts_table[1]

# Simulation at species level for 8 OMUs per species
Sim_clim_2pPCA <- sim.traits(tree = phylo.Ithomiini,                                  # Species tree
                             v = best_model_2pPCA$brownianVariance,                   # Evolutionary covariance matrix of pPCA-transformed climatic traits
                             anc = best_model_2pPCA$root.state,                       # Ancestral state of pPCA-transformed climatic traits
                             model = "lambda",                                        # Pagel's lambda model type
                             parameters = list(lambda = lamdba_best_model_2pPCA),     # Lambda = 0.409
                             nsim = 999,                                              # 999 simulations to get a p-value ranging from 0.001 to 1
                             nreps = max_nb_OMUs,                                     # nreps = 8, to be able to attribute randomly a value for each OMUs of a species
                             return.type = "matrix")
str(Sim_clim_2pPCA_2)
str(Sim_clim_2pPCA)
plot(Sim_clim_2pPCA$tree) # Original tree
plot(Sim_clim_2pPCA$sim_tree) # Transformed tree with lambda = 0.408
Sim_clim_2pPCA$original_X # List of simulation results for all species
Sim_clim_2pPCA$trait_data # List of simulation results for all OMUs (8 per species)

save(Sim_clim_2pPCA, file = paste0("./outputs/Niche_evolution/Evol_simul/Sim_clim_2pPCA.RData"))

# # Example for species 1 (Aeria eurimedia)
# Sim_clim_2pPCA$original_X[[1]][1,] # Values for species
# Sim_clim_2pPCA$trait_data[[1]][Sim_clim_2pPCA$trait_data[[1]]$species == row.names(Sim_clim_2pPCA$original_X[[1]])[1], ] # 8 values for potential OMUs of this species
# 
# # Example for species 1 (Aeria olena)
# Sim_clim_2pPCA$original_X[[1]][2,] # Values for species
# Sim_clim_2pPCA$trait_data[[1]][Sim_clim_2pPCA$trait_data[[1]]$species == row.names(Sim_clim_2pPCA$original_X[[1]])[2], ] # 8 values for potential OMUs of this species

### 3.2/ Store results for the 719 OMUs ####

Sim_clim_2pPCA_OMUs_719 <- list() # Create final list to store results
for (i in 1:length(Sim_clim_2pPCA$trait_data)) # Per simulation
{
  # i <- 1
  
  # Extract matrix of simulated data for species (8 values per species)
  Sim_clim_2pPCA_sp <- Sim_clim_2pPCA$trait_data[[i]]
  
  # Generate new template matrix to store simulated data per OMUs ordered as in list.unit_phyl_order_719
  pPCA_simul_template <- matrix(data = NA, nrow = nrow(list.unit_phyl_order_719), ncol = 2)
  
  # Extract data per species to fill matrix for OMUs
  for (j in 1:nrow(sp_env_table_339)) # Per species
  {
    # j <- 1
    
    sp <- sp_env_table_339$Sp_full[j] # Get species name
    sp_pPCA_8 <- Sim_clim_2pPCA_sp[Sim_clim_2pPCA_sp$species == sp, 2:3] # Extract the 8 values simulated for this species
    
    OMUs_indices <- which(list.unit_phyl_order_719$Sp_full == sp) # Get indices of OMUs for this species
    selection_indices <- sample(x = 1:8, size = length(OMUs_indices), replace = F) # Extract randomly a set of within-species simulation to use as OMUs values
    OMUs_pPCA <- sp_pPCA_8[selection_indices, ] # Extract values for OMUs
    pPCA_simul_template[OMUs_indices, ] <- as.matrix(OMUs_pPCA)  # Store OMUs value in the template for simulated pPCA climatic variables for OMUs
  }
  
  # Store matrix of simulated data in the list of simulations
  Sim_clim_2pPCA_OMUs_719[[i]] <- pPCA_simul_template
  
  save(Sim_clim_2pPCA_OMUs_719, file = paste0("./outputs/Niche_evolution/Evol_simul/Sim_clim_2pPCA_OMUs_719.RData"))
  
  if(i %% 10 == 0)
  {
    cat(paste0("Simulation n° ",i, "/999\n"))
  }

}


### 3.3/ Extract the 619 OMUs included in the perMANOVA and phlyoMANOVA analyses ####

OMUs_619_indices <- which(list.unit_phyl_order_719$Tag.model %in% reduced.list.unit_phyl_order$Tag.model)

Sim_clim_2pPCA_OMUs_619 <- list() # Create final list to store results
for (i in 1:length(Sim_clim_2pPCA_OMUs_719)) # Per simulation
{
  # i <- 1
  
  # Get the data for the 719 OMUs
  pPCA_simul_temp <- Sim_clim_2pPCA_OMUs_719[[i]]
  
  # Extract only the 619 OMUs and store them in the final list
  Sim_clim_2pPCA_OMUs_619 [[i]] <- pPCA_simul_temp[OMUs_619_indices, ]
  
  if(i %% 10 == 0)
  {
    cat(paste0("Simulation n° ",i, "/999\n"))
  }
  
}

save(Sim_clim_2pPCA_OMUs_619, file = paste0("./outputs/Niche_evolution/Evol_simul/Sim_clim_2pPCA_OMUs_619.RData"))


### 3.4/ Compute co-mimicry matrix for the 619 OMUs in the perMANOVA and phlyoMANOVA analyses ####

comimicry_matrix_units_619 <- matrix(nrow = nrow(reduced.list.unit_phyl_order), ncol = nrow(reduced.list.unit_phyl_order), data = 0)
for (i in 1:nrow(reduced.list.unit_phyl_order))
{
  ring_1 <- as.character(reduced.list.unit_phyl_order$Mimicry.model)[i]
  
  for (j in 1:nrow(reduced.list.unit_phyl_order)) 
  {
    ring_2 <- as.character(reduced.list.unit_phyl_order$Mimicry.model)[j]
    
    if (ring_1 == ring_2) 
    {
      comimicry_matrix_units_619[i,j] <- 1
    }
  }
  if (i %% 10 == 0) {print(i)}
}
save(comimicry_matrix_units_619, file = paste0("./outputs/Niche_evolution/comimicry_matrix_units_619.RData"))


##### 4/ Test for convergence in the climatic niche via mean climatic distance of comimics #####

### 4.1/ Project the 719 OMUs climatic data in the pPCA climatic space ####

load(file = paste0("./outputs/Niche_evolution/permMANOVA/PCA.Revell.RData"))
load(file = paste0("./outputs/Niche_evolution/unit.env.table_719.RData"))

# Standardize env traits following Revell's procedure
X <- as.matrix(unit.env.table_719) # Env data for OMU
V <- PCA.Revell$evolVCV        # Evolutionary variance/covariance matrix of environmental traits based on sp phenetic distances and env data
a <- PCA.Revell$ancestral      # Ancestral state of env traits based on sp phenetic distances and env data
eigenV <- PCA.Revell$Evec      # Eigenvectors to project in new space
one <- matrix(1,nrow(X),1)     # Vector to get means
X <- X/(one %*% t(sqrt(diag(V))))  # Standardized env data
pPC.env_units_719 <- (X-one%*%t(a)) %*% eigenV # Project OMu data in new space

save(pPC.env_units_719, file = paste0("./outputs/Niche_evolution/Evol_simul/pPC.env_units_719.RData"))


### 4.2/ Compute MCD and null distri for all comimics ####

load(file = paste0("./outputs/Niche_evolution/comimicry_matrix_units.RData"))
load(file = paste0("./outputs/Niche_evolution/Evol_simul/pPC.env_units_719.RData"))

# Compute the pairwise climatic distance matrix between the 719 OMUs based on pPCA variables
pairwise_climdist_pPCA_719 <- dist(x = pPC.env_units_719, method = "euclidian")
pairwise_climdist_pPCA_719_mat <- as.matrix(pairwise_climdist_pPCA_719)
save(pairwise_climdist_pPCA_719, pairwise_climdist_pPCA_719_mat, file = paste0("./outputs/Niche_evolution/Evol_simul/pairwise_climdist_pPCA_719.RData"))

# Compute mean climatic distances
Global_MCD_obs <- mean(pairwise_climdist_pPCA_719) # Global MCD for all pairs of OMUs = 4.41
Comimic_MCD_obs <- weighted.mean(x = pairwise_climdist_pPCA_719, w = as.dist(comimicry_matrix_units)) # Global MCD only for pairs of comimics = 3.45
# Need to standardized the comimic MCD because some simulations may create a more dilated or reduced climatic space which make comparison of mean distances biased
Comimic_MCD_obs_std <- Comimic_MCD_obs/Global_MCD_obs # Standardized MCD obs = 0.782

save(Global_MCD_obs, Comimic_MCD_obs, Comimic_MCD_obs_std, file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_obs_stats_719.RData"))

# Compute these stats for all the simulations
Global_MCD_null <- Comimic_MCD_null <- Comimic_MCD_std_null <-  NA # Initiate vectors to store results
for (i in 1:length(Sim_clim_2pPCA_OMUs_719)) # Per simulation
{
  # i <- 1
  
  # Compute the pairwise climatic distance matrix between the 719 OMUs based on simulated pPCA variables
  pairwise_climdist_pPCA_719_simul <- dist(x = Sim_clim_2pPCA_OMUs_719[[i]], method = "euclidian")
  
  # Compute mean climatic distances
  Global_MCD_simul <- mean(pairwise_climdist_pPCA_719_simul) # For all pairs of OMUs
  Comimic_MCD_simul <- weighted.mean(x = pairwise_climdist_pPCA_719_simul, w = as.dist(comimicry_matrix_units)) # For comimics only
  # Need to standardized the comimic MCD because some simulations may create a more dilated or reduced climatic space which make comparison of mean distances biased
  Comimic_MCD_std_simul <- Comimic_MCD_simul/Global_MCD_simul # Standardized MCD
                                     
  # Store results                             
  Global_MCD_null[i] <- Global_MCD_simul
  Comimic_MCD_null[i] <- Comimic_MCD_simul
  Comimic_MCD_std_null[i] <- Comimic_MCD_std_simul
  
  save(Global_MCD_null, Comimic_MCD_null, Comimic_MCD_std_null, file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_null_stats_719.RData"))
  
  if(i %% 10 == 0) {cat(paste0("Simulation n° ",i, "/999\n"))}
}
save(Global_MCD_null, Comimic_MCD_null, Comimic_MCD_std_null, file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_null_stats_719.RData"))

summary(Comimic_MCD_std_null)


### 4.3/ Plot the distri of the stats for all comimics ####

load(file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_obs_stats_719.RData"))
load(file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_null_stats_719.RData"))


pdf(file = "./graphs/Niche_evolution/Evol_simul/Comimics_MCD_null_719.pdf", height = 6, width = 7)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

hist(x = c(Comimic_MCD_std_null, Comimic_MCD_obs_std), breaks = 40, freq = TRUE, col = "gray",
     xlim = c(0.75, 1.05),
     # ylim = c(0, 200),
     main = "Distribution of the Mean Climatic Distance \n of co-mimetic OMUs \n under neutral evolution",
     # main = "",
     xlab = "Standardized Mean pairwise Climatic Distance",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)

arrows(Comimic_MCD_obs_std - 0.0003, 53, Comimic_MCD_obs_std - 0.0003, 5, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
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

legend(legend = as.expression(bquote(bold("A"))), 
       x = "topright", inset = c(0.05, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()


### 4.4/ Compute MCD obs and null distri per mimicry ring ####

load(file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_obs_stats_719.RData"))
load(file = paste0("./outputs/Niche_evolution/Evol_simul/pairwise_climdist_pPCA_719.RData"))

# Check if the order is the same
identical(colnames(pairwise_climdist_pPCA_719_mat), list.unit_phyl_order_719$Tag.model)

# Generate list of mimicry rings 
ring_list <- as.character(unique(list.unit_phyl_order_719$Mimicry.model))

# Compute MCD obs per mimicry ring
MCD_per_ring_obs <- MCD_std_per_ring_obs <- rep(NA, length(ring_list)) # Generate empty vector to store MPD obs and its standardized version

for (i in 1:length(ring_list)) # Per mimicry ring
{
  # i <- 10
  
  ring <- ring_list[i]
  
  # Get indices of rows/columns associated of OMUs for this ring
  ring_indices <- which(list.unit_phyl_order_719$Mimicry.model == ring)
  
  if(length(ring_indices) != 1) # Need computation only if there are more than 1 OMU, otherwise, no MCD possible
  {
    # Extract only the pairwise distances for this ring
    pairwise_dist_ring <- as.dist(pairwise_climdist_pPCA_719_mat[ring_indices, ring_indices])
    
    MCD_per_ring_obs[i] <- mean(pairwise_dist_ring) # Compute MCD for this ring
    MCD_std_per_ring_obs[i] <- MCD_per_ring_obs[i]/Global_MCD_obs # Compute standardized MCD for this ring
      
  }
  print(i)
}
names(MCD_per_ring_obs) <- names(MCD_std_per_ring_obs) <- ring_list
MCD_std_per_ring_obs
save(MCD_per_ring_obs, MCD_std_per_ring_obs, file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_per_ring_obs_719.RData"))


# Compute MCD per ring for all the simulations

# Load MCD null stats to extract global MCD for each simulation. Used for standardization.
load(file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_null_stats_719.RData"))

MCD_per_ring_simul <- NA # Initiate vectors to store results
MCD_per_ring_null <- MCD_std_per_ring_null <- matrix(data = NA, nrow = length(Sim_clim_2pPCA_OMUs_719), ncol = length(ring_list)) # Initiate matrix to store results
for (i in 1:length(Sim_clim_2pPCA_OMUs_719)) # Per simulation
{
  # i <- 1
  
  # Compute the pairwise climatic distance matrix between the 619 OMUs based on simulated pPCA variables
  pairwise_climdist_pPCA_719_simul <- as.matrix(dist(x = Sim_clim_2pPCA_OMUs_719[[i]], method = "euclidian"))
  
  for (j in 1:length(ring_list)) # Per mimicry ring
  {
    # j <- 1
    
    ring <- ring_list[j]
    
    # Get indices of rows/columns associated of OMUs for this ring
    ring_indices <- which(list.unit_phyl_order_719$Mimicry.model == ring)
    
    if(length(ring_indices) != 1) # Need computation only if there are more than 1 OMU, otherwise, no MCD possible
    {
      # Extract only the pairwise distances for this ring
      pairwise_dist_ring <- as.dist(pairwise_climdist_pPCA_719_simul[ring_indices, ring_indices])
      
      # Compute MCD for this ring
      MCD_per_ring_simul[j] <- mean(pairwise_dist_ring) 

    }
  }
  
  # Store results for each simulation                        
  MCD_per_ring_null[i, ] <- MCD_per_ring_simul
  MCD_std_per_ring_null[i, ] <- MCD_per_ring_simul/Global_MCD_null[i] # Compute standardized MCD for this simulation
  
  save(MCD_per_ring_null, MCD_std_per_ring_null, file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_per_ring_null_stats_719.RData"))
  
  if(i %% 10 == 0) {cat(paste0("Simulation n° ",i, "/999\n"))}
}
save(MCD_per_ring_null, MCD_std_per_ring_null, file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_per_ring_null_stats_719.RData"))

summary(MCD_std_per_ring_null)


### 4.5/ Plot MCD null distri per mimicry ring and generate summary table ####

load(file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_per_ring_obs_719.RData"))
load(file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_per_ring_null_stats_719.RData")) 

MCD_std_per_ring_obs
MCD_std_per_ring_null

# Reorder in alphabetic order

alphabetic_order <- order(names(MCD_std_per_ring_obs))
MCD_std_per_ring_obs <- MCD_std_per_ring_obs[alphabetic_order]
MCD_std_per_ring_null <- MCD_std_per_ring_null[ , alphabetic_order]

ring_list <- names(MCD_std_per_ring_obs)


MCD_ring_summary_table_719 <- as.data.frame(matrix(ncol = 9, nrow = 44, data = NA))
names(MCD_ring_summary_table_719) <- c("ring", "N_units", "N_pairs", "MCD_obs", "mean_MCD", "MCD_2.5", "MCD_97.5", "p_value", "pattern")


for (i in 1:length(ring_list)) # Per mimicry ring
{
  # i <- 2
  
  MCD_ring_summary_table_719$ring[i] <- ring <- ring_list[i]
  MCD_ring_summary_table_719$N_units[i] <- N_units <- sum(list.unit_phyl_order_719$Mimicry.model == ring)
  MCD_ring_summary_table_719$N_pairs[i] <- N_pairs <- N_units*(N_units-1)/2
  
  if (is.na(MCD_std_per_ring_obs[i])) # Case for ring with only one OMUs. No pairs. No MCD.
  {
    pdf(file = paste0("./graphs/Niche_evolution/Evol_simul/Per_ring/MCD_null_",ring,".pdf"), height = 6, width = 7)
    
    plot(1:100,1:100, type = "n", xlab = "Standardized Mean pairwise Climatic Distance",
         main = paste0("Distribution of Mean Climatic Distance \n of ", ring, " OMUs \n under neutral evolution"))
    text(x = 50, y = 50, labels = "Only one OMU for this mimicry ring \n No pair available for index computation")
    
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
    
    MCD_ring_summary_table_719$MCD_obs[i] <- round(MCD_std_per_ring_obs[i],3)
    MCD_ring_summary_table_719$mean_MCD[i]  <- mean_val
    MCD_ring_summary_table_719$MCD_2.5[i] <- round(quantile(c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i]), 0.025),3)
    MCD_ring_summary_table_719$MCD_97.5[i] <- round(quantile(c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i]), 0.975),3)
    MCD_ring_summary_table_719$p_value[i] <- p_value
    MCD_ring_summary_table_719$pattern[i] <- pattern    
    
    
    histo.save <- hist(MCD_std_per_ring_null[,i],
                       breaks = 30,
                       plot = F)
    
    pdf(file = paste0("./graphs/Niche_evolution/Evol_simul/Per_ring/MCD_null_",ring,".pdf"), height = 6, width = 7)
    
    hist(c(MCD_std_per_ring_null[,i], MCD_std_per_ring_obs[i]), 
         breaks = 30,
         col = "gray", xlab = "Standardized Mean pairwise Climatic Distance", 
         main = paste0("Distribution of the Mean Climatic Distance \n of co-mimetic OMUs \n under neutral evolution"),
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
  
  save(MCD_ring_summary_table_719, file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_ring_summary_table_719.Rdata"))
  
  cat(paste0("N° ",i, "/",length(ring_list)," - ",ring, " - Done \n"))
}

View(MCD_ring_summary_table_719)
write.csv2(MCD_ring_summary_table_719, file = "./tables/MCD_ring_summary_table_719.csv")



# Same but for the 619 OMUs included in the perMANOVA and phlyoMANOVA analyses 

# load(file = paste0("./outputs/Niche_evolution/comimicry_matrix_units_619.RData"))
# load(file = paste0("./outputs/Niche_evolution/permMANOVA/pPC.env_units.RData"))
# 
# # Compute the pairwise climatic distance matrix between the 619 OMUs based on pPCA variables
# pairwise_climdist_pPCA_619 <- dist(x = pPC.env_units, method = "euclidian")
# pairwise_climdist_pPCA_619_mat <- as.matrix(pairwise_climdist_pPCA_619)
# save(pairwise_climdist_pPCA_619, pairwise_climdist_pPCA_619_mat, file = paste0("./outputs/Niche_evolution/Evol_simul/pairwise_climdist_pPCA_619.RData"))
# 
# # Compute mean climatic distances
# Global_MCD_obs <- mean(pairwise_climdist_pPCA_619) # Global MCD for all pairs of OMUs = 4.34
# Comimic_MCD_obs <- weighted.mean(x = pairwise_climdist_pPCA_619, w = as.dist(comimicry_matrix_units_619)) # Global MCD only for pairs of comimics = 3.50
# # Need to standardized the comimic MCD because some simulations may create a more dilated or reduced climatic space which make comparison of mean distances biased
# Comimic_MCD_obs_std <- Comimic_MCD_obs/Global_MCD_obs # Standardized MCD obs = 0.805
# 
# save(Global_MCD_obs, Comimic_MCD_obs, Comimic_MCD_obs_std, file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_obs_stats.RData"))
# 
# # Compute these stats for all the simulations
# Global_MCD_null <- Comimic_MCD_null <- Comimic_MCD_std_null <-  NA # Initiate vectors to store results
# for (i in 1:length(Sim_clim_2pPCA_OMUs)) # Per simulation
# {
#   # i <- 1
#   
#   # Compute the pairwise climatic distance matrix between the 619 OMUs based on simulated pPCA variables
#   pairwise_climdist_pPCA_619_simul <- dist(x = Sim_clim_2pPCA_OMUs[[i]], method = "euclidian")
#   
#   # Compute mean climatic distances
#   Global_MCD_simul <- mean(pairwise_climdist_pPCA_619_simul) # For all pairs of OMUs
#   Comimic_MCD_simul <- weighted.mean(x = pairwise_climdist_pPCA_619_simul, w = as.dist(comimicry_matrix_units_619)) # For comimics only
#   # Need to standardized the comimic MCD because some simulations may create a more dilated or reduced climatic space which make comparison of mean distances biased
#   Comimic_MCD_std_simul <- Comimic_MCD_simul/Global_MCD_simul # Standardized MCD
#   
#   # Store results                             
#   Global_MCD_null[i] <- Global_MCD_simul
#   Comimic_MCD_null[i] <- Comimic_MCD_simul
#   Comimic_MCD_std_null[i] <- Comimic_MCD_std_simul
#   
#   save(Global_MCD_null, Comimic_MCD_null, Comimic_MCD_std_null, file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_null_stats.RData"))
#   
#   if(i %% 10 == 0) {cat(paste0("Simulation n° ",i, "/999\n"))}
# }
# save(Global_MCD_null, Comimic_MCD_null, Comimic_MCD_std_null, file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_null_stats.RData"))
# 
# summary(Comimic_MCD_std_null)
# 
# 
# ### 4.3/ Plot the distri of the stats for all comimics 
# 
# load(file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_obs_stats.RData"))
# load(file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_null_stats.RData"))
# 
# Comimic_MCD_obs_std
# Comimic_MCD_std_null
# 
# pdf(file = "./graphs/Niche_evolution/Evol_simul/Comimics_MCD_null.pdf", height = 6, width = 7)
# 
# original_ext_margins <- par()$oma
# original_int_margins <- par()$mar
# par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))
# 
# hist(x = c(Comimic_MCD_std_null, Comimic_MCD_obs_std), breaks = 40, freq = TRUE, col = "gray",
#      xlim = c(0.80, 1.05),
#      # ylim = c(0, 200),
#      main = "Distribution of the Mean Climatic Distance \n of co-mimetic OMUs \n under neutral evolution",
#      # main = "",
#      xlab = "Standardized Mean pairwise Climatic Distance",
#      cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)
# 
# arrows(Comimic_MCD_obs_std - 0.0031, 60, Comimic_MCD_obs_std - 0.0031, 5, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
# abline(v = mean(c(Comimic_MCD_std_null, Comimic_MCD_obs_std)), lwd = 2, lty = 2) # Add vertical line for mean value
# abline(v = quantile(c(Comimic_MCD_std_null, Comimic_MCD_obs_std), 0.05), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value
# 
# legend(legend = c(paste0("Mean = ", format(round(mean(c(Comimic_MCD_std_null, Comimic_MCD_obs_std)),3), nsmall = 3)), 
#                   paste0("CI 5% = ", round(quantile(c(Comimic_MCD_std_null, Comimic_MCD_obs_std), 0.05),3))), 
#        x = "topleft", inset = c(0.02, 0.05), 
#        lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty = "n")
# legend(legend = c(paste0("MCD obs = ", round(Comimic_MCD_obs_std, 3)),
#                   paste0("p = 0.001")),
#        x = "left", inset = c(-0.045, -0.05),
#        cex = 1.2, bty ="n", xjust = 1)
# 
# legend(legend = as.expression(bquote(bold("A"))), 
#        x = "topright", inset = c(0.05, 0.001), xjust = 0.5,
#        cex = 1.3, bty ="n")
# 
# par(oma = original_ext_margins, mar = original_int_margins)
# 
# dev.off()
# 
# 
# ### 4.4/ Compute MCD obs and null distri per mimicry ring 
# 
# load(file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_obs_stats.RData"))
# load(file = paste0("./outputs/Niche_evolution/Evol_simul/pairwise_climdist_pPCA_619.RData"))
# 
# # Check if the order is the same
# identical(colnames(pairwise_climdist_pPCA_619_mat), reduced.list.unit_phyl_order$Tag.model)
# 
# # Generate list of mimicry rings 
# ring_list <- as.character(unique(reduced.list.unit_phyl_order$Mimicry.model))
# 
# MCD_per_ring_obs <- MCD_std_per_ring_obs <- rep(NA, length(ring_list)) # Generate empty vector to store MPD obs and its standardized version
# 
# for (i in 1:length(ring_list)) # Per mimicry ring
# {
#   # i <-  3
#   
#   ring <- ring_list[i]
#   
#   # Get indices of rows/columns associated of OMUs for this ring
#   ring_indices <- which(reduced.list.unit_phyl_order$Mimicry.model == ring)
#   
#   if(length(ring_indices) != 1) # Need computation only if there are more than 1 OMU, otherwise, no MCD possible
#   {
#     # Extract only the pairwise distances for this ring
#     pairwise_dist_ring <- as.dist(pairwise_climdist_pPCA_619_mat[ring_indices, ring_indices])
#     
#     MCD_per_ring_obs[i] <- mean(pairwise_dist_ring) # Compute MCD for this ring
#     MCD_std_per_ring_obs[i] <- MCD_per_ring_obs[i]/Global_MCD_obs # Compute standardized MCD for this ring
#     
#   }
#   print(i)
# }
# names(MPD_per_ring_obs) <- names(MCD_std_per_ring_obs) <- ring_list
# MCD_std_per_ring_obs
# save(MPD_per_ring_obs, MCD_std_per_ring_obs, file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_per_ring_obs.RData"))



##### 5/ phlyoMANOVA to test for divergence of climatic niche among rings #####

# Need to use only the 619 OMUs included in rings with N >= 10
load(file = paste0("./outputs/Niche_evolution/reduced.list.unit_phyl_order.RData"))
load(file = paste0("./outputs/Niche_evolution/Evol_simul/Sim_clim_2pPCA_OMUs_619.RData"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/pPC.env_units.RData"))

reduced.list.unit_phyl_order
pPC.env_units
Sim_clim_2pPCA_OMUs_619

identical(rownames(pPC.env_units), reduced.list.unit_phyl_order$Tag.model)

### 5.1/ Compute MANOVA on observed data and extract the Wilks' summary stat and the pseudo-F ####

# No use of geiger::aov.phlyo because it does not allow to use data that have already been simulated
# library(geiger)
# ?aov.phylo

MANOVA_obs_pPCA <- manova(pPC.env_units[,1:2] ~ reduced.list.unit_phyl_order$Mimicry.model)
summary(MANOVA_obs_pPCA, test="Wilks") # Global test
summary.aov(MANOVA_obs_pPCA) # Test per response variable (the two pPCA-axis)
Wilks_lambda_obs_pPCA <- summary(MANOVA_obs_pPCA, test="Wilks")$stats[,2][1] # Save the wilk's lambda
Pseudo_F_obs_pPCA <- summary(MANOVA_obs_pPCA, test="Wilks")$stats[,3][1] # Save the associated pseudo-F used for the test (not the phylo test, the regular MANOVA test)
save(Wilks_lambda_obs_pPCA, file = paste0("./outputs/Niche_evolution/Evol_simul/Wilks_lambda_obs_pPCA.RData"))
save(Pseudo_F_obs_pPCA, file = paste0("./outputs/Niche_evolution/Evol_simul/Pseudo_F_obs_pPCA.RData"))


### 5.2/ Compute MANOVA on each simulation and extract the Wilks' summary stat and the pseudo-F ####
Wilks_lambda_null_pPCA <- Pseudo_F_null_pPCA <- NA
for (i in 1:length(Sim_clim_2pPCA_OMUs_619)) # Per simulation
{
  MANOVA_simul_pPCA <- manova(as.matrix(Sim_clim_2pPCA_OMUs_619[[i]]) ~ reduced.list.unit_phyl_order$Mimicry.model)
  Wilks_lambda_null_pPCA[i] <- summary(MANOVA_simul_pPCA, test="Wilks")$stats[,2][1] # Save the wilk's lambda
  Pseudo_F_null_pPCA[i] <- summary(MANOVA_simul_pPCA, test="Wilks")$stats[,3][1] # Save the associated pseudo-F used for the test (not the phylo test, the regular MANOVA test)
  
  if(i %% 100 == 0) {cat(paste0("Simulation n° ",i, "/",length(Sim_clim_2pPCA_OMUs_619),"\n"))}
}
summary(Wilks_lambda_null_pPCA)
summary(Pseudo_F_null_pPCA)
save(Wilks_lambda_null_pPCA, file = paste0("./outputs/Niche_evolution/Evol_simul/Wilks_lambda_null_pPCA.RData"))
save(Pseudo_F_null_pPCA, file = paste0("./outputs/Niche_evolution/Evol_simul/Pseudo_F_null_pPCA.RData"))

### 5.3/ Plot the distri of the Wilks lambda ####

load(file = paste0("./outputs/Niche_evolution/Evol_simul/Wilks_lambda_obs_pPCA.RData"))
load(file = paste0("./outputs/Niche_evolution/Evol_simul/Wilks_lambda_null_pPCA.RData"))


pdf(file = "./graphs/Niche_evolution/Evol_simul/MANOVA_phylo_pPCA_Wilks.pdf", height = 6, width = 7)

hist(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA), 30, freq = TRUE, col = "gray", 
     main = bquote("Phylogenetic MANOVA \n Distribution of Wilks" ~lambda~ "\n under neutral evolution"),
     xlab = bquote("Wilks'" ~lambda),
     xlim = c(0.2, 1.0),
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5, lwd = 2)

arrows(Wilks_lambda_obs_pPCA - 0.002, 120, Wilks_lambda_obs_pPCA - 0.002, 5, length = 0.1, lwd = 2)
abline(v = mean(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA)), lty = 2, lwd = 2)
abline(v = quantile(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA), 0.05), lty = 2, lwd = 2, col = "red")

legend(legend = c(paste0("Mean = ", round(mean(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA)),3)), 
                  paste0("CI 5% = ", round(quantile(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA), 0.05),3))), 
       x = "topleft", inset = c(0.02, 0.00), lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")

legend(legend = bquote(lambda ~ "obs =" ~ .(format(round(Wilks_lambda_obs_pPCA, 3), nsmall = 3))),
       x = "bottomleft", cex = 1.2, bty ="n",
       inset = c(0.03, 0.52))
legend(legend = paste0("p = ", round(ecdf(x = c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA))(Wilks_lambda_obs_pPCA),3)),
       x = "bottomleft", cex = 1.2, bty ="n", y.intersp = 2.1,
       inset = c(0.03, 0.43))

legend(legend = as.expression(bquote(bold("B"))), 
       x = "topright", inset = c(0.05, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

dev.off()

### 5.4/ Plot the distri of the pseudo-F ####

load(file = paste0("./outputs/Niche_evolution/Evol_simul/Pseudo_F_obs_pPCA.RData"))
load(file = paste0("./outputs/Niche_evolution/Evol_simul/Pseudo_F_null_pPCA.RData"))


pdf(file = "./graphs/Niche_evolution/Evol_simul/MANOVA_phylo_pPCA_pseudo_F.pdf", height = 6, width = 7)

hist(c(Pseudo_F_null_pPCA, Pseudo_F_obs_pPCA), 50, freq = TRUE, col = "gray", 
     main = "Phylogenetic MANOVA\n Distribution of pseudo-F\n under neutral evolution",
     xlab = "pseudo-F",
     # xlim = c(0.2, 1.0),
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5, lwd = 2)

arrows(Pseudo_F_obs_pPCA - 0.13, 150, Pseudo_F_obs_pPCA - 0.13, 15, length = 0.1, lwd = 2)
abline(v = mean(c(Pseudo_F_null_pPCA, Pseudo_F_obs_pPCA)), lty = 2, lwd = 2)
abline(v = quantile(c(Pseudo_F_null_pPCA, Pseudo_F_obs_pPCA), 0.95), lty = 2, lwd = 2, col = "red")

legend(legend = c(paste0("Mean = ", round(mean(c(Pseudo_F_null_pPCA, Pseudo_F_obs_pPCA)),2)), 
                  paste0("CI 95% = ", round(quantile(c(Pseudo_F_null_pPCA, Pseudo_F_obs_pPCA), 0.95),2))), 
       x = "topright", inset = c(0.02, 0.15), lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")

legend(legend = paste0("pseudo-F obs = ", round(Pseudo_F_obs_pPCA, 2)),
       x = "bottomright", cex = 1.2, bty ="n",
       inset = c(0.00, 0.44))
legend(legend = paste0("p = 0.001"),
       x = "bottomright", cex = 1.2, bty ="n", y.intersp = 2.1,
       inset = c(0.00, 0.35))

legend(legend = as.expression(bquote(bold("B"))), 
       x = "topright", inset = c(0.05, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

dev.off()


##### 6/ Plot both MCD and Wilks' lambda #####

load(file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_obs_stats_719.RData"))
load(file = paste0("./outputs/Niche_evolution/Evol_simul/MCD_null_stats_719.RData"))

load(file = paste0("./outputs/Niche_evolution/Evol_simul/Wilks_lambda_obs_pPCA.RData"))
load(file = paste0("./outputs/Niche_evolution/Evol_simul/Wilks_lambda_null_pPCA.RData"))

round(Comimic_MCD_obs_std, 3) # MCD obs = 0.782
format(round(Wilks_lambda_obs_pPCA, 3), nsmall = 3) # Lambda obs = 0.271

pdf(file = "./graphs/Niche_evolution/Evol_simul/Wilks_lambda_&_MCD_both_plots.pdf", height = 6.4, width = 14)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,5,1.5,1.5), mfrow = c(1,2))

# Panel A: Wilk's lambda from phylogenetic MANOVA

hist(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA),
     breaks = seq(0.2, 1.0, 0.2/6),
     freq = TRUE, col = "gray", 
     # main = bquote("Phylogenetic MANOVA \n Distribution of Wilks" ~lambda~ "\n under neutral evolution"),
     main = "",
     xlab = bquote("Wilks'" ~lambda),
     xlim = c(0.2, 1.0),
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.2, lwd = 2)

arrows(Wilks_lambda_obs_pPCA + 0.012, 160, Wilks_lambda_obs_pPCA + 0.012, 10, length = 0.1, lwd = 2)
abline(v = mean(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA)), lty = 2, lwd = 2)
abline(v = quantile(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA), 0.05), lty = 2, lwd = 2, col = "red")

legend(legend = c(paste0("Mean = ", round(mean(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA)),3)), 
                  paste0("CI 5% = ", round(quantile(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA), 0.05),3))), 
       x = "topleft", inset = c(0.02, 0.15), lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.6, bty ="n")

# legend(legend = bquote(lambda ~ "obs =" ~ .(format(round(Wilks_lambda_obs_pPCA, 3), nsmall = 3))),
#        x = "bottomleft", cex = 1.6, bty ="n",
#        inset = c(0.03, 0.45))
# legend(legend = paste0("p = ", round(ecdf(x = c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA))(Wilks_lambda_obs_pPCA),3)),
#        x = "bottomleft", cex = 1.6, bty ="n", y.intersp = 2.1,
#        inset = c(0.03, 0.365))

legend(legend = c(as.expression(bquote(bold(paste(lambda ~ "obs = 0.271")))),
                  as.expression(bquote(bold(paste("p = 0.001"))))),
       x = "bottomleft", inset = c(0.03, 0.40), xjust = 1, # Use inset to manually adjust position
       cex = 1.6, bty ="n", bg = "white", box.col = NA)

legend(legend = as.expression(bquote(bold("A"))), 
       x = "topleft", inset = c(-0.03, -0.03), xjust = 0.5,
       cex = 1.8, bty ="n")


# Panel B: Mean Climatic Distances (MCD)

hist(x = c(Comimic_MCD_std_null, Comimic_MCD_obs_std),
     breaks = seq(0.76, 1.05, 0.05/3),
     freq = TRUE, col = "gray",
     xlim = c(0.75, 1.05),
     # ylim = c(0, 200),
     # main = "Distribution of the Mean Climatic Distance \n of co-mimetic OMUs \n under neutral evolution",
     main = "",
     xlab = "Standardized Mean pairwise Climatic Distance",
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.2, lwd = 2)

arrows(Comimic_MCD_obs_std + 0.003, 160, Comimic_MCD_obs_std + 0.003, 10, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = mean(c(Comimic_MCD_std_null, Comimic_MCD_obs_std)), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(c(Comimic_MCD_std_null, Comimic_MCD_obs_std), 0.05), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("Mean = ", format(round(mean(c(Comimic_MCD_std_null, Comimic_MCD_obs_std)),3), nsmall = 3)), 
                  paste0("CI 5% = ", round(quantile(c(Comimic_MCD_std_null, Comimic_MCD_obs_std), 0.05),3))), 
       x = "topleft", inset = c(0.02, 0.15), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.6, bty = "n")

legend(legend = c(as.expression(bquote(bold("MCD obs = 0.782"))),
                  as.expression(bquote(bold("p = 0.001")))),
       x = "left", inset = c(0.045, -0.05),
       cex = 1.6, bty ="n", xjust = 1)

legend(legend = as.expression(bquote(bold("B"))), 
       x = "topleft", inset = c(-0.03, -0.03), xjust = 0.5,
       cex = 1.8, bty ="n")



par(oma = original_ext_margins, mar = original_int_margins, mfrow = c(1,1))

dev.off()


### Same but with MCD as panel A and Wilk's lambda as panel B

# pdf(file = "./graphs/Niche_evolution/Evol_simul/MCD_&_Wilks_lambda_both_plots.pdf", height = 6.4, width = 14)
# 
# original_ext_margins <- par()$oma
# original_int_margins <- par()$mar
# par(oma = c(0,0,0,0), mar = c(4.5,5,1.5,1.3), mfrow = c(1,2))
# 
# hist(x = c(Comimic_MCD_std_null, Comimic_MCD_obs_std),
#      breaks = seq(0.76, 1.05, 0.05/3),
#      freq = TRUE, col = "gray",
#      xlim = c(0.75, 1.05),
#      # ylim = c(0, 200),
#      # main = "Distribution of the Mean Climatic Distance \n of co-mimetic OMUs \n under neutral evolution",
#      main = "",
#      xlab = "Standardized Mean pairwise Climatic Distance",
#      cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.2, lwd = 2)
# 
# arrows(Comimic_MCD_obs_std + 0.003, 160, Comimic_MCD_obs_std + 0.003, 10, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
# abline(v = mean(c(Comimic_MCD_std_null, Comimic_MCD_obs_std)), lwd = 2, lty = 2) # Add vertical line for mean value
# abline(v = quantile(c(Comimic_MCD_std_null, Comimic_MCD_obs_std), 0.05), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value
# 
# legend(legend = c(paste0("Mean = ", format(round(mean(c(Comimic_MCD_std_null, Comimic_MCD_obs_std)),3), nsmall = 3)), 
#                   paste0("CI 5% = ", round(quantile(c(Comimic_MCD_std_null, Comimic_MCD_obs_std), 0.05),3))), 
#        x = "topleft", inset = c(0.02, 0.15), 
#        lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.6, bty = "n")
# 
# legend(legend = c(as.expression(bquote(bold("MCD obs = 0.782"))),
#                   as.expression(bquote(bold("p = 0.001")))),
#        x = "left", inset = c(0.045, -0.05),
#        cex = 1.6, bty ="n", xjust = 1)
# 
# legend(legend = as.expression(bquote(bold("A"))), 
#        x = "topleft", inset = c(-0.03, -0.03), xjust = 0.5,
#        cex = 1.8, bty ="n")
# 
# 
# hist(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA),
#      breaks = seq(0.2, 1.0, 0.2/6),
#      freq = TRUE, col = "gray", 
#      # main = bquote("Phylogenetic MANOVA \n Distribution of Wilks" ~lambda~ "\n under neutral evolution"),
#      main = "",
#      xlab = bquote("Wilks'" ~lambda),
#      xlim = c(0.2, 1.0),
#      cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.2, lwd = 2)
# 
# arrows(Wilks_lambda_obs_pPCA + 0.012, 160, Wilks_lambda_obs_pPCA + 0.012, 10, length = 0.1, lwd = 2)
# abline(v = mean(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA)), lty = 2, lwd = 2)
# abline(v = quantile(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA), 0.05), lty = 2, lwd = 2, col = "red")
# 
# legend(legend = c(paste0("Mean = ", round(mean(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA)),3)), 
#                   paste0("CI 5% = ", round(quantile(c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA), 0.05),3))), 
#        x = "topleft", inset = c(0.02, 0.15), lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.6, bty ="n")
# 
# # legend(legend = bquote(lambda ~ "obs =" ~ .(format(round(Wilks_lambda_obs_pPCA, 3), nsmall = 3))),
# #        x = "bottomleft", cex = 1.6, bty ="n",
# #        inset = c(0.03, 0.45))
# # legend(legend = paste0("p = ", round(ecdf(x = c(Wilks_lambda_null_pPCA, Wilks_lambda_obs_pPCA))(Wilks_lambda_obs_pPCA),3)),
# #        x = "bottomleft", cex = 1.6, bty ="n", y.intersp = 2.1,
# #        inset = c(0.03, 0.365))
# 
# legend(legend = c(as.expression(bquote(bold(paste(lambda ~ "obs = 0.271")))),
#                   as.expression(bquote(bold(paste("p = 0.001"))))),
#        x = "bottomleft", inset = c(0.03, 0.40), xjust = 1, # Use inset to manually adjust position
#        cex = 1.6, bty ="n", bg = "white", box.col = NA)
# 
# legend(legend = as.expression(bquote(bold("B"))), 
#        x = "topleft", inset = c(-0.03, -0.03), xjust = 0.5,
#        cex = 1.8, bty ="n")
# 
# 
# par(oma = original_ext_margins, mar = original_int_margins, mfrow = c(1,1))
# 
# dev.off()
