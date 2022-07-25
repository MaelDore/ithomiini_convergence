##### Script 07: Test significance of convergence using C1 metrics #####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

# Test for convergence in the climatic niche via C1 metrics (Stayton, 2015)
# Null model based on simulation of trait evolution

### Input files

# Comimcry matrix for the 719 OMUs included in the phylogeny
# Phylogeny of the Ithomiini (Chazot et al., 2019), with OMUs
# pPCA climatic values of OMUS
# Simulated pPCA climatic values from the best fitted evolutionary model  


### Output files

# Global mean C1 for comimics, and per mimicry rings + null distribution
# Plot null distribution for test of C1, global and per mimicry ring


# Clean environment
rm(list = ls())

##### 1/ Load stuff  #####

### Load library
library(geiger)
library(gdata)
library(convevol)

### Load phylogenies

# 719 OMUs in the phylogeny
load(file = paste0("./input_data/Phylogenies/Final_units_phylogeny.RData"))
Ithomiini_phylo_OMU <- phylo.Ithomiini.units

### Load the observed data in the pPCA-space
load(file = paste0("./outputs/Niche_evolution/Evol_simul/pPC.env_units_719.RData"))
pPCA_obs_719 <- pPC.env_units_719[, 1:2]

### Load simulated data under best fitted evolutionary model (Lambda = 0.41)

# Chose the way to model intraspecific variance within species (between OMUs from the same species)
# intraspecific_choice <- ""
# intraspecific_choice <- "_null"
intraspecific_choice <- "_from_obs"

# Load the associated simulated data
load(paste0("./outputs/Niche_evolution/Evol_simul/Sim_clim_2pPCA_OMUs_719",intraspecific_choice,".RData"))

### Load the comimicry matrix for the 719 OMUs
load(file = paste0("./outputs/Niche_evolution/comimicry_matrix_units.RData"))

### Load the best fitted evolutionary model
load(file = paste0("./outputs/Niche_evolution/Phylo_signal/best_model_2pPCA.RData"))
lambda_best_model_2pPCA # Get the lambda parameter

### Check order of units

identical(Ithomiini_phylo_OMU$tip.label, row.names(comimicry_matrix_units))


##### 2/ Compute pairwise C1 for all pairs #####

?convevol::convrat
?convevol::convratsig

# Ideally, should reconstruct ancestral state using the best fitted model, and not a BM !
# Solution = provide the transformed phylogeny applying the Lambda from the best fitted model

### 2.1/ Transform the phylogeny as in the best fitted model ####

Ithomiini_phylo_OMU_Lambda <- geiger::rescale(x = Ithomiini_phylo_OMU, model = "lambda", lambda = lambda_best_model_2pPCA)

plot(Ithomiini_phylo_OMU)
plot(Ithomiini_phylo_OMU_Lambda)

### 2.2/ Compute observed C1 values ####

### Function to compute all pairwise C1 values from a dataset

all_pairwise_C1 <- function (df, phylo)
{
  source(file = "./functions/match_df_and_phylo.R")
  clean_df_and_phylo <- match_df_and_phylo(df = df, phylo = phylo)
  
  clean_df <- clean_df_and_phylo[[1]]
  clean_phylo <- clean_df_and_phylo[[2]]
  
  # Initiate the matrix
  pairwise_C1_obs <- matrix(data = NA, nrow = nrow(df), ncol = nrow(df))
  
  # Loop per pair of OMU
  for (i in 2:nrow(df))
  {
    for (j in 1:(i-1))
    {
      pairwise_C1 <- convevol::convrat(phyl = phylo, phendata = df, convtips = row.names(df)[c(i,j)])[1]
      # Fill the lower triangle
      pairwise_C1_obs[i,j] <- pairwise_C1
    }
    
    # Print progress
    if(i %% 1 == 0)
    {
      cat(paste0(Sys.time(), " - C1 computed for species ", i, " / ", nrow(df), "\n"))
    }
  }
  # Duplicate the lower triangle in the upper triangle
  gdata::upperTriangle(pairwise_C1_obs) <- gdata::lowerTriangle(pairwise_C1_obs)
  
  # Names rows and columns
  row.names(pairwise_C1_obs) <- colnames(pairwise_C1_obs) <- row.names(df)
  
  # Export
  return(pairwise_C1_obs)
}

### Compute pairwise observed C1 values
pairwise_C1_obs <- all_pairwise_C1(df = pPCA_obs_719, phylo = Ithomiini_phylo_OMU_Lambda)

# Save C1 matrix of observed data
saveRDS(object = pairwise_C1_obs, file = "./outputs/Niche_evolution/C1_metrics/pairwise_C1_obs.rds")

# # Loop per pair of OMU
# 
# df <- pPCA_obs_719
# phylo <-  Ithomiini_phylo_OMU_Lambda
# 
# # for (i in 2:nrow(df))
# for (i in 701:nrow(df))
# {
#   for (j in 1:(i-1))
#   {
#     pairwise_C1 <- convevol::convrat(phyl = phylo, phendata = df, convtips = row.names(df)[c(i,j)])[1]
#     # Fill the lower triangle
#     pairwise_C1_obs[i,j] <- pairwise_C1
#   }
#   
#   # Print progress
#   if(i %% 1 == 0)
#   {
#     saveRDS(object = pairwise_C1_obs, file = "./outputs/Niche_evolution/C1_metrics/pairwise_C1_obs.rds")
#     cat(paste0(Sys.time(), " - C1 computed for species ", i, " / ", nrow(df), "\n"))
#   }
# }

### 2.3/ Summarize observed values for all comimics ####

pairwise_C1_obs <- readRDS(file = "./outputs/Niche_evolution/C1_metrics/pairwise_C1_obs.rds")
pairwise_C1_obs_dist <- as.dist(pairwise_C1_obs)

# Load comimicry matrix
load(file = paste0("./outputs/Niche_evolution/comimicry_matrix_units.RData"))

identical(row.names(pairwise_C1_obs), row.names(comimicry_matrix_units))

# Compute mean C1 for all pairs
Global_C1_obs <- mean(pairwise_C1_obs_dist) 
Global_C1_obs # Global C1 for all pairs of OMUs = 0.187

# Compute mean C1 for all comimics
Comimic_C1_obs <- weighted.mean(x = pairwise_C1_obs_dist, w = as.dist(comimicry_matrix_units)) 
Comimic_C1_obs # Global C1 only for pairs of comimics = 0.301

# Compute mean C1 for all non-comimics
No_comimic_C1_obs <- weighted.mean(x = pairwise_C1_obs_dist, w = as.dist((comimicry_matrix_units*-1)+1)) 
No_comimic_C1_obs # Global C1 only for pairs of non-comimics = 0.182

save(Global_C1_obs, Comimic_C1_obs, No_comimic_C1_obs, file = paste0("./outputs/Niche_evolution/C1_metrics/C1_obs_global_stats_719.RData"))

### 2.4/ Quick test and plot of comimics vs. non-comimics C1 values ####

## Build df of C1 values
All_comimic_C1_obs <- pairwise_C1_obs_dist[as.dist(comimicry_matrix_units) == 1]
All_not_comimic_C1_obs <- pairwise_C1_obs_dist[as.dist((comimicry_matrix_units*-1)+1) == 1]

C1_values <- c(All_comimic_C1_obs, All_not_comimic_C1_obs)
Mimicry_type <- c(rep("Comimic", length(All_comimic_C1_obs)), rep("Not_mimic", length(All_not_comimic_C1_obs)))
C1_obs_values_df <- as.data.frame(cbind(C1_values, Mimicry_type))

saveRDS(C1_obs_values_df, file = paste0("./outputs/Niche_evolution/C1_metrics/C1_obs_values_df.rds"))

## Kruskal-Wallis test

kruskal_global_C1 <- kruskal.test(C1_values ~ Mimicry_type, data = C1_obs_values_df)
kruskal_global_C1

## ggplot

ggplot(data = C1_obs_values_df, aes(x = Mimicry_type, y = C1_values)) +
  geom_boxplot()
  # Add result of test on the plot


### 2.5/ Summarize observed values per mimicry ring ####

pairwise_C1_obs <- readRDS(file = "./outputs/Niche_evolution/C1_metrics/pairwise_C1_obs.rds")

# Load summary table for the 719 OMUs
list.unit_phyl_order_719 <- readRDS(file = "./outputs/Niche_evolution/list.unit_phyl_order_719.rds")

identical(row.names(pairwise_C1_obs), list.unit_phyl_order_719$Tag.model)

# Generate list of mimicry rings 
ring_list <- as.character(unique(list.unit_phyl_order_719$Mimicry.model))
ring_list <- ring_list[order(ring_list)]

# Compute C1 obs per mimicry ring
mean_C1_per_ring_obs <- sd_C1_per_ring_obs <- rep(NA, length(ring_list)) # Generate empty vector to store C1 obs

for (i in 1:length(ring_list)) # Per mimicry ring
{
  # i <- 10
  
  ring <- ring_list[i]
  
  # Get indices of rows/columns associated of OMUs for this ring
  ring_indices <- which(list.unit_phyl_order_719$Mimicry.model == ring)
  
  if(length(ring_indices) != 1) # Need computation only if there are more than 1 OMU, otherwise, no MCD possible
  {
    # Extract only the pairwise distances for this ring
    pairwise_C1_ring <- as.dist(pairwise_C1_obs[ring_indices, ring_indices])
    
    mean_C1_per_ring_obs[i] <- mean(pairwise_C1_ring) # Compute mean C1 for this ring
    sd_C1_per_ring_obs[i] <- sd(pairwise_C1_ring) # Compute sd of C1 for this ring
  }
  print(i)
}
names(mean_C1_per_ring_obs) <- names(sd_C1_per_ring_obs) <- ring_list
C1_obs_ring_df <- as.data.frame(cbind(mean_C1_per_ring_obs, sd_C1_per_ring_obs))
names(C1_obs_ring_df) <- c("mean_C1", "sd_C1")

View(C1_obs_ring_df)

saveRDS(C1_obs_ring_df, file = paste0("./outputs/Niche_evolution/C1_metrics/C1_obs_ring_df.rds"))

### 2.6/ Test for significance based on Kruskal-Wallis (focal ring vs. non-comimics)

pairwise_C1_obs <- readRDS(file = "./outputs/Niche_evolution/C1_metrics/pairwise_C1_obs.rds")

C1_obs_ring_df <- readRDS(file = paste0("./outputs/Niche_evolution/C1_metrics/C1_obs_ring_df.rds"))
list.unit_phyl_order_719 <- readRDS(file = "./outputs/Niche_evolution/list.unit_phyl_order_719.rds")

C1_obs_values_df <- readRDS(file = paste0("./outputs/Niche_evolution/C1_metrics/C1_obs_values_df.rds"))
All_not_comimic_C1_obs <- C1_obs_values_df$C1_values[C1_obs_values_df$Mimicry_type == "Not_mimic"]

# Loop per ring
C1_obs_ring_df$W_score <- NA
C1_obs_ring_df$p_value <- NA
for (i in 1:length(ring_list))
{
  # i <- 1
  
  ring <- ring_list[i]
  
  ring_indices <- which((list.unit_phyl_order_719$Mimicry.model == ring))
  
  if (length(ring_indices) > 1)
  {
    C1_obs_ring <- as.dist(pairwise_C1_obs[ring_indices, ring_indices])
    
    # Generate temporary df for test
    C1_values <- c(C1_obs_ring, All_not_comimic_C1_obs)
    Mimicry_type <- c(rep("Focal_ring", length(C1_obs_ring)), rep("Not_mimic", length(All_not_comimic_C1_obs)))
    C1_obs_focal_ring_df <- as.data.frame(cbind(C1_values, Mimicry_type))
    
    # Test against not comimics values
    kruskal_ring_C1 <- kruskal.test(C1_values ~ Mimicry_type, data = C1_obs_focal_ring_df)
    
    # Extract results
    W_score <- kruskal_ring_C1$statistic
    p_value <- round(kruskal_ring_C1$p.value, 3)
    
  } else { # Case with no pairs of OMUs
    W_score <- p_value <- NA
  }
  
  # Inform summary table
  C1_obs_ring_df$W_score[i] <- W_score
  C1_obs_ring_df$p_value[i] <- p_value
  
  # Print progress
  cat(paste0(ring, " - Ring n°", i, "/", length(ring_list), "\n"))
  
}

View(C1_obs_ring_df)

saveRDS(C1_obs_ring_df, file = paste0("./outputs/Niche_evolution/C1_metrics/C1_obs_ring_df.rds"))

## More or less everything happens to be significant just because of the power of the test with many not-comimic values...

## Try by subsetting to equal sample size both groups ???


### 3/ Compute simulated C1 values ####

# Store in a list (or an array ?)

### 4/ Tests per pairs of OMU

# Compute summary stats and p-values per pairs

### 5/ Test for all comimics

# Compute summary stats and p-values for all comimics

# Plot histogram for all comimics

### 6/ Tests per mimicry rings

# Compute summary stats and p-values per mimicry ring

# Plot histogram for per mimicry ring