# Function to compute pairwise betadiversity matrix from a sites x species matrix ####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

# Based on betapart R package (Baselga et al., 2022)

# library(betapart)
# ?betapart::phylo.beta.pair

### Inputs
  # Sites/communities x biological units matrix of incidence, abundances, likelihood of presence, ...
  # Phylogeny with all the units

compute_pairwise_betadiversity <- function (regional_data, diversity_type, # "taxo" or "phylo"
                                            phylo, index_family, # "jaccard" or "sorensen"
                                            beta_part) # "total", "nestedness" or "turnover"
{
  ### Prepare data 
  
  # Initiate check for the presence of data
  no_regional_data <- F
  
  # Remove species not found in the phylogeny
  if (diversity_type == "phylo")
  {
    species_in_phylo <- phylo$tip.label
    species_in_communities <- colnames(regional_data)
    species_not_found <- species_in_communities[!(colnames(regional_data) %in% species_in_phylo)]
    regional_data <- regional_data[, colnames(regional_data) %in% species_in_phylo]
    
    nb_sp_not_found <- length(species_not_found)
    
    if(nb_sp_not_found > 0)
    {
      cat(paste0(nb_sp_not_found, " species were removed from community data because they were not found in the phylogeny:\n", paste(species_not_found, collapse = " ")))
    }
  }

  # Clean sites with NA
  regional_data_clean <- na.omit(regional_data) 
  
  # Check there are at least two communities with data in the region
  if (nrow(regional_data_clean) < 2)
  {
    no_regional_data <- T
    beta_value <- NA
  }
  
  ### Compute indices if needed 
  if (!no_regional_data)
  {
    # For taxonomic beta-diversity based on species identity
    if (diversity_type == "taxo")
    {
      # Compute pairwise taxonomic beta-diversity for all pairs of sites
      beta_results <- betapart::beta.pair(x = regional_data_clean, index.family = index_family)
    }
    
    # For phylogenetic beta-diversity based on Faith's Phylogenetic Diversity
    if (diversity_type == "phylo")
    {
      # Compute pairwise phylogenetic beta-diversity for all pairs of sites
      beta_results <- betapart::phylo.beta.pair(x = regional_data_clean, tree = phylo, index.family = index_family)
    }
    
    ### Extract results
    
    # Depending on the partition of beta-diversity
    if (beta_part == "turnover") {beta_mat <- as.matrix(beta_results[[1]])}
    if (beta_part == "nestedness") {beta_mat <- as.matrix(beta_results[[2]])}
    if (beta_part == "total") {beta_mat <- as.matrix(beta_results[[3]])}
    
  }  
  return(beta_mat)
}
