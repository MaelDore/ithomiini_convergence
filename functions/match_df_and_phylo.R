#### Function that match species list in a data frame and phylogeny and clean them ####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

### Inputs
 # Dataframe with species as rows. Typically, community data with sites as columns
 # Phylogeny of the group

### Outputs
 # Dataframe with only species in the phylogeny
 # Phylogeny with only species in the dataframe
 # Rows are ordered similarly as in the phylogeny
 
### Options
 # Remove species with no data (all values = NA) => argument clean_NA = T
 # Remove species with no presence recorded (all values = 0) => argument clean_null = T

### Function

match_df_and_phylo <- function (df, phylo,
                                species_var = "row.names", # To provide the location of the species list
                                clean_NA = T,   # To remove or not species with data (all values = NA)
                                clean_null = F) # To remove or not species with no presence recorded (all values = 0)
{
  
  # Extract species list from row.names
  if (species_var == "row.names")
  {
    species_list <- row.names(df)
  } else {
  # Extract species list from column
    species_list <- df[[species_var]]
  }
  
  # Remove species not in the phylogeny from the df
  df_for_phylo <- df[(species_list %in% phylo$tip.label), ]
  
  # Print species remove from stack because they are not in the phylogeny
  not_in_phylo <- species_list[!(species_list %in% phylo$tip.label)]
  species_list <- species_list[(species_list %in% phylo$tip.label)]
  if (length(not_in_phylo) > 0)
  {
    cat(paste0("\n", length(not_in_phylo), " species removed from df because they are absent from the phylogeny:\n\n"))
    print(not_in_phylo)
    cat("\n")
  }
  
  if (clean_NA)
  {
    # Clean species with no data if needed
    no_data_indices <- !(complete.cases(df_for_phylo))
    df_for_phylo <- df_for_phylo[!no_data_indices, ]
    
    # Print species remove from df because they have no data
    no_data <- species_list[no_data_indices]
    species_list <- species_list[!no_data_indices]
    if (length(no_data) > 0)
    {
      cat(paste0("\n", length(no_data), " species removed from df because they have no data:\n\n"))
      print(no_data)
      cat("\n")
    }
  }
  
  if (clean_null)
  {
    # Clean species with no presence recorded if needed
    no_presence_indices <- apply(X = df_for_phylo, MARGIN = 1, FUN = function(x) all(x == 0))
    df_for_phylo <- df_for_phylo[!no_presence_indices, ]
    
    # Print species remove from df because they have no data
    no_presence <- species_list[no_presence_indices]
    species_list <- species_list[!no_presence_indices]
    if (length(no_presence) > 0)
    {
      cat(paste0("\n", length(no_presence), " species removed from df because they have no presence recorded:\n\n"))
      print(no_presence)
      cat("\n")
    }
  }
  
  # Prune the phylogeny to keep only species in the df
  pruned_tree <- ape::keep.tip(phylo, species_list)
  
  # Print species remove from phylogeny because they are not in the df
  not_in_df <- phylo$tip.label[!(phylo$tip.label %in% species_list)]
  if (length(not_in_df) > 0)
  {
    cat(paste0("\n", length(not_in_df), " species removed from the phylogeny because they are absent from the curated df:\n\n"))
    print(not_in_df)
    cat("\n")
  }
  
  # Rearrange df so the species are in the same order than in the phylogeny
  df_for_phylo <- df_for_phylo[match(species_list, pruned_tree$tip.label), ]
  
  # Return cleaned df and phylogeny in a list
  cleaned_df_and_phylo <- list(df_for_phylo, pruned_tree)
  return(cleaned_df_and_phylo)
}

