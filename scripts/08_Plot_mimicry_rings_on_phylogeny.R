##### Script 08: Plot mimicry rings on the phylogeny #####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

# Plot species phylogeny with associated mimicry membership table

### Input files

# Summary table for species, including mimicry ring membership
# Phylogeny of the Ithomiini (Chazot et al., 2019), with and without OMUs

### Output files

# Species phylogeny with associated mimicry membership table

###

# Clean environment
rm(list = ls())

##### 1/ Load stuff  #####

### Load libraries

# install.packages("BiocManager")
# BiocManager::install("ggtree")

library(ggtree)
library(readxl)
library(tidyverse)

### Load phylogeny

# 339 species
# load(file = paste0("./input_data/Phylogenies/Final_phylogeny.RData"))
Ithomiini_phylo <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny_with_updated_names.rds")

# ### Load summary table species, including mimicry ring membership
# load(file = "./input_data/Summary_tables/list.sp.RData")
# 
# ## Update species names to match the phylogeny
# Taxonomy_update_2021 <- read_excel("./input_data/Taxonomy_update_2021.xlsx")
# 
# # Generate new df
# list_sp_updated_names <- list.sp
# 
# # Convert Genus, Species to character string
# list_sp_updated_names[, c("Genus", "Species")] <- apply(X = list_sp_updated_names[, c("Genus", "Species")], MARGIN = 2, FUN = as.character)
# 
# # Loop to update dataset
# for (i in 1:nrow(Taxonomy_update_2021))
# {
#   # i <- 12
# 
#   Old_sp_full <- paste(Taxonomy_update_2021$Genus_old[i], Taxonomy_update_2021$Species_old[i], sep = ".")
# 
#   # Check and update match for species
#   sp_update_indices <- which(list_sp_updated_names$Sp_full == Old_sp_full)
#   list_sp_updated_names$Genus[sp_update_indices] <- Taxonomy_update_2021$Genus_new[i]
#   list_sp_updated_names$Species[sp_update_indices] <- Taxonomy_update_2021$Species_new[i]
# 
# }
# 
# # Update Sp_full and Subsp_full
# list_sp_updated_names$Sp_full <- paste(list_sp_updated_names$Genus, list_sp_updated_names$Species, sep = ".")
# 
# # Check if some species in the phylogeny are still missing a match in the summary table
# setdiff(Ithomiini_phylo$tip.label, unique(list_sp_updated_names$Sp_full))
# 
# # Save updated species summary table
# saveRDS(list_sp_updated_names, file = "./input_data/Summary_tables/list_sp_updated_names.rds")
# 
# ###

# Load updated species summary table
list_sp <- readRDS(file = "./input_data/Summary_tables/list_sp_updated_names.rds")


##### 2/ Compute mimicry membership incidence matrix #####

# Species not included in the phylogeny
setdiff(list_sp$Sp_full, Ithomiini_phylo$tip.label)

# Species names to update
setdiff(Ithomiini_phylo$tip.label, list_sp$Sp_full)

# Extract ring membership per species
sp_mimicry_membership_list <- str_split(string = list_sp$Mimicry_full, pattern = "(\\.)|(_)") 
sp_mimicry_membership_list <- lapply(X = sp_mimicry_membership_list, FUN = unique)
names(sp_mimicry_membership_list) <- list_sp$Sp_full

# Save species mimicry membership list
saveRDS(object = sp_mimicry_membership_list, file = "./input_data/sp_mimicry_membership_list.rds")

# Get mimicry ring list in alphabetic order
mimicry_list <- sort(unique(unlist(sp_mimicry_membership_list)))

# Build an incidence table
mimicry_membership_incidence_matrix <- matrix(data = 0, nrow = length(Ithomiini_phylo$tip.label), ncol = length(mimicry_list),
                                              dimnames = list(Ithomiini_phylo$tip.label, mimicry_list))
for (i in 1:length(Ithomiini_phylo$tip.label))
{
  # i <- 4
  
 species_name <- Ithomiini_phylo$tip.label[i]
 species_index <- which(names(sp_mimicry_membership_list) == species_name)
 
 sp_ring_list <- sp_mimicry_membership_list[[species_index]]
 
 mimicry_membership_incidence_matrix[i, sp_ring_list] <- 1
}

# Save mimicry membership incidence matrix
saveRDS(mimicry_membership_incidence_matrix, file = "./input_data/mimicry_membership_incidence_matrix.rds")

##### 3/ Plot phylogeny with associated mimicry ring membership ####

# Load species mimicry membership list
sp_mimicry_membership_list <- readRDS(file = "./input_data/sp_mimicry_membership_list.rds")

# Load mimicry membership incidence matrix
mimicry_membership_incidence_matrix <- readRDS(file = "./input_data/mimicry_membership_incidence_matrix.rds")

### 3.1/ Create the color palette to use ####

# Find a sequences of 11 colors to repeat
heatmap_colors <- c("grey", "tan", "brown",
                     "red", "orange", "magenta",
                     "purple", "blue", "skyblue2",
                     "seagreen3","darkgreen")

# Repeat 4 times to have 44 colors = 44 mimicry rings
heatmap_colors_full <- rep(heatmap_colors, times = 4)

# Name colors with their own name to fit values in the df and make the manual discrete scale works
names(heatmap_colors) <- heatmap_colors

# Check palette results
plot(1:11, pch = 16, cex = 3, col = heatmap_colors)
plot(1:44, pch = 16, cex = 3, col = heatmap_colors_full)

### 3.2/ Build the heatmap dataframe with the values associated with the color palette

# Rename tip labels with space instead of "."
Ithomiini_phylo_clean_labels <- Ithomiini_phylo
Ithomiini_phylo_clean_labels$tip.label <- str_replace(string = Ithomiini_phylo$tip.label, pattern = "\\.", replacement = " ")

# Transform matrix into df
mimicry_membership_incidence_df <- as.data.frame(mimicry_membership_incidence_matrix)
row.names(mimicry_membership_incidence_df) <- Ithomiini_phylo_clean_labels$tip.label

# Reorder mimicry ring columns per occurrences ?
mimicry_richness <- apply(X = mimicry_membership_incidence_df , MARGIN = 2, FUN = sum)
mimicry_membership_incidence_df <- mimicry_membership_incidence_df[ , rev(order(mimicry_richness))]

# Replace values with the associated index in the color palette
for (i in 1:ncol(mimicry_membership_incidence_df))
{
  # Replace incidence with each color names associated to this column
  mimicry_membership_incidence_df[mimicry_membership_incidence_df[, i] == 1, i] <- heatmap_colours_full[i]
}


# Generate base plot
tree_plot <- ggtree(tr = Ithomiini_phylo_clean_labels, layout = "rectangular")
# tree_plot <- ggtree(tr = Ithomiini_phylo_clean_labels, layout="circular")

?ggtree::gheatmap

# Add the ring membership data
tree_plot_heatmap <- gheatmap(p = tree_plot,                      # A tree plot (not a data_tree_obj)
                              data = mimicry_membership_incidence_df,         # Must be a dataframe with row = tips. Columns = data
                              width = 1.5,
                              offset = 9,                 # Offset to allow to add multiple heatmaps with different offset aside. (Easier to add multiple variable in a single heatmap, but allow to apply a different color scheme for each matrix)
                              color = "grey80",               # Color of heatmap cell border
                              colnames_position = "top",  # Where to put the label for each heatmap column
                              colnames_angle = 90,        # Rotate labels
                              colnames_offset_y = 0.2,      # Vertical offset for column names
                              hjust = 0,                  # Horizontal adjustment for column names
                              font.size = 2,              # Font size for columns names  
                              legend_title = "Mimicry ring membership")  +        # Legend title for the entire heatmap    
  
  scale_fill_manual(name = "Mimicry ring", values = heatmap_colors, na.value = NA) + # To set manually the colors
  
  geom_tiplab(hjust = 0, size = 2, fontface = 3) + # Add species name on the tips
  # geom_tiplab2(hjust = 0, size = 2, , font = 3) + # Circular version

  # geom_tiplab(aes(label = paste0('italic(', Ithomiini_phylo_clean_labels$tip.label, ')'), parse = T), hjust = 0, size = 2, font = 3) +

  guides(fill="none") + # Remove legend
  
  ggtree::vexpand(ratio = 0.015, direction = 1) # Expand vertically the graph area to have space for heatmap column labels

print(tree_plot_heatmap)

# Save plot
ggsave(plot = tree_plot_heatmap, file = "./graphs/Phylo_with_mimicry_heatmap.pdf",
       height = 12000, width = 3500, units = "px", dpi = 300)




