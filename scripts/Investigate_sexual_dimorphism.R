
##### Script to investigate the frequency of sexual dimorphism and the number of subspecies per OMU #####


# Clean environment
rm(list = ls())

# Load Ithomiini records
Ithomiini_records <- readRDS(file = "./input_data/Ithomiini_records_with_ssp.rds")

# Check nb of sample sites from original coordinates
Ithomiini_records_true_coords <- readRDS(file = "../ithomiini_diversity/sensitive_input_data/Ithomiini_records_with_ssp_and_true_coordinates.rds")
nrow(unique(Ithomiini_records_true_coords[, c("Latitude", "Longitude")]))

# Extract unique combination of species and mimicry (for both sexes)

Ithomiini_units <- Ithomiini_records[, c("Genus", "Species", "Sub.species", "M.mimicry", "F.mimicry")]
Ithomiini_units <- Ithomiini_units[!duplicated(Ithomiini_units), ]

# Detect sexual dimorphism
Ithomiini_units$Dimorphism <- (Ithomiini_units$M.mimicry != Ithomiini_units$F.mimicry)

# Keep one entry per subspecies
Ithomiini_units$Taxa <- paste(Ithomiini_units$Genus, Ithomiini_units$Species, Ithomiini_units$Sub.species, sep = "_")
all_subspecies <- unique(Ithomiini_units$Taxa)

subspecies_dimorphism_df <- data.frame(Taxa = all_subspecies, Dimorphism = NA)
for(i in seq_along(all_subspecies))
{
  # i <- 1
  
  subspecies <- all_subspecies[i]
  
  # Extract entries for one subsepcies
  subsepcies_subset <- Ithomiini_units[Ithomiini_units$taxa == subspecies, ]
  
  # Aggregate information for dimorphism. One case = TRUE
  subspecies_dimorphism_df$Dimorphism[i] <- any(subsepcies_subset$Dimorphism)
}

# Summary of sexual dimorphism

table(Ithomiini_units$Dimorphism)
table(subspecies_dimorphism_df$Dimorphism)

# Summary in %

table(subspecies_dimorphism_df$Dimorphism)/sum(table(subspecies_dimorphism_df$Dimorphism))*100


### Count the number of subspecies per OMUs

load(file = "./input_data/list.models.RData")
all_OMU <- list.models$Tag.model

OMU_df <- data.frame(OMU = all_OMU, count_ssp = NA)
for(i in seq_along(all_OMU))
{
  # i <- 1
  
  OMU <- all_OMU[i]
  split_OMU <- unlist(strsplit(x = OMU, split = "\\."))

  # Extract entries for one OMU
  OMU_subset <- Ithomiini_units[(Ithomiini_units$Genus == split_OMU[1]) & (Ithomiini_units$Species == split_OMU[2]) & ((Ithomiini_units$M.mimicry == split_OMU[3]) | (Ithomiini_units$F.mimicry == split_OMU[3])), ]
  
  # Count nb of subspecies
  OMU_df$count_ssp[i] <- length(unique(OMU_subset$Sub.species))
}

# Summary

summary(OMU_df$count_ssp)
