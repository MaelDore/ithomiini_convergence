##### Script 10: Test for niche similarity: Niche overlap #####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

# Quantify niche overlap across pairs of OMUs
   # Apply Kernel-density smoothers on occurrences in the environmental space
   # Correct for environmental availability
   # Define hypervolume boundaries
   # Compute Jaccard overlap
# Investigate correlation between centroid distances and Jaccard's niche overlap
# Tests for significance of niche overlap
   # Global tests for all comimics
   # Test per mimicry rings

### Input files

# Occurrence database
# Stack of environmental variables in the study area (Neotropics)

### Output files

# Density maps of species/OMUs occurrence in the environmental space
# Density map of environmental availability
# Density maps of species/OMUs occupancy in the environmental space

# Hypervolume models per OMUs
# Pairwise niche overlap matrix per OMUs
# Correlation plot: centroid distances vs. Jaccard's niche overlap

# Tests for significance of niche overlap
   # Global tests for all comimics => Plot of global null distribution
   # Test per mimicry rings => Plot of null distribution per mimicry ring


# Clean environment
rm(list = ls())


##### 1/ Prepare data #####

### 1.1/ Load libraries

library(tidyverse)
library(hypervolume)
library(raster)
library(ecodist)
library(isocat)

### 1.2/ Load data ####

### Load occurrence dataset
Ithomiini_dataset <- readRDS(file = "./input_data/Occurrences/Ithomiini_dataset_updated.rds")

# Keep only useful variables
Ithomiini_dataset <- Ithomiini_dataset %>% 
  dplyr::select(ID_obs, Latitude, Longitude, Genus, Species, M.mimicry, F.mimicry) %>% 
  mutate(Dimorphism = !(M.mimicry == F.mimicry))

# Duplicate occurrences for dimorphism cases
Ithomiini_dataset_dimorph <- Ithomiini_dataset %>%
  filter(Dimorphism == TRUE) %>% 
  mutate(Mimicry = F.mimicry)
Ithomiini_dataset <- Ithomiini_dataset %>% 
  mutate(Mimicry = M.mimicry)
Ithomiini_dataset <- rbind(Ithomiini_dataset, Ithomiini_dataset_dimorph) %>% 
  dplyr::select(-c(F.mimicry, M.mimicry, Dimorphism)) %>% 
  mutate(OMU = paste(Genus, Species, Mimicry, sep = "_")) %>%
  dplyr::select(ID_obs, Latitude, Longitude, OMU, Genus, Species, Mimicry)

# Save curated occurrence dataset
saveRDS(object = Ithomiini_dataset, file = "./input_data/Occurrences/Ithomiini_dataset_curated.rds")


### Check compatibility with OMU in BC dataset ###

load(file = paste0("./input_data/Summary_tables/list.unit.RData"))

# New names
unique(Ithomiini_dataset$OMU)[!(str_replace_all(string = unique(Ithomiini_dataset$OMU), pattern = "_", replacement = ".") %in% list.unit$Tag.model)]
# Old names to change
list.unit$Tag.model[!( list.unit$Tag.model %in% str_replace_all(string = unique(Ithomiini_dataset$OMU), pattern = "_", replacement = "."))]

### Load environmental data
# Merra_clim_5_range <- readRDS(file = "./input_data/Env_data/Merra_clim_5_range.rds")
# Merra_clim_5_range_new_names <- readRDS(file = "./input_data/Env_data/Merra_clim_5_range_new_names.rds")
# Env_data <- readRDS(file = "./input_data/Env_data/Env_data_range_Merra_clim_5.rds")

### Load stack of environmental data
Env_stack <- readRDS(file = "./input_data/Env_data/Select_env_15.rds")

# Extract only climatic variables
names(Env_stack) <- c("Tmean", "Tvar", "Ptot", "Pvar", "Elevation", "Forests")
Env_stack <- Env_stack[[1:4]]

plot(Env_stack)

# Save
saveRDS(object = Env_stack, file = "./input_data/Env_data/Env_stack_curated.rds")

### 1.3/ Create environmental stack bounded by Ithomiini range ####

# Load richness map to define Ithomiini range as the study area boundaries
sp_richness <- readRDS("D:/Mael/R_projects/ithomiini_convergence/input_data/SDM_stacks/sp_richness.rds")

Ithomiini_range_raster <- sp_richness
Ithomiini_range_raster[!(Ithomiini_range_raster[] > 0)] <- NA

plot(Ithomiini_range_raster)

Env_stack_within_range <- raster::mask(x = Env_stack, mask = Ithomiini_range_raster)

plot(Env_stack_within_range)

saveRDS(object = Env_stack_within_range, file = paste0("./input_data/Env_data/Env_stack_within_range.rds"))


##### 2/ Extract environmental variables #####

library(raster)

Env_data_per_occ <- raster::extract(x = Env_stack, y = Ithomiini_dataset[, c("Longitude", "Latitude")])
Ithomiini_dataset_with_env_data <- cbind(Ithomiini_dataset, Env_data_per_occ)

saveRDS(object = Ithomiini_dataset_with_env_data, file = paste0("./input_data/Occurrences/Ithomiini_dataset_with_env_data.rds"))



##### 3/ Compute hypervolumes #####


library(hypervolume)

?hypervolume::hypervolume_gaussian


### By setting a weight parameter, the algorithm can instead take a weighted average the kernel functions centered on each observation

### Use the weighted option to build directly the occupancy maps ???
# Apply weights as the inverse of the environment availability ?
# Need to standardize between 0 and 1 if not already.
# Output need to be standardized between 0 and 1 by dividing by max


### 3.1/ Function to build and save all hypervolumes under Gaussian KDE models ####

build_hypervolume_gaussian_KDE_models <- function (data, # Dataframe with occurrences, species names and associated environmental variables
                                                   taxa_range_raster_stack = NULL, # To provide if PCA should be applied on environmental data from the whole clade range (recommended to avoid species sampling bias)
                                                   env_variables, # Names of environmental variables
                                                   unit_names_variable, # Name of the column that host the biological unit names
                                                   center = T, scale = T, # Apply Z-transformation (center and scale)
                                                   apply_PCA = T, # Apply PCA to work on independent axes and account for multicollinearity
                                                   nb_PCA_axes = NULL, # Select the number of PC-axes to keep. Default is all.
                                                   quantile_value = 0.95,
                                                   quantile_type = "probability",
                                                   seed = 1,  # Set seed for reproducibility
                                                   output_path = "./outputs/Niche_similarity/Hypervolumes/") # Output path to save the hypervolume models under .rds format
  
{ 
  # Extract list of biological units
  unit_list <- unique(data[ , unit_names_variable])
  
  # Clean the dataset from occurrences with no environmental data
  data <- data[complete.cases(data[, env_variables]), ]
  
  # Extract only requested environmental variables
  selected_env_data_raw <- data[, selected_variables]
  
  ### Do not apply transformation on the data independently per unit otherwise the axes of the hypervolumes are not comparable anymore !
  
  # Apply Z-transformation (scale and center data) on the whole dataset if requested
  selected_env_data_transfo <- scale(x = selected_env_data_raw, center = center, scale = scale)
  
  # Apply PCA if requested
  if (apply_PCA)
  {
    # If environmental data from the taxa range is provided as a Raster Stack, compute PCA on the clade range data
    if (!is.null(taxa_range_raster_stack))
    {
      # Extract environmental data for the whole study area
      range_env_data <- as.data.frame(raster::getValues(taxa_range_raster_stack))
      
      # Clean the dataset from occurrences with no environmental data
      range_env_data <- range_env_data[complete.cases(range_env_data[, env_variables]), env_variables]
      
      # Select the number of PC-axes to keep
      if (is.null(nb_PCA_axes)) {nb_PCA_axis <- ncol(range_env_data)}
      
      # Run PCA on all study area data
      PCA_results <- ade4::dudi.pca(df = range_env_data, center = center, scale = scale,
                                    scannf = F, nf = nb_PCA_axes)
      
      # Project coordinates of occurrences for all biological units
      selected_env_data_transfo <- ade4::suprow(x = PCA_results, Xsup = selected_env_data_raw)$lisup
      
    } else { # If no data is provided for the clade range, apply PCA on the raw dataset for all biological units at once
      
      PCA_results <- ade4::dudi.pca(df = selected_env_data_transfo, center = center, scale = scale,
                                    scannf = F, nf = ncol(selected_env_data_transfo))
      selected_env_data_transfo <- PCA_results$li
    }
  }
  
  ### Create directory to store model outputs
  if(!dir.exists(paths = output_path)) # Test if the folder exists already or not
  { 
    # Create folder if absent
    dir.create(path = output_path, recursive = T) 
  }
  
  # Set seed for reproducibility
  set.seed(seed = seed)
  
  ### Loop for each biological unit
  for (i in 1:length(unit_list))
  {
    # Get unit name
    unit_name <- unit_list[i]
    
    # Extract biological unit data form the transformed dataset
    unit_dataset <- selected_env_data_transfo[data[, unit_names_variable] == unit_name, , drop = F]
    
    # Remove environmental duplicates
    unit_dataset <- unit_dataset[!duplicated(unit_dataset), ]
    
    ### Estimate the bandwith for the Kernel
    # Silverman estimator =  (4/(n+2))^(1/(n+4)) * m^(-1/(n+4))*sd(X) with n = nb of observations, m = dimensionality, X = data in each dimensions
    # For biological unit with multiple occurrences, use the Silverman's estimate
    if(nrow(unit_dataset) > 1)
    {
      kde_bandwidth <- hypervolume::estimate_bandwidth(data = unit_dataset, method = "silverman")
    } else {
      # For biological unit with a single occurrence, assume a sd = 1 for the Silverman estimate
      kde_bandwidth <- hypervolume::estimate_bandwidth(data = unit_dataset,
                                                       method = "fixed",
                                                       value = (4/(1+2))^(1/(1+4)) * ncol(unit_dataset)^(-1/(1+4)) * rep(x = 1, times = ncol(unit_dataset)))
    }
    
    ### Build the hypervolume model following a Gaussian Kernel Estimate
    unit_hypervolume_model <- hypervolume::hypervolume_gaussian(data = unit_dataset,
                                                                name = unit_name,
                                                                kde.bandwidth = kde_bandwidth, # Default Silverman's estimate for Kernel bandwith
                                                                quantile.requested = quantile_value,
                                                                quantile.requested.type = quantile_type,
                                                                verbose = F)
    # str(unit_hypervolume_model)
    
    # Save unit hypervolume model
    saveRDS(object = unit_hypervolume_model, file = paste0(output_path, "hypervolume_",unit_name,"_gaussianKBE_",quantile_value,"_",quantile_type, ".rds"))
    
    cat(paste0(Sys.time(), " - Gaussian KBE Hypervolume generation - ", i, " out of ", length(unit_list), "\n"))
  }
}

### 3.2/ Build and save hypervolumes for all OMUs ####


# Provide selection of environmental variables
selected_variables <- c("Tmean", "Tvar", "Ptot", "Pvar")

### Run the function to build and save hypervolume Gaussian KDE models
build_hypervolume_gaussian_KDE_models(data = Ithomiini_dataset_with_env_data, # Dataframe with occurrences, biological unit names and associated environmental variables
                                      taxa_range_raster_stack = Env_stack_within_range, # Provide environmental data to apply PCA on the whole clade range (recommended to avoid sampling bias depending on the selection of biological units)
                                      env_variables = selected_variables, # Names of environmental variables
                                      unit_names_variable = "OMU", # Name of the column that host the biological unit names
                                      center = T, scale = T, # Apply Z-transformation (center and scale)
                                      apply_PCA = T, # Apply PCA to work on independent axes and account for multicollinearity
                                      nb_PCA_axes = 2, # Select the number of PCA axes to keep
                                      quantile_value = 1, # Set to one because we don't want any existing occurrence points to be left out of the hypervolume. Binarization to 95% will occur after standardization
                                      quantile_type = "probability",
                                      seed = 1,  # Set seed for reproducibility
                                      output_path = "./outputs/Niche_similarity/Hypervolumes/") # Output path to save the hypervolume models under .rds format

# Be careful when binarizing the hypervolume.
# Need to be rebinarized after weighting by smoothed environmental availability


### 3.3/ Explore hypervolume output ####

# Load
hypervolume_output <- readRDS(file = "./outputs/Niche_similarity/Hypervolumes/hypervolume_Aeria_eurimedia_EURIMEDIA_gaussianKBE_1_probability.rds")

plot(hypervolume_output)

hypervolume_output@RandomPoints # Coordinates of Random points
hypervolume_output@ValueAtRandomPoints # Initial smoothed density scores
hist(hypervolume_output@ValueAtRandomPoints) # Do not range from 0 to 1. Need to be standardized

hypervolume_output@Volume # Volume as computed applying the 100% probability threshold

class(hypervolume_output) # Not thresholded


### 3.4/ Compute hypervolume for the environmental availability ####

Env_stack_within_range <- readRDS(file = paste0("./input_data/Env_data/Env_stack_within_range.rds"))

# Get environmental data as a "species"
Env_brick_within_range <- Env_stack_within_range*1
Available_env_data <- Env_brick_within_range[]
dim(Available_env_data) # 93600 pixel in the raster Stack
Available_env_data <- as.data.frame(Available_env_data[complete.cases(Available_env_data), ])
dim(Available_env_data) # 21415 communities within Ithomiini range

# Add column to define the unit for hypervolume modeling
Available_env_data$Data_type <- "Env"

# Save Available environment dataframe
saveRDS(Available_env_data, file = "./input_data/Env_data/Available_env_data.rds")

# Provide selection of environmental variables
selected_variables <- c("Tmean", "Tvar", "Ptot", "Pvar")

### Do not apply threshold since even if the output is not a TresholdedHypervolume, the RandomPoints provided seems to be limited (i.e., threshold is applied)
# # Select quantile threshold options (will be rerun with higher precision afterwards)
# quantile_value <- 0.95
# quantile_type <- "probability"
# # quantile_type <- "volume"

### Run the function to build and save hypervolume Gaussian KDE models
build_hypervolume_gaussian_KDE_models(data = Available_env_data, # Dataframe with all available environmental data
                                      taxa_range_raster_stack = Env_stack_within_range, # Provide environmental data to apply PCA on the whole clade range (recommended to avoid sampling bias depending on the selection of biological units)
                                      env_variables = selected_variables, # Names of environmental variables
                                      unit_names_variable = "Data_type", # Name of the column that host the biological unit names
                                      center = T, scale = T, # Apply Z-transformation (center and scale)
                                      apply_PCA = T, # Apply PCA to work on independent axes and account for multicollinearity
                                      nb_PCA_axes = 2, # Select the number of PCA axis to keep
                                      quantile_value = 1, # Set to one because we don't want any existing environment to be left out of the hypervolume
                                      quantile_type = "probability",
                                      seed = 1,  # Set seed for reproducibility
                                      output_path = "./outputs/Niche_similarity/Hypervolumes/") # Output path to save the hypervolume models under .rds format

## If too long, set manually the number of point per sample to a lower value


### 3.5/ Rasterize into smoothed density of environmental availiability ####

# Create raster template based on the first 2D of the environment availiability hypervolume
?raster::raster(x)

# Load environment availability hypervolume
Env_availability_hypervolume <- readRDS(file = "./outputs/Niche_similarity/Hypervolumes/hypervolume_Env_gaussianKBE_1_probability.rds")

plot(Env_availability_hypervolume)

Env_random_points <- Env_availability_hypervolume@RandomPoints
Env_random_points <- cbind(Env_random_points, Env_availability_hypervolume@ValueAtRandomPoints)
dim(Env_random_points)
names(Env_random_points)[length(names(Env_random_points))] <- "Density"

# Extract extent of the environmental space
xmin <- floor(min(Env_random_points[,1])*100)/100
xmax <- ceiling(max(Env_random_points[,1])*100)/100
ymin <- floor(min(Env_random_points[,2])*100)/100
ymax <- ceiling(max(Env_random_points[,2])*100)/100

# Set resolution
resolution = c(0.05, 0.05)

Template_raster <- raster(xmn = xmin, xmx = xmax, ymn = ymin, ymx = ymax,
                          resolution = resolution,
                          crs = "", # provide no CRS
                          vals = NULL)

Env_availability_density_map <- rasterize(x = Env_random_points[, c(1,2)], # Matrix of datapoints (or Spatial Object)
                                     y = Template_raster, # Provide the grid to fill with CRS, bbox and resolution
                                     field = Env_random_points[, ncol(Env_random_points)], # How to fill non empty cells. With the value of a varibel in the df of the sp_obj, or directly with a fixed value ?
                                     background = NA, # Value to use to fill empty cells
                                     fun = mean) # What to do when multiple sp features overlap with a cell ?

plot(Env_availability_density_map)

# Resample x3 to smooth output
Env_availability_density_map <- raster::resample(x = Env_availability_density_map, y = Env_availability_density_map, method = "bilinear") 
Env_availability_density_map <- raster::resample(x = Env_availability_density_map, y = Env_availability_density_map, method = "bilinear") 
Env_availability_density_map <- raster::resample(x = Env_availability_density_map, y = Env_availability_density_map, method = "bilinear") 

plot(Env_availability_density_map)

# Standardize range within [0,1]
Env_availability_density_map <- Env_availability_density_map/max(Env_availability_density_map[], na.rm = T)

# Remove low values below the lowest available environment to avoid infinite effect when close to zero
PCA_Available_env_data <- ade4::dudi.pca(df = Available_env_data[, 1:4], center = T, scale = T, scannf = F, nf = 2)$li
Available_env_density_values <- raster::extract(x = Env_availability_density_map, y = PCA_Available_env_data)
Density_threshold <- min(Available_env_density_values)*0.95
Test_map <- Env_availability_density_map
Env_availability_density_map[Env_availability_density_map[] < Density_threshold] <- NA
Test_map[Test_map[] >= Density_threshold] <- NA

# Check all environments are included
plot(Test_map)
points(PCA_Available_env_data, cex = 0.1)
plot(Env_availability_density_map)
points(PCA_Available_env_data, cex = 0.1)

plot(Env_availability_density_map, main = "Environmental Availability")

# Save smoothed density of environmental availability
saveRDS(object = Env_availability_density_map, file = paste0("./outputs/Niche_similarity/Density_maps/Env_availability_density_map.rds"))

### Transform into binary map NA/0 to delineate density maps

Env_availability_binary_map <- Env_availability_density_map
Env_availability_binary_map[Env_availability_binary_map[] >= 0] <- 0
plot(Env_availability_binary_map, main = "Binary Environmental Availability")

# Save Environmental binary map
saveRDS(object = Env_availability_binary_map, file = paste0("./outputs/Niche_similarity/Density_maps/Env_availability_binary_map.rds"))


##### 4/ Build density maps of smoothed occupancy for all OMUs ####

# Load smoothed density of environmental availability
Env_availability_density_map <- readRDS(file = paste0("./outputs/Niche_similarity/Density_maps/Env_availability_density_map.rds"))
# Env_availability_density_map <- readRDS(file = paste0("./outputs/Niche_similarity/Density_maps/Env_availability_density_map_0.95.rds"))

# Load binary map of environmental availability used to delineate all unit density maps
Env_availability_binary_map <- readRDS(file = paste0("./outputs/Niche_similarity/Density_maps/Env_availability_binary_map.rds"))

# Load occurrence dataset
Ithomiini_dataset_with_env_data <- readRDS(file = paste0("./input_data/Occurrences/Ithomiini_dataset_with_env_data.rds"))

# Extract list of OMU
OMU_list <- unique(Ithomiini_dataset_with_env_data$OMU)[order(unique(Ithomiini_dataset_with_env_data$OMU))]

### Loop per OMU
for (i in 1:length(OMU_list))
{
  # i <- 1
  # i <- 6
  
  # Extract OMU name
  OMU <- OMU_list[i]
  
  ### Load OMU's hypervolume
  OMU_hypervolume <- readRDS(file = paste0("./outputs/Niche_similarity/Hypervolumes/hypervolume_",OMU,"_gaussianKBE_1_probability.rds"))
  
  # plot(OMU_hypervolume)
  
  ### Extract random points data
  OMU_random_points <- OMU_hypervolume@RandomPoints
  OMU_random_points <- cbind(OMU_random_points, OMU_hypervolume@ValueAtRandomPoints)
  dim(OMU_random_points)
  names(OMU_random_points)[length(names(OMU_random_points))] <- "Density"
  
  ### Rasterize
  OMU_occurrence_density_map <- rasterize(x = OMU_random_points[, c(1,2)], # Matrix of datapoints (or Spatial Object)
                                          y = Env_availability_density_map, # Provide the grid to fill with CRS, bbox and resolution
                                          field = OMU_random_points[, ncol(OMU_random_points)], # How to fill non empty cells. With the value of a varibel in the df of the sp_obj, or directly with a fixed value ?
                                          background = NA, # Value to use to fill empty cells
                                          fun = mean) # What to do when multiple sp features overlap with a cell ?
  
  # plot(OMU_occurrence_density_map)
  
  ### Resample x2 to smooth output
  OMU_occurrence_density_map <- raster::resample(x = OMU_occurrence_density_map, y = OMU_occurrence_density_map, method = "bilinear") 
  OMU_occurrence_density_map <- raster::resample(x = OMU_occurrence_density_map, y = OMU_occurrence_density_map, method = "bilinear") 
  # plot(OMU_occurrence_density_map)
  
  ### Standardize range within [0,1]
  OMU_occurrence_density_map_std <- OMU_occurrence_density_map/max(OMU_occurrence_density_map[], na.rm = T)
  
  ### Add null values in existing environment
  OMU_occurrence_density_map_std_with_background <- Env_availability_binary_map
  OMU_occurrence_density_map_std_with_background[!is.na(OMU_occurrence_density_map_std[])] <- OMU_occurrence_density_map_std[!is.na(OMU_occurrence_density_map_std[])]
  # plot(OMU_occurrence_density_map_std_with_background, main = "OMU occurrence density map")
  
  ### Save density map of OMU's occurrences ###
  saveRDS(OMU_occurrence_density_map_std_with_background, file = paste0("./outputs/Niche_similarity/Density_maps/Occurrence_density_map_",OMU,".rds"))
  
  ### Obtain density of species occupancy = density of occurrence/environmental availability
  OMU_occupancy_density_map <- OMU_occurrence_density_map_std/Env_availability_density_map
  # plot(OMU_occupancy_density_map)
  
  ### Standardize range within [0,1]
  OMU_occupancy_density_map_std <- OMU_occupancy_density_map/max(OMU_occupancy_density_map[], na.rm = T)

  ### Add null values in existing environment
  OMU_occupancy_density_map_std_with_background <- Env_availability_binary_map
  OMU_occupancy_density_map_std_with_background[!is.na(OMU_occupancy_density_map_std[])] <- OMU_occupancy_density_map_std[!is.na(OMU_occupancy_density_map_std[])]
  # plot(OMU_occupancy_density_map_std_with_background, main = "OMU occupancy density map")
  
  ### Save density map of OMU's occupancy ###
  saveRDS(OMU_occupancy_density_map_std_with_background, file = paste0("./outputs/Niche_similarity/Density_maps/Occupancy_density_map_",OMU,".rds"))
  
  # ### Plot for visual checking
  # par(mfrow = c(2,2))
  # 
  # plot(OMU_occurrence_density_map_std, main = "OMU occurrence density map")
  # points(OMU_hypervolume@Data, cex = 0.5, pch = 16)
  # plot(OMU_occupancy_density_map_std, main = "OMU occupancy density map")
  # points(OMU_hypervolume@Data, cex = 0.5, pch = 16)
  # plot(Env_availability_density_map, main = "Environmental Availability")
  # 
  # par(mfrow = c(1,1))
  
  ### Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Density maps created - ", i, " out of ", length(OMU_list), "\n"))
  }
}


### 5/ Build hypervolumes for occupancy density (i.e., weighted by environmental availability) ####

# Load occurrence dataset
Ithomiini_dataset_with_env_data <- readRDS(file = paste0("./input_data/Occurrences/Ithomiini_dataset_with_env_data.rds"))

# Extract list of OMU
OMU_list <- unique(Ithomiini_dataset_with_env_data$OMU)[order(unique(Ithomiini_dataset_with_env_data$OMU))]

### Loop per OMU
for (i in 1:length(OMU_list))
{
  # i <- 1
  # i <- 6
  
  # Extract OMU name
  OMU <- OMU_list[i]
  
  ### Load OMU's hypervolume
  OMU_hypervolume <- readRDS(file = paste0("./outputs/Niche_similarity/Hypervolumes/hypervolume_",OMU,"_gaussianKBE_1_probability.rds"))
  # plot(OMU_hypervolume)
  
  # Copy to create the hypervolume fo roccupancy
  OMU_hypervolume_occupancy <- OMU_hypervolume 
  
  ### Load OMU's density map of occcupancy
  OMU_occupancy_density_map <- readRDS(file = paste0("./outputs/Niche_similarity/Density_maps/Occupancy_density_map_",OMU,".rds"))
  # plot(OMU_occupancy_density_map)
  
  ### Obtain occupancy value for each RandomPoint
  occupancy_values <- raster::extract(x = OMU_occupancy_density_map, y = OMU_hypervolume@RandomPoints[, c(1,2)])
  occupancy_values <- replace_na(occupancy_values, replace = 0)
  
  ### Paste values in hypervolume
  OMU_hypervolume_occupancy@ValueAtRandomPoints <- occupancy_values
  # plot(OMU_hypervolume_occupancy)
  
  ### Save occupancy hypervolume
  saveRDS(OMU_hypervolume_occupancy, file = paste0("./outputs/Niche_similarity/Hypervolumes/Occupancy_hypervolume_",OMU,".rds"))
  
  ### Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Occupancy hypervolume created - ", i, " out of ", length(OMU_list), "\n"))
  }
  
}


### 6/ Apply 95% probability threshold ####

?hypervolume::hypervolume_threshold

# Select quantile threshold options
quantile_value <- 0.95
quantile_type <- "probability"
# quantile_type <- "volume"

# Load occurrence dataset
Ithomiini_dataset_with_env_data <- readRDS(file = paste0("./input_data/Occurrences/Ithomiini_dataset_with_env_data.rds"))

# Extract list of OMU
OMU_list <- unique(Ithomiini_dataset_with_env_data$OMU)[order(unique(Ithomiini_dataset_with_env_data$OMU))]

### Loop per OMU
for (i in 1:length(OMU_list))
{
  # i <- 1
  # i <- 6
  
  # Extract OMU name
  OMU <- OMU_list[i]
  
  ### Threshold occurrence hypervolume
  
  # Load occurrence hypervolume model
  OMU_hypervolume_occurrence <- readRDS(file = paste0("./outputs/Niche_similarity/Hypervolumes/hypervolume_",OMU,"_gaussianKBE_1_probability.rds"))
  # plot(OMU_hypervolume_occurrence)
  
  # Apply threshold to define the occurrence hypervolume boundaries
  OMU_hypervolume_occurrence_thresholded <- hypervolume::hypervolume_threshold(OMU_hypervolume_occurrence,
                                                                               num.thresholds = 100,
                                                                               quantile.requested = quantile_value,
                                                                               quantile.requested.type = quantile_type,
                                                                               uniform.density = F, plot = F)
    
  plot(OMU_hypervolume_occurrence_thresholded$HypervolumesThresholded)                                                                             
                                                                                                                                                           
  # Save the thresholded occurrence hypervolume 
  saveRDS(object = OMU_hypervolume_occurrence_thresholded,
          file = paste0("./outputs/Niche_similarity/Hypervolumes/hypervolume_occurrence_thresholded_",OMU,"_",quantile_value,"_",quantile_type,".rds"))
  
  ### Threshold occupancy hypervolume
  
  # Load occupancy hypervolume model
  OMU_hypervolume_occupancy <- readRDS(file = paste0("./outputs/Niche_similarity/Hypervolumes/Occupancy_hypervolume_",OMU,".rds"))
  # plot(OMU_hypervolume_occupancy)
  
  # Apply threshold to define the occupancy hypervolume boundaries
  OMU_hypervolume_occupancy_thresholded <- hypervolume::hypervolume_threshold(OMU_hypervolume_occupancy,
                                                                              num.thresholds = 100,
                                                                              quantile.requested = quantile_value,
                                                                              quantile.requested.type = quantile_type,
                                                                              uniform.density = F, plot = F)
  
  plot(OMU_hypervolume_occupancy_thresholded$HypervolumesThresholded)                                                                             
  
  # Save the thresholded occupancy hypervolume 
  saveRDS(object = OMU_hypervolume_occupancy_thresholded,
          file = paste0("./outputs/Niche_similarity/Hypervolumes/hypervolume_occupancy_thresholded_",OMU,"_",quantile_value,"_",quantile_type,".rds"))
  

  ### Print progress
  if (i %% 10 == 0)
  {
    cat(paste0(Sys.time(), " - Occupancy hypervolume created - ", i, " out of ", length(OMU_list), "\n"))
  }
  
}


### 7/ Compute Jaccard's overlap between pairs of hypervolume ####

# See statistics from	Hypervolume_overlap_statistics()
?hypervolume_overlap_statistics  # Favor Jaccard similarity for the meaningfulness and ease of interpretation
?hypervolume_overlap_test  # Could be used to test for niche similarity based on overlap for pairs of OMUs

### 7.1/ Create dissimilarity matrices to fill ####

# Load occurrence dataset
Ithomiini_dataset_with_env_data <- readRDS(file = paste0("./input_data/Occurrences/Ithomiini_dataset_with_env_data.rds"))

# Extract list of OMU
OMU_list <- unique(Ithomiini_dataset_with_env_data$OMU)[order(unique(Ithomiini_dataset_with_env_data$OMU))]

# Build template matrices 
jaccard_occurrence_niche_dissimilarity_matrix <- matrix(data = 0, nrow = length(OMU_list), ncol = length(OMU_list))
jaccard_occupancy_niche_dissimilarity_matrix <- matrix(data = 0, nrow = length(OMU_list), ncol = length(OMU_list))


### Double loop to compute all pairs

# Select quantile threshold options
quantile_value <- 0.95
quantile_type <- "probability"
# quantile_type <- "volume"

for (i in 1:(length(OMU_list)-1))
{
  # i <- 1
  
  # Extract OMU i name
  OMU_i <- OMU_list[i]
  
  # Load thresholded OMU occurrence hypervolume model
  OMU_hypervolume_occurrence_thresholded_i <- readRDS(file = paste0("./outputs/Niche_similarity/Hypervolumes/hypervolume_occurrence_thresholded_",OMU_i,"_",quantile_value,"_",quantile_type,".rds"))
  OMU_hypervolume_occurrence_thresholded_i <- OMU_hypervolume_occurrence_thresholded_i$HypervolumesThresholded

  # Load thresholded OMU occupancy hypervolume model
  OMU_hypervolume_occupancy_thresholded_i <- readRDS(file = paste0("./outputs/Niche_similarity/Hypervolumes/hypervolume_occupancy_thresholded_",OMU_i,"_",quantile_value,"_",quantile_type,".rds"))
  OMU_hypervolume_occupancy_thresholded_i <- OMU_hypervolume_occupancy_thresholded_i$HypervolumesThresholded
  
  for (j in (i+1):length(OMU_list))
  {
    # j <- 6
    
    # Extract species j name
    OMU_j <- OMU_list[j]
    
    # Load thresholded OMU occurrence hypervolume model
    OMU_hypervolume_occurrence_thresholded_j <- readRDS(file = paste0("./outputs/Niche_similarity/Hypervolumes/hypervolume_occurrence_thresholded_",OMU_j,"_",quantile_value,"_",quantile_type,".rds"))
    OMU_hypervolume_occurrence_thresholded_j <- OMU_hypervolume_occurrence_thresholded_j$HypervolumesThresholded
    
    # plot(OMU_hypervolume_occurrence_thresholded_i)
    # plot(OMU_hypervolume_occurrence_thresholded_j)
    
    # Load thresholded OMU occupancy hypervolume model
    OMU_hypervolume_occupancy_thresholded_j <- readRDS(file = paste0("./outputs/Niche_similarity/Hypervolumes/hypervolume_occupancy_thresholded_",OMU_j,"_",quantile_value,"_",quantile_type,".rds"))
    OMU_hypervolume_occupancy_thresholded_j <- OMU_hypervolume_occupancy_thresholded_j$HypervolumesThresholded
    
    # plot(OMU_hypervolume_occupancy_thresholded_i)
    # plot(OMU_hypervolume_occupancy_thresholded_j)
    
    ### 7.2/ Compute Jaccard Similarity index for Niches based on Occurrence density ####
    occurrence_density_hypervolume_set <- suppressMessages(hypervolume_set(hv1 = OMU_hypervolume_occurrence_thresholded_i,
                                                          hv2 = OMU_hypervolume_occurrence_thresholded_j,
                                                          num.points.max = 10000,
                                                          check.memory = FALSE,
                                                          verbose = F))
    jaccard_similarity_occurrence <- hypervolume::hypervolume_overlap_statistics(hvlist = occurrence_density_hypervolume_set)[1]
    jaccard_dissimilarity_occurrence <- round(1 - jaccard_similarity_occurrence, 3)
    
    # Store Jaccard distances for Niches based on Occurrence density
    jaccard_occurrence_niche_dissimilarity_matrix[i,j] <- jaccard_dissimilarity_occurrence

    ### 7.3/ Compute Jaccard Similarity index for Niches based on occupancy density ####
    occupancy_density_hypervolume_set <- suppressMessages(hypervolume_set(hv1 = OMU_hypervolume_occupancy_thresholded_i,
                                                         hv2 = OMU_hypervolume_occupancy_thresholded_j,
                                                         num.points.max = 10000,
                                                         check.memory = FALSE,
                                                         verbose = F))
    jaccard_similarity_occupancy <- hypervolume::hypervolume_overlap_statistics(hvlist = occupancy_density_hypervolume_set)[1]
    jaccard_dissimilarity_occupancy <- round(1 - jaccard_similarity_occupancy, 3)
    
    # Store Jaccard distances for Niches based on occupancy density
    jaccard_occupancy_niche_dissimilarity_matrix[i,j] <- jaccard_dissimilarity_occupancy
    
    
  }
  
  # Save temporary matrices
  saveRDS(jaccard_occurrence_niche_dissimilarity_matrix, file = paste0("./outputs/Niche_similarity/Jaccard_occurrence_niche_dissimilarity_matrix_",quantile_value,"_",quantile_type,".rds"))
  saveRDS(jaccard_occupancy_niche_dissimilarity_matrix, file = paste0("./outputs/Niche_similarity/Jaccard_occupancy_niche_dissimilarity_matrix_",quantile_value,"_",quantile_type,".rds"))
  
  cat(paste0(Sys.time(), " - Pairwise Jaccard niche overlap computed for OMU ", i, " out of ", length(OMU_list), "\n"))
}

### Fill the lower-triangle
fill_lower_triangle <- function (m)
{
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
}

jaccard_occupancy_niche_dissimilarity_matrix <- fill_lower_triangle(jaccard_occupancy_niche_dissimilarity_matrix)
jaccard_occurrence_niche_dissimilarity_matrix <- fill_lower_triangle(jaccard_occurrence_niche_dissimilarity_matrix)

# Add row and col names
row.names(jaccard_occupancy_niche_dissimilarity_matrix) <- colnames(jaccard_occupancy_niche_dissimilarity_matrix) <- OMU_list
row.names(jaccard_occurrence_niche_dissimilarity_matrix) <- colnames(jaccard_occurrence_niche_dissimilarity_matrix) <- OMU_list

# Save
saveRDS(jaccard_occurrence_niche_dissimilarity_matrix, file = paste0("./outputs/Niche_similarity/Jaccard_occurrence_niche_dissimilarity_matrix_",quantile_value,"_",quantile_type,".rds"))
saveRDS(jaccard_occupancy_niche_dissimilarity_matrix, file = paste0("./outputs/Niche_similarity/Jaccard_occupancy_niche_dissimilarity_matrix_",quantile_value,"_",quantile_type,".rds"))

### 7.4/ Explore Jaccard dissimilarities ####

# Load dissimilarity matrices
jaccard_occurrence_niche_dissimilarity_matrix <- readRDS(file = paste0("./outputs/Niche_similarity/Jaccard_occurrence_niche_dissimilarity_matrix_",quantile_value,"_",quantile_type,".rds"))
jaccard_occupancy_niche_dissimilarity_matrix <- readRDS(file = paste0("./outputs/Niche_similarity/Jaccard_occupancy_niche_dissimilarity_matrix_",quantile_value,"_",quantile_type,".rds"))

hist(jaccard_occurrence_niche_dissimilarity_matrix)
View(jaccard_occurrence_niche_dissimilarity_matrix)

hist(jaccard_occupancy_niche_dissimilarity_matrix)
View(jaccard_occupancy_niche_dissimilarity_matrix)


### Plot histogram of Jaccard's niche overlap based on occurrences

pdf(file = paste0("./graphs/Niche_similarity/Jaccard_overlap_occurrences_histo.pdf"), height = 6, width = 7)

hist(jaccard_occurrence_niche_dissimilarity_matrix, 
     col = "gray", xlab = "Jaccard's niche non-overlap", 
     ylab = "Frequency",
     main = paste0("Jaccard's niche non-overlap\nacross all pairs of OMUs"),
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)

dev.off()



### Mantel test

# Mantel_test_niche_overlaps <- ecodist::mantel(formula = jaccard_occupancy_niche_dissimilarity_matrix ~ jaccard_occurrence_niche_dissimilarity_matrix)

Mantel_test_niche_overlaps <- vegan::mantel(xdis = jaccard_occurrence_niche_dissimilarity_matrix,
                                            ydis = jaccard_occupancy_niche_dissimilarity_matrix,
                                            permutations = 999,
                                            method = "spearman", na.rm = T)
Mantel_test_niche_overlaps

# Save Mantel test
saveRDS(Mantel_test_niche_overlaps, file = "./outputs/Niche_similarity/Correlations/Mantel_test_niche_overlaps.rds")

### MRM test

library(ecodist)
?ecodist::MRM

# Load home-made function which can provide summary of permutation outputs
source("./functions/MRM_verbose.R")

# Prepare distances data = linearized in vectors
jaccard_occurrence_niche_dist <- as.vector(as.dist(jaccard_occurrence_niche_dissimilarity_matrix))
jaccard_occupancy_niche_dist <- as.vector(as.dist(jaccard_occupancy_niche_dissimilarity_matrix))

# # Standardized data to obtain beta-coefficents
# jaccard_occurrence_niche_dist <- scale(jaccard_occurrence_niche_dist, center = T, scale = T)
# jaccard_occupancy_niche_dist <- scale(jaccard_occupancy_niche_dist, center = T, scale = T)

## MRM: Occupancy ~ Occurrence
MRM_ranks_niche_overlaps <- MRM_verbose(formula = jaccard_occupancy_niche_dist ~ jaccard_occurrence_niche_dist,
                                        nperm = 1000,
                                        method = "linear", # To use Linear Model with normal distribution assumed for the conditional response variable
                                        mrank = TRUE,  # Transform data in ranks
                                        perm_output = TRUE)  # Provide permutation outputs

print(MRM_ranks_niche_overlaps)
print(MRM_ranks_niche_overlaps$coef)
print(MRM_ranks_niche_overlaps$r.squared)
print(MRM_ranks_niche_overlaps$F.test)


### Correlation plot

# Quick plot
plot(jaccard_occupancy_niche_dissimilarity_matrix ~ jaccard_occurrence_niche_dissimilarity_matrix,
     main = "Correlation between pairwise niche overlap metrics",
     xlab = "Occurrence-based niche overlap",
     ylab = "Occupancy-based niche overlap")


# Function to plot scatterplot of pairwise distances
plot_pairwise_distances_with_MRM <- function(title, cex_title = 1.3,
                                             y, x, y_lab, x_lab,
                                             y_lim = NULL,
                                             cex_axis = 1.5, cex_lab = 1.6, cex_legend = 1.4,
                                             beta_value, p_value, regression,
                                             legend_position = "topleft", legend_hjust = 1,
                                             inset_blank = c(0, 0.02), inset_beta = c(-0.01, 0.02), inset_p_value = c(-0.01, 0.02),
                                             panel_letter = "", cex_panel_letter = 2.0)
{
  plot(x = x, y = y, 
       ylim = y_lim,
       main = title, ylab = y_lab, xlab = x_lab,
       cex.axis = cex_axis, cex.lab = cex_lab, cex.main = cex_title, lwd = 1, type = "n", axes = F)
  
  points(x = x, y = y, pch = 16, col = "#00000070")
  
  # Insert blank
  legend(legend = c("              ",
                    "              "),
         x = legend_position, inset = inset_blank, xjust = legend_hjust, y.intersp = 1.1, # Use inset to manually adjust position
         cex = cex_legend, bty = "o", bg = "white", box.col = NA)
  
  # Insert beta value legend
  legend(legend = bquote(bold(beta) ~ bold('=') ~ bold(.(beta_value))), text.font = 2,
         x = legend_position, inset = inset_beta, xjust = legend_hjust, # Use inset to manually adjust position
         cex = cex_legend, bty = "n", bg = "white", box.col = NA)
  
  # Insert p-value legend
  legend(legend = c("",paste0("p = ",p_value)), text.font = 2,
         x = legend_position, inset = inset_p_value, xjust = legend_hjust, # Use inset to manually adjust position
         cex = cex_legend, bty = "n", bg = "white", box.col = NA)
  
  axis(side = 1, lwd = 2, cex.axis = cex_axis)
  axis(side = 2, lwd = 2, cex.axis = cex_axis)
  abline(regression, lwd = 3, col = "red")
  
  # Add panel legend
  legend(legend = panel_letter, text.font = 2,
         x = "topright", inset = c(0.00, 0.03), xjust = 0.5,
         cex = cex_panel_letter, bty ="n")
  
}


# Extract stats for legend
beta_value <- format(round(MRM_ranks_niche_overlaps$coef[2,1], 3), nsmall = 3)
p_value <- format(MRM_ranks_niche_overlaps$coef[2,4], nsmall = 3)

# Compute a linear regression to get coefficients to draw a predict line
reg_niche_overlaps <- lm(jaccard_occupancy_niche_dist ~ jaccard_occurrence_niche_dist)
summary(reg_niche_overlaps)

# Plot
pdf(file = paste0("./graphs/Niche_similarity/Correlation_niche_overlaps_occupancy_occurrence.pdf"), height = 6, width = 7)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,6.1,4.1,2.1)) # tlbr

plot_pairwise_distances_with_MRM(title = "Compare Jaccard's niche overlaps\n Occupancy vs. Occurrence", cex_title = 1.3, 
                                 y = jaccard_occupancy_niche_dist, x = jaccard_occurrence_niche_dist, 
                                 y_lab = "Jaccard's niche overlap for Occupancy", x_lab = "Jaccard's niche overlap for Occurrence", 
                                 cex_axis = 1.3, cex_lab = 1.4, panel_letter = "",
                                 beta_value = beta_value, p_value = p_value, regression = reg_niche_overlaps)

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()


##### 8/ Compute Schoener's D in environmental space between pairs of density maps ####

library(raster)

# Load occurrence dataset
Ithomiini_dataset_with_env_data <- readRDS(file = paste0("./input_data/Occurrences/Ithomiini_dataset_with_env_data.rds"))

# Extract list of OMU
OMU_list <- unique(Ithomiini_dataset_with_env_data$OMU)[order(unique(Ithomiini_dataset_with_env_data$OMU))]

# Build template matrices 
Schoener_occurrence_niche_dissimilarity_matrix <- matrix(data = 0, nrow = length(OMU_list), ncol = length(OMU_list))
Schoener_occupancy_niche_dissimilarity_matrix <- matrix(data = 0, nrow = length(OMU_list), ncol = length(OMU_list))

### Function to apply probability threshold based on cumsum of values

apply_probability_threshold_on_raster <- function(raster, quantile = 0.95)
{
  # Extract pixels with values
  pixel_values <- raster[!is.na(raster[])]
  # Order pixels by increasing values
  ordered_pixels <- pixel_values[order(pixel_values)]
  # Compute % cumsum to reach
  cumsum_threshold <- sum(ordered_pixels, na.rm = T) * (1 - quantile)
  # Obtain probability threshold needed to reach the % cumsum threshold
  probability_threshold <- ordered_pixels[which.max(cumsum(ordered_pixels) >= cumsum_threshold)]
  # Apply threshold on density map
  raster[raster[] < probability_threshold] <- 0
  # Return thresholded raster
  return(raster)
}

### Function to compute Schoener's from two raster layers

# See isocat::schoenersD

compute_SchoenerD_from_rasters <- function(raster1, raster2)
{
  ### Transform rasters into density probability by making all pixel values sum to 1
  
  # Raster 1
  sum1 <- sum(raster1[], na.rm = T)
  raster1_prob <- raster1/sum1
  
  # Raster 2
  sum2 <- sum(raster2[], na.rm = T)
  raster2_prob <- raster2/sum2
  
  # Compute Schoener's D from probability rasters
  SchoenerD <- 1 - (0.5 * raster::cellStats(abs(raster1_prob - raster2_prob), stat = "sum"))
  
  # Return index
  return(SchoenerD)
}

### Double loop to compute all pairs

# Select quantile threshold options
quantile_value <- 0.95

for (i in 1:(length(OMU_list)-1))
{
  # i <- 1
  
  # Extract OMU i name
  OMU_i <- OMU_list[i]
  
  # Load OMU occurrence density map
  OMU_occurrence_density_map_i <- readRDS(file = paste0("./outputs/Niche_similarity/Density_maps/occurrence_density_map_",OMU_i,".rds"))

  # Load OMU occurrence density map
  OMU_occupancy_density_map_i <- readRDS(file = paste0("./outputs/Niche_similarity/Density_maps/Occupancy_density_map_",OMU_i,".rds"))

  for (j in (i+1):length(OMU_list))
  {
    # j <- 6
    
    # Extract species j name
    OMU_j <- OMU_list[j]
    
    # Load OMU occurrence density map
    OMU_occurrence_density_map_j <- readRDS(file = paste0("./outputs/Niche_similarity/Density_maps/occurrence_density_map_",OMU_j,".rds"))
    # plot(OMU_occurrence_density_map_i)
    # plot(OMU_occurrence_density_map_j)
    
    # Load OMU occurrence density map
    OMU_occupancy_density_map_j <- readRDS(file = paste0("./outputs/Niche_similarity/Density_maps/Occupancy_density_map_",OMU_j,".rds"))
    # plot(OMU_occupancy_density_map_i)
    # plot(OMU_occupancy_density_map_j)
    
    ### 8.1/ Apply threshold ###
    
    ### Remove pixels below 5% of sum of values (apply 95% probability threshold)
    
    OMU_occurrence_density_map_thresholded_i <- apply_probability_threshold_on_raster(OMU_occurrence_density_map_i, quantile = quantile_value)
    # plot(OMU_occurrence_density_map_thresholded_i)
    
    OMU_occurrence_density_map_thresholded_j <- apply_probability_threshold_on_raster(OMU_occurrence_density_map_j, quantile = quantile_value)
    # plot(OMU_occurrence_density_map_thresholded_j)
    
    OMU_occupancy_density_map_thresholded_i <- apply_probability_threshold_on_raster(OMU_occupancy_density_map_i, quantile = quantile_value)
    # plot(OMU_occupancy_density_map_thresholded_i)
    
    OMU_occupancy_density_map_thresholded_j <- apply_probability_threshold_on_raster(OMU_occupancy_density_map_j, quantile = quantile_value)
    # plot(OMU_occupancy_density_map_thresholded_j)
    
    ### 8.2/ Compute Schoener's D index for Niches based on Occurrence density ####
    
    # plot(abs(OMU_occurrence_density_map_thresholded_i - OMU_occurrence_density_map_thresholded_j))
    
    SchoenerD_occurrence <- compute_SchoenerD_from_rasters(OMU_occurrence_density_map_thresholded_i, OMU_occurrence_density_map_thresholded_j)
    
    # Transform into dissimilarity metric
    SchoenerD_dist_occurrence <- round(1 - SchoenerD_occurrence, 3)
    
    # Store into final matrix
    Schoener_occurrence_niche_dissimilarity_matrix[i,j] <- SchoenerD_dist_occurrence
    
    ### 8.3/ Compute Schoener's D index for Niches based on Occurrence density ####
    
    # plot(abs(OMU_occupancy_density_map_thresholded_i - OMU_occupancy_density_map_thresholded_j))
    
    SchoenerD_occupancy <- compute_SchoenerD_from_rasters(OMU_occupancy_density_map_thresholded_i, OMU_occupancy_density_map_thresholded_j)
  
    # Transform into dissimilarity metric
    SchoenerD_dist_occupancy <- round(1 - SchoenerD_occupancy, 3)
    
    # Store into final matrix
    Schoener_occupancy_niche_dissimilarity_matrix[i,j] <- SchoenerD_dist_occupancy
    
  }
  
  # Save temporary matrices
  saveRDS(Schoener_occurrence_niche_dissimilarity_matrix, file = paste0("./outputs/Niche_similarity/Schoener_occurrence_niche_dissimilarity_matrix_",quantile_value,".rds"))
  saveRDS(Schoener_occupancy_niche_dissimilarity_matrix, file = paste0("./outputs/Niche_similarity/Schoener_occupancy_niche_dissimilarity_matrix_",quantile_value,".rds"))
  
  cat(paste0(Sys.time(), " - Pairwise Schoener's D niche overlap computed for OMU ", i, " out of ", length(OMU_list), "\n"))
  
}

### Fill the lower-triangle
fill_lower_triangle <- function (m)
{
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
}

Schoener_occurrence_niche_dissimilarity_matrix <- fill_lower_triangle(Schoener_occurrence_niche_dissimilarity_matrix)
Schoener_occupancy_niche_dissimilarity_matrix <- fill_lower_triangle(Schoener_occupancy_niche_dissimilarity_matrix)

# Add row and col names
row.names(Schoener_occurrence_niche_dissimilarity_matrix) <- colnames(Schoener_occurrence_niche_dissimilarity_matrix) <- OMU_list
row.names(Schoener_occupancy_niche_dissimilarity_matrix) <- colnames(Schoener_occupancy_niche_dissimilarity_matrix) <- OMU_list

# Save
saveRDS(Schoener_occurrence_niche_dissimilarity_matrix, file = paste0("./outputs/Niche_similarity/Schoener_occurrence_niche_dissimilarity_matrix_",quantile_value,".rds"))
saveRDS(Schoener_occupancy_niche_dissimilarity_matrix, file = paste0("./outputs/Niche_similarity/Schoener_occupancy_niche_dissimilarity_matrix_",quantile_value,".rds"))


##### 9/ Compute Euclidean centroid distance matrix #####

### 9.1/ Compute centroid coordinates in global PCA applied in available envrionment within Ithomiini range ###

# Load available environmental data
Available_env_data <- readRDS(file = "./input_data/Env_data/Available_env_data.rds")

# Load Ithomiini occurrence dataset
Ithomiini_dataset_with_env_data <- readRDS(file = paste0("./input_data/Occurrences/Ithomiini_dataset_with_env_data.rds"))

# Provide selection of environmental variables
selected_variables <- c("Tmean", "Tvar", "Ptot", "Pvar")

# Run PCA on all study area data
PCA_Available_env_data <- ade4::dudi.pca(df = Available_env_data[, selected_variables], center = T, scale = T, scannf = F, nf = 2)
# Project coordinates of occurrences for all biological units
OMU_occurrences_PCA_coords <- ade4::suprow(x = PCA_Available_env_data, Xsup = Ithomiini_dataset_with_env_data[, selected_variables])$lisup
names(OMU_occurrences_PCA_coords) <- c("PC1", "PC2")

# Add PCA coordinates to Ithomiini dataset
Ithomiini_dataset_with_env_data <- cbind(Ithomiini_dataset_with_env_data, OMU_occurrences_PCA_coords)

# Save Ithomiini dataset
saveRDS(object = Ithomiini_dataset_with_env_data, file = paste0("./input_data/Occurrences/Ithomiini_dataset_with_env_data.rds"))

### Create dataframe for centroids data

OMU_list <- unique(Ithomiini_dataset_with_env_data$OMU)[order(unique(Ithomiini_dataset_with_env_data$OMU))]

OMU_centroids_df <- data.frame(ID_OMU = 1:length(OMU_list), OMU = OMU_list)
OMU_centroids_df <- left_join(x = OMU_centroids_df, y = Ithomiini_dataset_with_env_data, by = "OMU") %>% 
  dplyr::select(ID_OMU, OMU, Genus, Species, Mimicry) %>% 
  dplyr::distinct() %>% 
  mutate(PC1 = NA) %>% 
  mutate(PC2 = NA)
  

### Compute centroids per OMU

# Loop per OMU
for (i in 1:length(OMU_list))
{
  # i <- 2
  
  OMU <- OMU_list[i]
  
  # Extract coordinates of OMU's occurrences in PCA
  OMU_coordinates <- Ithomiini_dataset_with_env_data[Ithomiini_dataset_with_env_data$OMU == OMU, c("PC1", "PC2"), drop = F]
  
  # Compute centroids
  OMU_centroids <- apply(X = OMU_coordinates, MARGIN = 2, FUN = mean, na.rm = T)
  
  # Fill summary df
  OMU_centroids_df[i, c("PC1", "PC2")] <- OMU_centroids
}

# Save OMU's centroids coordinates in PCA
saveRDS(object = OMU_centroids_df, file = paste0("./input_data/Occurrences/OMU_centroids_df.rds"))


### 9.2/ Compute pairwise Euclidean distances

# Load df with OMU's centroids coordinates in PCA
OMU_centroids_df <- readRDS(file = paste0("./input_data/Occurrences/OMU_centroids_df.rds"))

# Compute all pairwise distances
OMU_pairwise_Euclidean_controid_distances <- as.matrix(dist(x = OMU_centroids_df[, c("PC1", "PC2")], method = "euclidean"))

# Add row and col names
row.names(OMU_pairwise_Euclidean_controid_distances) <- colnames(OMU_pairwise_Euclidean_controid_distances) <- OMU_list

# Save
saveRDS(OMU_pairwise_Euclidean_controid_distances, file = paste0("./outputs/Niche_similarity/OMU_pairwise_Euclidean_controid_distances.rds"))


### 9.3/ Plot histogram of pairwise Euclidean centroid distances

pdf(file = paste0("./graphs/Niche_similarity/Euclidean_controid_distances_histo.pdf"), height = 6, width = 7)

hist(OMU_pairwise_Euclidean_controid_distances, 
     col = "gray", xlab = "Euclidean centroid distances", 
     ylab = "Frequency",
     main = paste0("Euclidean centroid distances\nacross all pairs of OMUs"),
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)

dev.off()


##### 10/ Compare all dissimilarity indices ####

### 10.1/ Load dissimilarity matrices ####

quantile_value <-  0.95
quantile_type <- "probability"

# Load Jaccard's dissimilarity matrices
jaccard_occurrence_niche_dissimilarity_matrix <- readRDS(file = paste0("./outputs/Niche_similarity/Jaccard_occurrence_niche_dissimilarity_matrix_",quantile_value,"_",quantile_type,".rds"))
jaccard_occupancy_niche_dissimilarity_matrix <- readRDS(file = paste0("./outputs/Niche_similarity/Jaccard_occupancy_niche_dissimilarity_matrix_",quantile_value,"_",quantile_type,".rds"))

# Load Schoener's D dissimilarity matrices
Schoener_occurrence_niche_dissimilarity_matrix <- readRDS(file = paste0("./outputs/Niche_similarity/Schoener_occurrence_niche_dissimilarity_matrix_",quantile_value,".rds"))
Schoener_occupancy_niche_dissimilarity_matrix <- readRDS(file = paste0("./outputs/Niche_similarity/Schoener_occupancy_niche_dissimilarity_matrix_",quantile_value,".rds"))

# Load Euclidean centroid distances matrix
OMU_pairwise_Euclidean_controid_distances <- readRDS(file = paste0("./outputs/Niche_similarity/OMU_pairwise_Euclidean_controid_distances.rds"))

### 10.2/ Build dataframe to compare dissimilarity metrics ####

niche_dissimilarity_metrics_df <- as.data.frame(cbind(as.dist(OMU_pairwise_Euclidean_controid_distances),
                                                as.dist(jaccard_occurrence_niche_dissimilarity_matrix),
                                                as.dist(jaccard_occupancy_niche_dissimilarity_matrix),
                                                as.dist(Schoener_occurrence_niche_dissimilarity_matrix),
                                                as.dist(Schoener_occupancy_niche_dissimilarity_matrix)))

names(niche_dissimilarity_metrics_df) <- c("Euclidean_distances_occurrence_centroids",
                                           "Jaccard_occurrence_niche_overlap",
                                           "Jaccard_occupancy_niche_overlap",
                                           "Schoener_occurrence_niche_overlap",
                                           "Schoener_occupancy_niche_overlap")

View(niche_dissimilarity_metrics_df)

# Save
saveRDS(niche_dissimilarity_metrics_df, file = "./outputs/Niche_similarity/Correlations/niche_dissimilarity_metrics_df.rds")

### 10.3/ Make correlation plots ####

# Load
niche_dissimilarity_metrics_df <- readRDS(file = "./outputs/Niche_similarity/Correlations/niche_dissimilarity_metrics_df.rds")

cor_mat <- cor(niche_dissimilarity_metrics_df, method = "spearman")
# cor_mat <- cor(niche_dissimilarity_metrics_df, method = "pearson")
col_ramp <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200)

# Plot
pdf(file = paste0("./graphs/Niche_similarity/Correlation_niche_dissimilarity_metrics.pdf"), height = 7, width = 7)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(6.1,6.5,4.1,2.1)) # tlbr

corrplot::corrplot.mixed(corr = cor_mat,
                         lower.col = col_ramp,
                         upper.col = col_ramp,
                         tl.pos = "lt",
                         tl.col = "black")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()


### 10.4/ Run MRM tests for correlation and LM regressions ####

library(ecodist)
?ecodist::MRM

# Load home-made function which can provide summary of permutation outputs
source("./functions/MRM_verbose.R")

# Prepare distances data = linearized in vectors
Euclidean_distances_occurrence_centroids_dist <- as.vector(niche_dissimilarity_metrics_df$Euclidean_distances_occurrence_centroids)
Jaccard_occurrence_niche_overlap_dist <- as.vector(niche_dissimilarity_metrics_df$Jaccard_occurrence_niche_overlap)
Schoener_occurrence_niche_overlap_dist <- as.vector(niche_dissimilarity_metrics_df$Schoener_occurrence_niche_overlap)


### 10.4.1/ Jaccard's occurrence niche overlap ~ Euclidean centroid distances

## MRM: Jaccard overlap ~ Euclidean distances
MRM_ranks_Jaccard_Euclidean_dist <- MRM_verbose(formula = Jaccard_occurrence_niche_overlap_dist ~ Euclidean_distances_occurrence_centroids_dist,
                                        nperm = 1000,
                                        method = "linear", # To use Linear Model with normal distribution assumed for the conditional response variable
                                        mrank = TRUE,  # Transform data in ranks
                                        perm_output = TRUE)  # Provide permutation outputs

print(MRM_ranks_Jaccard_Euclidean_dist)
print(MRM_ranks_Jaccard_Euclidean_dist$coef)
print(MRM_ranks_Jaccard_Euclidean_dist$r.squared)
print(MRM_ranks_Jaccard_Euclidean_dist$F.test)

# Save MRM
saveRDS(object = MRM_ranks_Jaccard_Euclidean_dist, file = "./outputs/Niche_similarity/Correlations/MRM_ranks_Jaccard_Euclidean_dist.rds")

## LM: Jaccard overlap ~ Euclidean distances

# Compute a linear regression to get coefficients to draw a predict line
LM_Jaccard_Euclidean_dist <- lm(Jaccard_occurrence_niche_overlap_dist ~ Euclidean_distances_occurrence_centroids_dist)
summary(LM_Jaccard_Euclidean_dist)

# Save LM
saveRDS(object = LM_Jaccard_Euclidean_dist, file = "./outputs/Niche_similarity/Correlations/LM_Jaccard_Euclidean_dist.rds")


### 10.4.2/ Schoener's D occurrence niche overlap ~ Euclidean centroid distances

## MRM: Schoener's D overlap ~ Euclidean distances
MRM_ranks_Schoener_Euclidean_dist <- MRM_verbose(formula = Schoener_occurrence_niche_overlap_dist ~ Euclidean_distances_occurrence_centroids_dist,
                                                 nperm = 1000,
                                                 method = "linear", # To use Linear Model with normal distribution assumed for the conditional response variable
                                                 mrank = TRUE,  # Transform data in ranks
                                                 perm_output = TRUE)  # Provide permutation outputs

print(MRM_ranks_Schoener_Euclidean_dist)
print(MRM_ranks_Schoener_Euclidean_dist$coef)
print(MRM_ranks_Schoener_Euclidean_dist$r.squared)
print(MRM_ranks_Schoener_Euclidean_dist$F.test)

# Save MRM
saveRDS(object = MRM_ranks_Schoener_Euclidean_dist, file = "./outputs/Niche_similarity/Correlations/MRM_ranks_Schoener_Euclidean_dist.rds")

## LM: Schoener's D overlap ~ Euclidean distances

# Compute a linear regression to get coefficients to draw a predict line
LM_Schoener_Euclidean_dist <- lm(Schoener_occurrence_niche_overlap_dist ~ Euclidean_distances_occurrence_centroids_dist)
summary(LM_Schoener_Euclidean_dist)

# Save LM
saveRDS(object = LM_Schoener_Euclidean_dist, file = "./outputs/Niche_similarity/Correlations/LM_Schoener_Euclidean_dist.rds")


### 10.4.3/ Jaccard's occurrence niche overlap ~ Schoener's D occurrence niche overlap

## MRM: Jaccard overlap ~ Schoener's D overlap
MRM_ranks_Jaccard_Schoener <- MRM_verbose(formula = Jaccard_occurrence_niche_overlap_dist ~ Schoener_occurrence_niche_overlap_dist,
                                                nperm = 1000,
                                                method = "linear", # To use Linear Model with normal distribution assumed for the conditional response variable
                                                mrank = TRUE,  # Transform data in ranks
                                                perm_output = TRUE)  # Provide permutation outputs

print(MRM_ranks_Jaccard_Schoener)
print(MRM_ranks_Jaccard_Schoener$coef)
print(MRM_ranks_Jaccard_Schoener$r.squared)
print(MRM_ranks_Jaccard_Schoener$F.test)

# Save MRM
saveRDS(object = MRM_ranks_Jaccard_Schoener, file = "./outputs/Niche_similarity/Correlations/MRM_ranks_Jaccard_Schoener.rds")

## LM: Jaccard overlap ~ Schoener's D overlap

# Compute a linear regression to get coefficients to draw a predict line
LM_Jaccard_Schoener <- lm(Jaccard_occurrence_niche_overlap_dist ~ Schoener_occurrence_niche_overlap_dist)
summary(LM_Jaccard_Schoener)

# Save LM
saveRDS(object = LM_Jaccard_Schoener, file = "./outputs/Niche_similarity/Correlations/LM_Jaccard_Schoener.rds")


### 10.5/ Plot pairwise distances scatterplot ####

### 10.5.1/ Jaccard's occurrence niche overlap ~ Euclidean centroid distances

# Load MRM output
MRM_ranks_Jaccard_Euclidean_dist <- readRDS(file = "./outputs/Niche_similarity/Correlations/MRM_ranks_Jaccard_Euclidean_dist.rds")

# Load LM output
LM_Jaccard_Euclidean_dist <- readRDS(file = "./outputs/Niche_similarity/Correlations/LM_Jaccard_Euclidean_dist.rds")

# Extract stats for legend
beta_value <- format(round(MRM_ranks_Jaccard_Euclidean_dist$coef[2,1], 3), nsmall = 3)
p_value <- format(MRM_ranks_Jaccard_Euclidean_dist$coef[2,4], nsmall = 3)

# Plot
pdf(file = paste0("./graphs/Niche_similarity/Correlation_plot_Jaccard_Euclidean_dist.pdf"), height = 6, width = 7)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,6.1,4.1,2.1)) # tlbr

plot_pairwise_distances_with_MRM(title = "Jaccard's niche non-overlaps\n vs. Euclidean centroid distances", cex_title = 1.3, 
                                 y = Jaccard_occurrence_niche_overlap_dist, x = Euclidean_distances_occurrence_centroids_dist, 
                                 y_lab = "Jaccard's niche non-overlap for Occurrence", x_lab = "Euclidean centroid distances",
                                 legend_position = "bottomright", legend_hjust = 0,
                                 inset_blank = c(-0.03, 0.06), inset_beta = c(0.00,0.14), inset_p_value = c(0,0.07),
                                 cex_axis = 1.3, cex_lab = 1.4, panel_letter = "",
                                 beta_value = beta_value, p_value = p_value, regression = LM_Jaccard_Euclidean_dist)

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()



### 10.5.2/ Schoener's D occurrence niche overlap ~ Euclidean centroid distances

# Load MRM output
MRM_ranks_Schoener_Euclidean_dist <- readRDS(file = "./outputs/Niche_similarity/Correlations/MRM_ranks_Schoener_Euclidean_dist.rds")

# Load LM output
LM_Schoener_Euclidean_dist <- readRDS(file = "./outputs/Niche_similarity/Correlations/LM_Schoener_Euclidean_dist.rds")

# Extract stats for legend
beta_value <- format(round(MRM_ranks_Schoener_Euclidean_dist$coef[2,1], 3), nsmall = 3)
p_value <- format(MRM_ranks_Schoener_Euclidean_dist$coef[2,4], nsmall = 3)

# Plot
pdf(file = paste0("./graphs/Niche_similarity/Correlation_plot_Schoener_Euclidean_dist.pdf"), height = 6, width = 7)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,6.1,4.1,2.1)) # tlbr

plot_pairwise_distances_with_MRM(title = "Schoener's D niche dissimilarity\n vs. Euclidean centroid distances", cex_title = 1.3, 
                                 y = Schoener_occurrence_niche_overlap_dist, x = Euclidean_distances_occurrence_centroids_dist, 
                                 y_lab = "Schoener's D niche dissimilarity for Occurrence", x_lab = "Euclidean centroid distances",
                                 legend_position = "bottomright", legend_hjust = 0,
                                 inset_blank = c(-0.03, 0.06), inset_beta = c(0.00,0.14), inset_p_value = c(0,0.07),
                                 cex_axis = 1.3, cex_lab = 1.4, panel_letter = "",
                                 beta_value = beta_value, p_value = p_value, regression = LM_Schoener_Euclidean_dist)

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()


### 10.5.3/ Jaccard's occurrence niche overlap ~ Schoener's D occurrence niche overlap

# Load MRM output
MRM_ranks_Jaccard_Schoener <- readRDS(file = "./outputs/Niche_similarity/Correlations/MRM_ranks_Jaccard_Schoener.rds")

# Load LM output
LM_Jaccard_Schoener <- readRDS(file = "./outputs/Niche_similarity/Correlations/LM_Jaccard_Schoener.rds")

# Extract stats for legend
beta_value <- format(round(MRM_ranks_Jaccard_Schoener$coef[2,1], 3), nsmall = 3)
p_value <- format(MRM_ranks_Jaccard_Schoener$coef[2,4], nsmall = 3)

# Plot
pdf(file = paste0("./graphs/Niche_similarity/Correlation_plot_Jaccard_Schoener.pdf"), height = 6, width = 7)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,6.1,4.1,2.1)) # tlbr

plot_pairwise_distances_with_MRM(title = "Jaccard's niche non-overlaps\n vs. Schoener's D niche dissimilarity", cex_title = 1.3, 
                                 y = Jaccard_occurrence_niche_overlap_dist, x = Schoener_occurrence_niche_overlap_dist, 
                                 y_lab = "Jaccard's niche non-overlap for Occurrence", x_lab = "Schoener's D niche dissimilarity for Occurrence",
                                 legend_position = "bottomright", legend_hjust = 0,
                                 inset_blank = c(-0.03, 0.06), inset_beta = c(0.00,0.14), inset_p_value = c(0,0.07),
                                 cex_axis = 1.3, cex_lab = 1.4, panel_letter = "",
                                 beta_value = beta_value, p_value = p_value, regression = LM_Jaccard_Schoener)

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()


##### 11/ Test for niche similarity based on Jaccard's niche overlap ####

### Permutation test: 
 # Null hypothesis = No relationship between pairwise niche overlap and mimicry membership
 # Null model = permutation of niche overlap across OMUs
 # Statistic = mean niche overlap between comimics and within each mimicry rings
 # Output = null distribution (plot) + stats, quantiles and p_value (summary table)
   # Global test = for all comimics
   # Per mimicry ring

### Load pairwise Jaccard's niche non-overlap based on Occurrence density
jaccard_occurrence_niche_dissimilarity_matrix <- readRDS(file = paste0("./outputs/Niche_similarity/Jaccard_occurrence_niche_dissimilarity_matrix_",quantile_value,"_",quantile_type,".rds"))

# Use Jaccard's niche overlap (similarity) metric because it is more intuitive
jaccard_occurrence_niche_similarity_matrix <- (1 - jaccard_occurrence_niche_dissimilarity_matrix)

# Load df with OMU's centroids coordinates in PCA and mimicry membership
OMU_centroids_df <- readRDS(file = paste0("./input_data/Occurrences/OMU_centroids_df.rds"))

### 11.1/ Compute co-mimicry matrix for OMUs ####

OMU_comimicry_matrix <- matrix(nrow = nrow(OMU_centroids_df), ncol = nrow(OMU_centroids_df), data = 0)
for (i in 1:nrow(OMU_centroids_df))
{
  ring_1 <- as.character(OMU_centroids_df$Mimicry)[i]
  
  for (j in 1:nrow(OMU_centroids_df)) 
  {
    ring_2 <- as.character(OMU_centroids_df$Mimicry)[j]
    
    if (ring_1 == ring_2) 
    {
      OMU_comimicry_matrix[i,j] <- 1
    }
  }
  if (i %% 10 == 0) {print(i)}
}
row.names(OMU_comimicry_matrix) <- colnames(OMU_comimicry_matrix) <- OMU_centroids_df$OMU
saveRDS(OMU_comimicry_matrix, file = paste0("./outputs/Niche_similarity/OMU_comimicry_matrix.rds"))

load(file = paste0("./outputs/Niche_evolution/OMU_comimicry_matrix.RData"))


### 11.2/ Compute the global observed Jaccard's niche overlap weighted by mimicry similarity ####

mean(as.dist(jaccard_occurrence_niche_similarity_matrix)) # Global Jaccard's niche overlap for all pairs of OMUs = 14.8 %
Global_Jaccard_overlap_obs <- weighted.mean(x = as.dist(jaccard_occurrence_niche_similarity_matrix), w = as.dist(OMU_comimicry_matrix)) # Global Jaccard's niche overlap only for pairs of comimics = 22.8%

saveRDS(Global_Jaccard_overlap_obs, file = paste0("./outputs/Niche_similarity/Similarity_tests/Global_Jaccard_overlap_obs.rds"))

### 11.3/ Compute the observed Jaccard's niche overlap per mimicry ring ####

mimicry_rings_list <- as.character(unique(OMU_centroids_df$Mimicry))
mimicry_rings_list <- mimicry_rings_list[order(mimicry_rings_list)]

Jaccard_overlap_per_mimic_obs <- rep(NA, length(mimicry_rings_list)) # Generate empty vector to store Jaccard's niche overlap obs

for (i in 1:length(mimicry_rings_list)) # Per mimic ring
{
  # i <-  3
  
  ring <- mimicry_rings_list[i]
  
  # Get indices of rows/columns associated of OMUs for this ring
  ring_indices <- which(OMU_centroids_df$Mimicry == ring)
  
  if(length(ring_indices) != 1) # Need computation only if there are more than 1 OMU, otherwise, no Jaccard_overlap possible
  {
    # Extract only the pairwise Jaccard's niche overlaps for this ring
    pairwise_Jaccard_overlap_ring <- as.dist(jaccard_occurrence_niche_similarity_matrix[ring_indices, ring_indices])
    
    Jaccard_overlap_per_mimic_obs[i] <- mean(pairwise_Jaccard_overlap_ring)
  }
  print(i)
}
names(Jaccard_overlap_per_mimic_obs) <- mimicry_rings_list
Jaccard_overlap_per_mimic_obs

# Save the observed Jaccard's niche overlap per mimicry ring
saveRDS(Jaccard_overlap_per_mimic_obs, file = paste0("./outputs/Niche_similarity/Similarity_tests/Jaccard_overlap_per_mimic_obs.rds"))


### 11.4/ Simulate Jaccard's niche overlap under null hypothesis by permutation ####

## Start the loop
Global_Jaccard_overlap_null <-  NA
Jaccard_overlap_per_mimic_null <- data.frame(matrix(nrow = 999, ncol = length(mimicry_rings_list))) # Generate df to store results
for (l in 1:999) # 999 permutations
{
  # l <- 1
  
  ## Shuffle mimicry ring among units in the comimicry matrix
  shuffle.indices <- sample(x = 1:nrow(OMU_centroids_df), size = nrow(OMU_centroids_df), replace = F)
  OMU_comimicry_matrix_shuffle <- OMU_comimicry_matrix[shuffle.indices, shuffle.indices]
  
  # Compute the global Jaccard's niche overlap
  Global_Jaccard_overlap_simul <- weighted.mean(x = as.dist(jaccard_occurrence_niche_similarity_matrix), w = as.dist(OMU_comimicry_matrix_shuffle))
  Global_Jaccard_overlap_null[l] <- Global_Jaccard_overlap_simul
  
  saveRDS(Global_Jaccard_overlap_null, file = paste0("./outputs/Niche_similarity/Similarity_tests/Global_Jaccard_overlap_null.rds"))
  
  
  # Compute the Jaccard's niche overlap per ring
  Jaccard_overlap_per_mimic_simul <- rep(NA, length(mimicry_rings_list)) # Generate empty vector to store Jaccard_overlap obs
  
  for (i in 1:length(mimicry_rings_list)) # Per mimic ring
  {
    # i <-  3
    
    ring <- mimicry_rings_list[i]
    
    # Get indices of rows/columns associated of OMUs for this ring in the real dataset
    ring_indices <- which(OMU_centroids_df$Mimicry == ring)
    
    if (length(ring_indices) != 1) # Need computation only if there are more than 1 OMU, otherwise, no Jaccard's niche overlap possible
    {
      shuffle_ring_indices <- shuffle.indices[ring_indices] # Get the new indices of the OMU in the shuffled comimicry matrix
      
      # Extract only the pairwise Jaccard's niche overlaps for this ring
      pairwise_dist_ring <- as.dist(jaccard_occurrence_niche_similarity_matrix[shuffle_ring_indices, shuffle_ring_indices])
      
      Jaccard_overlap_per_mimic_simul[i] <- mean(pairwise_dist_ring)
    }
  }
  
  Jaccard_overlap_per_mimic_null[l, ] <- Jaccard_overlap_per_mimic_simul
  
  # Temporary save of the simulated indices
  saveRDS(Jaccard_overlap_per_mimic_null, file = paste0("./outputs/Niche_similarity/Similarity_tests/Jaccard_overlap_per_mimic_null.rds")) 
  
  # Print progress
  if (l %% 10 == 0) { cat(paste0(Sys.time(), " - Simul n°", l, " out of 999\n")) }
}
names(Jaccard_overlap_per_mimic_null) <- mimicry_rings_list
saveRDS(Jaccard_overlap_per_mimic_null, file = paste0("./outputs/Niche_similarity/Similarity_tests/Jaccard_overlap_per_mimic_null.rds")) 


### 11.5/ Plot null distri for global comimics Jaccard's niche overlap ####

Global_Jaccard_overlap_obs <- readRDS(file = paste0("./outputs/Niche_similarity/Similarity_tests/Global_Jaccard_overlap_obs.rds"))
Global_Jaccard_overlap_null <- readRDS(file = paste0("./outputs/Niche_similarity/Similarity_tests/Global_Jaccard_overlap_null.rds"))

Null_distri_Jaccard_overlap <- c(Global_Jaccard_overlap_null, Global_Jaccard_overlap_obs)

pdf(file = "./graphs/Niche_similarity/Global_Jaccard_overlap_null.pdf", height = 6, width = 7)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

hist(x = Null_distri_Jaccard_overlap, breaks = 40, freq = TRUE, col = "gray",
     # xlim = c(34.5, 38.0),
     # ylim = c(0, 200),
     main = "Niches are more similar between co-mimetic OMUs\nbased on Jaccard's niche overlaps",
     # main = "",
     ylab = "Frequency",
     xlab = "Mean pairwise Jaccard's niche overlap",
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)

arrows(Global_Jaccard_overlap_obs + 0.0007, 95, Global_Jaccard_overlap_obs + 0.0007, 5, length = 0.1, lwd = 2)  # Draw arrow above mean Jaccard's overlap obs
abline(v = mean(Null_distri_Jaccard_overlap), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(Null_distri_Jaccard_overlap, 0.95), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("Mean = ", round(mean(Null_distri_Jaccard_overlap),3)), 
                  paste0("CI 95% = ", round(quantile(Null_distri_Jaccard_overlap, 0.95),3))), 
       x = "topright", inset = c(0.02, 0.08), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty = "n")
legend(legend = c(paste0("J obs = ", round(Global_Jaccard_overlap_obs, 2)),
                  paste0("p = 0.001")),
       x = "right", inset = c(0.00, 0.00),
       cex = 1.2, bty ="n", xjust = 0)

legend(legend = as.expression(bquote(bold(""))), 
       x = "topright", inset = c(0.05, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()



### 11.6/ Plot null distri Jaccard's niche overlap per mimic ring ####

Jaccard_overlap_per_mimic_obs <- readRDS(file = paste0("./outputs/Niche_similarity/Similarity_tests/Jaccard_overlap_per_mimic_obs.rds"))
Jaccard_overlap_per_mimic_null <- readRDS(file = paste0("./outputs/Niche_similarity/Similarity_tests/Jaccard_overlap_per_mimic_null.rds")) 

Jaccard_overlap_per_mimic_obs
Jaccard_overlap_per_mimic_null


# Create summary table for phylogenetic signal of mimicry rings
Jaccard_overlap_ring_summary_table <- as.data.frame(matrix(ncol = 9, nrow = 44, data = NA))
names(Jaccard_overlap_ring_summary_table) <- c("ring", "N_units", "N_pairs", "Jaccard_overlap_obs", "mean_Jaccard_overlap", "Jaccard_overlap_5", "Jaccard_overlap_95", "p_value", "pattern")

for (i in 1:length(mimicry_rings_list)) 
{
  # i <- 2
  
  Jaccard_overlap_ring_summary_table$ring[i] <- ring <- mimicry_rings_list[i]
  Jaccard_overlap_ring_summary_table$N_units[i] <- N_units <- sum(OMU_centroids_df$Mimicry == ring)
  Jaccard_overlap_ring_summary_table$N_pairs[i] <- N_pairs <- N_units*(N_units-1)/2
  
  if (is.na(Jaccard_overlap_per_mimic_obs[i]))  # Case for ring with only one OMUs. No pairs. No Jaccard's niche overlap.
  {
    pdf(file = paste0("./graphs/Niche_similarity/Per_ring/Jaccard_overlap_null_",ring,".pdf"), height = 6, width = 7)
    plot(1:100,1:100, type = "n", xlab = "Mean pairwise Jaccard's niche overlap", ylab = "Frequency",
         main = paste0("Distribution of Mean pairwise Jaccard's niche overlap\n of ", ring, " OMUs \n under null Hypothesis"))
    text(x = 50, y = 50, labels = "Only one OMU for this mimicry ring \n No pair available for index computation")
    dev.off()
    
  } else {
    
    mean_val <- round(mean(Jaccard_overlap_per_mimic_null[,i]),3)
    
    if (Jaccard_overlap_per_mimic_obs[i] < mean_val)  # Overlap observed is lower than null statistics. Signal for niche dissimilarity
    {
      pattern <- "dissimilarity"
      p_value <- round(ecdf(x = Jaccard_overlap_per_mimic_null[,i])(Jaccard_overlap_per_mimic_obs[i]),3)
    } else {   # Overlap observed is higher than null statistics. Signal for niche similarity
      pattern <- "similarity"
      p_value <- round(1 - ecdf(x = Jaccard_overlap_per_mimic_null[,i])(Jaccard_overlap_per_mimic_obs[i]),3)
    }
    
    Jaccard_overlap_ring_summary_table$Jaccard_overlap_obs[i] <- round(Jaccard_overlap_per_mimic_obs[i],3)
    Jaccard_overlap_ring_summary_table$mean_Jaccard_overlap[i] <- mean_val
    Jaccard_overlap_ring_summary_table$Jaccard_overlap_5[i] <- round(quantile(Jaccard_overlap_per_mimic_null[,i], 0.05, na.rm = T),3)
    Jaccard_overlap_ring_summary_table$Jaccard_overlap_95[i] <- round(quantile(Jaccard_overlap_per_mimic_null[,i], 0.95, na.rm = T),3)
    Jaccard_overlap_ring_summary_table$p_value[i] <- p_value
    Jaccard_overlap_ring_summary_table$pattern[i] <- pattern
    
    histo.save <- hist(Jaccard_overlap_per_mimic_null[,i],
                       breaks = seq(from = floor(min(min(Jaccard_overlap_per_mimic_null[,i], na.rm = T), Jaccard_overlap_per_mimic_obs[i])*100)/100 - 0.01, to = ceiling(max(max(Jaccard_overlap_per_mimic_null[,i], na.rm = T), Jaccard_overlap_per_mimic_obs[i])*100)/100 + 0.01, by = 0.005),
                       plot = F)
    
    pdf(file = paste0("./graphs/Niche_similarity/Per_ring/Jaccard_overlap_null_",ring,".pdf"), height = 6, width = 7)
    
    hist(Jaccard_overlap_per_mimic_null[,i], 
         breaks = seq(from = floor(min(min(Jaccard_overlap_per_mimic_null[,i], na.rm = T), Jaccard_overlap_per_mimic_obs[i])*100)/100 - 0.01, to = ceiling(max(max(Jaccard_overlap_per_mimic_null[,i], na.rm = T), Jaccard_overlap_per_mimic_obs[i])*100)/100 + 0.01, by = 0.005),
         col = "gray", xlab = "Mean pairwise Jaccard's niche overlap", 
         ylab = "Frequency",
         main = paste0("Distribution of the Mean Jaccard's niche overlaps\n of ",ring," OMUs \n under the null Hypothesis"),
         cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.2, lwd = 2)
    arrows(Jaccard_overlap_per_mimic_obs[i], max(histo.save$counts)/3, Jaccard_overlap_per_mimic_obs[i], max(histo.save$counts)/30, length = 0.1, lwd = 2)
    abline(v = mean(Jaccard_overlap_per_mimic_null[,i]), lty = 2, lwd = 2)
    abline(v = quantile(Jaccard_overlap_per_mimic_null[,i], 0.025, na.rm = T), lty = 2, lwd = 2, col = "red")
    abline(v = quantile(Jaccard_overlap_per_mimic_null[,i], 0.975, na.rm = T), lty = 2, lwd = 2, col = "red")
    legend(legend = c(paste0("N units = ", N_units), 
                      paste0("N pairs = ", N_pairs)),
           x = "topleft", cex = 1, bty ="n")
    
    legend(legend = c(paste0("Mean = ", mean_val), 
                      paste0("Q 5% = ", round(quantile(Jaccard_overlap_per_mimic_null[,i], 0.05, na.rm = T),3)),
                      paste0("Q 95% = ", round(quantile(Jaccard_overlap_per_mimic_null[,i], 0.95, na.rm = T),3))), 
           x = "topright", inset = c(0.00, 0.00),
           lty = 2 , lwd = 2, col = c("black", "red"),
           cex = 1, bty ="n")
    
    legend(legend = c(paste0("J obs = ", round(Jaccard_overlap_per_mimic_obs[i], 3)),
                      paste0("p = ", p_value)),
           x = "topright", inset = c(0.00, 0.50),
           cex = 1, bty ="n")
    
    dev.off()
  }
  saveRDS(Jaccard_overlap_ring_summary_table, file = paste0("./outputs/Niche_similarity/Similarity_tests/Jaccard_overlap_ring_summary_table.rds"))
  
  cat(paste0("N° ",i, "/",length(mimicry_rings_list)," - ",ring, " - Done \n"))
}

View(Jaccard_overlap_ring_summary_table)
write.csv2(Jaccard_overlap_ring_summary_table, file = "./tables/Jaccard_overlap_ring_summary_table.csv")




