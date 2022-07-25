
##### Script to investigate the effect of the spatial resolution on Bray-Curtis measurements #####

# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

### Load libraries
library(raster)
library(GGally)

### Load raster Stack of Ithomiini OMUs

# Initial resolution = 0.25°
OMU_proba_stack_0.25 <- readRDS(file = "./input_data/SDM_stacks/All_OMU_stack_Jaccard.80.rds")

##### 2/ Aggregate data to higher spatial scale = lower resolution ####

# Function to compute probability of presence with multiple OMUs 
aggreg_prob = function(x, na.rm)
{ 
  # Case with all NA
  if (all(is.na(x)))
  { 
    y <- NA 
  } else {
    # Case with some NA but not all
    x <- x[!is.na(x)]
    
    # General case
    y <- 1-prod(1-x) # Probability of presence of species = probability of presence of at least one OMU = opposite of probability of absence of all OMU
  }

  # Output
  return(y) 
}

# Aggregate by increasing factors
OMU_proba_stack_0.5 <- aggregate(x = OMU_proba_stack_0.25, fact = 2, fun = aggreg_prob)
OMU_proba_stack_1 <- aggregate(x = OMU_proba_stack_0.25, fact = 4, fun = aggreg_prob)
OMU_proba_stack_2 <- aggregate(x = OMU_proba_stack_0.25, fact = 8, fun = aggreg_prob)

plot(OMU_proba_stack_0.25)
plot(OMU_proba_stack_0.5)
plot(OMU_proba_stack_1)
plot(OMU_proba_stack_2)

# Save
saveRDS(object = OMU_proba_stack_0.5, file = "./input_data/SDM_stacks/OMU_proba_stack_0.5.rds")
saveRDS(object = OMU_proba_stack_1, file = "./input_data/SDM_stacks/OMU_proba_stack_1.rds")
saveRDS(object = OMU_proba_stack_2, file = "./input_data/SDM_stacks/OMU_proba_stack_2.rds")

# Load
OMU_proba_stack_0.25 <- readRDS(file = "./input_data/SDM_stacks/All_OMU_stack_Jaccard.80.rds")
OMU_proba_stack_0.5 <- readRDS(file = "./input_data/SDM_stacks/OMU_proba_stack_0.5.rds")
OMU_proba_stack_1 <- readRDS(file = "./input_data/SDM_stacks/OMU_proba_stack_1.rds")
OMU_proba_stack_2 <- readRDS(file = "./input_data/SDM_stacks/OMU_proba_stack_2.rds")

##### 3/ Compute pairwise Bray-Curtis values ####

### 3.1/ Generate matrix communities * OMU for probabilities of presence ####

# Function to build community matrix from raster Stack
build_community_matrix <- function(raster_stack, filter_NA = T, filter_null = T)
{
  raster_brick <- raster_stack*1
  community_matrix <- NA
  for (i in 1:raster::nlayers(raster_brick)) 
  {
    if (i==1){
      community_matrix <- raster_brick[[i]]@data@values
    }else{
      unit_values <- raster_brick[[i]]@data@values
      community_matrix <- cbind(community_matrix, unit_values)
    }
    if (i %% 100 == 0) {cat(paste0("\nUnits n°", i, " out of ", nlayers(raster_brick)))}
  }
  colnames(community_matrix) <- names(raster_brick)
  
  # Remove communities with NA if requested
  if (filter_NA)
  {
    # Filter out communities with NA
    real_com_index <- NA
    for (i in 1:nrow(community_matrix)) 
    {
      com_row <- community_matrix[i,]
      real_com_index[i] <- (!any(is.na(com_row)))&(!any(is.nan(com_row)))
    }
    filtered_community_matrix <- community_matrix[real_com_index, ]
    
    # Print results of filtering
    cat(paste0("\n\n", nrow(community_matrix), " communities before filtering of NA"))
    cat(paste0("\n", nrow(filtered_community_matrix), " communities after filtering of NA"))

    community_matrix <- filtered_community_matrix
  }
  
  # Remove communities with no records if requested
  if (filter_null)
  {
    # Filter out communities with no records (less than one predicted)
    empty_com_index <- apply(X = community_matrix, MARGIN = 1, FUN = sum) < 1
    filtered_community_matrix <- community_matrix[!empty_com_index, ]

    # Print results of filtering
    cat(paste0("\n\n", nrow(community_matrix), " communities before filtering of empty communities"))
    cat(paste0("\n", nrow(filtered_community_matrix), " communities after filtering of empty communities"))
    
    community_matrix <- filtered_community_matrix
  }
  
  return(community_matrix)

}

# Build community matrices for all spatial scales
OMU_community_matrix_0.25 <- build_community_matrix(OMU_proba_stack_0.25)
OMU_community_matrix_0.5 <- build_community_matrix(OMU_proba_stack_0.5)
OMU_community_matrix_1 <- build_community_matrix(OMU_proba_stack_1)
OMU_community_matrix_2 <- build_community_matrix(OMU_proba_stack_2)

# Save
saveRDS(OMU_community_matrix_0.25, file = "./outputs/Community_Structure/Spatial_scale/OMU_community_matrix_0.25.rds")
saveRDS(OMU_community_matrix_0.5, file = "./outputs/Community_Structure/Spatial_scale/OMU_community_matrix_0.5.rds")
saveRDS(OMU_community_matrix_1, file = "./outputs/Community_Structure/Spatial_scale/OMU_community_matrix_1.rds")
saveRDS(OMU_community_matrix_2, file = "./outputs/Community_Structure/Spatial_scale/OMU_community_matrix_2.rds")

# Load 
OMU_community_matrix_0.25 <- readRDS(file = "./outputs/Community_Structure/Spatial_scale/OMU_community_matrix_0.25.rds")
OMU_community_matrix_0.5 <- readRDS(file = "./outputs/Community_Structure/Spatial_scale/OMU_community_matrix_0.5.rds")
OMU_community_matrix_1 <- readRDS(file = "./outputs/Community_Structure/Spatial_scale/OMU_community_matrix_1.rds")
OMU_community_matrix_2 <- readRDS(file = "./outputs/Community_Structure/Spatial_scale/OMU_community_matrix_2.rds")

### 3/ Compute Bray-Curtis indices ####

### 3.1/ Compute Bray-Curtis index for all pairs of units ####

library(vegan)

dim(OMU_community_matrix_0.25)

# Compute dissimilarities between rows = between communities
OMU_BC_dist_0.25 <- vegdist(x = t(OMU_community_matrix_0.25), method = "bray")
OMU_BC_dist_0.5 <- vegdist(x = t(OMU_community_matrix_0.5), method = "bray")
OMU_BC_dist_1 <- vegdist(x = t(OMU_community_matrix_1), method = "bray")
OMU_BC_dist_2 <- vegdist(x = t(OMU_community_matrix_2), method = "bray")

# Save
saveRDS(OMU_BC_dist_0.25, file = "./outputs/Community_Structure/Spatial_scale/OMU_BC_dist_0.25.rds")
saveRDS(OMU_BC_dist_0.5, file = "./outputs/Community_Structure/Spatial_scale/OMU_BC_dist_0.5.rds")
saveRDS(OMU_BC_dist_1, file = "./outputs/Community_Structure/Spatial_scale/OMU_BC_dist_1.rds")
saveRDS(OMU_BC_dist_2, file = "./outputs/Community_Structure/Spatial_scale/OMU_BC_dist_2.rds")

# Load
OMU_BC_dist_0.25 <- readRDS(file = "./outputs/Community_Structure/Spatial_scale/OMU_BC_dist_0.25.rds")
OMU_BC_dist_0.5 <- saveRDS(file = "./outputs/Community_Structure/Spatial_scale/OMU_BC_dist_0.5.rds")
OMU_BC_dist_1 <- saveRDS(file = "./outputs/Community_Structure/Spatial_scale/OMU_BC_dist_1.rds")
OMU_BC_dist_2 <- saveRDS(file = "./outputs/Community_Structure/Spatial_scale/OMU_BC_dist_2.rds")


### Plot correlogram and paired scatterplot

# Build df
BC_spatial_scales_df <- data.frame(quarter_degree = as.vector(OMU_BC_dist_0.25),
                                   half_degree = as.vector(OMU_BC_dist_0.5),
                                   one_degree = as.vector(OMU_BC_dist_1),
                                   two_degrees = as.vector(OMU_BC_dist_2))

# Save

subsample <- sample(x = 1:nrow(BC_spatial_scales_df), size = 1000)
BC_spatial_scales_df_subsample <- BC_spatial_scales_df[subsample, ]

# Plot correlogram and scatterplots

?GGally::ggpairs()

GGally::ggpairs(BC_spatial_scales_df_subsample,
                title = "Correlation between BC values for several spatial resolution",
                diag = list(continuous = "barDiag"))

library(corrgram)
?corrgram::corrgram

pdf(file = "./graphs/Community_Structure/Spatial_scales_BC.pdf", height = 6, width = 7)
corrgram::corrgram(x = BC_spatial_scales_df_subsample, type = "data", order = F,
                   labels = c("0.25°", "0.5°", "1°", "2°"),
                   lower.panel = panel.pts,
                   upper.panel = panel.cor,
                   diag.panel = NULL,
                   cor.method = "spearman",
                   col.regions = colorRampPalette(c("navy", "royalblue", "white", "salmon", "red"))) # pour plotter les correlations entre variables de manière stylée
dev.off()

