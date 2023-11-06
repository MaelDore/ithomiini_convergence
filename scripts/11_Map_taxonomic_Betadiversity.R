##### Script: 35: Map Taxonomic betadiversity based on NMDS from pairwise distances #####

# Author: Maël Doré
# Contact: mael.dore@gmail.com

### Goals = 
   # Plot RGB maps to represent taxonomic betadiversity in Ithomiini
###


### Inputs
  # RGB map already computed from 'Compute_indices_with_functions.R'
###

### Outputs
  # RGB map of taxonomic betadiversity
###

# Clean environment
rm(list = ls())

##### 1/ Prepare data #####

### 1.1/ Load libraries

library(tidyverse)
library(raster)
library(sf)

source(file = "../ithomiini_diversity/functions/compute_betadiversity.R")

### 1.2/ Load/Compute RBG map from pairwise taxonomic Betadiversity ####

### Load the complete stack of species SDM outputs
sp_proba_stack <- readRDS(file = paste0("./input_data/SDM_stacks/All_sp_proba_stack_Jaccard.80.rds"))

# Binarize following ranks (needed for Betadiversity based on Baselga's index)
sp_binary_stack <- binarize_output_with_ranks(sp_proba_stack)

betadiversity_RGB_stack <- map_pairwise_betadiversity(sp_binary_stack, diversity_type = "taxo",
                                                    phylo = NULL, index_family = "sorensen",
                                                    beta_part = "turnover",
                                                    color_space = "3D",
                                                    min_sp = 5,
                                                    subsample_size = 5000,
                                                    RGB_plot = T, NMDS_plot = T)

betadiversity_RGB_stack <- betadiversity_RGB_map
plotRGB(betadiversity_RGB_stack)

#### Add smoothing in the initial resolution 

template_raster <- sp_binary_stack[[1]]

# Resample at initial resolution
betadiversity_RGB_stack <- raster::resample(x = betadiversity_RGB_stack, y = template_raster, method = "bilinear")
  
# Correct for values outside of the [0, 255] range
betadiversity_RGB_stack[betadiversity_RGB_stack[] < 0] <- 0
betadiversity_RGB_stack[betadiversity_RGB_stack[] > 255] <- 255

# Add zero values in the initial range
Red_layer <- Green_layer <- Blue_layer <- template_raster
Red_layer[!is.na(betadiversity_RGB_stack[[1]]@data@values)] <- betadiversity_RGB_stack[[1]]@data@values[!is.na(betadiversity_RGB_stack[[1]]@data@values)]
Green_layer[!is.na(betadiversity_RGB_stack[[2]]@data@values)] <- betadiversity_RGB_stack[[2]]@data@values[!is.na(betadiversity_RGB_stack[[2]]@data@values)]
Blue_layer[!is.na(betadiversity_RGB_stack[[3]]@data@values)] <- betadiversity_RGB_stack[[3]]@data@values[!is.na(betadiversity_RGB_stack[[3]]@data@values)]
betadiversity_RGB_stack <- stack(Red_layer, Green_layer, Blue_layer)
names(betadiversity_RGB_stack) <- c("Red", "Green", "Blue")

# Mask values outside of the initial range
betadiversity_RGB_stack <- mask(betadiversity_RGB_stack, template_raster)

plot(betadiversity_RGB_stack[[1]])
plot(betadiversity_RGB_stack[[2]])
plot(betadiversity_RGB_stack[[3]])

plotRGB(betadiversity_RGB_stack)

saveRDS(betadiversity_RGB_stack, file = paste0("./outputs/Community_structure/Betadiv/betadiversity_RGB_stack.rds"))

# Load already computed RGB map of Taxonomic beta-diversity based on Sorensen's turnover

betadiversity_RGB_stack <- readRDS(file = paste0("./outputs/Community_structure/Betadiv/betadiversity_RGB_stack.rds"))


##### 2/ Plot RBG Map #####

### 2.1/ Load useful stuff ####

# Load country borders
country_borders_sf <- readRDS(file = "./input_data/Map_stuff/country_borders_sf.rds")
# country_borders <- as(country_borders_sf, "Spatial")
plot(country_borders_sf)

# country_borders_Mollweide <- spTransform(x = country_borders, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
country_borders_sf_Mollweide <- st_transform(x = country_borders_sf, crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")

# Load external_borders
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
plot(large_bg_mask_Mollweide)
large_bg_mask_sf_Mollweide <- st_as_sf(large_bg_mask_Mollweide)

### Color cube legend

# Load PNG
color_cube <- png::readPNG(source = "./maps/Community_structure/RGB_cube.png", native = TRUE)
# Convert to raster image
color_cube <- rasterGrob(image = color_cube, interpolate = TRUE)


### 2.2/ Project to Mollweide ####

# Load
betadiversity_RGB_stack <- readRDS(file = paste0("./outputs/Community_structure/Betadiv/betadiversity_RGB_stack.rds"))

plotRGB(betadiversity_RGB_stack)

# Project
RGB_stack_TBD_Mollweide <- projectRaster(from = betadiversity_RGB_stack,
                                   method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                   crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                   alignOnly = F)

# Map
plotRGB(RGB_stack_TBD_Mollweide, axes = T,
        stretch = "lin",
        main = "Taxonomic Betadiversity", cex.main = 1.5,
        colNA = "aliceblue")

# Save
saveRDS(object = RGB_stack_TBD_Mollweide, file = "./maps/Community_structure/RGB_stack_TBD_Mollweide.rds")


### 2.3/ Chose your favorite color scheme ####

# Extract the RGB stack
# RGB_TBD_stack <- betadiversity_RGB_stack[[1]]
RGB_TBD_stack <- RGB_stack_TBD_Mollweide

# Try other options for the RGB bands
RBG_TBD_stack <- subset(x = RGB_TBD_stack, subset = c(1,3,2))
GRB_TBD_stack <- subset(x = RGB_TBD_stack, subset = c(2,1,3))
GBR_TBD_stack <- subset(x = RGB_TBD_stack, subset = c(2,3,1))
BRG_TBD_stack <- subset(x = RGB_TBD_stack, subset = c(3,1,2))
BGR_TBD_stack <- subset(x = RGB_TBD_stack, subset = c(3,2,1))

# Plot all options

par(mfrow = c(2, 3))

plotRGB(RGB_TBD_stack, axes = T,
        stretch = "lin",
        main = "RGB Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(RBG_TBD_stack, axes = T,
        stretch = "lin",
        main = "RBG Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(GRB_TBD_stack, axes = T,
        stretch = "lin",
        main = "GRB Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(GBR_TBD_stack, axes = T,
        stretch = "lin",
        main = "GBR Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(BRG_TBD_stack, axes = T,
        stretch = "lin",
        main = "BRG Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(BGR_TBD_stack, axes = T,
        stretch = "lin",
        main = "BGR Plot", cex.main = 1.5,
        colNA = "aliceblue")

par(mfrow = c(1,1))

# Choose the final map
RGB_TBD_stack_chosen <- BRG_TBD_stack

### 5.5.3/ Correct background color ####

### Extract color scheme per pixel once linear stretching has been applied
new_col_TBD_df <- RStoolbox::ggRGB(img = RGB_TBD_stack_chosen,
                                       r = 1, g = 2, b = 3,
                                       # Can play around with type of stretch to obtain a nicer result
                                       # stretch = "none" 
                                       stretch = "lin",
                                       # stretch = "log"
                                       # stretch = "sqrt"
                                       # stretch = "hist"
                                       ggObj = F)

### Convert Hexadecimal into RGB
new_RGB_TBD_df <- as.data.frame(t(col2rgb(col = new_col_TBD_df$fill)))
colnames(new_RGB_TBD_df) <- colnames(RGB_TBD_stack_chosen[])
# Correct NA values
new_RGB_TBD_df[is.na(new_col_TBD_df$fill), ] <- NA

### Get the cell indices of the background color (= the most frequent one)
background_color <- names(table(new_col_TBD_df$fill)[order(table(new_col_TBD_df$fill), decreasing = T)])[1]
background_coordinates <- new_col_TBD_df[new_col_TBD_df$fill == background_color & !is.na(new_col_TBD_df$fill), c("x", "y")]
background_cells <- raster::cellFromXY(object = RGB_TBD_stack_chosen, xy = background_coordinates)

### Correct background cells values to be grey
new_RGB_TBD_df[background_cells, ] <- t(col2rgb(col = "grey90"))

### Assign new colors to raster
RGB_TBD_background <- RGB_TBD_stack_chosen
RGB_TBD_background@data@values <- as.matrix(new_RGB_TBD_df)

# View(RGB_TBD_background[])

### Check plot
plotRGB(RGB_TBD_background, axes = T, scale = 255,
        stretch = NULL,
        main = "Taxonomic BetaDiversity", cex.main = 1.5,
        colNA = "aliceblue")

# Save
saveRDS(object = RGB_TBD_background, file = "./maps/Community_structure/RGB_TBD_background.rds")

# Load
readRDS(file = "./maps/Community_structure/RGB_TBD_background.rds")


### 5.5.4/ Plot with ggplot2 ####

?RStoolbox::ggRGB()
?ggplot2::theme

ggplot_RGB_TBD_map <- RStoolbox::ggRGB(img = RGB_TBD_background,
                                           r = 1, g = 2, b = 3,
                                           scale = 255,
                                           # Can play around with type of stretch to obtain a nicer result
                                           stretch = "none" 
                                           # stretch = "lin"
                                           # stretch = "log"
                                           # stretch = "sqrt"
                                           # stretch = "hist"
) +
  
  ### Add country borders
  geom_sf(data = country_borders_sf_Mollweide, fill = NA, color = "grey40") +
  
  ### Add external borders
  geom_sf(data = large_bg_mask_sf_Mollweide, fill = NA, color = "black", size = 0.8) +
  
  ### Change projection
  coord_sf(crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs",
           xlim = c(-4660, 4560), ylim = c(-4350, 3400), expand = FALSE) +
  
  ### Add legend color scheme as an image
  annotation_custom(grob = color_cube, xmin = -3750, xmax = -2250, ymin = -3000, ymax = -1500) +
  annotate(geom = "rect", xmin = -4450, xmax = -1550, ymin = -1400, ymax = -550,
           alpha = 1, fill = "aliceblue") +
  annotate(geom = "text", x = -3000, y = -900, hjust = 0.5,
           label = "Taxonomic diversity\n[RGB color space]", fontface = 2, size = 7) +
  
  # ### Add butterfly/mimicry symbol (to oppose to climate symbol)
  # annotation_custom(grob = Climate_symbol, xmin = -3950, xmax = -2150, ymin = -050, ymax = 1900) +
  
  ### Add scale
  ggspatial::annotation_scale(width_hint = 0.2,
                              location = "bl", height = unit(0.010, "npc"),
                              pad_x = unit(0.08, "npc"), pad_y = unit(0.05, "npc"),
                              bar_cols = c("black"),
                              text_cex = 0, text_pad = unit(-0.05, "npc")) +
  
  annotate(geom = "text", x = -3000, y = -3550, hjust = 0.5,
           label = "2000 km", fontface = 2, size = 6.5) +
  
  ### Add north indication
  ggspatial::annotation_north_arrow(
    location = "tr", which_north = "grid",
    height = unit(0.08, "npc"), width = unit(0.06, "npc"),
    pad_x = unit(0.05, "npc"), pad_y = unit(0.05, "npc"),
    style = ggspatial::north_arrow_orienteering(text_col = NA)
  ) +
  
  ### Deal with axes title labels
  # xlab(label = "Longitude  [°]") +
  # ylab(label = "Latitude  [°]") +
  xlab(label = "") +
  ylab(label = "") +
  
  # ### Add panel label
  # annotate(geom = "text", x = 3900, y = -3800, hjust = 0.5,
  #          label = "B", fontface = 2, size = 10) +
  
  ### Deal with general aesthetics
  theme_classic() +
  
  theme(panel.background = element_rect(fill = "aliceblue"),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey70", linetype = "dashed", size = 0.5), # Plot graticules
        axis.ticks = element_line(size = 1.0),
        axis.ticks.length = unit(10, "pt"),
        axis.text = element_text(size = 21, color = "black", face = "bold"),
        axis.text.y = element_text(angle = 90, hjust = 0.5, margin = margin(l = 5, r = 10)),
        axis.text.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title = element_blank(),
        # axis.title = element_text(size = 14, face = "bold"),
        # axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        # axis.title.y = element_text(margin = margin(l = 5, r = 10))
  )


### Export
pdf(file = "./supplementaries/Betadiv_maps/RGB_map_TBD_diversity.pdf", width = 8, height = 8)
print(ggplot_RGB_TBD_map)
dev.off()
