
##### Script 09: Map mimicry turnover and climate beta-diversity with NMDS approach #####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

# Goal: Map mimicry turnover and climate beta-diversity with NMDS approach 

### Input files

# Binary range map
# Ring richness stack
# Borders shp files
# RGB cube, mimicry turnover and climate PNG

### Output files

# Pairwise Ist matrix
# Mimicry turnover map
# Pairwise climatic distances
# Climate beta-diversity map

###

# Clean environment
rm(list = ls())

##### 1/ Load stuff ####

# Load libraries
library(raster)
library(vegan)
library(plot3D)
library(tidyverse)
library(sf)
library(ggimage)
library(ggspatial)
library(fields) # For Thin-plates spline surface model (TPS) for interpolation
library(gstat) # For Inverse Distance Weighting model (IDW) for interpolation

# Load richness maps for mimicry rings

ring_richness_stack <- readRDS(file = "./input_data/SDM_stacks/All_ring_rich_stack_Jaccard.80.rds")
plot(ring_richness_stack)

# Create binary mask of ithomiini range
ring_richness <- raster::calc(x = ring_richness_stack, fun = sum)
plot(ring_richness)
binary_mask <- ring_richness
binary_mask[binary_mask[] > 0] <- 1
plot(binary_mask)

# Save binary mask of ithomiini range
saveRDS(object = binary_mask, file = "./input_data/Map_stuff/binary_mask.rds")

# Load binary mask of ithomiini range
binary_mask <- readRDS(file = "./input_data/Map_stuff/binary_mask.rds")


#  Function to compute simpson index = probability to find species from a different ring when two species are drawn at random
simpson <- function (x, na.rm) {
  if (any(is.na(x))) { # Case with some NA values
    y <- NA
  }else{ # Case with all layers with data
    y <- x/sum(x)
    y <- 1-(sum(y*y))
  }
  return(y)
}

### 2/ Compute and map pairwise mimetic Ist (based on species abundances across mimicry rings) ####

### 2.1/ Function to compute and map pairwise mimetic Ist ####
map_pairwise_Ist_from_stack <- function (ring_richness_stack,
                                         subsample_size = 1000, # Provide number of pixels to subsample globally and regularly in order to save computing time
                                         disaggregate = T, # To resample in the intial resolution for final rasters
                                         interpolate = T, # To interpolate values for the missing communities
                                         # interpolation_method = "IDW", # Chose interpolation model between Thin-Plate Spline (TPS) and Inverse Distance Weighting (IDW)
                                         min_sp = 1, # Minimum number of species to include a community
                                         maxit = 20, # Maximum numbe rof iterations for the NMDS optimization
                                         NMDS_plot = T, # To plot community coordinates in 2D and 3D scatterplots 
                                         RGB_plot = T)  # To plot rasters of pairwise community betadiversity in RGB color bands and a composite RGB layer
{
  ### Prepare data
  
  # Keep a template of initial raster resolution and range
  template_raster <- ((ring_richness_stack >= 0)[[1]] + 219) # Add 219 to get a grey color in the RGB plot (220, 220, 220) in terrestrial communities with no data 
  
  # Subsample sites in a standardized fashion such as density of sampling remains the same for all regions
  if (!is.na(subsample_size))
  {
    ring_richness_stack <- raster::sampleRegular(x = ring_richness_stack, size = subsample_size, asRaster = T)
  }
  
  # Extract community data in matrix
  ring_richness_brick <- ring_richness_stack*1
  community_matrix <- ring_richness_brick@data@values
  
  # Remove communities with NA
  community_matrix_clean_NA <- na.omit(community_matrix)
  NA_com_indices <- as.numeric(attr(community_matrix_clean_NA, which = "na.action"))
  
  # Remove communities with no or less than the minimum species threshold
  empty_com_indices <- which(!(apply(X = community_matrix, MARGIN = 1, FUN = sum) > min_sp))
  
  # Remove both NA and empty communities
  removed_com_indices <- c(NA_com_indices, empty_com_indices)
  community_matrix_clean <- community_matrix[-removed_com_indices, ]
  
  ### Compute pairwise Ist matrix
  
  # Generate squared matrix to store pairwise Ist between communities
  pairwise_mimetic_Ist <- matrix(ncol = nrow(community_matrix_clean), nrow = nrow(community_matrix_clean)) 
  
  for (i in 1:nrow(pairwise_mimetic_Ist)) # For each community in row
  {
    for (j in 1:ncol(pairwise_mimetic_Ist)) # For each community in column
    {
      com1 <- community_matrix_clean[i,] # Get mimicry richnesses of all rings for community 1
      com2 <- community_matrix_clean[j,] # Get mimicry richnesses of all rings for community 2
      D1 <- simpson(com1) # Compute local simpson index for community 1
      D2 <- simpson(com2) # Compute local simpson index for community 2
      Ds <- mean(c(D1,D2)) # Compute mean simpson index among the pair
      full_com <- apply(X = community_matrix_clean[c(i,j),], MARGIN = 2, FUN = sum) # Get mimicry richnesses of all rings for the two communities combined
      Dt <- simpson(full_com) # Compute Simpson for both communities merged
      pairwise_mimetic_Ist[i,j] <- 1-(Ds/Dt) # Compute pairwise Ist
    }
    
    # Show i every 100 iterations
    if (i %% 10 == 0) { cat(paste0("\nCompute Ist for community n°",i,"/",nrow(pairwise_mimetic_Ist)))}
  }
  
  ### Convert pairwise Ist matrix into 3D NMDS
  
  NMDS <- vegan::metaMDS(comm = pairwise_mimetic_Ist, k = 3, maxit = maxit)
  
  # Scale NMDS output into 0 to 255 range
  color_data <- apply(X = NMDS$points, MARGIN = 2, FUN = scales::rescale, to = c(0, 255))
  color_data <- apply(X = color_data, MARGIN = 2, FUN = round)
  
  # Add empty sites
  color_raster_data <- matrix(data = NA, nrow = nrow(community_matrix), ncol = 3)
  color_raster_data[-removed_com_indices, ] <- color_data
  colnames(color_raster_data) <- c("Red", "Green", "Blue")
  
  # Create template raster Brick
  RGB_brick <- brick(stack(ring_richness_stack[[1]], ring_richness_stack[[1]], ring_richness_stack[[1]]))
  names(RGB_brick) <- c("Red", "Green", "Blue")
  # Fill with data
  RGB_brick@data@values <- color_raster_data
  # Convert to raster Stack
  RGB_stack <- stack(RGB_brick)
  
  ### Plot community coordinates in 2D and 3D scatterplots if requested
  
  if(NMDS_plot)
  {
    par(mfrow = c(2,2))
    
    plot(color_data[, 1:2], type = "n", xlab = "Red Band (NMDS1)", ylab = "Green Band (NMDS2)")
    points(x = color_data[, 1:2], pch = 16, 
           col = rgb(red = color_data[, 1], green = color_data[, 2], blue = 0, maxColorValue = 255))
    
    plot(color_data[, c(1,3)], type = "n", xlab = "Red Band (NMDS1)", ylab = "Blue Band (NMDS3)")
    points(x = color_data[, c(1,3)], pch = 16, 
           col = rgb(red = color_data[, 1], green = 0, blue = color_data[, 3], maxColorValue = 255))
    
    plot(color_data[, c(3,2)], type = "n", xlab = "Blue Band (NMDS1)", ylab = "Green Band (NMDS3)")
    points(x = color_data[, c(3,2)], pch = 16, 
           col = rgb(red = 0, green = color_data[, 2], blue = color_data[, 3], maxColorValue = 255))
    
    ### 3D plot in RGB
    
    plot3D::scatter3D(x = color_data[, 1], y =  color_data[, 2], z =  color_data[, 3],
                      bty = "b2", colkey = FALSE, theta = 45, phi = 25,
                      ticktype = "detailed",
                      pch = 16, alpha = 0.7,
                      colvar = rep(NA, nrow(color_data)),
                      NAcol = rgb(red = color_data[, 1], green = color_data[, 2], blue = color_data[, 3], maxColorValue = 255),
                      main = "RGB plot from NMDS", xlab = "\nRed Band",
                      ylab = "\nGreen Band", zlab = "\nBlue Band")
    par(mfrow = c(1,1))
  } 
 
  #### Add interpolation for values in communities excluded from the computation (low species richness) but still within the groupe range
  if (interpolate)
  {
    # Loop per RGB band
    for (RGB_band in 1:3)
    {
      # RGB_band <- 1
      
      # Extract coordinates of raster cells
      xy <- data.frame(xyFromCell(RGB_stack[[RGB_band]], 1:ncell(RGB_stack[[RGB_band]])))
      
      # Extract values
      v <- getValues(RGB_stack[[RGB_band]])
      # Remove NAs
      i <- !is.na(v)
      # Filter coordinates
      xy <- xy[i,]
      # Filter values
      v <- v[i]
      
      #### Inverse Distance Weighting (IDW)
      
      # ?gstat::gstat
      df <- cbind(v, xy) # Build df of data and coordinates
      
      # Calibrate model
      IDW <- gstat(id = "v", formula = v ~ 1, locations = ~ x + y, data = df, 
                   nmax = 8, set=list(idp = .5))
      interpolation_model <- IDW
      
      # #### Thin-Plate Spline model (TPS)
      # if (interpolation_method == "TPS")
      # {
      #   # ?fields::Tps 
      #   tps <- fields::Tps(xy, v)
      #   interpolation_model <- tps
      # }
      
      ### Use model to predict values at all locations
      
      # ?raster::interpolate
      RGB_stack[[RGB_band]] <- raster::interpolate(RGB_stack[[RGB_band]], model = interpolation_model)
    }
    
    ### Get mask of group range to proper resolution
    richness_map <- calc(x = ring_richness_stack, fun = sum)
    range_mask <- richness_map
    range_mask[range_mask[] == 0] <- NA
    
    ### Apply mask
    # plot(RGB_stack)
    RGB_stack <- mask(RGB_stack, range_mask)
    # plot(RGB_stack)
  }
  
  
  #### Add smoothing in the initial resolution if requested ####  
  if (disaggregate)
  {
    # Resample at initial resolution
    RGB_stack <- raster::resample(x = RGB_stack, y = template_raster, method = "bilinear")
    
    # Correct for values outside of the [0, 255] range
    RGB_stack[RGB_stack[] < 0] <- 0
    RGB_stack[RGB_stack[] > 255] <- 255
    
    # Add zero values in the initial range
    Red_layer <- Green_layer <- Blue_layer <- template_raster
    Red_layer[!is.na(RGB_stack[[1]]@data@values)] <- RGB_stack[[1]]@data@values[!is.na(RGB_stack[[1]]@data@values)]
    Green_layer[!is.na(RGB_stack[[2]]@data@values)] <- RGB_stack[[2]]@data@values[!is.na(RGB_stack[[2]]@data@values)]
    Blue_layer[!is.na(RGB_stack[[3]]@data@values)] <- RGB_stack[[3]]@data@values[!is.na(RGB_stack[[3]]@data@values)]
    RGB_stack <- stack(Red_layer, Green_layer, Blue_layer)
    names(RGB_stack) <- c("Red", "Green", "Blue")
    
    # Mask values outside of the initial range
    RGB_stack <- mask(RGB_stack, template_raster)
  }

  #### Plot RGB maps ####

  if (RGB_plot)
  {
    red_palette <- colorRampPalette(c("white","red"))(256)
    green_palette <- colorRampPalette(c("white","green"))(256)
    blue_palette <- colorRampPalette(c("white","blue"))(256)
    
    par(mfrow = c(2,2))
    
    image(RGB_stack[[1]], col = red_palette)
    title(main = "Red band (NMDS1)", cex.main = 1.5, line = 1.5)
    
    image(RGB_stack[[2]], col = green_palette, colNA = "aliceblue")
    title(main = "Green band (NMDS2)", cex.main = 1.5, line = 1.5)
    
    image(RGB_stack[[3]], col = blue_palette, colNA = "aliceblue")
    title(main = "Blue band (NMDS3)", cex.main = 1.5, line = 1.5)
    
    plotRGB(RGB_stack, axes = T, main = "RGB Plot", cex.main = 1.5)
    # title(main = "RGB Plot", cex.main = 1.5, line = 1.5)
    
    par(mfrow = c(1,1))
    
  }

  # Export the RGB stack to map and the pairwise mimetic Ist matrix (on subsampled communities)
  output <- list(RGB_stack, pairwise_mimetic_Ist)
  return(output)
}

### 2.2/ Compute and map pairwise mimetic Ist ####

pairwise_Ist_RGB_map <- map_pairwise_Ist_from_stack(ring_richness_stack,
                                                    subsample_size = 10000, min_sp = 5,
                                                    disaggregate = T, interpolate = T)

# Save
saveRDS(object = pairwise_Ist_RGB_map, file = "./supplementaries/Betadiv_maps/pairwise_Ist_RGB_map.rds")

# Load
pairwise_Ist_RGB_map <- readRDS(file = "./supplementaries/Betadiv_maps/pairwise_Ist_RGB_map.rds")

### 3/ Reproject RGB raster in Mollweide ####

Mollweide_projection <- function(x) # Raster to project
{
  new_map <- projectRaster(from = x,
                           method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                           crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                           alignOnly = F)
  return(new_map)
}

RGB_mimicry_stack_Mollweide <- Mollweide_projection(pairwise_Ist_RGB_map[[1]])

plotRGB(RGB_mimicry_stack_Mollweide, axes = T,
        stretch = "lin",
        main = "RGB Plot", cex.main = 1.5,
        colNA = "aliceblue")

# Save
saveRDS(object = RGB_mimicry_stack_Mollweide, file = "./supplementaries/Betadiv_maps/RGB_mimicry_stack_Mollweide.rds")

# Load
RGB_mimicry_stack_Mollweide <- readRDS(file = "./supplementaries/Betadiv_maps/RGB_mimicry_stack_Mollweide.rds")


### 4/ Plot mimicry turnover ####

### 4.0/ Load useful stuff ####

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

### Grid = graticules. No need for ggplot2 that can plot the graticules itself

# grid_Mollweide <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
# grid_Mollweide_sf <- st_as_sf(grid_Mollweide) # Convert to sf

# Conversion WGS84 to Mollweide
# -120 = -4660
# -87.5 = -2330
# -75 = 0
# -52.5 = 2330
# -30 = 4660

# 3560
# -4595

### Color cube legend

# Load PNG
color_cube <- png::readPNG(source = "./maps/Community_structure/RGB_cube.png", native = TRUE)
# Convert to raster image
color_cube <- rasterGrob(image = color_cube, interpolate = TRUE)

### Mimicry legend

# Load PNG
Mimicry_symbol <- png::readPNG(source = "./maps/Community_structure/Mimicry_symbol.png", native = TRUE)
# Convert to raster image
Mimicry_symbol <- rasterGrob(image = Mimicry_symbol, interpolate = TRUE)


### 4.1/ Plot with raster ####

?raster::plotRGB

# Extract the RGB stack
# RGB_mimicry_stack <- pairwise_Ist_RGB_map[[1]]
RGB_mimicry_stack <- RGB_mimicry_stack_Mollweide

# Try other options for the RGB bands
RBG_mimicry_stack <- subset(x = RGB_mimicry_stack, subset = c(1,3,2))
GRB_mimicry_stack <- subset(x = RGB_mimicry_stack, subset = c(2,1,3))
GBR_mimicry_stack <- subset(x = RGB_mimicry_stack, subset = c(2,3,1))
BRG_mimicry_stack <- subset(x = RGB_mimicry_stack, subset = c(3,1,2))
BGR_mimicry_stack <- subset(x = RGB_mimicry_stack, subset = c(3,2,1))

# Plot all options

par(mfrow = c(2, 3))

plotRGB(RGB_mimicry_stack, axes = T,
        stretch = "lin",
        main = "RGB Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(RBG_mimicry_stack, axes = T,
        stretch = "lin",
        main = "RBG Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(GRB_mimicry_stack, axes = T,
        stretch = "lin",
        main = "GRB Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(GBR_mimicry_stack, axes = T,
        stretch = "lin",
        main = "GBR Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(BRG_mimicry_stack, axes = T,
        stretch = "lin",
        main = "BRG Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(BGR_mimicry_stack, axes = T,
        stretch = "lin",
        main = "BGR Plot", cex.main = 1.5,
        colNA = "aliceblue")

par(mfrow = c(1,1))

# Choose the final map
# RGB_mimicry_stack_chosen <- RGB_mimicry_stack
RGB_mimicry_stack_chosen <- GRB_mimicry_stack
# RGB_mimicry_stack_chosen <- BGR_mimicry_stack

### 4.2/ Correct background color ####

### Extract color scheme per pixel once linear stretching has been applied
new_col_mimicry_df <- RStoolbox::ggRGB(img = RGB_mimicry_stack_chosen,
                               r = 1, g = 2, b = 3,
                               # Can play around with type of stretch to obtain a nicer result
                               # stretch = "none" 
                               stretch = "lin",
                               # stretch = "log"
                               # stretch = "sqrt"
                               # stretch = "hist"
                               ggObj = F)

### Convert Hexadecimal into RGB
new_RGB_mimicry_df <- as.data.frame(t(col2rgb(col = new_col_mimicry_df$fill)))
colnames(new_RGB_mimicry_df) <- colnames(RGB_mimicry_stack_chosen[])
# Correct NA values
new_RGB_mimicry_df[is.na(new_col_mimicry_df$fill), ] <- NA

### Get the cell indices of the background color (= the most frequent one)
background_color <- names(table(new_col_mimicry_df$fill)[order(table(new_col_mimicry_df$fill), decreasing = T)])[1]
background_coordinates <- new_col_mimicry_df[new_col_mimicry_df$fill == background_color & !is.na(new_col_mimicry_df$fill), c("x", "y")]
background_cells <- raster::cellFromXY(object = RGB_mimicry_stack_chosen, xy = background_coordinates)

### Correct background cells values to be grey
new_RGB_mimicry_df[background_cells, ] <- t(col2rgb(col = "grey90"))

### Assign new colors to raster
RGB_mimicry_background <- RGB_mimicry_stack_chosen
RGB_mimicry_background@data@values <- as.matrix(new_RGB_mimicry_df)

# View(RGB_mimicry_background[])

### Check plot
plotRGB(RGB_mimicry_background, axes = T, scale = 255,
        stretch = NULL,
        main = "Mimicry turnover", cex.main = 1.5,
        colNA = "aliceblue")

# Save
saveRDS(object = RGB_mimicry_background, file = "./supplementaries/Betadiv_maps/RGB_mimicry_background.rds")

# Load
RGB_mimicry_background <- readRDS(file = "./supplementaries/Betadiv_maps/RGB_mimicry_background.rds")



### 4.3/ Plot with ggplot2 ####

?RStoolbox::ggRGB()
?ggplot2::theme

ggplot_RGB_mimicry_map <- RStoolbox::ggRGB(img = RGB_mimicry_background,
                                           r = 1, g = 2, b = 3,
                                           scale = 255,
                                           # Can play around with type of stretch to obtain a nicer result
                                           stretch = "none" 
                                           # stretch = "lin"
                                           # stretch = "log"
                                           # stretch = "sqrt"
                                           # stretch = "hist"
                                           ) +
  
  ### Find a way to obtain grey for background values 
  
  # geom_sf(mapping = aes(), data = grid_Mollweide_sf) +
  
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
           label = "Mimicry turnover\n[RGB color space]", fontface = 2, size = 7) +

  
  ### Add butterfly/mimicry symbol (to oppose to climate symbol)
  annotation_custom(grob = Mimicry_symbol, xmin = -4200, xmax = -1900, ymin = -250, ymax = 1900) +
  
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
  
  ### Add panel label
  annotate(geom = "text", x = 3900, y = -3800, hjust = 0.5,
           label = "A", fontface = 2, size = 10) +
  
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
pdf(file = "./supplementaries/Betadiv_maps/RGB_map_Mimicry_turnover.pdf", width = 8, height = 8)
print(ggplot_RGB_mimicry_map)
dev.off()



### 4.4/ Plot with base plot ####

# Need to convert to a Raster Layer with RAT table that provide color name to each pixel value (not efficient at all...)

r <- raster(ncol=10, nrow=10)
values(r) <- sample(0:255, ncell(r), replace=TRUE)
ctab <- sample(rainbow(256))
colortable(r) <- ctab
plot(r)
head(colortable(r))

colortable(RGB_mimicry_stack_chosen)

?image

image(x, col = color_palette,
      xlim = xlim, ylim = ylim, axes = F,
      xlab = xlab, ylab = ylab)


##### 5/ Map climatic diversity #####

### 5.0/ Load climate data

# Load SDM input data
Env_stack <- readRDS("D:/Mael/R_projects/ithomiini_convergence/input_data/Select_env_15.rds")
names(Env_stack)

# Select only climate variables
Climate_stack <- raster::subset(x = Env_stack, subset = c("bio1", "bio4", "bio12", "bio15"))
names(Climate_stack) <- c("Tmean", "Tvar", "Ptot", "Pvar")

mask <- binary_mask

### 5.1/ Function to compute and map pairwise distances from stack ####
map_pairwise_distances_from_stack <- function (variables_stack,
                                               distance_method = "euclidean", # Set the type of dissimilarity metric to use
                                               mask = NULL, # Provide a binary raster to use as a mask to limit the range of the grid cells involved in the analyses
                                               subsample_size = 1000, # Provide number of pixels to subsample globally and regularly in order to save computing time
                                               disaggregate = T, # To resample in the intial resolution for final rasters
                                               center = T, # Apply centering of variables before computing distances
                                               scale = T, # Apply scaling of variables before computing distances
                                               maxit = 20, # Maximum number of iterations for the NMDS optimization
                                               NMDS_plot = T, # To plot community coordinates in 2D and 3D scatterplots 
                                               RGB_plot = T)  # To plot rasters of pairwise community dissimilarity in RGB color bands and a composite RGB layer
{
  ### Prepare data
  
  # Keep a template of initial raster resolution and range
  template_raster <- ((variables_stack >= 0)[[1]] + 219) # Add 219 to get a grey color in the RGB plot (220, 220, 220) in terrestrial communities with no data 
  
  # Set to NA community outside of the mask if required
  if (!is.null(mask))
  {
    # Transform to NA the 0 values in the mask
    mask[mask[] == 0] <- NA
    # Apply mask
    variables_stack <- raster::mask(x = variables_stack, mask = mask)
  }
  
  # Subsample sites in a standardized fashion such as density of sampling remains the same for all regions
  if (!is.na(subsample_size))
  {
    variables_stack <- raster::sampleRegular(x = variables_stack, size = subsample_size, asRaster = T)
  }
  
  # Extract community data in matrix
  variables_brick <- variables_stack*1
  community_matrix <- variables_brick@data@values
  
  # Remove communities with NA
  community_matrix_clean <- na.omit(community_matrix)
  NA_com_indices <- as.numeric(attr(community_matrix_clean, which = "na.action"))

  ### Standardized data if required
  community_matrix_clean <- scale(x = community_matrix_clean, center = center, scale = scale)
  
  ### Compute pairwise distances
  
  cat(paste0("\n", Sys.time(), " - Compute pairwise distances\n"))
  
  com_pairwise_distances <- vegan::vegdist(x = community_matrix_clean, method = distance_method)
  
  ### Convert pairwise distances into 3D NMDS
  
  cat(paste0("\n", Sys.time(), " - Convert pairwise distances into 3D NMDS\n"))
  
  NMDS <- vegan::metaMDS(comm = com_pairwise_distances, k = 3, maxit = maxit)
  
  # Scale NMDS output into 0 to 255 range
  color_data <- apply(X = NMDS$points, MARGIN = 2, FUN = scales::rescale, to = c(0, 255))
  color_data <- apply(X = color_data, MARGIN = 2, FUN = round)
  
  # Add empty sites
  color_raster_data <- matrix(data = NA, nrow = nrow(community_matrix), ncol = 3)
  color_raster_data[-NA_com_indices, ] <- color_data
  colnames(color_raster_data) <- c("Red", "Green", "Blue")
  
  # Create template raster Brick
  RGB_brick <- brick(stack(variables_stack[[1]], variables_stack[[1]], variables_stack[[1]]))
  names(RGB_brick) <- c("Red", "Green", "Blue")
  # Fill with data
  RGB_brick@data@values <- color_raster_data
  # Convert to raster Stack
  RGB_stack <- stack(RGB_brick)
  
  ### Plot community coordinates in 2D and 3D scatterplots if requested
  
  if(NMDS_plot)
  {
    cat(paste0("\n", Sys.time(), " - Plot 2D and 3D scatterplots\n"))
    
    par(mfrow = c(2,2))
    
    plot(color_data[, 1:2], type = "n", xlab = "Red Band (NMDS1)", ylab = "Green Band (NMDS2)")
    points(x = color_data[, 1:2], pch = 16, 
           col = rgb(red = color_data[, 1], green = color_data[, 2], blue = 0, maxColorValue = 255))
    
    plot(color_data[, c(1,3)], type = "n", xlab = "Red Band (NMDS1)", ylab = "Blue Band (NMDS3)")
    points(x = color_data[, c(1,3)], pch = 16, 
           col = rgb(red = color_data[, 1], green = 0, blue = color_data[, 3], maxColorValue = 255))
    
    plot(color_data[, c(3,2)], type = "n", xlab = "Blue Band (NMDS1)", ylab = "Green Band (NMDS3)")
    points(x = color_data[, c(3,2)], pch = 16, 
           col = rgb(red = 0, green = color_data[, 2], blue = color_data[, 3], maxColorValue = 255))
    
    ### 3D plot in RGB
    
    plot3D::scatter3D(x = color_data[, 1], y =  color_data[, 2], z =  color_data[, 3],
                      bty = "b2", colkey = FALSE, theta = 45, phi = 25,
                      ticktype = "detailed",
                      pch = 16, alpha = 0.7,
                      colvar = rep(NA, nrow(color_data)),
                      NAcol = rgb(red = color_data[, 1], green = color_data[, 2], blue = color_data[, 3], maxColorValue = 255),
                      main = "RGB plot from NMDS", xlab = "\nRed Band",
                      ylab = "\nGreen Band", zlab = "\nBlue Band")
    par(mfrow = c(1,1))
  } 
  
  #### Add smoothing in the initial resolution if requested ####  
  if (disaggregate)
  {
    # Resample at initial resolution
    RGB_stack <- raster::resample(x = RGB_stack, y = template_raster, method = "bilinear")
    
    # Correct for values outside of the [0, 255] range
    RGB_stack[RGB_stack[] < 0] <- 0
    RGB_stack[RGB_stack[] > 255] <- 255
    
    # Add zero values in the initial range
    Red_layer <- Green_layer <- Blue_layer <- template_raster
    Red_layer[!is.na(RGB_stack[[1]]@data@values)] <- RGB_stack[[1]]@data@values[!is.na(RGB_stack[[1]]@data@values)]
    Green_layer[!is.na(RGB_stack[[2]]@data@values)] <- RGB_stack[[2]]@data@values[!is.na(RGB_stack[[2]]@data@values)]
    Blue_layer[!is.na(RGB_stack[[3]]@data@values)] <- RGB_stack[[3]]@data@values[!is.na(RGB_stack[[3]]@data@values)]
    RGB_stack <- stack(Red_layer, Green_layer, Blue_layer)
    names(RGB_stack) <- c("Red", "Green", "Blue")
    
    # Mask values outside of the initial range
    RGB_stack <- mask(RGB_stack, template_raster)
  }
  
  #### Plot RGB maps ####
  
  if (RGB_plot)
  {
    cat(paste0("\n", Sys.time(), " - Plot RGB bands and RGB plot\n"))
    
    red_palette <- colorRampPalette(c("white","red"))(256)
    green_palette <- colorRampPalette(c("white","green"))(256)
    blue_palette <- colorRampPalette(c("white","blue"))(256)
    
    par(mfrow = c(2,2))
    
    image(RGB_stack[[1]], col = red_palette)
    title(main = "Red band (NMDS1)", cex.main = 1.5, line = 1.5)
    
    image(RGB_stack[[2]], col = green_palette, colNA = "aliceblue")
    title(main = "Green band (NMDS2)", cex.main = 1.5, line = 1.5)
    
    image(RGB_stack[[3]], col = blue_palette, colNA = "aliceblue")
    title(main = "Blue band (NMDS3)", cex.main = 1.5, line = 1.5)
    
    plotRGB(RGB_stack, axes = T, scale = 255, main = "RGB Plot", cex.main = 1.5)
    # title(main = "RGB Plot", cex.main = 1.5, line = 1.5)
    
    par(mfrow = c(1,1))
    
  }
  
  # Export the RGB stack to map and the pairwise distances (on subsampled communities)
  output <- list(RGB_stack, com_pairwise_distances)
  return(output)
}


### 5.2/ Compute and map climatic distances across communities ####

pairwise_climate_RGB_map <- map_pairwise_distances_from_stack(variables_stack = Climate_stack,
                                                              distance_method = "euclidean", # Set the type of dissimilarity metric to use
                                                              mask = binary_mask,
                                                              subsample_size = 10000, # Provide number of pixels to subsample globally and regularly in order to save computing time
                                                              disaggregate = T, # To resample in the intial resolution for final rasters
                                                              center = T, # Apply centering of variables before computing distances
                                                              scale = T, # Apply scaling of variables before computing distances
                                                              maxit = 20, # Maximum number of iterations for the NMDS optimization
                                                              NMDS_plot = T, # To plot community coordinates in 2D and 3D scatterplots 
                                                              RGB_plot = T) 

# Save
saveRDS(object = pairwise_climate_RGB_map, file = "./supplementaries/Betadiv_maps/pairwise_climate_RGB_map.rds")

# Load
pairwise_climate_RGB_map <- readRDS(file = "./supplementaries/Betadiv_maps/pairwise_climate_RGB_map.rds")


### 5.3/ Procrustes adjustment to match with mimicry turnover color gradient ####

# Rationale = only distances matter in the NMDS space, so we can use Procrustes adjustment with translation and rotation (and eventually scaling),
# to make coordinates in the RGB space match between mimicry turnover and climate turnover

# Load mimicry turnover output
pairwise_Ist_RGB_map <- readRDS(file = "./supplementaries/Betadiv_maps/pairwise_Ist_RGB_map.rds")

RGB_mimicry_stack <- pairwise_Ist_RGB_map[[1]]
RGB_climate_stack <- pairwise_climate_RGB_map[[1]]

# Remove data out of the range
binary_mask <- readRDS(file = "./input_data/Map_stuff/binary_mask.rds")
NA_mask <- binary_mask
NA_mask[NA_mask[] == 0] <- NA

RGB_mimicry_stack_masked <- raster::mask(x = RGB_mimicry_stack, mask = NA_mask) 
RGB_climate_stack_masked <- raster::mask(x = RGB_climate_stack, mask = NA_mask) 

# Extract NMDS coordinates for communities without NA
RGB_mimicry_coords <- na.omit(RGB_mimicry_stack_masked[])
RGB_climate_coords <- na.omit(RGB_climate_stack_masked[])
NA_com_indices <- as.numeric(attr(RGB_mimicry_coords, which = "na.action"))

RGB_climate_coords_Procrutes <- vegan::procrustes(X = RGB_mimicry_coords, Y = RGB_climate_coords, scale = TRUE, symmetric = FALSE)

# Visualize residuals after Procrustes adjustment = Rotated Y vs target X
plot(RGB_climate_coords_Procrutes)

# Visualize the Procrustes adjustment = Rotated Y vs initial Y
par(mfrow = c(2,2))
plot(RGB_mimicry_coords[, c(1,2)], main = "Mimicry turnover") ; plot(RGB_mimicry_coords[, c(1,3)], main = "Mimicry turnover") ; plot(RGB_mimicry_coords[, c(2,3)], main = "Mimicry turnover") ; plot.new()
plot(RGB_climate_coords[, c(1,2)], main = "Climate turnover") ; plot(RGB_climate_coords[, c(1,3)], main = "Climate turnover") ; plot(RGB_climate_coords[, c(2,3)], main = "Climate turnover") ; plot.new()
plot(RGB_climate_coords_Procrutes$Yrot[, c(1,2)], main = "Climate (Procrustes)") ; plot(RGB_climate_coords_Procrutes$Yrot[, c(1,3)], main = "Climate (Procrustes)") ; plot(RGB_climate_coords_Procrutes$Yrot[, c(2,3)], main = "Climate (Procrustes)") ; plot.new()
par(mfrow = c(1,1))
plot(RGB_climate_coords_Procrutes$Yrot)
arrows(x0 = RGB_climate_coords[,1], y0 = RGB_climate_coords[,2], x1 = RGB_climate_coords_Procrutes$Yrot[,1], y1 = RGB_climate_coords_Procrutes$Yrot[,2], length = 0.02, col = "red")

# Need to reajust into [0, 255] dimensions...
RGB_climate_coords_Procrutes_rescale <- apply(X = RGB_climate_coords_Procrutes$Yrot, MARGIN = 2, FUN = scales::rescale, to = c(0, 255))
RGB_climate_coords_Procrutes_rescale <- apply(X = RGB_climate_coords_Procrutes_rescale, MARGIN = 2, FUN = round)

# Add empty sites
RGB_climate_Procrutes_data <- matrix(data = NA, nrow = nrow(RGB_mimicry_stack[]), ncol = 3)
RGB_climate_Procrutes_data[-NA_com_indices, ] <- RGB_climate_coords_Procrutes_rescale
colnames(RGB_climate_Procrutes_data) <- c("Red", "Green", "Blue")

# Create template raster Brick
RGB_climate_Procrustes_brick <- brick(RGB_climate_stack)
names(RGB_climate_Procrustes_brick) <- c("Red", "Green", "Blue")
# Fill with data
RGB_climate_Procrustes_brick@data@values <- RGB_climate_Procrutes_data
RGB_climate_Procrustes_brick@data@inmemory <- T
# Add background grey data
RGB_climate_Procrustes_brick[binary_mask[] == 0, ] <- col2rgb(col = "grey90")
# Convert to raster Stack
RGB_climate_Procrustes_stack <- stack(RGB_climate_Procrustes_brick)


### Compare maps

par(mfrow = c(2,2))

# A/ Plot mimicry turnover
plotRGB(RGB_mimicry_stack, axes = T, scale = 255, main = "Mimicry turnover", cex.main = 1.5)
# B/ Plot climate turnover
plotRGB(RGB_climate_stack, axes = T, scale = 255, main = "Climate turnover", cex.main = 1.5)
# C/ Plot climate turnover
plotRGB(RGB_climate_Procrustes_stack, axes = T, scale = 255, main = "Climate turnover (Procrustes)", cex.main = 1.5)

par(mfrow = c(1,1))
# Results is not really convincing...


### 5.4/ Project to Mollweide ####

# Load
pairwise_climate_RGB_map <- readRDS(file = "./supplementaries/Betadiv_maps/pairwise_climate_RGB_map.rds")

plotRGB(pairwise_climate_RGB_map[[1]])

# Project
RGB_climate_stack_Mollweide <- Mollweide_projection(pairwise_climate_RGB_map[[1]])

# Map
plotRGB(RGB_climate_stack_Mollweide, axes = T,
        stretch = "lin",
        main = "RGB Climate Plot", cex.main = 1.5,
        colNA = "aliceblue")

# Save
saveRDS(object = RGB_climate_stack_Mollweide, file = "./supplementaries/Betadiv_maps/RGB_climate_stack_Mollweide.rds")

# Load
RGB_climate_stack_Mollweide <- readRDS(file = "./supplementaries/Betadiv_maps/RGB_climate_stack_Mollweide.rds")

### 5.5/ Plot climate diversity ####

### 5.5.1/ Load useful stuff ####

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

### Climate symbol

# Load PNG
Climate_symbol <- png::readPNG(source = "./maps/Community_structure/Climate_symbol.png", native = TRUE)
# Convert to raster image
Climate_symbol <- rasterGrob(image = Climate_symbol, interpolate = TRUE)

### 5.5.2/ Chose your favorite color scheme ####

# Extract the RGB stack
# RGB_climate_stack <- pairwise_climate_RGB_map[[1]]
RGB_climate_stack <- RGB_climate_stack_Mollweide

# Try other options for the RGB bands
RBG_climate_stack <- subset(x = RGB_climate_stack, subset = c(1,3,2))
GRB_climate_stack <- subset(x = RGB_climate_stack, subset = c(2,1,3))
GBR_climate_stack <- subset(x = RGB_climate_stack, subset = c(2,3,1))
BRG_climate_stack <- subset(x = RGB_climate_stack, subset = c(3,1,2))
BGR_climate_stack <- subset(x = RGB_climate_stack, subset = c(3,2,1))

# Plot all options

par(mfrow = c(2, 3))

plotRGB(RGB_climate_stack, axes = T,
        stretch = "lin",
        main = "RGB Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(RBG_climate_stack, axes = T,
        stretch = "lin",
        main = "RBG Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(GRB_climate_stack, axes = T,
        stretch = "lin",
        main = "GRB Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(GBR_climate_stack, axes = T,
        stretch = "lin",
        main = "GBR Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(BRG_climate_stack, axes = T,
        stretch = "lin",
        main = "BRG Plot", cex.main = 1.5,
        colNA = "aliceblue")

plotRGB(BGR_climate_stack, axes = T,
        stretch = "lin",
        main = "BGR Plot", cex.main = 1.5,
        colNA = "aliceblue")

par(mfrow = c(1,1))

# Choose the final map
RGB_climate_stack_chosen <- GBR_climate_stack

### 5.5.3/ Correct background color ####

### Extract color scheme per pixel once linear stretching has been applied
new_col_climate_df <- RStoolbox::ggRGB(img = RGB_climate_stack_chosen,
                               r = 1, g = 2, b = 3,
                               # Can play around with type of stretch to obtain a nicer result
                               # stretch = "none" 
                               stretch = "lin",
                               # stretch = "log"
                               # stretch = "sqrt"
                               # stretch = "hist"
                               ggObj = F)

### Convert Hexadecimal into RGB
new_RGB_climate_df <- as.data.frame(t(col2rgb(col = new_col_climate_df$fill)))
colnames(new_RGB_climate_df) <- colnames(RGB_climate_stack_chosen[])
# Correct NA values
new_RGB_climate_df[is.na(new_col_climate_df$fill), ] <- NA

### Get the cell indices of the background color (= the most frequent one)
background_color <- names(table(new_col_climate_df$fill)[order(table(new_col_climate_df$fill), decreasing = T)])[1]
background_coordinates <- new_col_climate_df[new_col_climate_df$fill == background_color & !is.na(new_col_climate_df$fill), c("x", "y")]
background_cells <- raster::cellFromXY(object = RGB_climate_stack_chosen, xy = background_coordinates)

### Correct background cells values to be grey
new_RGB_climate_df[background_cells, ] <- t(col2rgb(col = "grey90"))

### Assign new colors to raster
RGB_climate_background <- RGB_climate_stack_chosen
RGB_climate_background@data@values <- as.matrix(new_RGB_climate_df)

# View(RGB_climate_background[])

### Check plot
plotRGB(RGB_climate_background, axes = T, scale = 255,
        stretch = NULL,
        main = "Climate diversity", cex.main = 1.5,
        colNA = "aliceblue")

# Save
saveRDS(object = RGB_climate_background, file = "./supplementaries/Betadiv_maps/RGB_climate_background.rds")

# Load
readRDS(file = "./supplementaries/Betadiv_maps/RGB_climate_background.rds")


### 5.5.4/ Plot with ggplot2 ####

?RStoolbox::ggRGB()
?ggplot2::theme

ggplot_RGB_climate_map <- RStoolbox::ggRGB(img = RGB_climate_background,
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
           label = "Climate diversity\n[RGB color space]", fontface = 2, size = 7) +
  
  ### Add butterfly/mimicry symbol (to oppose to climate symbol)
  annotation_custom(grob = Climate_symbol, xmin = -3950, xmax = -2150, ymin = -050, ymax = 1900) +
  
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
  
  ### Add panel label
  annotate(geom = "text", x = 3900, y = -3800, hjust = 0.5,
           label = "B", fontface = 2, size = 10) +
  
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
pdf(file = "./supplementaries/Betadiv_maps/RGB_map_Climate_diversity.pdf", width = 8, height = 8)
print(ggplot_RGB_climate_map)
dev.off()

##### 6/ Plot both Mimicry turnover & Climate diversity ####

library(ggpubr)

pdf(file = "./supplementaries/Betadiv_maps/RGB_map_Mimicry_turnover_&_Climate_diversity.pdf", width = 16, height = 7)

ggarrange(ggplot_RGB_mimicry_map, ggplot_RGB_climate_map,                   # List of plots
          ncol = 2, nrow = 1)

dev.off()

###### Remove labels and set margins in individual and double plot !

plot.margin = margin(0.1,0.1,2,0.1, "cm")
