
##### Script 1: Plot Ithomiini species richness map #####

### Plot Ithomiini species richness map based on species distribution maps from Doré et al., 2021


#####

# Inputs
    # Species distribution maps as raster stack from Doré et al., 2021

# Outputs
    # Ithomiini species richness map 

#####


# Clean environment
rm(list = ls())

### 1/ Load stuff ####

# Packages
library(raster)
library(prettymapr)
library(rangeBuilder)

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./supplementaries/Richness_map/Map_stuff/pal_bl_red_Mannion.rds")

# Load mask for continent borders, plot border, and grid
grid_Mollweide_out <- readRDS(file = "./supplementaries/Richness_map/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./supplementaries/Richness_map/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./supplementaries/Richness_map/Map_stuff/bbox_sp_Mollweide.rds")

load(file = "./supplementaries/Richness_map/Map_stuff/country_borders.RData")
country_borders <- as(country_borders, "Spatial")

### Load map file
sp_richness <- readRDS(file = paste0("./supplementaries/Richness_map/Map_stuff/sp_richness.rds"))

### 2/ Project map stuff in Mollweide projection ####

# 2.1/ Projection functions
Mollweide_projection <- function(x) # Raster to project
{
  new_map <- projectRaster(from = x, 
                           method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                           crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                           alignOnly = F)
  return(new_map)
}

Mollweide_shp_projection <-  function(x) # Shp to project
{
  x_name <- deparse(substitute(x)) # Get the name of the initial shp as a character string
  
  new_shp <- spTransform(x, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
  
  return(new_shp)
}

# 2.2/ Project all map stuff

sp_richness_Mollweide <- Mollweide_projection(sp_richness)
country_borders_Mollweide <- Mollweide_shp_projection(country_borders)

### 3/ Plot final figure ####

# 3.1/ Plotting function ####

{
  map_indices_Mollweide <- function(x,                                    # Raster to map
                                  color_palette = pal_bl_red_Mannion,   # Color palette
                                  main_title,                           # Main title
                                  main_title_cex = 1.5,                 # Main title size
                                  
                                  xlim = c(-4600, 4600),   # Limit of plot on x-axis (Longitude)
                                  ylim = c(-4450, 3400),    # Limit of plot on y-axis (Latitude)
                                  axis_cex = 1.4,             # Axes size
                                  
                                  xlab = "",                # X-axis label
                                  ylab = "",                # Y-axis label
                                  x_axis_breaks = c(-3930, -2170, -420, 1500, 3050),            # X-axis tick breaks
                                  y_axis_breaks = c(-3650, -2450, -1220, 0, 1230, 2445),        # Y-axis tick breaks
                                  x_axis_labels = c("120°E", "100°E", "80°E", "60°E", "40°E"),      # X-axis tick labels
                                  y_axis_labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"),  # Y-axis tick labels
                                  
                                  legend_title,             # Legend title
                                  legend_title_cex = 1.4,   # Legend title size
                                  legend_title_x = -3550,   # Legend title x position
                                  legend_title_y = 430,     # Legend title y position
                                  legend_cex = 1.4,         # Legend size
                                  legend_breaks,            # Legend tick positions
                                  legend_location = c(-4100, -3800, -3950, 0),  # Legend position
                                  
                                  scale_bar_position = c(-2600, -4000),  # Scale bar position
                                  
                                  arrow_scale = 0.45,           # North arrow size
                                  arrow_padin = c(0.15, 0.15),  # North arrow position adjustement
                                  
                                  facet_letter = "",                  # Small case letter for facet
                                  facet_letter_col = "black",         # Color of case letter for facet
                                  facet_letter_cex = 2.2,             # Size of small case letter for facet
                                  facet_letter_inset = c(0, -0.008))  # Position adjustment of small case letter for facet

{
  # Plot raster background without axis
  image(x, col = color_palette,
        xlim = xlim, ylim = ylim, axes = F,
        xlab = xlab, ylab = ylab)
  title(main = main_title, cex.main = main_title_cex, line = 1)
  
  # Generate axes with manual positioning of ticks
  axis(side = 1, at = x_axis_breaks, labels = x_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1, gap.axis = 0, padj = 0.5)
  axis(side = 2, at = y_axis_breaks, labels = y_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1, gap.axis = 0)
  
  # Add background, borders and graticules
  plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
  plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
  plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
  plot(country_borders_Mollweide, lwd = 1, border = "#00000030", col = NA, add = T)
  
  # Add scale bar in legend
  scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = scale_bar_position, label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1.2)
  prettymapr::addnortharrow(scale = arrow_scale, padin = arrow_padin, text.col = "#00000000")
  rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
  rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
  graphics::text(x = legend_title_x, y = legend_title_y, font = 2, cex = legend_title_cex, label = legend_title)
  
  # Add facet letter
  legend(legend = facet_letter, x = "bottomright", bty = "n", text.col = facet_letter_col,
         text.font = 2, cex = facet_letter_cex, inset = facet_letter_inset)
  
  }
}


### 3.2/ Design custom color palette ####

# tmaptools::palette_explorer()

custom_pal <- rev(tmaptools::get_brewer_pal("RdYlBu", n = 200))
custom_pal <- custom_pal[c(seq(2, 120, 2), 122:200)]
custom_pal[1] <- pal_bl_red_Mannion[1]

### 3.3/ Exporting plot ####

pdf(file = paste0("./supplementaries/Richness_map/Species_richness_map.pdf"), height = 6, width = 6)

internal_margins <- par()$mar
par(mar = c(3.1, 3.1, 2.7, 1.6))

map_indices_Mollweide(x = sp_richness_Mollweide,
                      main_title = "Species richness",
                      color_palette = custom_pal,
                      legend_title = "Species",
                      legend_title_x = -3550,
                      legend_title_y = 430,
                      legend_breaks = seq(0, 120, 20), 
                      facet_letter = "")

par(mar = internal_margins)

dev.off()


