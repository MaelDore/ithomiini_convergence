##### Script 01: Community composition analysis with Ist #####

# Compute Ist for all communities
# Test signifiance through randomisation of mimicry pattern
# Plot test results

### Input files

# Summary table of OMUs
# Stack of OMUs probailities of presence
# Mimicry richness stack for all rings

### Output files

# Simpson's Mimicry Diversity Map
# Ist values for observed data and randomizations
# Plot of Ist null distribution and test



### Preparation ###

# Effacer l'environnement
rm(list = ls())

library(raster)
library(rangeBuilder)

### 1/ Load functions and files ####

# Function to aggregate probabilities among several OMUs of the same species
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x)
  return(y) 
}

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

# From Dore et al., 2021 - Anthropogenic pressures coincide with Neotropical biodiversity hotspots in a flagship butterfly group ###

# Load summary table for OMU list
load(file = paste0("./input_data/Summary_tables/list.unit.RData"))

# Load the OMU probability stack
OMU_proba_stack <- readRDS(file = paste0("./input_data/SDM_stacks/All_OMU_stack_Jaccard.80.rds"))
# Load the mimicry ring richness stack
ring_richness_stack <- readRDS(file = paste0("./input_data/SDM_stacks/All_ring_rich_stack_Jaccard.80.rds"))

# Color palet for plot
pal_bl_red_Mannion <- readRDS(file = "./input_data/Map_stuff/pal_bl_red_Mannion.rds")

# Load mask for continent borders
continent_mask <- readRDS(file = "./input_data/Map_stuff/continent_mask_15.rds")
crop_mask_shp <- readRDS(file = "./input_data/Map_stuff/crop_mask_shp_15.rds")
# bg_mask_pixel <- readRDS(file = "./input_data/Map_stuff/bg_mask_pixel.rds")
bg_mask <- readRDS(file = "./input_data/Map_stuff/bg_mask.rds")

##### 2/ Compute observed Ist value for whole communities #####

### 2.1/ Shift to raster brick to extract data under a single df for all layers
ring_richness_brick <- ring_richness_stack*1

### 2.2/ Compute Community diversity index/map = Dk
Dk <- calc(ring_richness_stack, fun = simpson, na.rm = T)*1
saveRDS(Dk, file = "./outputs/Community_Structure/Dk_map.rds")

### Plot Simpson's Mimicry Diversity Map
Dk <- readRDS(file = "./outputs/Community_Structure/Dk_map.rds")
Dk_map <- merge(Dk, continent_mask)
# Add 0 values in terrestrial pixel with no Ithomiini

pdf(file = paste0("./maps/Community_structure/Dk.pdf"), height = 6.3, width = 6.5)

image(Dk_map, col = pal_bl_red_Mannion, main = "Simpson's Mimicry Diversity Map", 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
addRasterLegend(Dk_map, locs = seq(0, 0.9, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.1, label = "Simpson's\n diversity")

dev.off()


### 2.3/ Compute Mean community diversity = Ds = 0.762
Global_Ds <- mean(Dk@data@values, na.rm = T)

### 2.4/ Compute Total diversity index = Dt = for all communities merged in one = 0.911
mimicry.abundances <- NA # To store each mimicry global expected richness (i.e., nb of species of this ring among all communities, duplicates included)
for (i in 1:nlayers(ring_richness_stack)) {
  mimicry.abundances[i] <- sum(ring_richness_stack[[i]]@data@values, na.rm = T)
}
# hist(mimicry.abundances)
Global_Dt <- simpson(mimicry.abundances)


### 2.5/ Compute Global Ist = 0.164
Global_Ist <- 1-(Global_Ds/Global_Dt)

# Save stuff
save(Global_Ds, Global_Dt, Global_Ist, file = "./outputs/Community_Structure/Global.Ist.RData")


##### 3/ Generate new virtual community matrices under null hypothesis of no effect of mimicry ring on species presence ####
# analogous to DeVries et al.'s [1999] and Hill's [2010] test of mimicry structure across microhabitats

mimicry.list <- as.character(unique(list.unit$Mimicry.model)) # 44 Mimicry rings

## Start the loop to shuffle mimicry patterns among OMUs
Ist_null <- NA
for (k in 1:999) # 999 permutations
{ 
  # k <- 1
  
  ## Shuffle mimicry ring among units
  shuffle.list.unit <- list.unit
  shuffle.list.unit$Mimicry.model <- sample(as.character(shuffle.list.unit$Mimicry.model))
  
  ## Generate the new Mimicry rings Richness Stack with random attribution of species to mimicry ring
  
  simul.ring_richness_stack <- Dk # 1st layer of the stack to remove later
  for (i in 1:length(mimicry.list)) # Per mimicry ring
  { 
    # i <- 1
    
    ring <- as.character(mimicry.list[i]) # Load mimicry ring name
    
    temp_ring_stack <- subset(x = OMU_proba_stack, subset = which(shuffle.list.unit$Mimicry.model == ring), drop = F)
    
    simul_ring_richness <- calc(temp_ring_stack, fun = sum)*1
    
    simul.ring_richness_stack <- addLayer(simul.ring_richness_stack, simul_ring_richness)
    
    # print(i)
    
  }
  simul.ring_richness_stack <- dropLayer(simul.ring_richness_stack, i = 1)*1
  names(simul.ring_richness_stack) <- mimicry.list
  
  save(simul.ring_richness_stack, file = paste0("./outputs/Community_Structure/Permutations/Ist/Ring_Richness_Stack_simul_",k,".RData"))
  
  # plot(ring_richness_stack)        # Original  
  # plot(simul.ring_richness_stack)  # Simulated
  
  
  ## Compute Community diversity index/map = Dk
  simul.Dk <- calc(simul.ring_richness_stack, fun = simpson, na.rm = T)*1
  # Plot Simulated Simpson's Mimicry Diversity Map
  # plot(simul.Dk, col = pal_bl_red_Mannion, main = paste0("Simulated Simpson's Mimicry Diversity Map n°", k))
  
  ## Compute index for all communities
  
  # Mean community diversity = Ds
  Simul_global_Ds <- mean(simul.Dk@data@values, na.rm = T)
  
  # Total diversity index = Dt
  mimicry.abundances <- NA # To store each mimicry global expected richness
  for (i in 1:nlayers(simul.ring_richness_stack)) 
  {
    mimicry.abundances[i] <- sum(simul.ring_richness_stack[[i]]@data@values, na.rm = T)
  }
  # hist(mimicry.abundances)
  Simul_global_Dt <- simpson(mimicry.abundances)
  
  # Global Ist
  Simul_global_Ist <- 1-(Simul_global_Ds/Simul_global_Dt)
  
  ## Save stuff
  Ist_null[k] <- Simul_global_Ist
  
  save(Ist_null, file = "./outputs/Community_Structure/Ist_null.RData")
  saveRDS(Ist_null, file = "./outputs/Community_Structure/Ist_null.rds")
  
  cat(paste0(Sys.time(), " - Simul n°", k,"\n"))
}

summary(Ist_null)

##### 4/ Plot distri of global null Ist ####

load(file = "./outputs/Community_Structure/Ist_null.RData")
load(file = "./outputs//Community_Structure/Global.Ist.RData")

# Final plot

pdf(file = "./graphs/Community_Structure/Global_Ist_null.pdf", height = 5.3, width = 6.5)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))

hist(c(Ist_null, Global_Ist), 50, freq = TRUE, col = "gray", xlab = bquote(~I[ST]), 
     main = bquote('Distribution of' ~I[ST]~ '\nunder null Hypothesis'),
     cex.axis = 1.3, cex.lab = 1.4, cex.main = 1.5, lwd = 2)
arrows(Global_Ist - 0.00055, 75, Global_Ist- 0.00055, 10, length = 0.1, lwd = 2) # Draw arrow above Ist obs
abline(v = mean(c(Ist_null, Global_Ist)), lty = 2, lwd = 2) # Add vertical line for mean value
abline(v = quantile(c(Ist_null, Global_Ist), 0.95), lty = 2, lwd = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("   Mean = ", format(round(mean(c(Ist_null, Global_Ist)),3), nsmall = 3)), 
                  paste0("CI 95% = ", round(quantile(c(Ist_null, Global_Ist), 0.95),3))), 
       x = "topright", inset = c(0, 0.17), xjust = 1, # Use inset to manually adjust position
       cex = 1.2, # text.width = c(0.9, 1),
       lty = 2, lwd = 2, col = c("black", "red"), bty ="n")
# legend(legend = c(paste0("          Ist obs = ", round(Global_Ist, 3)),
#                   paste0("                   p = 0.001")),
#        x = "right", cex = 1.2, bty ="n", xjust = 1)
legend(legend = bquote("           " ~ I[ST] ~ 'obs =' ~ .(round(Global_Ist, 3))),
       x = "right", cex = 1.2, bty ="n", xjust = 1)
legend(legend = c("",
                  paste0("                   p = 0.001")),
       x = "right", cex = 1.2, bty ="n", xjust = 1, y.intersp = 1.8)

legend(legend = as.expression(bquote(bold("A"))), 
       x = "topright", inset = c(0.001, 0.001), xjust = 0.5,
       cex = 1.3, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()



