##### Community composition analysis #####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

### Preparation ###

# Effacer l'environnement
rm(list = ls())

library(raster)

# Set Wd
internal.wd <- "F:/Documents/Etudes/Fac/Cours/TROPIMUNDO/Stage_S4/Projet_Papillons/R_codes"
initial.wd <-  "F:/Documents/Etudes/Fac/Cours/TROPIMUNDO/Modeling/"

# Chose the resolution
res <-  "5"

modeling.wd <- paste0(initial.wd,res,"/")
setwd(modeling.wd)
load(file = paste0(initial.wd,"/Models_Resume_table.",res,".RData"))
load(file = paste0(internal.wd,"/list.sp.RData"))
load(file = paste0(internal.wd,"/list.models.RData"))

# Load the mimicry ring richness stack
load(file = paste0(internal.wd,"/Com.ring.richness.stack.RData"))
# Load the species proba. stack
load(file = paste0(internal.wd,"/Com.sp.stack.RData"))

# Chargement d'un stack pour les N?otropiques ? la r?solution des Final Maps
load(file = paste0(internal.wd,"/envCom.RData"))

# Color palet for plot
cool = rainbow(50, s = 1, v = 1, start=rgb2hsv(col2rgb('yellow'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, s = 1, v= 1, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
pal_bl_red  = c(rev(cool), rev(warm))



##### Compute observed Ist value

### Load directly the mimicry ring richness stack (each layer = expected local number of species for a mimicry ring)
load(file = paste0(internal.wd,"/Com.ring.richness.stack.RData"))
unique.df.ring.richness.stack <- ring.richness.stack*1
save(unique.df.ring.richness.stack, file = paste0(internal.wd,"/Unique.df.com.ring.richness.stack.RData"))

# Fonction pour calculer l'indice de simpson car vegan fait de la merde !!!
simpson <- function (x, na.rm) {
  if (any(is.na(x))) { # Cas avec une valeur ? NA
    y <- NA
  }else{ # Cas avec toutes valeurs renseign?es
    y <- x/sum(x)
    y <- 1-(sum(y*y))
  }
  return(y)
}


### Index computation

## Community diversity index/map = Dk
library(vegan)
# Dk <- calc(ring.richness.stack, fun = vegan::diversity, index = "simpson", MARGIN = 1)*1
Dk <- calc(ring.richness.stack, fun = simpson, na.rm = T)*1
save(Dk, file = paste0(internal.wd,"/Community_Structure/Dk_map.RData"))


# Plot Simpson's Mimicry Diversity Map
load(file = paste0(internal.wd,"/Community_Structure/Dk_map.RData"))
jpeg(filename = paste0(internal.wd,"/Maps/Indices/Simpson_Mimicry_Diversity.jpeg"), quality =100)
plot(Dk, col = pal_bl_red, main = "Simpson's Mimicry Diversity Map")
plot(crop_mask_shp, add = T)
dev.off()

## Compute index for all communities

# Mean community diversity = Ds
Global_Ds <- mean(Dk@data@values, na.rm = T)

# Total diversity index = Dt  => ? recalculer sur la somme des communaut?s ?chantillonn?es en cas de sub-sampling
mimicry.abundances <- NA # To store each mimicry global expected richness
for (i in 1:nlayers(ring.richness.stack)) {
  mimicry.abundances[i] <- sum(ring.richness.stack[[i]]@data@values, na.rm = T)
}
# hist(mimicry.abundances)
# Dt <- vegan::diversity(x = as.matrix(mimicry.abundances), index = "simpson", MARGIN = 2)
Global_Dt <- simpson(mimicry.abundances)

# Global Ist
Global_Ist <- 1-(Global_Ds/Global_Dt)

save(Global_Ds, Global_Dt, Global_Ist, file = paste0(internal.wd,"/Community_Structure/Global.Ist.RData"))

## Subsampling global
# Generate distribution of "observed Ist"
Ist_obs <- NA
subsampling <- 500
for (i in 1:1000) {
  # Subsample X communities to look at spatially independant communities
  sample.index <- sample(x = which(!is.na(Dk_values)), size = subsampling, replace = F)
  sampled.Dk <- Dk_values[sample.index]
  Ds <- mean(sampled.Dk)
  
  sampled.Com <- unique.df.ring.richness.stack@data@values[sample.index,]
  mimicry.abundances <- apply(X = sampled.Com, MARGIN = 2, FUN = sum)
  Dt <- vegan::diversity(x = as.matrix(mimicry.abundances), index = "simpson", MARGIN = 2)
  
  ?diversity
  
  Ist <- 1-(Ds/Dt)
  Ist_obs[i] <- Ist
  # Show i every 10 iterations
  if (i%%10==0) {
    print(i)
  }
}

summary(Ist_obs)


# Plot distri of observed Ist
jpeg(filename = paste0(internal.wd,"/Community_Structure/Plots/Ist_obs_for_Sub_",subsampling,".jpeg"), quality =100)
hist(Ist_obs, xlab = "Ist obs", main = paste0("Distribution of observed Ist \n following subsampling to n = ", subsampling))
abline(v = mean(Ist_obs, na.rm = T), col = "red", lty = 2, lwd = 2)
legend(legend = c(paste0("Mean = ", round(mean(Ist_obs, na.rm = T),3)), 
                  paste0("CI 5% = ", round(quantile(Ist_obs, 0.05),3)),
                  paste0("CI 95% = ", round(quantile(Ist_obs, 0.95),3))),
       x = "topright", cex = 1, bty ="n") 
dev.off()

# Save
save(Ist_obs, file = paste0(internal.wd,"/Community_Structure/Ist_obs_",subsampling,".RData"))

### Load directly 
load(file = paste0(internal.wd,"/Community_Structure/Ist_obs_",subsampling,".RData"))


##### Generate new virtual community matrices under null hypothesis of no effect of mimicry ring on species presence
# analogous to DeVries et al.'s [1999] and Hill's [2010] test of mimicry structure across microhabitats

# Retrieve the community matrix for species
sp.stack.single.df <- sp.stack/1000
com.mat.sp.388 <- sp.stack.single.df@data@values

# Filter out communities with NA
real.com.index <- NA
for (i in 1:nrow(com.mat.sp.388)) {
  com.row <- com.mat.sp.388[i,]
  real.com.index[i] <- !any(is.na(com.row))
}
com.mat.sp.388 <- com.mat.sp.388[real.com.index,]
nrow(com.mat.sp.388) # 24986 community with no NA

save(com.mat.sp.388, file = paste0(internal.wd,"/Filtered.Com.mat.sp.388.RData"))
# load(file = paste0(internal.wd,"/Filtered.Com.mat.sp.388.RData"))

## Start the loop
Ist_null <- NA
for (k in 1:1000) {
  
  ## Suffle mimicry ring among units
  shuffle.list.unit <- list.models
  shuffle.list.unit$Mimicry.model <- sample(as.character(shuffle.list.unit$Mimicry.model))
  
  ## Generate the new Mimicry rings Richness Stack with random attribution of species to mimicry ring
  
  # Mimicry list
  mimicry.list <- as.character(unique(shuffle.list.unit$Mimicry.model)) # 44 Mimicry rings
  
  # Units avec Niche Model ou pas
  modeled.shuffle.list.unit <- shuffle.list.unit[which(Models_Resume_table$Models_select_0.5 >= 5),]
  occ.shuffle.list.unit <- shuffle.list.unit[-(which(Models_Resume_table$Models_select_0.5 >= 5)),]
  
  # Fonction pour calculer les probabilit?s ? plus grande ?chelle
  aggreg_prob = function(x, na.rm) { # D?claration des arguments en input
    y <- 1-prod(1-x) # D?finition du calcul interne
    return(y) # Output
  }
  
  # i <- 1
  simul.ring.richness.stack <- stack(envCom[[1]]) # 1e couche du stack final ? supprimer par la suite
  for (i in 1:length(mimicry.list)) { # Par mimicry ring
    
    ring.stack <- stack(envCom[[1]]) # 1e couche du stack temporaire des units ? supprimer par la suite
    
    ring <- as.character(mimicry.list[i]) # Charge le nom de sp
    
    Models <- modeled.shuffle.list.unit[modeled.shuffle.list.unit$Mimicry.model == ring,]
    if (nrow(Models)>0) {
      for (j in 1:nrow(Models)) {
        unit <-  as.character(Models$Tag.model[j])
        load(file = paste0(modeling.wd,unit,"/Com_map_quanti.RData"))
        ring.stack <- addLayer(ring.stack, EM.raster.com)
      }
    }
    
    Binaries <- occ.shuffle.list.unit[occ.shuffle.list.unit$Mimicry.model == ring,]
    if (nrow(Binaries)>0) {
      for (j in 1:nrow(Binaries)) {
        unit <-  as.character(Binaries$Tag.model[j])
        load(file = paste0(modeling.wd,unit,"/Com_map_Occ.RData"))
        ring.stack <- addLayer(ring.stack, Com_map_Occ)
      }
    }
    
    ring.stack <- dropLayer(ring.stack, i = 1)/1000  
    
    ring.richness.Com.map <- calc(ring.stack, fun = sum)*1
    
    simul.ring.richness.stack <- addLayer(simul.ring.richness.stack, ring.richness.Com.map)
    
    # print(i)
    
  }
  simul.ring.richness.stack <- dropLayer(simul.ring.richness.stack, i = 1)*1
  names(simul.ring.richness.stack) <- mimicry.list
  
  save(simul.ring.richness.stack, file = paste0(internal.wd,"/Community_Structure/Ring_Richness_Stack_simul_",j,".RData"))
  
  # plot(simul.ring.richness.stack)
  # plot(ring.richness.stack)
  
  ## Community diversity index/map = Dk
  # Dk <- calc(ring.richness.stack, fun = vegan::diversity, index = "simpson", MARGIN = 1)*1
  simul.Dk <- calc(simul.ring.richness.stack, fun = simpson, na.rm = T)*1
  # Plot Simulated Simpson's Mimicry Diversity Map
  plot(simul.Dk, col = pal_bl_red, main = paste0("Simulated Simpson's Mimicry Diversity Map n?", k))
  
  ## Compute index for all communities
  
  # Mean community diversity = Ds
  Simul_global_Ds <- mean(simul.Dk@data@values, na.rm = T)
  
  # Total diversity index = Dt  => ? recalculer sur la somme des communaut?s ?chantillonn?es en cas de sub-sampling
  mimicry.abundances <- NA # To store each mimicry global expected richness
  for (i in 1:nlayers(simul.ring.richness.stack)) {
    mimicry.abundances[i] <- sum(simul.ring.richness.stack[[i]]@data@values, na.rm = T)
  }
  # hist(mimicry.abundances)
  # Dt <- vegan::diversity(x = as.matrix(mimicry.abundances), index = "simpson", MARGIN = 2)
  Simul_global_Dt <- simpson(mimicry.abundances)
  
  # Global Ist
  Simul_global_Ist <- 1-(Simul_global_Ds/Simul_global_Dt)
  
  Ist_null[k] <- Simul_global_Ist
  
  save(Ist_null, file = paste0(internal.wd,"/Community_Structure/Ist_null.RData"))
  
  # Subsampling ?
  
  cat(paste0(Sys.time(), " - Simul n?", k,"\n"))
}

summary(Ist_null)

# Plot distri of global null Ist
load(file = paste0(internal.wd,"/Community_Structure/Ist_null.RData"))
load(file = paste0(internal.wd,"/Community_Structure/Global.Ist.RData"))
jpeg(filename = paste0(internal.wd,"/Community_Structure/Plots/Global_Ist_null.jpeg"), quality =100)
hist(Ist_null, xlab = "Ist under H0", main = "Distribution of Ist \n under null Hypothesis", breaks = seq(from = 0.07, to = 0.16, by = 0.005))
abline(v = mean(Ist_null, na.rm = T), col = "grey20", lty = 2, lwd = 2)
legend(legend = c(paste0("Mean = ", round(mean(Ist_null, na.rm = T),3)), 
                  paste0("CI 5% = ", round(quantile(Ist_null, 0.05),3)),
                  paste0("CI 95% = ", round(quantile(Ist_null, 0.95),3))),
       x = "top", cex = 1, bty ="n") 
abline(v = Global_Ist, col = "red", lty = 2, lwd = 2)
text(x = Global_Ist-0.015, y = 70, labels = paste0("Ist obs = ", round(Global_Ist, 3)), col = "red")
text(x = Global_Ist-0.015, y = 50, labels = "p < 0.001", col = "red")
dev.off()


# Plot
load(file = paste0(internal.wd,"/Maps_exploitation/MPD_mean_values.RData"))
tiff(filename = paste0(internal.wd,"/Community_Structure/Plots/Global_Ist_null.tif"))
original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))
hist(c(Ist_null, Global_Ist), 30, freq = TRUE, col = "gray", xlab = "Ist", 
     main = "Distribution of Ist \n under null Hypothesis",
     cex.axis = 1.7, cex.lab = 1.7, cex.main = 1, axis.args=list(cex.axis=1.7))
arrows(Global_Ist + 0.0001, 60, Global_Ist+ 0.0001, 10, length = 0.1, lwd = 2)
abline(v = mean(c(Ist_null, Global_Ist)), lty = 2)
abline(v = quantile(c(Ist_null, Global_Ist), 0.95), lty = 2, col = "red")
legend(legend = c(paste0("Mean = ", round(mean(c(Ist_null, Global_Ist)),3)), 
                  paste0("CI 95% = ", round(quantile(c(Ist_null, Global_Ist), 0.95),3))), 
       x = "topright", lty = 2 , lwd = 1.5, col = c("black", "red"), cex = 1.2, bty ="n")
legend(legend = c(paste0("          Ist obs = ", round(Global_Ist, 3)),
                  paste0("                   p = 0.001")),
       x = "right", cex = 1.2, bty ="n", xjust = 1)
par(oma = original_ext_margins, mar = original_int_margins)
dev.off()


# Plot distri of null Ist with subsampling
jpeg(filename = paste0(internal.wd,"/Community_Structure/Ist_null_for_Sub_",subsampling,".jpeg"), quality =100)
hist(Ist_null, xlab = "Ist under H0", main = paste0("Distribution of Ist \n under null Hypothesis \n following subsampling to n = ", subsampling))
abline(v = mean(Ist_null, na.rm = T), col = "red", lty = 2, lwd = 2)
legend(legend = c(paste0("Mean = ", round(mean(Ist_null, na.rm = T),3)), 
                  paste0("CI 5% = ", round(quantile(Ist_null, 0.05),3)),
                  paste0("CI 95% = ", round(quantile(Ist_null, 0.95),3))),
       x = "topright", cex = 1, bty ="n") 
dev.off()


### Final test

# Build df for t-test and boxplot
Ist <- c(Ist_obs, Ist_null)
Distri <- as.factor(c(rep("Obs",length(Ist_obs)),rep("Null",length(Ist_null))))
Ist.df <- data.frame(Ist, Distri)

## Hypoth?ses T-test impossible. On passe en non-param?trique
# Hypoth?se de normalit?
shapiro.test(Ist_obs)
shapiro.test(Ist_null)
# Hypoth?se d'homosc?dascticit?
var.test(Ist_obs,Ist_null,alt="two.sided") 
##  T-test impossible. On passe en non-param?trique.

# Wilcox test
Ist_test <- wilcox.test(Ist_obs, Ist_null, alt = "greater", paired = F) ; Ist_test
# Significatif !

# Boxplot
jpeg(filename = paste0(internal.wd,"/Community_Structure/Boxplot_sub_",subsampling,".jpeg"), quality =100)
boxplot(formula=Ist.df$Ist~Ist.df$Distri, col = rainbow(2), ylim = c(-1,1), main = paste0("Comparison of observed and null distribution of Ist \n for subsampling = ", subsampling))
legend(legend = c("W = 840020","p-value < 0.001"), x = "bottomright", cex = 1, bty ="o")
dev.off()


##### Approche par Bray-Curtis sur les paires d'units co-mim?tiques

### Load directly the unit probability stack
load(file = paste0(internal.wd,"/Com.unit.stack.RData"))
# unique.df.unit.stack <- unit.stack*1
com.unit.mat.783 <- NA
for (i in 1:nlayers(unit.stack)) {
  if (i==1){
    com.unit.mat.783 <- unit.stack[[i]]@data@values
  }else{
    unit.values <- unit.stack[[i]]@data@values
    com.unit.mat.783 <- cbind(com.unit.mat.783, unit.values)
  }
  print(i)
}
names(com.unit.mat.783) <- names(unit.stack)
rm(unit.stack)

save(com.unit.mat.783, file = paste0(internal.wd,"/full.unit.mat.RData"))

# Filter out communities with NA
real.com.index <- NA
for (i in 1:nrow(com.unit.mat.783)) {
  com.row <- com.unit.mat.783[i,]
  real.com.index[i] <- (!any(is.na(com.row)))&(!any(is.nan(com.row)))
}
com.unit.mat.783 <- com.unit.mat.783[real.com.index,]
nrow(com.unit.mat.783) # 24986 community with no NA

save(com.unit.mat.783, file = paste0(internal.wd,"/filtered.unit.mat.RData"))

# Load directly the matrix of community * units probability presence
load(file = paste0(internal.wd,"/filtered.unit.mat.RData"))

# Compute Bray-Curtis index for all pairs of units
library(vegan)

unit.BC.dist <- vegdist(x = t(com.unit.mat.783), method = "bray") # Calcul de dissimilarit?s entre lignes

save(unit.BC.dist, file = paste0(internal.wd,"/Community_Structure/unit.BC.dist.RData"))

# Load directly the vector of BD distances among all units
load(file = paste0(internal.wd,"/Community_Structure/unit.BC.dist.RData"))
load(file = paste0(internal.wd,"/Com.unit.stack.RData"))

hist(unit.BC.dist)

unit.BC.dist.mat <- as.matrix(unit.BC.dist)

mimicry.list <- as.character(unique(list.models$Mimicry.model))

mean_BC <- NA
for (i in 1:length(mimicry.list)) { # Par cercles mim?tiques
  mimic <- mimicry.list[i]
  for (j in 1:nrow(list.models)) { # Recherche de l'index des units
    tags <- as.character(list.models$Tag.model[list.models$Mimicry.model == mimic])
    index <- which(names(unit.stack) %in% tags)
  }
  if (length(tags)==1) { # Cas avec une seule unit dans le cercle. Impossible de calculer une distance moyenne entre paires.
    mat.index <- data.frame(matrix(nrow = 2, ncol = 0)) # Empty df of pairs coordinates
    save(mat.index, file = paste0(internal.wd,"/Community_Structure/mat.index.",mimic,".RData"))
    BC <- c() # Empty vector of Bray-Curtis values
    save(BC, file = paste0(internal.wd,"/Community_Structure/BC.",mimic,".RData"))
    mean_BC[i] <- NA
  }else{ # Cas normal avec au moins 2 units.
    mat.index <- combn(x = index, m = 2, FUN = c) # R?cup?ration de toutes les combinaisons de paires de units possibles
    save(mat.index, file = paste0(internal.wd,"/Community_Structure/mat.index.",mimic,".RData"))
    BC <- NA
    for(j in 1:ncol(mat.index)){
      BC[j] <- unit.BC.dist.mat[mat.index[1,j],mat.index[2,j]] # Extraction des Bray-Curtis index
    }
    save(BC, file = paste0(internal.wd,"/Community_Structure/BC.",mimic,".RData"))
    mean_BC[i]  <- mean(BC) # Calcule la moyenne des BC pour toutes les paires de co-mimics pour ce cercle
  }
  
  print(i)
}
names(mean_BC) <- mimicry.list
save(mean_BC, file  = paste0(internal.wd,"/Community_Structure/mean_BC.RData"))

# Load mean BC of mimicry rings
load(file  = paste0(internal.wd,"/Community_Structure/mean_BC.RData"))
mean_BC

# Retrive all mimic coordinates
all_mimic_mat.index <- data.frame(matrix(ncol = 0, nrow=2))
for (i in 1:length(mimicry.list)) { # Par cercles mim?tiques
  mimic <- mimicry.list[i]
  load(file = paste0(internal.wd,"/Community_Structure/mat.index.",mimic,".RData"))
  all_mimic_mat.index <- cbind(all_mimic_mat.index,mat.index)
}
dim(all_mimic_mat.index) # 14 659 pairs of mimic units
# Retrive all mimic values
BC_mimic <- NA ; unit.BC.dist.mat.no.mimic <- unit.BC.dist.mat
for(j in 1:ncol(all_mimic_mat.index)){
  BC_mimic[j] <- unit.BC.dist.mat[all_mimic_mat.index[1,j],all_mimic_mat.index[2,j]] # Extraction des Bray-Curtis index des mimics
  unit.BC.dist.mat.no.mimic[all_mimic_mat.index[1,j],all_mimic_mat.index[2,j]] <- NA # Ecrasement des BC des mimics
}
length(BC_mimic) # 14 659 pairs of mimic units
save(BC_mimic, file = paste0(internal.wd,"/Community_Structure/BC_mimic.RData"))

# Retrieve all non-mimic values
BC_no.mimic <- unit.BC.dist.mat.no.mimic[upper.tri(unit.BC.dist.mat.no.mimic)]
BC_no.mimic <- na.omit(BC_no.mimic)
length(BC_no.mimic) # 291 494 pairs of non-mimic units
save(BC_no.mimic, file = paste0(internal.wd,"/Community_Structure/BC_no.mimic.RData"))

# Global mean computation
load(file = paste0(internal.wd,"/Community_Structure/unit.BC.dist.RData")) # All BC index
load(file = paste0(internal.wd,"/Community_Structure/BC_mimic.RData")) # Only mimic pairs
load(file = paste0(internal.wd,"/Community_Structure/BC_no.mimic.RData")) # Only non mimic pairs

Global_mean_BC <- mean(unit.BC.dist) ; Global_mean_BC # 0.961
Global_mimic_mean_BC <- mean(BC_mimic) ; Global_mimic_mean_BC # 0.913
Global_no.mimic_mean_BC <- mean(BC_no.mimic) ; Global_no.mimic_mean_BC # 0.963

save(Global_mean_BC, Global_mimic_mean_BC, Global_no.mimic_mean_BC, file = paste0(internal.wd,"/Community_Structure/All_Global_BC.RData"))

# Boxplots plut?t moches car beaucoup d'outliers
All_BC <- c(BC_mimic, BC_no.mimic)
Status <- as.factor(c(rep("Mimic", length(BC_mimic)), rep("Non-Mimic", length(BC_no.mimic))))
boxplot(All_BC~Status)

# Histogrammes
# Plot distri of BC for mimic pairs
jpeg(filename = paste0(internal.wd,"/Community_Structure/hist_BC_mimic.jpeg"), quality =100)
hist(BC_mimic, xlab = "Bray-Curtis index", main = "Bray-Curtis index of mimic pairs")
abline(v = mean(BC_mimic), col = "red", lty = 2, lwd = 2)
legend(legend = c(paste0("Mean = ", round(mean(BC_mimic, na.rm = T),3)), 
                  paste0("CI 5% = ", round(quantile(BC_mimic, 0.05),3)),
                  paste0("CI 95% = ", round(quantile(BC_mimic, 0.95),3))),
       x = "topleft", cex = 1, bty ="n") 
dev.off()


##### Generate new virtual community matrices under null hypothesis of no effect of mimicry ring on species presence
# analogous to DeVries et al.'s [1999] and Hill's [2010] test of mimicry structure across microhabitats

# Load directly the vector of BD distances among all units
load(file = paste0(internal.wd,"/Community_Structure/unit.BC.dist.RData"))
load(file = paste0(internal.wd,"/Com.unit.stack.RData"))

unit.BC.dist.mat <- as.matrix(unit.BC.dist)

mimicry.list <- as.character(unique(list.models$Mimicry.model))


## Start the loop
BC_mimic_null <- BC_no.mimic_null <- NA
mean_BC_null <- matrix(ncol=44, nrow =0)
for (k in 36:1000) { # 1000 simulations
  
  ## Suffle mimicry ring among units
  shuffle.list.unit <- list.models
  shuffle.list.unit$Mimicry.model <- sample(as.character(shuffle.list.unit$Mimicry.model))
  
  table(shuffle.list.unit$Mimicry.model)
  
  ## Generate the new Mimicry rings Richness Stack with random attribution of species to mimicry ring
  
  # Mimicry list
  mimicry.list <- as.character(unique(list.models$Mimicry.model)) # 44 Mimicry rings
  
  mean_BC <- NA
  for (i in 1:length(mimicry.list)) { # Par cercles mim?tiques
    mimic <- mimicry.list[i]
    for (j in 1:nrow(list.models)) { # Recherche de l'index des units originelles
      tags <- as.character(list.models$Tag.model[shuffle.list.unit$Mimicry.model == mimic])
      index <- which(names(unit.stack) %in% tags)
    }
    if (length(tags)==1) { # Cas avec une seule unit dans le cercle. Impossible de calculer une distance moyenne entre paires.
      mat.index <- data.frame(matrix(nrow = 2, ncol = 0)) # Empty df of pairs coordinates
      save(mat.index, file = paste0(internal.wd,"/Community_Structure/Simuls/mat.index.",mimic,".RData"))
      mean_BC[i] <- NA
    }else{ # Cas normal avec au moins 2 units.
      mat.index <- combn(x = index, m = 2, FUN = c) # R?cup?ration de toutes les combinaisons de paires de units possibles
      save(mat.index, file = paste0(internal.wd,"/Community_Structure/Simuls/mat.index.",mimic,".RData"))
      BC <- NA
      for(j in 1:ncol(mat.index)){
        BC[j] <- unit.BC.dist.mat[mat.index[1,j],mat.index[2,j]] # Extraction des Bray-Curtis index
      }
      mean_BC[i]  <- mean(BC) # Calcule la moyenne des BC pour toutes les paires de co-mimics pour ce cercle
    }
    
    print(i)
  }
  names(mean_BC) <- mimicry.list
  mean_BC_null <- rbind(mean_BC_null,mean_BC)
  
  # Retrive all mimic coordinates
  all_mimic_mat.index <- data.frame(matrix(ncol = 0, nrow=2))
  for (i in 1:length(mimicry.list)) { # Par cercles mim?tiques
    mimic <- mimicry.list[i]
    load(file = paste0(internal.wd,"/Community_Structure/Simuls/mat.index.",mimic,".RData"))
    all_mimic_mat.index <- cbind(all_mimic_mat.index,mat.index)
  }
  # dim(all_mimic_mat.index) # 14 659 pairs of mimic units
  # Retrive all mimic values
  BC_mimic <- NA ; unit.BC.dist.mat.no.mimic <- unit.BC.dist.mat
  for(j in 1:ncol(all_mimic_mat.index)){
    BC_mimic[j] <- unit.BC.dist.mat[all_mimic_mat.index[1,j],all_mimic_mat.index[2,j]] # Extraction des Bray-Curtis index des mimics
    unit.BC.dist.mat.no.mimic[all_mimic_mat.index[1,j],all_mimic_mat.index[2,j]] <- NA # Ecrasement des BC des mimics
  }
  # length(BC_mimic) # 14 659 pairs of mimic units
  
  # Retrieve all non-mimic values
  BC_no.mimic <- unit.BC.dist.mat.no.mimic[upper.tri(unit.BC.dist.mat.no.mimic)]
  BC_no.mimic <- na.omit(BC_no.mimic)
  # length(BC_no.mimic) # 291 494 pairs of non-mimic units
  
  # Global mean computation
  BC_mimic_null[k] <- mean(BC_mimic)
  BC_no.mimic_null[k] <- mean(BC_no.mimic)
  
  cat(paste0(Sys.time(), " - Simul n?", k,"\n"))
  save(mean_BC_null, BC_mimic_null, BC_no.mimic_null, file = paste0(internal.wd,"/Community_Structure/Simuls/All_simul_BC.RData"))
}

summary(mean_BC_null)
summary(BC_mimic_null)
summary(BC_no.mimic_null)

load(file = paste0(internal.wd,"/Community_Structure/All_Global_BC.RData"))
load(file  = paste0(internal.wd,"/Community_Structure/mean_BC.RData"))

# Plot distri of null BC for co-mimics
load(file = paste0(internal.wd,"/Community_Structure/All_Global_BC.RData"))
load(file  = paste0(internal.wd,"/Community_Structure/mean_BC.RData"))
load(file = paste0(internal.wd,"/Community_Structure/Simuls/All_simul_BC.RData"))

# Plot
load(file = paste0(internal.wd,"/Community_Structure/All_Global_BC.RData"))
load(file  = paste0(internal.wd,"/Community_Structure/mean_BC.RData"))
tiff(filename = paste0(internal.wd,"/Community_Structure/Plots/BC_mimic_null.tif"))
original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(5.1,8,4.1,4))
hist(c(BC_mimic_null, Global_mimic_mean_BC), 40, freq = TRUE, col = "gray", xlab = "Ist", xlim = c(0.91, 0.97), ylim = c(0, 250),
     main = "Distribution of Mean pairwise BC Index \n of co-mimetic species \n under null Hypothesis",
     xlab = "Mean pairwise BC Index of co-mimetic species",
     cex.axis = 1.7, cex.lab = 1.7, cex.main = 1, axis.args=list(cex.axis=1.7))
arrows(Global_mimic_mean_BC + 0.0001, 60, Global_mimic_mean_BC+ 0.0001, 10, length = 0.1, lwd = 2)
abline(v = mean(c(BC_mimic_null, Global_mimic_mean_BC)), lty = 2)
abline(v = quantile(c(BC_mimic_null, Global_mimic_mean_BC), 0.05), lty = 2, col = "red")
legend(legend = c(paste0("Mean = ", round(mean(c(BC_mimic_null, Global_mimic_mean_BC)),3)), 
                  paste0("CI 5% = ", round(quantile(c(BC_mimic_null, Global_mimic_mean_BC), 0.95),3))), 
       x = "topleft", lty = 2 , lwd = 1.5, col = c("black", "red"), cex = 1.2, bty ="n")
legend(legend = c(paste0("Mean BC obs = ", round(Global_mimic_mean_BC, 3)),
                  paste0("p = 0.001")),
       x = "left", cex = 1.2, bty ="n", xjust = 1)
par(oma = original_ext_margins, mar = original_int_margins)
dev.off()

BC_ring_summary_table <- as.data.frame(matrix(ncol = 8, nrow = 44, data = NA))
names(BC_ring_summary_table) <- c("ring", "N_units", "N_pairs", "BC_obs", "mean_BC", "BC_2.5", "BC_97.5", "p_value")

# Plot distri of null BC per mimicry ring co-mimics
for (i in 1:length(mimicry.list)) {
  BC_ring_summary_table$ring[i] <- mimic <- mimicry.list[i]
  BC_ring_summary_table$N_units[i] <- N_units <- sum(unit.list$Mimicry.model==mimic)
  BC_ring_summary_table$N_pairs[i] <- N_pairs <- N_units*(N_units-1)/2
  if (is.na(mean_BC[i])) {
    jpeg(filename = paste0(internal.wd,"/Community_Structure/Plots/BC_mimic_null_",mimic,".jpeg"), quality =100)
    plot(1:100,1:100, type = "n", main = paste0("Distribution of mean pairwise Bray-Curtis index \n of ", mimic, " species \n under null Hypothesis"))
    text(x = 50, y = 50, labels = "Only one unit for this mimicry ring \n No pair available index for computation")
    dev.off()
  }else{
    BC_ring_summary_table$BC_obs[i] <- round(mean_BC[i],3)
    BC_ring_summary_table$mean_BC[i] <- round(mean(mean_BC_null[,i], na.rm = T),3)
    BC_ring_summary_table$BC_2.5[i] <- round(quantile(mean_BC_null[,i], 0.025),3)
    BC_ring_summary_table$BC_97.5[i] <- round(quantile(mean_BC_null[,i], 0.975),3)
    BC_ring_summary_table$p_value[i] <- round(ecdf(x = mean_BC_null[,i])(mean_BC[i]),3)
    jpeg(filename = paste0(internal.wd,"/Community_Structure/Plots/BC_mimic_null_",mimic,".jpeg"), quality =100)
    hist(mean_BC_null[,i], breaks = seq(from = floor(floor(min(min(mean_BC_null[,i]), mean_BC[i])*1000)/2)/500, to = ceiling(ceiling(max(max(mean_BC_null[,i]), mean_BC[i])*1000)/2)/500, by = 0.002), xlab = "Mean pairwise Bray-Curtis index", main = paste0("Distribution of mean pairwise Bray-Curtis index \n of ", mimic, " species \n under null Hypothesis"))
    abline(v = mean(mean_BC_null[,i], na.rm = T), col = "grey20", lty = 2, lwd = 2)
    legend(legend = c(paste0("N units = ", N_units), 
                      paste0("N pairs = ", N_pairs)),
           x = "topright", cex = 1, bty ="n")
    legend(legend = c(paste0("Mean = ", round(mean(mean_BC_null[,i], na.rm = T),3)), 
                      paste0("CI 2.5% = ", round(quantile(mean_BC_null[,i], 0.025),3)),
                      paste0("CI 97.5% = ", round(quantile(mean_BC_null[,i], 0.975),3))),
           x = "topleft", cex = 1, bty ="n") 
    abline(v = mean_BC[i], col = "red", lty = 2, lwd = 2)
    legend(legend = paste0("Mean BC obs = ", round(mean_BC[i], 3), "\n p = ", round(ecdf(x = mean_BC_null[,i])(mean_BC[i]),3), "\n \n \n \n \n \n \n \n"), text.col = "red", x = "bottomleft", cex = 1, bty ="n")
    dev.off()
    save(BC_ring_summary_table, file = paste0(internal.wd,"/Community_Structure/BC_ring_summary_table.Rdata"))
  }
  cat(paste0(i, " - ",mimic, " - Done \n"))
}

View(BC_ring_summary_table)

##### Approche Mantel pour pairwise Ist ~ Dclim

### Load directly the mimicry ring richness brick (each layer = expected local number of species for a mimicry ring)
load(file = paste0(internal.wd,"/Unique.df.com.ring.richness.stack.RData"))
load(file = paste0(internal.wd,"/Com.sp.richness.RData"))
load(file = paste0(internal.wd,"/envCom.RData"))

# Fonction pour calculer l'indice de simpson car vegan fait de la merde !!!
simpson <- function (x, na.rm) {
  if (any(is.na(x))) { # Cas avec une valeur ? NA
    y <- NA
  }else{ # Cas avec toutes valeurs renseign?es
    y <- x/sum(x)
    y <- 1-(sum(y*y))
  }
  return(y)
}

# Sous-?chantillonnage de N communaut?s
subsampling <- 1000
sample.index <- sample(x = which(tot.sp.richness@data@values>0), size = subsampling, replace = F)

# Subsample the richness data
sampled.Com <- unique.df.ring.richness.stack@data@values[sample.index,]

# Compute pairwise_Ist
pairwise_Ist <- matrix(ncol = subsampling, nrow = subsampling)
for (i in 1:subsampling) {
  for (j in 1:subsampling) {
    com1 <- sampled.Com[i,]
    com2 <- sampled.Com[j,]
    D1 <- simpson(com1)
    D2 <- simpson(com2)
    Ds <- mean(c(D1,D2))
    full_com <- apply(X = sampled.Com[c(i,j),], MARGIN = 2, FUN = sum)
    Dt <- simpson(full_com)
    pairwise_Ist[i,j] <- 1-(Ds/Dt)
  }
  # Show i every 10 iterations
  if (i%%10==0) {
    print(i)
  }
}

# Compute pairwise geographical distances
coord <- coordinates(tot.sp.richness)[sample.index,]
# plot(coord)
Pairwise_geodist = geosphere::distm(x = coord)/1000 # Dist Geo sur l'ellipsoid de r?f?rence WGS84, en km

# Compute pairwise euclidian climatic distances on standardized climatic variables
sampled_com_env <- envCom@data@values[sample.index,]
sampled_com_env <- scale(x = sampled_com_env, center = T, scale = T)
Pairwise_climdist <- as.matrix(dist(x = sampled_com_env, method = "euclidian"))

# Save all
save(pairwise_Ist, Pairwise_climdist, Pairwise_geodist, file = paste0(internal.wd, "/Community_Structure/Pairwise_values.RData"))
load(file = paste0(internal.wd, "/Community_Structure/Pairwise_values.RData"))

## Mantel classique
library("vegan")

# Ist ~ clim avec rho de Spearman
Mantel.Spearman_Ist_clim <- mantel(xdis = pairwise_Ist, ydis = Pairwise_climdist, method = "spearman", na.rm = T)
Mantel.Spearman_Ist_clim
save(Mantel.Spearman_Ist_clim , file = paste0(internal.wd, "/Community_Structure/Resultats_Mantel_Ist_clim_Spearman.RData"))
load(file = paste0(internal.wd, "/Community_Structure/Resultats_Mantel_Ist_clim_Spearman.RData"))

# Ist ~ geo avec rho de Spearman
Mantel.Spearman_Ist_geo <- mantel(xdis = pairwise_Ist, ydis = Pairwise_geodist, method = "spearman", na.rm = T)
Mantel.Spearman_Ist_geo
save(Mantel.Spearman_Ist_geo , file = paste0(internal.wd, "/Community_Structure/Resultats_Mantel_Ist_geo_Spearman.RData"))
load(file = paste0(internal.wd, "/Community_Structure/Resultats_Mantel_Ist_geo_Spearman.RData"))


## Mantel Partial Test : Ist ~ Dclim + cov(Dgeo)
library("vegan")

# V?rifie la normalit?
hist(Pairwise_climdist)
hist(Pairwise_geodist)
hist(pairwise_Ist)

# R?alise le test avec r de Pearson
Mantel.Pearson <- mantel.partial(xdis = pairwise_Ist, ydis = Pairwise_climdist, zdis = Pairwise_geodist, method = "pearson", na.rm = T) # Z = effet confondant ? retirer lors de la mesure de correlation partielle entre X et Y.
Mantel.Pearson
save(Mantel.Pearson, file = paste0(internal.wd, "/Community_Structure/Resultats_Mantel_partiel_Ist_Pearson.RData"))

# Plot de la distri de la stat
jpeg(filename = paste0(internal.wd,"/Community_Structure/Plots/Distri_Mantel_Pearson.jpeg"), quality =100)
hist(Mantel.Pearson$perm, xlab = "Pearson's r", breaks = seq(from = -0.06, to = 0.24, by = 0.02 ), main = "Mantel test statistic distribution\n under 999 permutations")
abline(v = Mantel.Pearson$statistic, lty = 2, lwd = 2, col = "red")
text(x = 0.18, y = 50, labels = paste0("r = ",round(Mantel.Pearson$statistic, 3), "\n p <0.001"), col = "red")
dev.off()

# R?alise le test avec rho de Spearman
Mantel.Spearman <- mantel.partial(xdis = pairwise_Ist, ydis = Pairwise_climdist, zdis = Pairwise_geodist, method = "spearman", na.rm = T) # Z = effet confondant ? retirer lors de la mesure de correlation partielle entre X et Y.
Mantel.Spearman
save(Mantel.Spearman, file = paste0(internal.wd, "/Community_Structure/Resultats_Mantel_partiel_Ist_Spearman.RData"))
load(file = paste0(internal.wd, "/Community_Structure/Resultats_Mantel_partiel_Ist_Spearman.RData"))


# Plot de la distri de la stat
load(file = paste0(internal.wd, "/Community_Structure/Resultats_Mantel_partiel_Ist_Spearman.RData"))
tiff(filename = paste0(internal.wd,"/Community_Structure/Plots/Distri_Mantel_Spearman.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
hist(c(Mantel.Spearman$perm, Mantel.Spearman$statistic), 30, freq = TRUE, col = "gray", 
     main = "Mantel test statistic distribution\n under 999 permutations", xlab = "Spearman's rho",
     cex.axis = 1.7, cex.main = 1, cex.lab = 1.7)
arrows(Mantel.Spearman$statistic + 0.0051, 80, Mantel.Spearman$statistic+ 0.0051, 10, length = 0.1, lwd = 2)
abline(v = mean(c(Mantel.Spearman$perm, Mantel.Spearman$statistic)), lty = 2, lwd = 2)
abline(v = quantile(c(Mantel.Spearman$perm, Mantel.Spearman$statistic), 0.95), lty = 2, lwd = 2, col = "red")
legend(legend = c(paste0("Mean = ", round(mean(c(Mantel.Spearman$perm, Mantel.Spearman$statistic)),3)), 
                  paste0("CI 95% = ", round(quantile(c(Mantel.Spearman$perm, Mantel.Spearman$statistic), 0.95),3))), 
       x = "topright", lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")
legend(legend = c(paste0("Rho obs = ", round(Mantel.Spearman$statistic, 3)),
                  paste0("                  p = ", 1-(ecdf(x = c(Mantel.Spearman$perm, Mantel.Spearman$statistic))(Mantel.Spearman$statistic)))),
       x = "right", cex = 1.2, bty ="n", xjust = 1)
par(mar = original_int_margins)
dev.off()


## Plot de la relation avec 1000 points tir?s au hasard parmi les n(n-1)/2 valeurs de distances

# Extraction des vecteurs
pairwise_Ist_vec <- as.dist(pairwise_Ist)[which(!is.na(as.dist(pairwise_Ist)))]
pairwise_climdist_vec <- as.dist(Pairwise_climdist)[which(!is.na(as.dist(pairwise_Ist)))]
pairwise_geodist_vec <- as.dist(Pairwise_geodist)[which(!is.na(as.dist(pairwise_Ist)))]

# Ist ~ Dclim
reg.clim <- lm(pairwise_Ist_vec~pairwise_climdist_vec)
summary(reg.clim)
random.sample <- sample(x= 1: length(pairwise_Ist_vec), size = 1000, replace = F)
tiff(filename = paste0(internal.wd,"/Community_Structure/Plots/Ist_Dclim.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
plot(pairwise_climdist_vec[random.sample], pairwise_Ist_vec[random.sample], 
     main = "Climate structures \n mimicry community composition", ylab = "Pairwise Ist", xlab = "Standardized climatic distance",
     cex.axis = 1.7, cex.main = 1, axis.args=list(cex.axis=1.7), cex.lab = 1.7)
abline(reg.clim, lwd = 2, col = "red")
# text(x = 7, y = 0.4, labels = paste0("Rho = ", round(Mantel.Spearman_Ist_clim$statistic,3), "\n p < 0.001 "), col = "red")
par(mar = original_int_margins)
dev.off()

# Ist ~ Dgeo
reg.geo <- lm(pairwise_Ist_vec~pairwise_geodist_vec)
summary(reg.geo)
random.sample <- sample(x= 1: length(pairwise_Ist_vec), size = 1000, replace = F)
tiff(filename = paste0(internal.wd,"/Community_Structure/Plots/Ist_Dgeo.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
plot(pairwise_geodist_vec[random.sample], pairwise_Ist_vec[random.sample],
     main = "Space structures \n mimicry community composition", ylab = "Pairwise Ist", xlab = "Geographic distance",
     cex.axis = 1.7, cex.main = 1, axis.args=list(cex.axis=1.7), cex.lab = 1.7)
abline(reg.geo, lwd = 2, col = "red")
# text(x = 2000, y = 0.65, labels = paste0("Rho = ", round(Mantel.Spearman$statistic,3), "\n p < 0.001 "), col = "red")
par(mar = original_int_margins)
dev.off()

# Ist ~ Dclim + covar(Dgeo) (r?gression sur les r?sidus)
Ist_resid <- residuals(reg.geo, type = "response")
reg.clim.geo <- lm(Ist_resid~pairwise_climdist_vec)
summary(reg.clim.geo)
random.sample <- sample(x= 1: length(pairwise_Ist_vec), size = 1000, replace = F)
tiff(filename = paste0(internal.wd,"/Community_Structure/Plots/Ist_Dclim_geo.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
plot(pairwise_climdist_vec[random.sample], Ist_resid[random.sample], 
     main = "Climate structures \n mimicry community composition \n even when space is taking into account", ylab = "Residual Pairwise Ist", xlab = "Standardized climatic distance",
     cex.axis = 1.7, cex.main = 1, axis.args=list(cex.axis=1.7), cex.lab = 1.7)
abline(reg.clim.geo, lwd = 2, col = "red")
# text(x = 8, y = 0.3, labels = paste0("Rho = ", round(Mantel.Spearman_Ist_geo$statistic,3), "\n p < 0.001 "), col = "red")
par(mar = original_int_margins)
dev.off()
