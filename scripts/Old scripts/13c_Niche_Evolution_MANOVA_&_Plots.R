##### Niche Evolution : MANOVA tests et plots des esp?ces/units dans l'espace climatique #####

### Preparation ###

# Effacer l'environnement
rm(list = ls())

library(raster)

# Set Wd
internal.wd <- "F:/Documents/Etudes/Fac/Cours/TROPIMUNDO/Stage_S4/Projet_Papillons/R_codes"

# Load phylogeny, sp list and unit list
load(file = paste0(internal.wd, "/Phylogenies/Final_phylogeny.RData"))
load(file = paste0(internal.wd,"/list.sp.RData"))
load(file = paste0(internal.wd,"/unit.env.table.RData"))

# Modeling resolution
res <-  "5"

# Extract only the 339 species included in the phylogeny from list.sp, list.models and the sp.stack
list.sp <- list.sp[list.sp$Sp_full %in% phylo.Ithomiini$tip.label,] # 339 species
list.models <- list.models[list.models$Sp_full %in% phylo.Ithomiini$tip.label,] # 719 units


### Test sur les units donc besoin d'ajouter les units en tant que polytomies de longueurs nulles sur l'arbre !
write.nexus(phylo.Ithomiini, file = paste0(internal.wd, "/Phylogenies/Species_phylogeny.nex"), translate = F)  

plot(phylo.Ithomiini)

library(phytools)

# phylo.Ithomiini.units <- phylo.Ithomiini

# Add units tips
?bind.tip
for (i in 1:nrow(list.models)) {
  unit <- list.models[i,]
  sp.nod <- which(phylo.Ithomiini.units$tip.label==unit$Sp_full)
  sp.edge <- which(phylo.Ithomiini.units$edge[,2]==sp.nod)
  phylo.Ithomiini.units <- bind.tip(tree = phylo.Ithomiini.units, tip.label = as.character(unit$Tag.model),
                                    edge.length = 0,
                                    where=sp.nod,
                                    position=0)
  print(i)
}
# Remove species tips
?drop.tip
for (i in 1:nrow(list.sp)) {
  sp <- list.sp$Sp_full[i]
  phylo.Ithomiini.units <- drop.tip(phy = phylo.Ithomiini.units, tip = sp)
  print(i)
}

# phylo.Ithomiini.units$tip.label

save(phylo.Ithomiini.units, file = paste0(internal.wd, "/Phylogenies/Final_units_phylogeny.RData"))
write.nexus(phylo.Ithomiini.units, file = paste0(internal.wd, "/Phylogenies/Units_phylogeny.nex"), translate = F)  

# Reorder list.models following Phylogeny
order.index <- NA
for (i in 1:length(phylo.Ithomiini.units$tip.label)) {
  unit <- as.character(phylo.Ithomiini.units$tip.label[i]) 
  order.index[i] <- which(list.models$Tag.model==unit)
}

list.models <- list.models[order.index,]
row.names(list.models) <- as.character(list.models$Tag.model)

save(list.models, file = paste0(internal.wd,"/unit.env.table.RData"))


# Load the final units phylogeny
load(file = paste0(internal.wd, "/Phylogenies/Final_units_phylogeny.RData"))

# Load the final env. variables df
load(file = paste0(internal.wd,"/unit.env.table.RData"))


### Reduce the number of circles !!! (N min = 10 per circle)
table(list.models$Mimicry.model)
N.sp.per.ring <- table(list.models$Mimicry.model)[order(decreasing = T,table(list.models$Mimicry.model))]
barplot(N.sp.per.ring, axes = T)
?barplot
sum(table(list.models$Mimicry.model)>=10) # 23 circles instead of 44
abline(v = 27.8, lty = 2)
big.rings <- names(table(list.models$Mimicry.model))[table(list.models$Mimicry.model)>=10] 

# Reduce dataset
reduced.list.models <- list.models[list.models$Mimicry.model %in% big.rings,] # 619 units instead of 719
save(reduced.list.models, file = paste0(internal.wd, "/Niche_evolution/reduced.list.models.RData"))

# Reduce phylogeny
remove.units <- as.character(list.models$Tag.model[!(list.models$Mimicry.model %in% big.rings)])
reduced.phylo.Ithomiini.units <- drop.tip(phy = phylo.Ithomiini.units, tip = remove.units)



##### permMANOVA without phylogenetic correction ####

### On 5 bioclim variables ####

unit.env.table <- reduced.list.models[,c("bio1","bio3","bio4","bio12","bio15")]

library("vegan")

?adonis
?adonis2

permMANOVA.5 <- adonis(formula = as.matrix(scale(unit.env.table)) ~ reduced.list.models$Mimicry.model, permutations = 999, method = "euclidian", strata = NULL, by = "margin") 
# Strata = pour contraindre les permutations au sein des groupes # Can be used to constrain under phylogeny ?
# By = "terms" pour Anova type I ; By = "margin" pour Anova type II. Ici, une seule variable explicative donc ?a ne change rien !
print(permMANOVA.5) # Sortie de la PERMANOVA
str(permMANOVA.5) # Pour explorer en d?tails les coeffs, les r?sultats des F des permutations, la model.matrix, ...

# Post-hoc test if significant to detect between each modalities of the factor the difference is significant
library("RVAideMemoire")

?pairwise.perm.manova
# Pour les diff?rents tests, voir ?anova.mlm
# pour les diff?rentes corrections de p-value due aux multiples tests, voir ?p.adjust
Pairwise.permMANOVA.5 <- pairwise.perm.manova(resp = as.matrix(scale(unit.env.table)), fact = reduced.list.models$Mimicry.model, nperm = 999, progress = T, p.method = "none", test = "Wilks")
Pairwise.permMANOVA.5
save(Pairwise.permMANOVA.5, file = paste0(internal.wd, "/Niche_evolution/Pairwise.permMANOVA.5.RData"))

### On 2 Revell's pPC axis ####

unit.env.table <- reduced.list.models[,c("bio1","bio3","bio4","bio12","bio15")]

# Revell's pPCA

source(paste0(internal.wd,"/Niche_evolution/Revell_phyl_pca.R"))
library(ape)

C <- vcv(reduced.phylo.Ithomiini.units)
X <- as.matrix(unit.env.table)

PCA.Revell <- Revell_phyl_pca(C, X, mode = "corr")
# Not working because the matrix C is not invertible...

# Classic PCA instead...
library(ade4) 



ACP <-  dudi.pca(df = unit.env.table, center = T, scale = T, scannf = F, nf = 2) # G?n?re l'ACP et le graph des ?boulis directement depuis les donn?es brutes
round(ACP$eig/sum(ACP$eig)*100,3) # % Variance expliqu?e
cumsum(round(ACP$eig/sum(ACP$eig)*100,3)) # % Variance expliqu?e cumulative
# 92% de la variance expliqu?e avec les 2 premiers axes, on ne garde que ces deux axes !

s.corcircle(ACP$co,xax=1,yax=2) # Pour tracer le cercle des correlations entre variables d'origine et CP. On peut choisir quelles CP sont repr?sent?es sur les axes x et y

reduced.PC.env <- ACP$li # Extraction des nouvelles variables

save(reduced.PC.env, file = paste0(internal.wd, "/Niche_evolution/reduced.PC.env.RData"))

load(file = paste0(internal.wd, "/Niche_evolution/reduced.PC.env.RData"))

library("vegan")

?adonis
?adonis2

permMANOVA.pPC <- adonis(formula = as.matrix(reduced.PC.env) ~ reduced.list.models$Mimicry.model, permutations = 999, method = "euclidian", strata = NULL, by = "margin") 
# Strata = pour contraindre les permutations au sein des groupes # Can be used to constrain under phylogeny ?
# By = "terms" pour Anova type I ; By = "margin" pour Anova type II. Ici, une seule variable explicative donc ?a ne change rien !
print(permMANOVA.pPC) # Sortie de la PERMANOVA
str(permMANOVA.pPC) # Pour explorer en d?tails les coeffs, les r?sultats des F des permutations, la model.matrix, ...
save(permMANOVA.pPC, file = paste0(internal.wd, "/Niche_evolution/permMANOVA.pPC.RData"))
load(file = paste0(internal.wd, "/Niche_evolution/permMANOVA.pPC.RData"))

# Post-hoc test if significant to detect between each modalities of the factor the difference is significant ####
library("RVAideMemoire")

?pairwise.perm.manova
# Pour les diff?rents tests, voir ?anova.mlm
# pour les diff?rentes corrections de p-value due aux multiples tests, voir ?p.adjust
Pairwise.permMANOVA.pPC <- pairwise.perm.manova(resp = as.matrix(reduced.PC.env), fact = reduced.list.models$Mimicry.model, nperm = 999, progress = T, p.method = "none", test = "Wilks")
Pairwise.permMANOVA.pPC
save(Pairwise.permMANOVA.pPC, file = paste0(internal.wd, "/Niche_evolution/Pairwise.permMANOVA.pPC.RData"))
load(file = paste0(internal.wd, "/Niche_evolution/Pairwise.permMANOVA.pPC.RData"))

Pairwise.permMANOVA.pPC
Pairwise.permMANOVA.pPC$p.value

sum(as.dist(Pairwise.permMANOVA.pPC$p.value)<0.002) # 168 out of 253
sum(as.dist(Pairwise.permMANOVA.pPC$p.value)<0.01) # 184 out of 253
sum(as.dist(Pairwise.permMANOVA.pPC$p.value)<0.05) # 205 out of 253

p.adjust(p = as.dist(Pairwise.permMANOVA.pPC$p.value), method = "holm")

Pseudo_Cor_mat_perMANOVA <- as.matrix(as.dist(Pairwise.permMANOVA.pPC$p.value))
diag(Pseudo_Cor_mat_perMANOVA) <- 1

hist(Pairwise.permMANOVA.pPC$p.value)
hist(log1p(Pairwise.permMANOVA.pPC$p.value))
hist(sqrt(sqrt(Pairwise.permMANOVA.pPC$p.value)))

# Plot de la heatmap
tiff(filename = paste0(internal.wd, "/Niche_evolution/Plots/Heatmap_PerMANOVA.tiff"))
heatmap.res <- heatmap(x = Pseudo_Cor_mat_perMANOVA^(1/5.8), hclustfun = hclust, symm = T, col = pal_bl_red)
dev.off()

?heatmap

ordered.Pseudo_Cor_mat_perMANOVA <- Pseudo_Cor_mat_perMANOVA[heatmap.res$rowInd, heatmap.res$colInd]

# Plot de la scale
tiff(filename = paste0(internal.wd, "/Niche_evolution/Plots/Heatmap_PerMANOVA_scale.tiff"))
original_ext_margins <- par()$oma
par(mar = c(0,0,0,5))
plot(raster(x = ordered.Pseudo_Cor_mat_perMANOVA), col = pal_bl_red,
     axis.args=list(at = c(0,0.2,0.4,0.6,0.8,1), cex.axis = 1.7,
                    labels=round((c(0,0.2,0.4,0.6,0.8,1)^5.8), digits = 3)))
dev.off()



##### Phylogenetic MANOVA ####

### Can only be done on PC axis because we don't have a good evolutionary model for the 5 original dimensions

# Simulate variables following the best evolutionary model (but calibrated at species level...)

load(file = paste0(internal.wd, "/Niche_evolution/best.model.for.pPC_evol"))

library(phylocurve)

Sim.Revell.res <- sim.traits(tree = phylo.Ithomiini.units, v = best.model$brownianVariance, anc = best.model$root.state, model = "lambda", parameters = list(lambda = best.model$Lambda[1,1]), nsim = 999, return.type = "matrix")

str(Sim.Revell.res)
plot(Sim.Revell.res$tree) # Arbre original
plot(Sim.Revell.res$sim_tree) # Arbre apr?s tranformation via lambda = 0.47

Sim.Revell.PC <- Sim.Revell.res$trait_data # Liste des simulations

# Compute MANOVA on observed data and extract the Wilks' summary stat

load(file = paste0(internal.wd, "/Niche_evolution/reduced.PC.env.RData"))
MANOVA_obs.pPC <- manova(as.matrix(reduced.PC.env) ~ reduced.list.models$Mimicry.model, test="Wilks")
MANOVA_obs.pPC 
summary(MANOVA_obs.pPC)$stats[,2][1]
Wilks_lambda_obs.pPC <- summary(MANOVA_obs.pPC)$stats[,2][1]
save(Wilks_lambda_obs.pPC, file = paste0(internal.wd, "/Niche_evolution/Wilks_lambda_obs.pPC.RData"))
summary.aov(MANOVA_obs.pPC)

# Compute MANOVA on each simulation and extract the Wilks' summary stat
Wilks_lambda_null.pPC <- NA
for (i in 1:length(Sim.Revell.PC)){
  MANOVA_sim.pPC <- manova(as.matrix(Sim.Revell.PC[[i]]) ~ reduced.list.models$Mimicry.model, test="Wilks")
  Wilks_lambda_sim.pPC <- summary(MANOVA_sim.pPC)$stats[,2][1]
  Wilks_lambda_null.pPC[i] <- Wilks_lambda_sim.pPC 
  print(i)
}
summary(Wilks_lambda_null.pPC)
save(Wilks_lambda_null.pPC, file = paste0(internal.wd, "/Niche_evolution/Wilks_lambda_null.pPC.RData"))
load(file = paste0(internal.wd, "/Niche_evolution/Wilks_lambda_null.pPC.RData"))

# Plot the distri


jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/MANOVA_phylo_PC.jpeg"), quality =100)
hist(c(Wilks_lambda_null.pPC, Wilks_lambda_obs.pPC), 30, freq = TRUE, col = "gray", main = "Phylogenetic MANOVA \n Distribution of Wilks' lambda \n under the null hypothesis", xlab = "Wilks' lambda")
arrows(Wilks_lambda_obs.pPC+0.004, 50, Wilks_lambda_obs.pPC+0.004, 5, length = 0.1, lwd = 2)
abline(v = mean(c(Wilks_lambda_null.pPC, Wilks_lambda_obs.pPC)), lty = 2)
abline(v = quantile(c(Wilks_lambda_null.pPC, Wilks_lambda_obs.pPC), 0.95), lty = 2, col = "red")
legend(legend = c(paste0("Mean = ", round(mean(c(Wilks_lambda_null.pPC, Wilks_lambda_obs.pPC)),3)), 
                  paste0("CI 95% = ", round(quantile(c(Wilks_lambda_null.pPC, Wilks_lambda_obs.pPC), 0.95),3))), 
       x = "topright", cex = 1, bty ="n")
legend(legend = c(paste0("Lambda obs = ", round(Wilks_lambda_obs.pPC, 3)),
                  "          p < 0.001"),
       x = "right", cex = 1, bty ="n")
dev.off()

load(file = paste0(internal.wd, "/Niche_evolution/Wilks_lambda_obs.pPC.RData"))
load(file = paste0(internal.wd, "/Niche_evolution/Wilks_lambda_null.pPC.RData"))
tiff(filename = paste0(internal.wd,"/Niche_evolution/Plots/MANOVA_phylo_PC.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
hist(c(Wilks_lambda_null.pPC, Wilks_lambda_obs.pPC), 30, freq = TRUE, col = "gray", 
     main = "Phylogenetic MANOVA \n Distribution of Wilks' lambda \n under the null hypothesis", xlab = "Wilks' lambda ( )",
     cex.axis = 1.7, cex.main = 1, cex.lab = 1.7)
arrows(Wilks_lambda_obs.pPC + 0.0001, 80, Wilks_lambda_obs.pPC+ 0.0001, 10, length = 0.1, lwd = 2)
abline(v = mean(c(Wilks_lambda_null.pPC, Wilks_lambda_obs.pPC)), lty = 2)
abline(v = quantile(c(Wilks_lambda_null.pPC, Wilks_lambda_obs.pPC), 0.95), lty = 2, col = "red")
legend(legend = c(paste0("Mean = ", round(mean(c(Wilks_lambda_null.pPC, Wilks_lambda_obs.pPC)),3)), 
                  paste0("CI 95% = ", round(quantile(c(Wilks_lambda_null.pPC, Wilks_lambda_obs.pPC), 0.95),3))), 
       x = "topright", lty = 2, lwd = 2, col = c("black", "red"), cex = 1.2, bty ="n")
legend(legend = c(paste0("l obs = ", round(Wilks_lambda_obs.pPC, 3)),
                  paste0("      p = 0.001")),
       x = "right", cex = 1, bty ="n", xjust = 1)
par(mar = original_int_margins)
dev.off()







##### Plot mimicry ring in reduced env. space for species

# Compute standardized euclidian distances in the climatic space for species
load(file = paste0(internal.wd,"/sp.env.table.RData"))
sp.env.table <- list.sp[,c("bio1","bio3","bio4","bio12","bio15")]

std_env <- scale(x = sp.env.table, center = T, scale = T)
clim_dist <- as.matrix(dist(x = std_env, method = "euclidian"))

save(clim_dist, file = paste0(internal.wd, "/Niche_evolution/sp.clim.dist.RData"))

### NMDS Version
# 2D NMDS
library(vegan)
index <- (list.sp$Mimicry=="Mono")&(list.sp$N.mimicry==1) # Plot only monomorph species
load(file = paste0(internal.wd, "/Community_Structure/Pairwise_values.RData"))
NMDS <- vegan::metaMDS(comm = clim_dist, k = 2, plot = F) # Compute the NMDS
plot(NMDS, type = "p", choices = c(1,2), col = rainbow(180)[as.factor(list.sp$Mimicry_full)]) # Plot type "t" = text/labels. # Choices = axes to show
points(NMDS$points[index,], pch = 16 , col = rainbow(44)[as.factor(list.sp$Mimicry_full[index])])
legend(legend = paste0("Stress = ", round(NMDS$stress,3)), # If stress > 0.15, the NMDS graph is not relevant.
       x = "topright", cex = 1, bty ="n") 

# Enhanced 3D NMDS plot
NMDS_3D <- vegan::metaMDS(comm = clim_dist, k = 3, plot = F, maxit = 1000) # Compute the NMDS

rgl::plot3d(x = NMDS$points[index,], type = "n") # To generate the windows using the axis
rgl::points3d(x = NMDS$points[index,], col= rainbow(44)[as.factor(list.sp$Mimicry_full[index])]) # To add points in the desired color in the 3D plot

### Classic PCA version

library(ade4) 

ACP <-  dudi.pca(df = sp.env.table, center = T, scale = T, scannf = F, nf = 2) # G?n?re l'ACP et le graph des ?boulis directement depuis les donn?es brutes
round(ACP$eig/sum(ACP$eig)*100,3) # % Variance expliqu?e
cumsum(round(ACP$eig/sum(ACP$eig)*100,3)) # % Variance expliqu?e cumulative

s.corcircle(ACP$co,xax=1,yax=2) # Pour tracer le cercle des correlations entre variables d'origine et CP. On peut choisir quelles CP sont repr?sent?es sur les axes x et y

index <- (list.sp$Mimicry=="Mono")&(list.sp$N.mimicry==1) # Plot only monomorph species

plot(ACP$li, col = "grey80", main = "PCA of species in climatic space", xlab = paste0("PC1 (", round(ACP$eig/sum(ACP$eig)*100,2)[1] ," %)"), ylab = paste0("PC2 (", round(ACP$eig/sum(ACP$eig)*100,2)[2] ," %)"))
points(ACP$li[index,], pch = 16, col = rainbow(44)[as.factor(list.sp$Mimicry_full[index])])
arrows(x0 = rep(0, nrow(ACP$co)), y0 = rep(0, nrow(ACP$co)), x1 = ACP$co[,1]*3, y1 =  ACP$co[,2]*3, length = 0.2, code = 2, lwd = 1)
text(ACP$co[-2,]*3.3, labels = row.names(ACP$co)[-2])
text(ACP$co[2,]*3.3, labels = row.names(ACP$co)[2], adj = c(1, 2))

table(list.sp$Mimicry_full[index])
# AGNOSIA 25 # LERIDA 21 # EURIMEDIA 17 # CONFUSA 10 # PANTHYALE 9 # HERMIAS 8

# Plot with only the 6 main circles in colors

top6.index <- (list.sp$Mimicry=="Mono")&(list.sp$N.mimicry==1)&(list.sp$Mimicry_full %in% c("AGNOSIA.AGNOSIA", "LERIDA.LERIDA", "EURIMEDIA.EURIMEDIA", "CONFUSA.CONFUSA", "PANTHYALE.PANTHYALE", "HERMIAS.HERMIAS")) # Plot only monomorph species of the top 6 mimicry rings

jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/PCA_Top6.jpeg"), quality =100)
plot(ACP$li, col = "grey80", main = "PCA of species in climatic space", xlab = paste0("PC1 (", round(ACP$eig/sum(ACP$eig)*100,2)[1] ," %)"), ylab = paste0("PC2 (", round(ACP$eig/sum(ACP$eig)*100,2)[2] ," %)"))
points(ACP$li[top6.index,], pch = 16, col = c("black","red","deepskyblue","limegreen","orange","blue")[as.factor(list.sp$Mimicry_full[top6.index])])
arrows(x0 = rep(0, nrow(ACP$co)), y0 = rep(0, nrow(ACP$co)), x1 = ACP$co[,1]*3, y1 =  ACP$co[,2]*3, length = 0.2, code = 2, lwd = 1)
text(ACP$co[-2,]*3.3, labels = row.names(ACP$co)[-2])
text(ACP$co[2,]*3.3, labels = row.names(ACP$co)[2], adj = c(1, 2))
legend(x = "topleft", pch = 16, cex = 0.75, legend = unique(unlist(strsplit(x = levels(as.factor(list.sp$Mimicry_full[top6.index])), split = ".", fixed = T))), col = c("black","red","deepskyblue","limegreen","orange","blue"))
# s.class(ACP$li[top6.index,], label = NULL, fac = as.factor(list.sp$Mimicry_full[top6.index]), col=c("black","red","deepskyblue","limegreen","orange","blue"), add.p=TRUE, xax=1, yax=2) # Permet d'associer les individus en classes sur le plot de s.label. # fac = vecteur de facteurs qui partitionne les lignes de ACP$li en classes d'individus # col = vecteur de couleur repr?senter pour chaque classe # add.plot=T pour ajouter le graph' des classes au graph' d?j? trac?
dev.off()

### Revell's pPCA version 

# Compute PCA to extract new axis

source(paste0(internal.wd,"/Niche_evolution/Revell_phyl_pca.R"))
library(ape)

C <- vcv(phylo.Ithomiini)
X <- as.matrix(sp.env.table)

PCA.Revell <- Revell_phyl_pca(C, X, mode = "corr")
diag(PCA.Revell$Eval) # eigenvalues
PCA.Revell$Evec # eigenvectors
PCA.Revell$S # scores
PCA.Revell$L # PC loadings

PCvar <- round(diag(PCA.Revell$Eval)/sum(PCA.Revell$Eval),3) ; PCvar # % Variance expliqu?e
cumVar <- round(cumsum(diag(PCA.Revell$Eval))/sum(PCA.Revell$Eval),3) ; cumVar # % Variance expliqu?e cumulative
# 91% de la variance expliqu?e avec les 2 premiers axes, on ne garde que ces deux axes !

library(ade4)
s.corcircle(PCA.Revell$L, xax=1, yax=2, label = names(sp.env.table)) # Pour tracer le cercle des correlations entre variables d'origine et CP. On peut choisir quelles CP sont repr?sent?es sur les axes x et y

pPC.env <- PCA.Revell$S[,1:2]

index <- (list.sp$Mimicry=="Mono")&(list.sp$N.mimicry==1) # Plot only monomorph species

plot(pPC.env, col = "grey80", main = "pPCA of species in climatic space", xlab = paste0("PC1 (", PCvar[1] ," %)"), ylab = paste0("PC2 (", PCvar[2] ," %)"))
points(pPC.env[index,], pch = 16, col = rainbow(44)[as.factor(list.sp$Mimicry_full[index])])
arrows(x0 = rep(0, nrow(ACP$co)), y0 = rep(0, nrow(ACP$co)), x1 = ACP$co[,1]*4.5, y1 =  ACP$co[,2]*4.5, length = 0.2, code = 2, lwd = 1)
text(ACP$co[-2,]*4.9, labels = row.names(ACP$co)[-2])
text(ACP$co[2,]*5, labels = row.names(ACP$co)[2], adj = c(0.5, -1))

table(list.sp$Mimicry_full[index])
# AGNOSIA 25 # LERIDA 21 # EURIMEDIA 17 # CONFUSA 10 # PANTHYALE 9 # HERMIAS 8

# Plot with only the 6 main circles in colors

top6.index <- (list.sp$Mimicry=="Mono")&(list.sp$N.mimicry==1)&(list.sp$Mimicry_full %in% c("AGNOSIA.AGNOSIA", "LERIDA.LERIDA", "EURIMEDIA.EURIMEDIA", "CONFUSA.CONFUSA", "PANTHYALE.PANTHYALE", "HERMIAS.HERMIAS")) # Plot only monomorph species of the top 6 mimicry rings

# Build minimum convex polygon for the 6 main circles

# Que les points, avec variables d'origine
jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/pPCA_Top6.jpeg"), quality =100)
plot(pPC.env, col = "grey80", main = "pPCA of species in climatic space", xlab = paste0("PC1 (", round(PCvar[1],2)*100 ," %)"), ylab = paste0("PC2 (", round(PCvar[2],2)*100 ," %)"))
points(pPC.env[top6.index,], pch = 16, col = c("black","red","deepskyblue","limegreen","orange","blue")[as.factor(list.sp$Mimicry_full[top6.index])])
legend(x = "bottomleft", pch = 16, cex = 0.75, legend = unique(unlist(strsplit(x = levels(as.factor(list.sp$Mimicry_full[top6.index])), split = ".", fixed = T))), col = c("black","red","deepskyblue","limegreen","orange","blue"))
arrows(x0 = rep(0, nrow(PCA.Revell$L)), y0 = rep(0, nrow(PCA.Revell$L)), x1 = PCA.Revell$L[,1]*4.5, y1 =  PCA.Revell$L[,2]*4.5, length = 0.2, code = 2, lwd = 1)
text(PCA.Revell$L[-2,]*4.9, labels = c("Tmean", "Tiso", "Tvar", "Hmean", "Hvar")[-2])
text(x = PCA.Revell$L[2,1]*5, y = PCA.Revell$L[2,2]*5, labels = c("Tmean", "Tiso", "Tvar", "Hmean", "Hvar")[2], adj = c(0.5, -1))
dev.off()

?text

# Avec ellipses, sans les variables d'origine
jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/pPCA_Top6_ellipse.jpeg"), quality =100)
plot(pPC.env, col = "grey80", main = "pPCA of species in climatic space", xlab = paste0("PC1 (", round(ACP$eig/sum(ACP$eig)*100,2)[1] ," %)"), ylab = paste0("PC2 (", round(ACP$eig/sum(ACP$eig)*100,2)[2] ," %)"))
points(pPC.env[top6.index,], pch = 16, col = c("black","red","deepskyblue","limegreen","orange","blue")[as.factor(list.sp$Mimicry_full[top6.index])])
legend(x = "bottomleft", pch = 16, cex = 0.75, legend = unique(unlist(strsplit(x = levels(as.factor(list.sp$Mimicry_full[top6.index])), split = ".", fixed = T))), col = c("black","red","deepskyblue","limegreen","orange","blue"))
s.class(pPC.env[top6.index,], label = NULL, fac = as.factor(list.sp$Mimicry_full[top6.index]), col=c("black","red","deepskyblue","limegreen","orange","blue"), add.p=TRUE, xax=1, yax=2) # Permet d'associer les individus en classes sur le plot de s.label. # fac = vecteur de facteurs qui partitionne les lignes de pPC.env en classes d'individus # col = vecteur de couleur repr?senter pour chaque classe # add.plot=T pour ajouter le graph' des classes au graph' d?j? trac?
dev.off()


##### Plot mimicry ring in reduced env. space for units

# Compute standardized euclidian distances in the climatic space for units
load(file = paste0(internal.wd,"/unit.env.table.RData"))
unit.env.table <- list.models[,c("bio1","bio3","bio4","bio12","bio15")]
names(unit.env.table) <- c("Tmean", "Tiso", "Tvar", "Hmean", "Hvar")

unit_std_env <- scale(x = unit.env.table, center = T, scale = T)
unit_clim_dist <- as.matrix(dist(x = unit_std_env, method = "euclidian"))

save(unit_clim_dist, file = paste0(internal.wd, "/Niche_evolution/unit.clim.dist.RData"))

### Revell's pPCA version (not working because the matrix cannot be inversed :-())

# Compute pPCA to extract new axis

source(paste0(internal.wd,"/Niche_evolution/Revell_phyl_pca.R"))
library(ape)

C <- vcv(phylo.Ithomiini.units)
X <- as.matrix(unit.env.table)

PCA.Revell.unit <- Revell_phyl_pca(C, X, mode = "corr")

# Classic PCA by default

library(ade4) 

ACP.unit <-  dudi.pca(df = unit.env.table, center = T, scale = T, scannf = F, nf = 2) # G?n?re l'ACP et le graph des ?boulis directement depuis les donn?es brutes
PCvar <- round(ACP.unit$eig/sum(ACP.unit$eig)*100,1) ; PCvar # % Variance expliqu?e
cumsum(round(ACP.unit$eig/sum(ACP.unit$eig)*100,1)) # 92 % Variance expliqu?e cumulative

s.corcircle(ACP.unit$co,xax=1,yax=2) # Pour tracer le cercle des correlations entre variables d'origine et CP. On peut choisir quelles CP sont repr?sent?es sur les axes x et y

# On utilise abusivement le terme pPCA... (mais le r?sultat est quand m?me hyper similaire dans tous les cas)

plot(ACP.unit$li, col = "grey80", main = "pPCA of OMUs in climatic space", xlab = paste0("pPC1 (", PCvar[1] ," %)"), ylab = paste0("pPC2 (", PCvar[2] ," %)"))
points(ACP.unit$li, pch = 16, col = rainbow(44)[as.factor(list.models$Mimicry.model)])
arrows(x0 = rep(0, nrow(ACP$co)), y0 = rep(0, nrow(ACP$co)), x1 = ACP$co[,1]*3, y1 =  ACP$co[,2]*3, length = 0.2, code = 2, lwd = 1)
text(ACP$co[-2,]*3.5, labels = row.names(ACP.unit$co)[-2])
text(ACP$co[2,]*3.5, labels = row.names(ACP.unit$co)[2], adj = c(1.1, 1.2))

Mimicry_counts_table <- table(list.models$Mimicry.model)[order(table(list.models$Mimicry.model), decreasing = T)] ; Mimicry_counts_table
tiff(filename = paste0(internal.wd,"/Niche_evolution/Plots/Barplot_OMU_mimic_count.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
barplot(height = Mimicry_counts_table, ylab = "Number of OMUs", xlab = "Mimicry rings", xaxt='n',
        cex.axis = 1.7, cex.main = 1, cex.lab = 1.7)
abline(v=27.8, lty = 2, lwd = 2, col = "red")
abline(h=9, lty = 2, lwd = 2, col = "red")
par(mar = original_int_margins)
dev.off()

# AGNOSIA 74 # LERIDA 63 # MAMERCUS 56 # HERMIAS 47 # BANJANA-M 43 # PANTHYALE 38 # DULICIDA 35 # EURIMEDIA 33 # HEWITSONI 27 

# Plot with only the 9 main circles in colors

top5.ring.list <- c("AGNOSIA", "LERIDA", "MAMERCUS", "HERMIAS", "BANJANAM")
top9.ring.list <- c("AGNOSIA", "LERIDA", "MAMERCUS", "HERMIAS", "BANJANAM", "PANTHYALE", "DULCIDA", "EURIMEDIA", "HEWITSONI")

top5.index <- (list.models$Mimicry.model %in% top5.ring.list)
top9.index <- (list.models$Mimicry.model %in% top9.ring.list) 

# Build minimum convex polygon for the 9 main circles

library(adehabitatHR)
library(sp)

?mcp

for (i in length(top9.ring.list)) {
  ring <- top9.ring.list[i]
  ring.index <- list.models$Mimicry.model == ring
  SPobj.ring <- SpatialPoints(coords = ACP.unit$li[ring.index,]) # Create a Spobj from the coordinates
  mcp.ring <- mcp(xy = SPobj.ring, percent = 100) # Generate the minimum convex polygon
  save(mcp.ring, file = paste0(internal.wd,"/Niche_evolution/MANOVAs/mcp.ring_", ring, ".RData"))
}


# Que le top 5
mimic.list <- levels(factor(list.models$Mimicry.model[top5.index]))
tiff(filename = paste0(internal.wd,"/Niche_evolution/Plots/pPCA_Top5.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))
plot(ACP.unit$li, col = "grey80", main = "pPCA of OMUs in climatic space", 
     xlab = paste0("pPC1 (", round(PCvar[1],1) ," %)"), ylab = paste0("pPC2 (", round(PCvar[2],1) ," %)"),
     cex.axis = 1.7, cex.lab = 1.7, cex.main = 1)
points(ACP.unit$li[top5.index,], pch = 16, col = c("black","red","deepskyblue","limegreen","orange","blue")[factor(list.models$Mimicry.model[top5.index])])
legend(x = "topleft", pch = 16, cex = 1.1, legend = mimic.list, col = c("black","red","deepskyblue","limegreen","orange","blue"))
arrows(x0 = rep(0, nrow(ACP$co)), y0 = rep(0, nrow(ACP$co)), x1 = ACP$co[,1]*3, y1 =  ACP$co[,2]*3, length = 0.2, code = 2, lwd = 1)
# for (i in 1:length(mimic.list)) {
#   ring <- mimic.list[i]
#   load(file = paste0(internal.wd,"/Niche_evolution/MANOVAs/mcp.ring_", ring, ".RData"))
#   plot(mcp.ring, add = T, border = c("black","red","deepskyblue","limegreen","orange","blue")[i], lwd = 2)
#  }
text(ACP$co[-2,]*3.5, labels = row.names(ACP.unit$co)[-2], font = 2)
text(ACP$co[2,]*3.5, labels = row.names(ACP.unit$co)[2], font = 2, adj = c(1.2, 1.3))
par(mar = original_int_margins)
dev.off()

?legend

?text
