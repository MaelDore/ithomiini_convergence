##### Approche brute de l'estimation du meilleur modèle Kappa - Lambda #####

### Preparation ###

# Effacer l'environnement
rm(list = ls())

library(raster)

# Set Wd
internal.wd <- "F:/Documents/Etudes/Fac/Cours/TROPIMUNDO/Stage_S4/Projet_Papillons/R_codes"
initial.wd <-  "F:/Documents/Etudes/Fac/Cours/TROPIMUNDO/Modeling/"

# Load phylogeny, sp list and unit list
load(file = paste0(internal.wd, "/Phylogenies/Final_phylogeny.RData"))
load(file = paste0(internal.wd,"/list.sp.RData"))
load(file = paste0(internal.wd,"/list.models.RData"))

# Modeling resolution
res <-  "5"

# Extract only the 339 species included in the phylogeny from list.sp, list.models and the sp.stack
list.sp <- list.sp[list.sp$Sp_full %in% phylo.Ithomiini$tip.label,] # 339 species
list.models <- list.models[list.models$Sp_full %in% phylo.Ithomiini$tip.label,] # 719 units

# Load the env. variables df
load(file = paste0(internal.wd,"/sp.env.table.RData"))


##### Approche bivariée après pPCA (car Revell dit que l'on ne peut pas faire de PCA sur des données phylo. sans prendre en compte la phylogénie) #####

load(file = paste0(internal.wd,"/sp.env.table.RData"))

sp.env.table <- list.sp[,c("bio1","bio3","bio4","bio12","bio15")]
row.names(sp.env.table) <- list.sp$Sp_full

# Vérifier que les espèces sont bien ordonnées de la même manière
list.sp$Sp_full
phylo.Ithomiini$tip.label

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

round(diag(PCA.Revell$Eval)/sum(PCA.Revell$Eval),3) # % Variance expliquée
round(cumsum(diag(PCA.Revell$Eval))/sum(PCA.Revell$Eval),3) # % Variance expliquée cumulative
# 91% de la variance expliquée avec les 2 premiers axes, on ne garde que ces deux axes !

library(ade4)
s.corcircle(PCA.Revell$L, xax=1, yax=2, label = names(sp.env.table)) # Pour tracer le cercle des correlations entre variables d'origine et CP. On peut choisir quelles CP sont repr?sent?es sur les axes x et y

PC.Revell.env <- PCA.Revell$S[,1:2]

hist(PC.Revell.env[,1])
hist(PC.Revell.env[,2])

library(motmot.2.0)

?transformPhylo.ML

par(mfcol=(c(1,1)))

# Create summary table
Models.PC.Revell.df <- data.frame(Model = character(), Likelihood = numeric(), AICc = numeric(), param1 = numeric(), param2 = numeric())
Models.PC.Revell.df$Model <- as.character(Models.PC.Revell.df$Model)
# Modèle Brownien
Revell.bm.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="bm")
Revell.bm.ml
Models.PC.Revell.df[1,] <-  c("BM", round(Revell.bm.ml$logLikelihood, 3), round(Revell.bm.ml$AICc, 3), NA, NA)

# Modèle avec Pagel's lambda
Revell.lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="lambda", profilePlot=T)
Revell.lambda.ml
Models.PC.Revell.df[2,] <-  c("Lambda", round(Revell.lambda.ml$MaximumLikelihood, 3), round(Revell.lambda.ml$AICc, 3), round(Revell.lambda.ml$Lambda[1,1],3), NA)
p.value <- 1 - pchisq(Revell.lambda.ml$MaximumLikelihood - Revell.bm.ml$logLikelihood, 1) ; p.value
Revell.bm.ml$AICc - Revell.lambda.ml$AICc

# Pagel's lambda univarié pour voir
var <- 2
test <- as.matrix(PC.Revell.env[,var]) ; row.names(test) <- list.sp$Sp_full 
test.lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=test, model="lambda", profilePlot=T)
test.lambda.ml

# Test du signal phylogénétique en multivariée avec lambda de Pagel
tree0 <- transformPhylo(phy = phylo.Ithomiini, model = "lambda", y = as.matrix(PC.Revell.env), lambda = 0) # Transformation de l'arbre sous Pagel's lambda = 0
Revell.lambda0.ml <- transformPhylo.ML(phy=tree0, y=as.matrix(PC.Revell.env), model="bm") # Estimation du meilleur modèle en BM sur l'arbre transformé
p.value <- 1 - pchisq(Revell.lambda.ml$MaximumLikelihood - Revell.lambda0.ml$logLikelihood, 1) ; p.value # LRT via modèles emboités

# Modèle avec Pagel's Kappa
Revell.kappa.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="kappa", profilePlot=T)
p.value <- 1 - pchisq(Revell.kappa.ml$MaximumLikelihood - Revell.bm.ml$logLikelihood, 1) ; p.value
Revell.bm.ml$AICc - Revell.kappa.ml$AICc
Models.PC.Revell.df[3,] <-  c("Kappa", round(Revell.kappa.ml$MaximumLikelihood, 3), round(Revell.kappa.ml$AICc, 3), round(Revell.kappa.ml$Kappa[1,1],3), NA)


# Modèle avec Pagel's lambda et Pagel's Kappa
Revell.kappa.lambda.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="kappa", lambdaEst = T, profilePlot=T)
Revell.kappa.lambda.ml

### Approche brute par transformation de l'arbre et comparaison des vraisemblances

# Generate the vectors of parameter
lambda <- seq(0, 1, 0.01) ; kappa <- rev(seq(0, 1, 0.01))
# Generate the matrix of results
Lk.mat <- matrix(data = NA, nrow = length(kappa), ncol = length(lambda))

# Loop
for (i in 1:nrow(Lk.mat)) {
  for (j in 1:ncol(Lk.mat)) {
    tree_transfo <- transformPhylo(phy = phylo.Ithomiini, model = "lambda", lambda = lambda[j]) # Applique d'abord le lambda
    tree_transfo <- transformPhylo(phy = tree_transfo, model = "kappa", kappa = kappa[i]) # Puis le kappa
    # plot(tree_transfo, show.tip.label = F)
    Revell.transfo <- transformPhylo.ML(phy=tree_transfo, y=as.matrix(PC.Revell.env), model="bm") # Estimation du meilleur modèle en BM sur l'arbre transformé
    Lk.mat[i,j] <- Revell.transfo$logLikelihood 
  }
  print(i)
}
save(Lk.mat, file = paste0(internal.wd,"/Niche_evolution/Lk.mat.RData"))
max(Lk.mat)
which.max(Lk.mat)

pos <- arrayInd(which.max(Lk.mat), dim(Lk.mat))
lambda.max <- lambda[pos[2]]
kappa.max <- kappa[pos[1]]

tiff(filename = paste0(internal.wd,"/Niche_evolution/Plots/Heatmap_Kappa_Lambda.tiff"))
original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,4))
image(z = t(Lk.mat[nrow(Lk.mat):1,]), xlab = "Lambda", ylab = "Kappa", zlim = c(-1600,-1550), col = hcl.colors(25, "YlOrRd", rev = T),
      main = "Heatmap of Likelihood of Lambda-Kappa models",
      cex.axis = 1.7, cex.lab = 1.8, cex.main = 1)
points(y = kappa.max, x = lambda.max, cex = 2, col = "black", pch = 16)
points(y = kappa.max, x = lambda.max, cex = 1.5, col = "red", pch = 16)
par(mar = original_int_margins)
dev.off()

tiff(filename = paste0(internal.wd,"/Niche_evolution/Plots/Scale_Heatmap_Kappa_Lambda.tiff"))
original_ext_margins <- par()$oma
par(oma = c(0,0,0,5))
plot(raster(Lk.mat), zlim = c(-1600,-1550), col = hcl.colors(25, "YlOrRd", rev = T), axis.args=list(cex.axis=1.7))
par(oma = original_ext_margins)
dev.off()

tree_maxML <- transformPhylo(phy = phylo.Ithomiini, model = "lambda", lambda = lambda.max) # Applique d'abord le lambda
tree_maxML <- transformPhylo(phy = tree_maxML, model = "kappa", kappa = kappa.max) # Puis le kappa
Revell.kappa.lambda <- transformPhylo.ML(phy=tree_maxML, y=as.matrix(PC.Revell.env), model="bm") # Estimation du meilleur modèle en BM sur l'arbre transformé
Revell.kappa.lambda.ml$MaximumLikelihood ; Revell.kappa.lambda$logLikelihood

Revell.kappa.lambda$AICc <- Revell.kappa.lambda$AICc+4 # on ajoute la pénalisation pour les 2 paramètres à estimer en plus dans l'AICc

Models.PC.Revell.df[4,] <-  c("Kappa - Lambda", round(Revell.kappa.lambda$logLikelihood, 3), round(Revell.kappa.lambda$AICc, 3), kappa.max, lambda.max)
p.value <- 1 - pchisq(Revell.kappa.lambda$logLikelihood - Revell.bm.ml$logLikelihood, 2) ; p.value
Revell.bm.ml$AICc - Revell.kappa.lambda$AICc
p.value <- 1 - pchisq(Revell.kappa.lambda$logLikelihood - Revell.kappa.ml$MaximumLikelihood, 1) ; p.value
Revell.kappa.ml$AICc - Revell.kappa.lambda$AICc
p.value <- 1 - pchisq(Revell.kappa.lambda$logLikelihood - Revell.lambda.ml$MaximumLikelihood, 1) ; p.value
Revell.lambda.ml$AICc - Revell.kappa.lambda$AICc
# Différence non-significative, on garde le modèle le plus simple soit le modèle avec simplement un lambda.

save(Models.PC.Revell.df, file = paste0(internal.wd,"/Niche_evolution/Models.PC.Revell.df.RData"))

# Modèles bonus
delta.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="delta", profilePlot=T)
ou.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="OU", profilePlot=T)
acdc.ml <- transformPhylo.ML(phy=phylo.Ithomiini, y=as.matrix(PC.Revell.env), model="ACDC", profilePlot=T)


## Traits simulation under null hypothesis

library(phylocurve)

best.model <- Revell.lambda.ml

?sim.traits # Besoin d'une matrice de variance-covariance entre les traits (obtenue depuis motmot)

Sim.Revell.res <- sim.traits(tree = phylo.Ithomiini, v = best.model$brownianVariance, anc = best.model$root.state, model = "lambda", parameters = list(lambda = best.model$Lambda[1,1]), nsim = 999, return.type = "matrix")

str(Sim.Revell.res)
plot(Sim.Revell.res$tree) # Arbre original
plot(Sim.Revell.res$sim_tree) # Arbre après tranformation via lambda = 0.47

Sim.Revell.PC <- Sim.Revell.res$trait_data # Liste des simulations

## Stat computation and distribution plot

# Compute the climatic distance matrix
Revell_Pairwise_climdist <- dist(x = PC.Revell.env, method = "euclidian")

# Load the Jaccard index similarities to use as weights
load(file = paste0(internal.wd,"/Niche_evolution/Jaccard_index.RData"))

# Compute the observed mean climatic distance (MCD) weighted by mimicry similarity
MCD.Revell_obs <- weighted.mean(x = Revell_Pairwise_climdist, w = mimicry.similarity.weights) # 3.544
# Standardize the MCD on full species mean pairwise climatic distance
MCD.Revell_fullsp <- mean(Revell_Pairwise_climdist) # 4.374
MCD.Revell_obs <- MCD.Revell_obs/MCD.Revell_fullsp # 0.810

# Compute the stat for all the simulation
MCD.Revell_null <-  NA
for (i in 1:length(Sim.Revell.PC)) {
  Revell_Pairwise_climdist_sim <- dist(x = Sim.Revell.PC[[i]], method = "euclidian")
  MCD.Revell_sim <- weighted.mean(x = Revell_Pairwise_climdist_sim, w = mimicry.similarity.weights)/mean(Revell_Pairwise_climdist_sim)
  MCD.Revell_null[i] <- MCD.Revell_sim
}
save(MCD.Revell_obs, MCD.Revell_null, file = paste0(internal.wd,"/Niche_evolution/MCD.Revell_values.RData"))

summary(MCD.Revell_null)

# Plot the distri of the stats

jpeg(filename = paste0(internal.wd,"/Niche_evolution/Plots/MCD.pPCA_test.jpeg"), quality =100)
hist(c(MCD.Revell_null, MCD.Revell_obs), 30, freq = TRUE, col = "gray", main = "Distribution of the Mean Climatic Distance \n of co-mimetic species \n under the null hypothesis", xlab = "Standardized Mean Climatic Distance")
arrows(MCD.Revell_obs+0.004, 50, MCD.Revell_obs+0.004, 5, length = 0.1, lwd = 2)
abline(v = mean(c(MCD.Revell_null, MCD.Revell_obs)), lty = 2)
abline(v = quantile(c(MCD.Revell_null, MCD.Revell_obs), 0.05), lty = 2, col = "red")
legend(legend = c(paste0("Mean = ", round(mean(c(MCD.Revell_null, MCD.Revell_obs)),3)), 
                  paste0("CI 5% = ", round(quantile(c(MCD.Revell_null, MCD.Revell_obs), 0.05),3))), 
       x = "topleft", cex = 1, bty ="n")
legend(legend = c(paste0("MCD obs = ", round(MCD.Revell_obs, 3)),
                  paste0("p = ", ecdf(x = c(MCD.Revell_null, MCD.Revell_obs))(MCD.Revell_obs))),
       x = "left", cex = 1, bty ="n")
dev.off()