##### perMANOVA to test association between climatic niches and mimicry rings without accounting for phylogeny #####

# Test for association between climatic niches and mimicry rings without accounting for phylogeny


### Input files

# Summary tables of OMUs and species
# Phylogeny of the Ithomiini (Chazot et al., 2019)

### Output files

# Barplot nb of OMUs per rings
# Phylogeny with OMUs branch added
# Summary table only for the 619 OMUs in the phylogeny and belonging to rings with N >= 10
# pPC-axis for species and OMUs in new climatic space
# Plot OMUs and rings in pPC climatic space
# permMANOVA on 4 climatic axis and 2 pPCA-axis
# Plot distri of pseudo-F from permMANOVA
# Post-hoc pairwise tests for the permMANOVA on 4 climatic axis and 2 pPCA-axis
# Plot heatmap of p-values from pairwise tests


# Effacer l'environnement
rm(list = ls())

library(raster)

##### 1/ Load stuff #####

# Load phylogeny, sp list and unit list
load(file = paste0("./input_data/Phylogenies/Final_phylogeny.RData"))
load(file = paste0("./input_data/Summary_tables/list.sp.RData"))
load(file = paste0("./input_data/Summary_tables/list.unit.RData"))

# Extract only the 339 species included in the phylogeny from list.sp and list.unit
# Need to discard OMUs not included in the phylogeny because next analysis is a phylogenetic MANOVA and we want to keep the same dataset for both.

list.sp <- list.sp[list.sp$Sp_full %in% phylo.Ithomiini$tip.label,] # 339 species
list.unit <- list.unit[list.unit$Sp_full %in% phylo.Ithomiini$tip.label,] # 719 units


##### 2/ Build a phylogeny including the OMUs on null distance terminal branches ####

# Needed to compute Revell's pPCA

library(phytools)

phylo.Ithomiini.units <- phylo.Ithomiini

# Add units tips
?bind.tip
for (i in 1:nrow(list.unit)) {
  unit <- list.unit[i,]
  sp.nod <- which(phylo.Ithomiini.units$tip.label == unit$Sp_full)
  sp.edge <- which(phylo.Ithomiini.units$edge[,2] == sp.nod)
  phylo.Ithomiini.units <- bind.tip(tree = phylo.Ithomiini.units, tip.label = as.character(unit$Tag.model),
                                    edge.length = 0,
                                    where = sp.nod,
                                    position = 0)
  if (i %% 10 == 0) {print(i)}
}
# Remove species tips
?drop.tip
for (i in 1:nrow(list.sp)) {
  sp <- list.sp$Sp_full[i]
  phylo.Ithomiini.units <- drop.tip(phy = phylo.Ithomiini.units, tip = sp)
  if (i %% 10 == 0) {print(i)}
}

# phylo.Ithomiini.units$tip.label

save(phylo.Ithomiini.units, file = paste0("./input_data/Phylogenies/Final_units_phylogeny.RData"))
write.nexus(phylo.Ithomiini.units, file = paste0("./input_data/Phylogenies/Units_phylogeny.nex"), translate = F)  

# Reorder list.unit following Phylogeny to match rows in the phenetic distance matrix with env data in the summary table
order.index <- NA
for (i in 1:length(phylo.Ithomiini.units$tip.label)) {
  unit <- as.character(phylo.Ithomiini.units$tip.label[i])
  order.index[i] <- which(list.unit$Tag.model==unit)
}

list.unit_phyl_order_719 <- list.unit[order.index,]
row.names(list.unit_phyl_order_719) <- as.character(list.unit_phyl_order_719$Tag.model)
save(list.unit_phyl_order_719, file = paste0("./outputs/Niche_evolution/list.unit_phyl_order_719.RData"))

##### 3/ Reduce the number of circles (N min = 10 per circle) #####

table(list.unit$Mimicry.model) # Table of nb of OMU per ring
N.sp.per.ring <- table(list.unit$Mimicry.model)[order(decreasing = T, table(list.unit$Mimicry.model))] # Reorganized in decreasing order
sum(table(list.unit$Mimicry.model)>=10) # 23 circles instead of 44

# Extract name of big mimicry rings (N >= 10) kept for analysis
big.rings <- names(table(list.unit$Mimicry.model))[table(list.unit$Mimicry.model)>=10] 

# Reduce dataset to keep only OMU from big mimicry rings (N >= 10)
reduced.list.unit_phyl_order <- list.unit[list.unit$Mimicry.model %in% big.rings,] 
nrow(reduced.list.unit_phyl_order) # 619 units instead of 719
save(reduced.list.unit_phyl_order, file = paste0("./outputs/Niche_evolution/reduced.list.unit_phyl_order.RData"))

# Remove OMUs tips not in the list of big rings
?drop.tip

OMU_to_keep <- reduced.list.unit_phyl_order$Tag.model
OMU_to_drop <- setdiff(list.unit$Tag.model, OMU_to_keep)

reduced.phylo.Ithomiini.units <- phylo.Ithomiini.units
for (i in 1:length(OMU_to_drop)) {
  unit <- OMU_to_drop[i]
  reduced.phylo.Ithomiini.units <- drop.tip(phy = reduced.phylo.Ithomiini.units, tip = unit)
  if (i %% 10 == 0) {print(i)}
}

save(reduced.phylo.Ithomiini.units, file = paste0("./input_data/Phylogenies/reduced.phylo.Ithomiini.units.RData"))
write.nexus(reduced.phylo.Ithomiini.units, file = paste0("./input_data/Phylogenies/reduced.phylo.Ithomiini.units.nex"), translate = F)  


##### 4/ Barplot of nb of OMUs per ring ####

pdf(file = "./graphs/Niche_evolution/Barplot_OMU_per_ring_count.pdf", height = 6, width = 8)

original_int_margins <- par()$mar
par(mar = c(3.5,5,1,1))

barplot(height = N.sp.per.ring, ylab = "Number of OMUs", xlab = "", xaxt='n',
        cex.axis = 1.7, cex.main = 1, cex.lab = 1.7, lwd = 2)
abline(v = 28.9, lty = 2, lwd = 2, col = "red")
abline(h = 9, lty = 2, lwd = 2, col = "red")
polygon(x = c(28.9, 28.9, 52.85, 52.85), y = c(0, 85, 85, 0), col = "#FF000050", border = NA)

title(xlab = "Mimicry rings", line = 1.3, cex.lab = 1.7)

par(mar = original_int_margins)

dev.off()


##### 5/ permMANOVA without phylogenetic correction #####

load(file = paste0("./outputs/Niche_evolution/reduced.list.unit_phyl_order.RData"))

# 5.1/ Extract mean data for the 4 bioclim variables for each OMU

unit.env.table_619 <- reduced.list.unit_phyl_order[,c("bio1","bio4","bio12","bio15")]
names(unit.env.table_619) <- c("Tmean", "Tvar", "Hmean", "Hvar")
save(unit.env.table_619, file = paste0("./outputs/Niche_evolution/unit.env.table_619.RData"))

# 5.2/ Run the perMANOVA

load(file = paste0("./outputs/Niche_evolution/unit.env.table_619.RData"))

library("vegan")

?adonis
?adonis2

permMANOVA.4 <- adonis(formula = as.matrix(scale(unit.env.table_619)) ~ reduced.list.unit_phyl_order$Mimicry.model, # Use standardized climatic space
                       permutations = 999,
                       method = "euclidian", # Use euclidian distances in standardized climatic space
                       strata = NULL, by = "margin")  # By = "terms" for Anova type I ; By = "margin" for Anova type II. No differences for a single predictive variable

print(permMANOVA.4) # output from the PERMANOVA
permMANOVA.4_stats <- permustats(permMANOVA.4) # Extract details about permutation test
summary(permMANOVA.4_stats)

save(permMANOVA.4, permMANOVA.4_stats, file = paste0("./outputs/Niche_evolution/permMANOVA/permMANOVA.4.RData"))


# 5.3/ Plot null distri of pseudo-F stat

load(file = paste0("./outputs/Niche_evolution/permMANOVA/permMANOVA.4.RData"))

pseudo_F_obs <- permMANOVA.4_stats$statistic
pseudo_F_null <- permMANOVA.4_stats$permutations[,1]

save(pseudo_F_obs, pseudo_F_null, file = paste0("./outputs/Niche_evolution/permMANOVA/permMANOVA.4_plot_stuff.RData"))

pdf(file = "./graphs/Niche_evolution/permMANOVA/perMANOVA_pseudo_F_null.pdf", height = 5.5, width = 5.5)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,5,3,1))

hist(x = log(c(pseudo_F_obs, pseudo_F_null)), 
     breaks = 20,
     freq = TRUE, col = "gray",
     # xlim = c(0, 20),
     # ylim = c(0, 400),
     main = "Distribution of pseudo-F from perMANOVA\nunder null Hypothesis",
     # main = "",
     xlab = "pseudo-F (log scale)",
     axes = F,
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.2, lwd = 2)

axis(side = 1, lwd = 2, cex.axis = 1.7, labels = c(0, 1, 3, 8, 20), at = seq(-1,3,1))
axis(side = 2, lwd = 2, cex.axis = 1.7)

arrows(log(pseudo_F_obs) - 0.065, 130, log(pseudo_F_obs) - 0.065, 10, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = log(mean(c(pseudo_F_null, pseudo_F_obs))), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = log(quantile(c(pseudo_F_null, pseudo_F_obs), 0.95)), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0(" Mean = ", format(round(mean(c(pseudo_F_null, pseudo_F_obs)),3), nsmall = 3)), 
                  paste0("CI 5% = ", round(quantile(c(pseudo_F_null, pseudo_F_obs), 0.95),3))), 
       x = "topright", inset = c(0, 0.15), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.6, bty = "n")
legend(legend = c(paste0("pseudo-F obs = ", round(pseudo_F_obs, 3)),
                  paste0("                     p = 0.001")),
       x = "bottomright", inset = c(0.00, 0.40),
       cex = 1.6, bty ="n", xjust = 1)

legend(legend = as.expression(bquote(bold("B"))), 
       x = "topright", inset = c(0.03, -0.03), xjust = 0.5,
       cex = 1.9, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()

##### 6/ Post-hoc test ####

# if significant to detect between each modalities of the factor the difference is significant

load(file = paste0("./outputs/Niche_evolution/unit.env.table_619.RData"))

library("RVAideMemoire")

?pairwise.perm.manova
# To look for the different statistical test available, see ?anova.mlm
# To see the different available correction of p-value for multiple testing, see ?p.adjust
Pairwise.permMANOVA.4 <- pairwise.perm.manova(resp = as.matrix(scale(unit.env.table_619)),
                                              fact = reduced.list.unit_phyl_order$Mimicry.model,
                                              nperm = 999, progress = T,
                                              p.method = "none", # No p-value correction because too many tests and p-value limited to 0.001
                                              test = "Wilks")
Pairwise.permMANOVA.4
save(Pairwise.permMANOVA.4, file = paste0("./outputs/Niche_evolution/permMANOVA/Pairwise.permMANOVA.4.RData"))

# 7/ Plot a heatmap of p-values for post-hoc tests #####

load(file = paste0("./outputs/Niche_evolution/permMANOVA/Pairwise.permMANOVA.4.RData"))

# 7.1/ Extract p-value to build a heatmap
extract_mask <- lower.tri(Pairwise.permMANOVA.4$p.value) ; diag(extract_mask) <- T
Pairwise.permMANOVA.4_pvalues <- Pairwise.permMANOVA.4$p.value[extract_mask]

length(Pairwise.permMANOVA.4_pvalues) # 253 pairwise p-values
sum(Pairwise.permMANOVA.4_pvalues == 0.001) # 186 out of 253
sum(Pairwise.permMANOVA.4_pvalues < 0.01) # 210 out of 253
sum(Pairwise.permMANOVA.4_pvalues < 0.05) # 226 out of 253

hist(Pairwise.permMANOVA.4_pvalues)
hist(log1p(Pairwise.permMANOVA.4_pvalues))
hist(sqrt(sqrt(Pairwise.permMANOVA.4_pvalues)))

# Try to adjust p-value for multiple testing with Holm correction. Non-sense because all p-value become non-significant due to the fact the smallest value possible is 0.001.
p.adjust(p = Pairwise.permMANOVA.4_pvalues, method = "holm")

# 7.2/ Set mimicry names properly
mimicry_ring_list <- c(attr(Pairwise.permMANOVA.4$p.value, "dimnames")[[2]], attr(Pairwise.permMANOVA.4$p.value, "dimnames")[[1]][22])
mimicry_ring_list[which(mimicry_ring_list == "THABENAF")] <- "THABENA-F"
mimicry_ring_list[which(mimicry_ring_list == "BANJANAM")] <- "BANJANA-M"

# 7.3/ Build a full matrix of pairwise p-values
library(gdata)
Pseudo_Cor_mat_perMANOVA <- matrix(nrow = 23, ncol = 23, dimnames = list(mimicry_ring_list, mimicry_ring_list))
Pseudo_Cor_mat_perMANOVA[2:23, 1:22] <- Pairwise.permMANOVA.4$p.value
upperTriangle(Pseudo_Cor_mat_perMANOVA) <- lowerTriangle(Pseudo_Cor_mat_perMANOVA, byrow=TRUE)
diag(Pseudo_Cor_mat_perMANOVA) <- 1

save(Pseudo_Cor_mat_perMANOVA, file = paste0("./outputs/Niche_evolution/permMANOVA/Pseudo_Cor_mat_perMANOVA.RData"))

# 7.4/ Compute the dendrogram outside of heatmap() from climatic values and set branch width manually
library(rlang)
library(tidyverse)
library(dendextend)

ring_env_table <- reduced.list.unit_phyl_order[,c("Mimicry.model","bio1","bio4","bio12","bio15")]
ring_env_table <- ring_env_table %>% 
  group_by(Mimicry.model) %>%
  summarise(Tmean = mean(bio1),
            Tvar = mean(bio4),
            Hmean = mean(bio12),
            Hvar = mean(bio15))
ring_env_mat <- as.matrix(ring_env_table[, c("Tmean","Tvar","Hmean","Hvar")])

ring_clim_dist <- dist((scale(ring_env_mat)))

dd <- set(as.dendrogram(hclust(ring_clim_dist, method = "complete")), "branches_lwd", 2)
plot(dd)

# Reorder leaves of the dendrogram manually
dd_reorder <- reorder(x = dd, wts = rowMeans(x = scale(ring_env_mat)))
plot(dd_reorder)

save(dd_reorder, file = paste0("./outputs/Niche_evolution/permMANOVA/dd_reorder.RData"))


# 7.5/ Plot de la heatmap

# Set color palette
heatmap_pal <- rev(tmaptools::get_brewer_pal("RdYlBu", n = 100))[10:95]

pdf(file = "./graphs/Niche_evolution/permMANOVA/Heatmap_PerMANOVA.pdf", height = 6, width = 8)

heatmap_posthoc_permMANOVA <- heatmap(x = Pseudo_Cor_mat_perMANOVA^(1/5.8),
                                      symm = T, Rowv = rev(dd_reorder), Colv = rev(dd_reorder), revC = T,
                                      col = heatmap_pal)
dev.off()

# Reorder matrix of p-values to fit the order on the heat map based on the dendrogram
ordered.Pseudo_Cor_mat_perMANOVA <- Pseudo_Cor_mat_perMANOVA[heatmap_posthoc_permMANOVA$rowInd, heatmap_posthoc_permMANOVA$colInd]

# Plot de la scale
library(raster)

pdf(file = "./graphs/Niche_evolution/permMANOVA/Heatmap_PerMANOVA_scale.pdf", height = 6, width = 8)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,5), mar = c(0,0,0,5))
# Set margin parameters: Bottom, left, top, right

plot(raster(x = ordered.Pseudo_Cor_mat_perMANOVA), col = heatmap_pal,
     axis.args=list(at = c(0,0.2,0.4,0.6,0.8,1), cex.axis = 1.6, font = 2,
                    labels = c("0.000", "0.001", "0.005", "0.05", "0.3", "1")),
     legend.args = list(text = "               p-value", side = 3,
                        font = 2, line = 1.8, cex = 1.6))

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()


##### 8/ Plot OMUs and rings in climatic space after PCA #####

load(file = paste0("./input_data/Phylogenies/reduced.phylo.Ithomiini.units.RData"))
load(file = paste0("./outputs/Niche_evolution/reduced.list.unit_phyl_order.RData"))
load(file = paste0("./outputs/Niche_evolution/unit.env.table_619.RData"))

# 8.1/ Try Revell's pPCA on OMUs to reduce number of axis to 2

source(paste0("./scripts/Revell_phyl_pca.R"))
library(ape)

C <- vcv(reduced.phylo.Ithomiini.units)
X <- as.matrix(unit.env.table_619)

PCA.Revell <- Revell_phyl_pca(C, X, mode = "corr")
# Not working because the matrix C is not invertible
# because OMUs from the same species have the same phylogenetic distances/covariance values towards other taxa,
# which make columns redondant, thus the matrix not invertible

# Solution = apply the pPCA on species, and retrieve the eigenvectors, ancestral states, and evolutionary VCV matrix to project environmental data into the pPCA space of species !

# 8.2/ Prepare environmental and phylogenetic data

# Get environmental data for species
sp_env_table <- reduced.list.unit_phyl_order[,c("Sp_full","bio1","bio4","bio12","bio15")]
sp_env_table <- sp_env_table %>% 
  group_by(Sp_full) %>%
  summarise(Tmean = mean(bio1),
            Tvar = mean(bio4),
            Hmean = mean(bio12),
            Hvar = mean(bio15))

# Remove species with no OMUs included in the permMANOVA analyses due to belonging to small mimicry rings
library(phytools)
?drop.tip

small_rings_sp <- setdiff(phylo.Ithomiini$tip.label, sp_env_table$Sp_full)
phylo.Ithomiini_permMANOVA <- phylo.Ithomiini
for (i in 1:length(small_rings_sp)) {
  sp <- small_rings_sp[i]
  phylo.Ithomiini_permMANOVA <- drop.tip(phy = phylo.Ithomiini_permMANOVA, tip = sp)
  print(i)
}

save(phylo.Ithomiini_permMANOVA, file = paste0("input_data/Phylogenies/phylo.Ithomiini_permMANOVA.RData"))
write.nexus(phylo.Ithomiini_permMANOVA, file = paste0("input_data/Phylogenies/phylo.Ithomiini_permMANOVA.nex"), translate = F)  

# Reorder sp_env_table following Phylogeny to match rows in the phenetic distance matrix with env data in the summary table
order.index <- NA
for (i in 1:length(phylo.Ithomiini_permMANOVA$tip.label)) {
  sp <- as.character(phylo.Ithomiini_permMANOVA$tip.label[i])
  order.index[i] <- which(sp_env_table$Sp_full == sp)
}

sp_env_table <- sp_env_table[order.index,]
row.names(sp_env_table) <- as.character(sp_env_table$Sp_full)

# Generate a sp env matrix from table
sp_env_mat <- as.matrix(sp_env_table[, c("Tmean","Tvar","Hmean","Hvar")])
row.names(sp_env_mat) <- sp_env_table$Sp_full
nrow(sp_env_mat) # 322 species in the phylogeny belonging to "big" rings (N >= 10)

save(sp_env_table, file = paste0("./outputs/Niche_evolution/permMANOVA/sp_env_table.RData"))
save(sp_env_mat, file = paste0("./outputs/Niche_evolution/permMANOVA/sp_env_mat.RData"))

# 8.3/ Run Revell's pPCA

load(file = paste0("input_data/Phylogenies/phylo.Ithomiini_permMANOVA.RData"))

load(file = paste0("./outputs/Niche_evolution/permMANOVA/sp_env_table.RData"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/sp_env_mat.RData"))

source(paste0("./scripts/Revell_phyl_pca.R"))

# Create input matrices
C <- vcv(phylo.Ithomiini_permMANOVA) # Variance-covariance matrix based on species phylogeny
X <- sp_env_mat # Environmental data matrix

# Check the rows in the two matrices match each other
identical(row.names(C), row.names(X))

PCA.Revell <- Revell_phyl_pca(C, X, mode = "corr")
diag(PCA.Revell$Eval) # eigenvalues
PCA.Revell$Evec # eigenvectors used to project new objects in the new pPCA space
PCA.Revell$S # scores = coordinates of objects in the new space
PCA.Revell$L # PC loadings
PCA.Revell$ancestral # Ancestral states of traits
PCA.Revell$evolVCV # Evolutionary variance/covariance matrix
PCA.Revell$evol_corr # Evolutionary correlation matrix
PCA.Revell$X # Standardized trait matrix

PCvar <- round(diag(PCA.Revell$Eval)/sum(PCA.Revell$Eval),3) ; PCvar # % Explained variance per axis
cumVar <- round(cumsum(diag(PCA.Revell$Eval))/sum(PCA.Revell$Eval),3) ; cumVar # % Cumulative explained variance
# 97% of variance explained with the 2 first axes!

save(PCA.Revell, PCvar, file = paste0("./outputs/Niche_evolution/permMANOVA/PCA.Revell.RData"))

library(ade4)
s.corcircle(PCA.Revell$L, xax=1, yax=2, label = colnames(sp_env_mat)) # Quick correlation circle

# Extract coordinates for species for comparison and later use
pPC.env_sp_322 <- PCA.Revell$S[,1:2]
save(pPC.env_sp_322, file = paste0("./outputs/Niche_evolution/permMANOVA/pPC.env_sp_322.RData"))

# Project OMUs in the new pPCA env space

# Standardize env traits following Revell's procedure
X <- as.matrix(unit.env.table_619) # Env data for OMU
V <- PCA.Revell$evolVCV        # Evolutionary variance/covariance matrix of environmental traits based on sp phenetic distances and env data
a <- PCA.Revell$ancestral      # Ancestral state of env traits based on sp phenetic distances and env data
eigenV <- PCA.Revell$Evec      # Eigenvectors to project in new space
one <- matrix(1,nrow(X),1)     # Vector to get means
X <- X/(one %*% t(sqrt(diag(V))))  # Standardized env data
pPC.env_units <- (X-one%*%t(a)) %*% eigenV # Project OMu data in new space

save(pPC.env_units, file = paste0("./outputs/Niche_evolution/permMANOVA/pPC.env_units.RData"))

# Get correlation values between original env variables and PC axis
pPCA_correl <- PCA.Revell$L[ ,1:2] ; row.names(pPCA_correl) <- row.names(a)
save(pPCA_correl, file = paste0("./outputs/Niche_evolution/permMANOVA/pPCA_correl.RData"))

# # Quick plot with all rings
# plot(pPC.env_units, col = "grey80", main = "pPCA of OMUs in climatic space", xlab = paste0("pPC1 (", PCvar[1]*100 ," %)"), ylab = paste0("pPC2 (", PCvar[2]*100 ," %)"))
# points(pPC.env_units, pch = 16, col = rainbow(44)[as.factor(reduced.list.unit_phyl_order$Mimicry.model)])
# arrows(x0 = rep(0, nrow(pPCA_correl)), y0 = rep(0, nrow(pPCA_correl)), x1 = pPCA_correl[,1]*3, y1 =  pPCA_correl[,2]*3, length = 0.2, code = 2, lwd = 1)
# text(pPCA_correl[,]*3.5, labels = row.names(pPCA_correl))

# Select only a reduced number of rings
Mimicry_counts_table <- table(reduced.list.unit_phyl_order$Mimicry.model)[order(table(reduced.list.unit_phyl_order$Mimicry.model), decreasing = T)] ; Mimicry_counts_table

# Extract only the 5 or 9 main circles in colors

top5.ring.list <- names(Mimicry_counts_table)[1:5]
# top9.ring.list <- names(Mimicry_counts_table)[1:9]

# AGNOSIA 74 # LERIDA 63 # MAMERCUS 56 # HERMIAS 47 # BANJANA-M 43 # PANTHYALE 38 # DULICIDA 35 # EURIMEDIA 33 # HEWITSONI 27 

top5.index <- (reduced.list.unit_phyl_order$Mimicry.model %in% top5.ring.list)
# top9.index <- (reduced.list.unit_phyl_order$Mimicry.model %in% top9.ring.list) 
save(top5.index, file = paste0("./outputs/Niche_evolution/permMANOVA/top5.index.RData"))


# # Build minimum convex polygon for the 9 main circles
# 
# library(adehabitatHR)
# library(sp)
# 
# ?mcp
# 
# for (i in 1:length(top9.ring.list)) 
# {
#   ring <- top9.ring.list[i]
#   ring.index <- reduced.list.unit_phyl_order$Mimicry.model == ring
#   SPobj.ring <- SpatialPoints(coords = pPC.env_units[ring.index, 1:2]) # Create a Spobj from the coordinates
#   mcp.ring <- mcp(xy = SPobj.ring, percent = 100) # Generate the minimum convex polygon
#   save(mcp.ring, file = paste0("./outputs/Niche_evolution/permMANOVA/mcp/mcp.ring_", ring, ".RData"))
# }

# Plot with only top 5 rings

mimic.list <- levels(factor(reduced.list.unit_phyl_order$Mimicry.model[top5.index]))
mimic.list_legend <- mimic.list ; mimic.list_legend[which(mimic.list_legend == "BANJANAM")] <- "BANJANA-M"
save(mimic.list, mimic.list_legend, file = paste0("./outputs/Niche_evolution/permMANOVA/mimic_lists_for_pPCA_plot.RData"))

load(file = paste0("./outputs/Niche_evolution/reduced.list.unit_phyl_order.RData"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/PCA.Revell.RData"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/pPC.env_units.RData"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/pPCA_correl.RData"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/mimic_lists_for_pPCA_plot.RData"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/top5.index.RData"))


pdf(file = "./graphs/Niche_evolution/permMANOVA/pPCA_Top5.pdf", height = 7.5, width = 8)

original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))

plot(pPC.env_units, col = "grey90", main = "pPCA of OMUs in climatic space", 
     xlab = paste0("pPC1 (", round(PCvar[1]*100,1) ," %)"), ylab = paste0("pPC2 (", round(PCvar[2]*100,1) ," %)"),
     pch = 16, type = "p", axes = F, ylim = c(-6, 6), xlim = c(-6, 6),
     cex.axis = 1.5, cex.lab = 1.6, cex.main = 1.3)
axis(side = 1, lwd = 2, cex.axis = 1.5)
axis(side = 2, lwd = 2, cex.axis = 1.5)
points(pPC.env_units[top5.index,], pch = 16, col = c("black","red","deepskyblue","limegreen","orange","blue")[factor(reduced.list.unit_phyl_order$Mimicry.model[top5.index])])
legend(x = "bottomright", pch = 16, cex = 1.1, legend = mimic.list_legend, 
       bty = "n", inset = c(0.02, 0.00), # horiz = T, 
       col = c("black","red","deepskyblue","limegreen","orange","blue"))
arrows(x0 = rep(0, nrow(pPCA_correl)), y0 = rep(0, nrow(pPCA_correl)), x1 = pPCA_correl[,1]*4, y1 =  pPCA_correl[,2]*4, length = 0.2, code = 2, lwd = 1)

# for (i in 1:length(mimic.list_legend)) 
# {
#   ring <- mimic.list[i]
#   ring_legend <- mimic.list_legend[i]
#   load(file = paste0("./outputs/Niche_evolution/permMANOVA/mcp/mcp.ring_", ring, ".RData"))
#   plot(mcp.ring, add = T, border = c("black","red","deepskyblue","limegreen","orange","blue")[i], lwd = 2)
# }

text(pPCA_correl[,1:2]*4.2 + cbind(c(0, -0.8, 0.0, -0.4), c(-0.2, -0.2, -0.2, -0.4)), labels = row.names(pPCA_correl), font = 2, adj = c(0, 0))

legend(legend = as.expression(bquote(bold("B"))), 
       x = "topright", inset = c(0.05, 0.00), xjust = 0.5,
       cex = 1.3, bty ="n")

par(mar = original_int_margins)

dev.off()


##### 9/ Plot both together #####

load(file = paste0("./outputs/Niche_evolution/permMANOVA/permMANOVA.4_plot_stuff.RData"))

load(file = paste0("./outputs/Niche_evolution/reduced.list.unit_phyl_order.RData"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/PCA.Revell.RData"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/pPC.env_units.RData"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/pPCA_correl.RData"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/mimic_lists_for_pPCA_plot.RData"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/top5.index.RData"))

round(pseudo_F_obs, 2) # pseudo-F obs = 19.35

pdf(file = "./graphs/Niche_evolution/permMANOVA/permMANOVA_both_plots.pdf", height = 7, width = 14)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,5,1.5,1), mfrow = c(1,2))

plot(pPC.env_units, col = "grey90", main = "", 
     xlab = paste0("pPC1 (", round(PCvar[1]*100,1) ," %)"), ylab = paste0("pPC2 (", round(PCvar[2]*100,1) ," %)"),
     pch = 16, type = "p", axes = F, ylim = c(-6, 6), xlim = c(-6, 6),
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.3)
axis(side = 1, lwd = 2, cex.axis = 1.7)
axis(side = 2, lwd = 2, cex.axis = 1.7)

points(pPC.env_units[top5.index,], pch = 16, col = c("black","red","deepskyblue","limegreen","orange","blue")[factor(reduced.list.unit_phyl_order$Mimicry.model[top5.index])])

legend(x = "bottomright", pch = 16, cex = 1.3, legend = mimic.list_legend[3:5], 
       bty = "n", inset = c(0.02, 0.00), # horiz = T, 
       col = c("black","red","deepskyblue","limegreen","orange","blue")[3:5])
legend(x = "bottomright", pch = 16, cex = 1.3, legend = mimic.list_legend[1:2], 
       bty = "n", inset = c(0.32, 0.00), # horiz = T, 
       col = c("black","red","deepskyblue","limegreen","orange","blue")[1:2])

arrows(x0 = rep(0, nrow(pPCA_correl)), y0 = rep(0, nrow(pPCA_correl)), x1 = pPCA_correl[,1]*4.5, y1 =  pPCA_correl[,2]*4.5, length = 0.2, code = 2, lwd = 1)

# for (i in 1:length(mimic.list_legend)) 
# {
#   ring <- mimic.list[i]
#   ring_legend <- mimic.list_legend[i]
#   load(file = paste0("./outputs/Niche_evolution/permMANOVA/mcp/mcp.ring_", ring, ".RData"))
#   plot(mcp.ring, add = T, border = c("black","red","deepskyblue","limegreen","orange","blue")[i], lwd = 2)
# }

text(pPCA_correl[,1:2]*4.7 + cbind(c(0, -0.9, 0.0, -0.6), c(-0.3, -0.3, -0.2, -0.4)),
     labels = row.names(pPCA_correl),
     cex = 1.2, font = 2, adj = c(0, 0))

legend(legend = as.expression(bquote(bold("A"))), 
       x = "topright", inset = c(0.05, -0.03), xjust = 0.5,
       cex = 1.9, bty ="n")

hist(x = log(c(pseudo_F_obs, pseudo_F_null)), 
     breaks = 20,
     freq = TRUE, col = "gray",
     # xlim = c(0, 20),
     # ylim = c(0, 400),
     main = "",
     xlab = "pseudo-F (log scale)",
     axes = F,
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.2, lwd = 2)

axis(side = 1, lwd = 2, cex.axis = 1.7, labels = c(0, 1, 3, 8, 20), at = seq(-1,3,1))
axis(side = 2, lwd = 2, cex.axis = 1.7)

arrows(log(pseudo_F_obs) - 0.065, 135, log(pseudo_F_obs) - 0.065, 10, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = log(mean(c(pseudo_F_null, pseudo_F_obs))), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = log(quantile(c(pseudo_F_null, pseudo_F_obs), 0.95)), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0("Mean = ", format(round(mean(c(pseudo_F_null, pseudo_F_obs)),2), nsmall = 2)), 
                  paste0("CI 95% = ", round(quantile(c(pseudo_F_null, pseudo_F_obs), 0.95), 2))), 
       x = "topright", inset = c(0, 0.15), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.6, bty = "n")

legend(legend = as.expression(bquote(bold("pseudo-F obs = 19.35"))),
       x = "bottomright", inset = c(0.00, 0.455),
       cex = 1.6, bty ="n", xjust = 1)

legend(legend = as.expression(bquote(bold("p = 0.001"))),
       x = "bottomright", inset = c(0.00, 0.40),
       cex = 1.6, bty ="n", xjust = 1)


legend(legend = as.expression(bquote(bold("B"))), 
       x = "topright", inset = c(0.03, -0.03), xjust = 0.5,
       cex = 1.9, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins, mfrow = c(1,1))

dev.off()


##### 10/ Bonus : permMANOVA on 2 Revell's pPC axis ####

# Load Revell's pPCA axis for units/OMUs
load(file = paste0("./outputs/Niche_evolution/permMANOVA/pPC.env_units.RData"))
# Load summary table for the 619 OMUs in the top 10 rings present in the phylogeny
load(file = paste0("./outputs/Niche_evolution/reduced.list.unit_phyl_order.RData"))

colnames(pPC.env_units) <- c("Tmean", "Tvar", "Htot", "Hvar")
reduced.list.unit_phyl_order

library("vegan")

?adonis
?adonis2

permMANOVA.pPC2 <- adonis(formula = pPC.env_units ~ reduced.list.unit_phyl_order$Mimicry.model,
                       permutations = 999,
                       method = "euclidian", # Use euclidian distances in standartized climatic space
                       strata = NULL, by = "margin")  # By = "terms" for Anova type I ; By = "margin" for Anova type II. No differences for a single predictive variable

print(permMANOVA.pPC2) # output from the PERMANOVA
permMANOVA.pPC2_stats <- permustats(permMANOVA.pPC2) # Extract details about permutation test
summary(permMANOVA.pPC2_stats)

save(permMANOVA.pPC2, permMANOVA.pPC2_stats, file = paste0("./outputs/Niche_evolution/permMANOVA/permMANOVA.pPC2.RData"))

load(file = paste0("./outputs/Niche_evolution/permMANOVA.pPC2.RData"))

# if significant to detect between each modalities of the factor the difference is significant

library("RVAideMemoire")

?pairwise.perm.manova
# To look for the different statistical test available, see ?anova.mlm
# To see the different available correction of p-value for multiple testing, see ?p.adjust
Pairwise.permMANOVA.pPC2 <- pairwise.perm.manova(resp = pPC.env_units,
                                              fact = reduced.list.unit_phyl_order$Mimicry.model,
                                              nperm = 999, progress = T,
                                              p.method = "none", # No p-value correction because too many tests and p-value limited to 0.001
                                              test = "Wilks")
Pairwise.permMANOVA.pPC2
save(Pairwise.permMANOVA.pPC2, file = paste0("./outputs/Niche_evolution/permMANOVA/Pairwise.permMANOVA.pPC2.RData"))

# Extract p-value (to build a heatmap)
extract_mask <- lower.tri(Pairwise.permMANOVA.pPC2$p.value) ; diag(extract_mask) <- T
Pairwise.permMANOVA.pPC2_pvalues <- Pairwise.permMANOVA.pPC2$p.value[extract_mask]

length(Pairwise.permMANOVA.pPC2_pvalues) # 253 pairwise p-values
sum(Pairwise.permMANOVA.pPC2_pvalues == 0.001) # 186 out of 253
sum(Pairwise.permMANOVA.pPC2_pvalues < 0.01) # 209 out of 253
sum(Pairwise.permMANOVA.pPC2_pvalues < 0.05) # 226 out of 253
