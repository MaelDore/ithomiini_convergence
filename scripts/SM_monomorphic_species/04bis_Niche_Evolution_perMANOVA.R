##### Script 04bis: perMANOVA to test association between climatic niches and mimicry rings without accounting for phylogeny #####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

# Test for association between climatic niches and mimicry rings without accounting for phylogeny

### Focus only on monomorphic/species => 1 OMU = 1 species

### Analyses limited to mimicry ring with N >= 5 => 106 monomorphic species/OMUs

### Input files

# Summary tables of OMUs and species
# Phylogeny of the Ithomiini (Chazot et al., 2019)

### Output files

# Barplot nb of OMUs per rings
# Phylogeny with the 106/139 monomorphic species included in the phylogeny
# Summary table with the 106/139 monomorphic species included in the phylogeny
# Extract pPC-axis for monomorphic species in new climatic space
# Plot OMUs and rings in pPC climatic space
# permMANOVA on 2 pPCA-axis
# Plot distri of pseudo-F from permMANOVA
# Post-hoc pairwise tests for the permMANOVA on 2 pPCA-axis
# Plot heatmap of p-values from pairwise tests


# Clean environment
rm(list = ls())

library(raster)

##### 1/ Load stuff #####

# Load phylogeny and OMU list for monomorphic species
load(file = paste0("./input_data/Phylogenies/Final_phylogeny.RData"))
list.unit_mono_178 <- readRDS(file = "./input_data/Summary_tables/list.unit_mono_178.rds")

#### Need to extract the pPCA-axis from the coordinates calibrated on the whole tribe data !!!

##### 2/ Extract data for the different subset of species of the phylogeny including only the 139 monomorphic species ####

# Extract only the 139 species included in the phylogeny from list.unit
# Need to discard monomorphic species/OMUs not included in the phylogeny because next analysis is a phylogenetic MANOVA and we want to keep the same dataset for both.

list.unit_mono_139 <- list.unit_mono_178[list.unit_mono_178$Sp_full %in% phylo.Ithomiini$tip.label, ] # 339 species

# Subset phylogeny is needed to simulated trait evolution 
# But in practice, we substract data from the simulations on complete dataset

library(ape)

# Keep only the monomorphic species
phylo.Ithomiini_mono_139 <- ape::keep.tip(phy = phylo.Ithomiini, tip = list.unit_mono_139$Sp_full)
length(phylo.Ithomiini_mono_139$tip.label)

# Tranform tip label into OMUs
phylo.Ithomiini_mono_139$tip.label <- list.unit_mono_139$Tag.model[match(phylo.Ithomiini_mono_139$tip.label, list.unit_mono_139$Sp_full)]

# Save
saveRDS(phylo.Ithomiini_mono_139, file = paste0("./input_data/Phylogenies/phylo.Ithomiini_mono_139.rds"))
write.nexus(phylo.Ithomiini_mono_139, file = paste0("./input_data/Phylogenies/phylo.Ithomiini_mono_139.nex"), translate = F)  

# Reorder list.unit following Phylogeny to match rows in the phenetic distance matrix with env data in the summary table

order_index <- NA
for (i in 1:length(phylo.Ithomiini_mono_139$tip.label)) {
  unit <- as.character(phylo.Ithomiini_mono_139$tip.label[i])
  order_index[i] <- which(list.unit_mono_139$Tag.model == unit)
}

list.unit_phyl_order_mono_139 <- list.unit_mono_139[order_index,]
row.names(list.unit_phyl_order_mono_139) <- as.character(list.unit_phyl_order_mono_139$Tag.model)
saveRDS(list.unit_phyl_order_mono_139, file = paste0("./outputs/Niche_evolution/Mono_sp/list.unit_phyl_order_mono_139.rds"))


##### 3/ Barplot of nb of monomorphic species/OMUs per ring ####

table(list.unit_phyl_order_mono_139$Mimicry.model) # Table of nb of OMU per ring
N.sp.per.ring <- table(list.unit_phyl_order_mono_139$Mimicry.model)[order(decreasing = T, table(list.unit_phyl_order_mono_139$Mimicry.model))]
N.sp.per.ring

pdf(file = "./graphs/Niche_evolution/Mono_139/Barplot_OMU_per_ring_count.pdf", height = 6, width = 8)

original_int_margins <- par()$mar
par(mar = c(3.5,5,1,1))

barplot(height = N.sp.per.ring, ylab = "Number of species", xlab = "", xaxt='n',
        main = "Monomorphic species\n(n = 139)\nper mimicry ring",
        cex.axis = 1.7, cex.main = 1.5, cex.lab = 1.7, lwd = 2)
abline(v = 10.95, lty = 2, lwd = 2, col = "red")
abline(h = 4, lty = 2, lwd = 2, col = "red")
polygon(x = c(11, 11, 53, 53), y = c(0, 4, 4, 0), col = "#FF000050", border = NA)

title(xlab = "Mimicry rings", line = 1.3, cex.lab = 1.7)

par(mar = original_int_margins)

dev.off()

##### 4/ Reduce the number of circles (N min = 5 per circle) #####

### 4.1/ Identify big rings (N >=5) ###

table(list.unit_phyl_order_mono_139$Mimicry.model) # Table of nb of OMU per ring
N.sp.per.ring <- table(list.unit_phyl_order_mono_139$Mimicry.model)[order(decreasing = T, table(list.unit_phyl_order_mono_139$Mimicry.model))]
N.sp.per.ring
sum(table(list.unit_phyl_order_mono_139$Mimicry.model) >= 5) # 9 circles instead of 30 (44)

# Extract names of big mimicry rings (N >= 5) kept for analysis
big_rings <- names(table(list.unit_phyl_order_mono_139$Mimicry.model))[table(list.unit_phyl_order_mono_139$Mimicry.model) >= 5] 

# Extract names of OMUs in big rings
OMU_mono_106 <- list.unit_phyl_order_mono_139$Tag.model[list.unit_phyl_order_mono_139$Mimicry.model %in% big_rings]

### 4.2/ Extract subtree for the 106 monomorphic species ###

phylo.Ithomiini_mono_106 <- ape::keep.tip(phy = phylo.Ithomiini_mono_139, tip = OMU_mono_106)
length(phylo.Ithomiini_mono_106$tip.label)

# Save
saveRDS(phylo.Ithomiini_mono_106, file = paste0("./input_data/Phylogenies/phylo.Ithomiini_mono_106.rds"))
write.nexus(phylo.Ithomiini_mono_106, file = paste0("./input_data/Phylogenies/phylo.Ithomiini_mono_106.nex"), translate = F)  

# Load the reduced phylogeny for monomorphic species in "big" rings (N >= 5)
phylo.Ithomiini_mono_106 <- readRDS(file = paste0("./input_data/Phylogenies/phylo.Ithomiini_mono_106.rds"))

### 4.3/ Reduce dataset to keep only OMU from big mimicry rings (N >= 5) ###

reduced.list.unit_phyl_order_mono_106 <- list.unit_phyl_order_mono_139[match(phylo.Ithomiini_mono_106$tip.label, list.unit_phyl_order_mono_139$Tag.model), ]
nrow(reduced.list.unit_phyl_order_mono_106) # 106 units instead of 139 (719 / 783)

identical(phylo.Ithomiini_mono_106$tip.label, reduced.list.unit_phyl_order_mono_106$Tag.model)

# Save
saveRDS(reduced.list.unit_phyl_order_mono_106, file = paste0("./outputs/Niche_evolution/Mono_sp/reduced.list.unit_phyl_order_mono_106.rds"))


##### 5/ permMANOVA without phylogenetic correction #####

# Load the reduced unit list for monomorphic species in the phylogeny, in "big" rings (N >= 5)
reduced.list.unit_phyl_order_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/reduced.list.unit_phyl_order_mono_106.rds"))

### 5.1/ Extract mean environmental data for monophonic species in big rings from the Revell's pPCA space

# Load PCA's space for the 619 OMUs in big rings
load(file = paste0("./outputs/Niche_evolution/permMANOVA/pPC.env_units.RData"))

# Extract data for the 106 monomorphic species in the phylogeny order
pPC.env_units_mono_106 <- pPC.env_units[match(reduced.list.unit_phyl_order_mono_106$Tag.model, row.names(pPC.env_units)), ]

# Save
saveRDS(pPC.env_units_mono_106, file = paste0("./outputs/Niche_evolution/Mono_sp/pPC.env_units_mono_106.rds"))

### 5.2/ Run the perMANOVA

# Load environmental pPCA data
pPC.env_units_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/pPC.env_units_mono_106.rds"))

library("vegan")

?adonis
?adonis2

permMANOVA_pPCA_mono_106 <- adonis(formula = as.matrix(scale(pPC.env_units_mono_106[, c(1,2)])) ~ reduced.list.unit_phyl_order_mono_106$Mimicry.model, # Use standardized climatic space
                                   permutations = 999,
                                   method = "euclidian", # Use euclidian distances in standardized climatic space
                                   strata = NULL, by = "margin")  # By = "terms" for Anova type I ; By = "margin" for Anova type II. No differences for a single predictive variable

print(permMANOVA_pPCA_mono_106) # output from the PERMANOVA
permMANOVA_pPCA_mono_106_stats <- permustats(permMANOVA_pPCA_mono_106) # Extract details about permutation test
summary(permMANOVA_pPCA_mono_106_stats)

save(permMANOVA_pPCA_mono_106, permMANOVA_pPCA_mono_106_stats, file = paste0("./outputs/Niche_evolution/Mono_sp/permMANOVA/permMANOVA_pPCA_mono_106.RData"))


### 5.3/ Plot null distri of pseudo-F stat ####

load(file = paste0("./outputs/Niche_evolution/Mono_sp/permMANOVA/permMANOVA_pPCA_mono_106.RData"))

pseudo_F_obs <- permMANOVA_pPCA_mono_106_stats$statistic
pseudo_F_null <- permMANOVA_pPCA_mono_106_stats$permutations[,1]

save(pseudo_F_obs, pseudo_F_null, file = paste0("./outputs/Niche_evolution/Mono_sp/permMANOVA/permMANOVA_pPCA_mono_106_plot_stuff.RData"))

pdf(file = "./graphs/Niche_evolution/Mono_sp/permMANOVA/perMANOVA_pseudo_F_null.pdf", height = 5.5, width = 5.5)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,5,3,1))

hist(x = c(pseudo_F_obs, pseudo_F_null), 
     breaks = 30,
     freq = TRUE, col = "gray",
     xlim = c(0, 6.5),
     ylim = c(0, 210),
     main = "Distribution of pseudo-F from perMANOVA\nunder null Hypothesis",
     # main = "",
     xlab = "pseudo-F",
     # axes = F,
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.2, lwd = 2)

# axis(side = 1, lwd = 2, cex.axis = 1.7, labels = c(0, 1.5, 3, 8, 20), at = seq(-1,1,0.5))
# axis(side = 2, lwd = 2, cex.axis = 1.7)

arrows(pseudo_F_obs + 0.08, 80, pseudo_F_obs + 0.08, 5, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = mean(c(pseudo_F_null, pseudo_F_obs)), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(c(pseudo_F_null, pseudo_F_obs), 0.95), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0(" Mean = ", format(round(mean(c(pseudo_F_null, pseudo_F_obs)),3), nsmall = 3)), 
                  paste0("CI 5% = ", round(quantile(c(pseudo_F_null, pseudo_F_obs), 0.95),3))), 
       x = "topright", inset = c(0, 0.15), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.4, bty = "n")
legend(legend = c(paste0("pseudo-F obs = ", round(pseudo_F_obs, 3)),
                  paste0("                     p = 0.001")),
       x = "bottomright", inset = c(0.00, 0.40),
       cex = 1.4, bty ="n", xjust = 1)

# legend(legend = as.expression(bquote(bold("B"))), 
#        x = "topright", inset = c(0.03, -0.03), xjust = 0.5,
#        cex = 1.9, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins)

dev.off()

##### 6/ Post-hoc test ####

# Test if differentce are significant between each modalities of the factor (i.e., each pair of mimicry rings)

# Load environmental pPCA data
pPC.env_units_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/pPC.env_units_mono_106.rds"))

library("RVAideMemoire")

?pairwise.perm.manova
# To look for the different statistical test available, see ?anova.mlm
# To see the different available correction of p-value for multiple testing, see ?p.adjust

Pairwise.permMANOVA_pPCA_mono_106 <- pairwise.perm.manova(resp = as.matrix(scale(pPC.env_units_mono_106[, c(1,2)])),
                                                          fact = reduced.list.unit_phyl_order_mono_106$Mimicry.model,
                                                          nperm = 999, progress = T,
                                                          p.method = "none", # No p-value correction because too many tests and p-value limited to 0.001
                                                          test = "Wilks")
Pairwise.permMANOVA_pPCA_mono_106
saveRDS(Pairwise.permMANOVA_pPCA_mono_106, file = paste0("./outputs/Niche_evolution/Mono_sp/permMANOVA/Pairwise.permMANOVA_pPCA_mono_106.rds"))

### 7/ Plot a heatmap of p-values for post-hoc tests #####

Pairwise.permMANOVA_pPCA_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/permMANOVA/Pairwise.permMANOVA_pPCA_mono_106.rds"))

# 7.1/ Extract p-value to build a heatmap

extract_mask <- lower.tri(Pairwise.permMANOVA_pPCA_mono_106$p.value) ; diag(extract_mask) <- T
Pairwise.permMANOVA_pPCA_mono_106_pvalues <- Pairwise.permMANOVA_pPCA_mono_106$p.value[extract_mask]

length(Pairwise.permMANOVA_pPCA_mono_106_pvalues) # 36 pairwise p-values
sum(Pairwise.permMANOVA_pPCA_mono_106_pvalues == 0.001) # 12 out of 36 (33.3%)
sum(Pairwise.permMANOVA_pPCA_mono_106_pvalues < 0.01) # 19 out of 36 (52.8%)
sum(Pairwise.permMANOVA_pPCA_mono_106_pvalues < 0.05) # 24 out of 36 (66.7%)

# Try to smooth the distribution for the heatmap
hist(Pairwise.permMANOVA_pPCA_mono_106_pvalues)
hist(log1p(Pairwise.permMANOVA_pPCA_mono_106_pvalues))
hist(sqrt(sqrt(Pairwise.permMANOVA_pPCA_mono_106_pvalues)))

# Try to adjust p-value for multiple testing with Holm correction. Most p-value become non-significant...
p.adjust(p = Pairwise.permMANOVA_pPCA_mono_106_pvalues, method = "holm")

# 7.2/ Build a full matrix of pairwise p-values

# Extract mimicry ring list
mimicry_ring_list <- big_rings
mimicry_ring_list <- c(attr(Pairwise.permMANOVA_pPCA_mono_106$p.value, "dimnames")[[2]], attr(Pairwise.permMANOVA_pPCA_mono_106$p.value, "dimnames")[[1]][nrow(Pairwise.permMANOVA_pPCA_mono_106$p.value)])

library(gdata)

Pseudo_Cor_mat_perMANOVA <- matrix(nrow = length(mimicry_ring_list), ncol = length(mimicry_ring_list), dimnames = list(mimicry_ring_list, mimicry_ring_list))
Pseudo_Cor_mat_perMANOVA[2:length(mimicry_ring_list), 1:(length(mimicry_ring_list)-1)] <- Pairwise.permMANOVA_pPCA_mono_106$p.value
upperTriangle(Pseudo_Cor_mat_perMANOVA) <- lowerTriangle(Pseudo_Cor_mat_perMANOVA, byrow=TRUE)
diag(Pseudo_Cor_mat_perMANOVA) <- 1

saveRDS(Pseudo_Cor_mat_perMANOVA, file = paste0("./outputs/Niche_evolution/Mono_sp/permMANOVA/Pseudo_Cor_mat_perMANOVA.rds"))

# 7.4/ Compute the dendrogram outside of heatmap() from climatic values and set branch width manually

library(rlang)
library(tidyverse)
library(dendextend)

ring_pPC_env_table_mono_106 <- as.data.frame(cbind(as.character(reduced.list.unit_phyl_order_mono_106$Mimicry.model), pPC.env_units_mono_106[, c(1,2)]))
names(ring_pPC_env_table_mono_106) <- c("Mimicry_ring", "PC1", "PC2")
ring_pPC_env_table_mono_106$PC1 <- as.numeric(ring_pPC_env_table_mono_106$PC1)
ring_pPC_env_table_mono_106$PC2 <- as.numeric(ring_pPC_env_table_mono_106$PC2)
ring_pPC_env_table_mono_106 <- ring_pPC_env_table_mono_106 %>% 
  group_by(Mimicry_ring) %>%
  summarise(PC1 = mean(PC1),
            PC2 = mean(PC2))

# Save
saveRDS(ring_pPC_env_table_mono_106, file = paste0("./outputs/Niche_evolution/Mono_sp/ring_pPC_env_table_mono_106.rds"))

# Convert to matrix and distance data
ring_pPC_env_mat_mono_106 <- as.matrix(ring_pPC_env_table_mono_106[, c("PC1","PC2")])

ring_pPC_clim_dist_mono_106 <- dist((scale(ring_pPC_env_mat_mono_106)))

# Set dendrogram for the heatmap
dd <- set(as.dendrogram(hclust(ring_pPC_clim_dist_mono_106, method = "complete")), "branches_lwd", 2)
plot(dd)

# Reorder leaves of the dendrogram manually
dd_reorder <- reorder(x = dd, wts = rowMeans(x = scale(ring_pPC_env_mat_mono_106)))
plot(dd_reorder)

saveRDS(dd_reorder, file = paste0("./outputs/Niche_evolution/Mono_sp/permMANOVA/dd_reorder.rds"))

# 7.5/ Plot the heatmap

# Set color palette
heatmap_pal <- rev(tmaptools::get_brewer_pal("RdYlBu", n = 100))[10:95]

pdf(file = "./graphs/Niche_evolution/Mono_sp/permMANOVA/Heatmap_PerMANOVA.pdf", height = 6, width = 5.5)

heatmap_posthoc_permMANOVA <- heatmap(# x = Pseudo_Cor_mat_perMANOVA,
                                      x = Pseudo_Cor_mat_perMANOVA^(1/5.8),
                                      symm = T, Rowv = rev(dd_reorder), Colv = rev(dd_reorder), revC = T,
                                      col = heatmap_pal,
                                      margins = c(6.5, 7.5))

dev.off()

# Reorder matrix of p-values to fit the order on the heat map based on the dendrogram
ordered.Pseudo_Cor_mat_perMANOVA <- Pseudo_Cor_mat_perMANOVA[heatmap_posthoc_permMANOVA$rowInd, heatmap_posthoc_permMANOVA$colInd]

# Plot the scale

library(raster)

pdf(file = "./graphs/Niche_evolution/Mono_sp/permMANOVA/Heatmap_PerMANOVA_scale.pdf", height = 6, width = 8)

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


##### 8/ Plot OMUs and rings in pPCA's climatic space #####

### 8.1/ Load stuff 

phylo.Ithomiini_mono_106 <- readRDS(file = paste0("./input_data/Phylogenies/phylo.Ithomiini_mono_106.rds"))
reduced.list.unit_phyl_order_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/reduced.list.unit_phyl_order_mono_106.rds"))
pPC.env_units_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/pPC.env_units_mono_106.rds"))
load(file = paste0("./outputs/Niche_evolution/permMANOVA/pPCA_correl.RData"))

### 8.2/ Select only a reduced number of rings

Mimicry_counts_table <- table(reduced.list.unit_phyl_order_mono_106$Mimicry.model)[order(table(reduced.list.unit_phyl_order_mono_106$Mimicry.model), decreasing = T)] ; Mimicry_counts_table

# Extract only the 3 or 6 main circles in colors

top3_ring_list <- names(Mimicry_counts_table)[1:3]
top5_ring_list <- names(Mimicry_counts_table)[1:5]

# AGNOSIA 25 # LERIDA 21 # EURIMEDIA 17 # CONFUSA 10 # PANTHYALE 9 # HERMIAS 8

top3_indices <- (reduced.list.unit_phyl_order_mono_106$Mimicry.model %in% top3_ring_list)
top5_indices <- (reduced.list.unit_phyl_order_mono_106$Mimicry.model %in% top5_ring_list)

saveRDS(top3_indices, file = paste0("./outputs/Niche_evolution/Mono_sp/permMANOVA/top3_indices.rds"))
saveRDS(top5_indices, file = paste0("./outputs/Niche_evolution/Mono_sp/permMANOVA/top5_indices.rds"))

top3_indices <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/permMANOVA/top3_indices.rds"))
top5_indices <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/permMANOVA/top5_indices.rds"))

### 8.3/ Plot for permMANOVA: top 3 rings

pdf(file = "./graphs/Niche_evolution/Mono_sp/permMANOVA/pPCA_Top3.pdf", height = 7.5, width = 8)

original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))

var_perc_pPCaxes <- diag(var(pPC.env_units_mono_106))/sum(diag(var(pPC.env_units_mono_106)))*100

plot(pPC.env_units_mono_106, col = "grey90", main = "", 
     xlab = paste0("pPC1 (", round(var_perc_pPCaxes[1],1) ," %)"), ylab = paste0("pPC2 (", round(var_perc_pPCaxes[2],1) ," %)"),
     pch = 16, type = "p", axes = F, ylim = c(-6, 6), xlim = c(-6, 6),
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.3)
axis(side = 1, lwd = 2, cex.axis = 1.7)
axis(side = 2, lwd = 2, cex.axis = 1.7)

points(pPC.env_units_mono_106[top3_indices, ], pch = 16,
       col = c("deepskyblue","limegreen","orange")[factor(reduced.list.unit_phyl_order_mono_106$Mimicry.model[top3_indices])])

legend(x = "bottomright", pch = 16, cex = 1.3, legend = top3_ring_list[1:3], 
       bty = "n", inset = c(0.26, 0.02), # horiz = T, 
       col = c("deepskyblue","limegreen","orange")[1:3])

car::dataEllipse(x = pPC.env_units_mono_106[top3_indices, 1],
                 y = pPC.env_units_mono_106[top3_indices, 2],
                 groups = factor(reduced.list.unit_phyl_order_mono_106$Mimicry.model[top3_indices]),
                 group.labels = "", ellipse.label = "",
                 levels = 0.80, center.pch = F,
                 center.cex = 1.5, draw=TRUE, plot.points = F, add = T,
                 col = c("deepskyblue","limegreen","orange"),
                 lwd = 2, fill = FALSE, fill.alpha = 0.3, grid = F)

arrows(x0 = rep(0, nrow(pPCA_correl)), y0 = rep(0, nrow(pPCA_correl)), x1 = pPCA_correl[,1]*4.5, y1 =  pPCA_correl[,2]*4.5, length = 0.2, code = 2, lwd = 1.5)

text(pPCA_correl[,1:2]*4.7 + cbind(c(0, -0.9, 0.0, -0.6), c(-0.3, -0.3, -0.2, -0.4)),
     labels = row.names(pPCA_correl),
     cex = 1.2, font = 2, adj = c(0, 0))

# legend(legend = as.expression(bquote(bold("D"))),
#        x = "topright", inset = c(0.05, 0.00), xjust = 0.5,
#        cex = 1.9, bty ="n")

par(mar = original_int_margins)

dev.off()


### 8.4/ Plot for permMANOVA: top 5 rings

pdf(file = "./graphs/Niche_evolution/Mono_sp/permMANOVA/pPCA_Top5.pdf", height = 7.5, width = 8)

original_int_margins <- par()$mar
par(mar = c(5.1,5,4.1,2.1))

var_perc_pPCaxes <- diag(var(pPC.env_units_mono_106))/sum(diag(var(pPC.env_units_mono_106)))*100

plot(pPC.env_units_mono_106, col = "grey90", main = "", 
     xlab = paste0("pPC1 (", round(var_perc_pPCaxes[1],1) ," %)"), ylab = paste0("pPC2 (", round(var_perc_pPCaxes[2],1) ," %)"),
     pch = 16, type = "p", axes = F, ylim = c(-6, 6), xlim = c(-6, 6),
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.3)
axis(side = 1, lwd = 2, cex.axis = 1.7)
axis(side = 2, lwd = 2, cex.axis = 1.7)

points(pPC.env_units_mono_106[top5_indices, ], pch = 16,
       col = c("deepskyblue","limegreen","orange", "red", "black")[factor(reduced.list.unit_phyl_order_mono_106$Mimicry.model[top5_indices])])

legend(x = "bottomright", pch = 16, cex = 1.3, legend = top5_ring_list[1:3], 
       bty = "n", inset = c(0.02, 0.02), # horiz = T, 
       col = c("deepskyblue","limegreen","orange", "red", "black")[1:3])

legend(x = "bottomright", pch = 16, cex = 1.3, legend = top5_ring_list[4:5], 
       bty = "n", inset = c(0.36, 0.02), # horiz = T, 
       col = c("deepskyblue","limegreen","orange", "red", "black")[4:5])

car::dataEllipse(x = pPC.env_units_mono_106[top5_indices, 1],
                 y = pPC.env_units_mono_106[top5_indices, 2],
                 groups = factor(reduced.list.unit_phyl_order_mono_106$Mimicry.model[top5_indices]),
                 group.labels = "", ellipse.label = "",
                 levels = 0.80, center.pch = F,
                 center.cex = 1.5, draw=TRUE, plot.points = F, add = T,
                 col = c("deepskyblue","limegreen","orange", "red", "black"),
                 lwd = 2, fill = FALSE, fill.alpha = 0.3, grid = F)

arrows(x0 = rep(0, nrow(pPCA_correl)), y0 = rep(0, nrow(pPCA_correl)), x1 = pPCA_correl[,1]*4.5, y1 =  pPCA_correl[,2]*4.5, length = 0.2, code = 2, lwd = 1.5)

text(pPCA_correl[,1:2]*4.7 + cbind(c(0, -0.9, 0.0, -0.6), c(-0.3, -0.3, -0.2, -0.4)),
     labels = row.names(pPCA_correl),
     cex = 1.2, font = 2, adj = c(0, 0))

# legend(legend = as.expression(bquote(bold("D"))),
#        x = "topright", inset = c(0.05, 0.00), xjust = 0.5,
#        cex = 1.9, bty ="n")

par(mar = original_int_margins)

dev.off()


##### 9/ Plot both together #####

# For Pseudo-F histogram
load(file = paste0("./outputs/Niche_evolution/Mono_sp/permMANOVA/permMANOVA_pPCA_mono_106_plot_stuff.RData"))
round(pseudo_F_obs, 2) # pseudo-F obs = 6.42

# For pPCA-space
reduced.list.unit_phyl_order_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/reduced.list.unit_phyl_order_mono_106.rds"))
pPC.env_units_mono_106 <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/pPC.env_units_mono_106.rds"))
var_perc_pPCaxes <- diag(var(pPC.env_units_mono_106))/sum(diag(var(pPC.env_units_mono_106)))*100
load(file = paste0("./outputs/Niche_evolution/permMANOVA/pPCA_correl.RData"))
top5_indices <- readRDS(file = paste0("./outputs/Niche_evolution/Mono_sp/permMANOVA/top5_indices.rds"))



pdf(file = "./graphs/Niche_evolution/Mono_sp/permMANOVA/permMANOVA_both_plots.pdf", height = 7, width = 14)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,5,1.5,1), mfrow = c(1,2))

### Panel A = PCA-space

plot(pPC.env_units_mono_106, col = "grey90", main = "", 
     xlab = paste0("pPC1 (", round(var_perc_pPCaxes[1],1) ," %)"), ylab = paste0("pPC2 (", round(var_perc_pPCaxes[2],1) ," %)"),
     pch = 16, type = "p", axes = F, ylim = c(-6, 6), xlim = c(-6, 6),
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.3)
axis(side = 1, lwd = 2, cex.axis = 1.7)
axis(side = 2, lwd = 2, cex.axis = 1.7)

points(pPC.env_units_mono_106[top5_indices, ], pch = 16,
       col = c("deepskyblue","limegreen","orange", "red", "black")[factor(reduced.list.unit_phyl_order_mono_106$Mimicry.model[top5_indices])])

legend(x = "bottomright", pch = 16, cex = 1.3, legend = top5_ring_list[1:3], 
       bty = "n", inset = c(0.02, 0.02), # horiz = T, 
       col = c("deepskyblue","limegreen","orange", "red", "black")[1:3])

legend(x = "bottomright", pch = 16, cex = 1.3, legend = top5_ring_list[4:5], 
       bty = "n", inset = c(0.36, 0.02), # horiz = T, 
       col = c("deepskyblue","limegreen","orange", "red", "black")[4:5])

car::dataEllipse(x = pPC.env_units_mono_106[top5_indices, 1],
                 y = pPC.env_units_mono_106[top5_indices, 2],
                 groups = factor(reduced.list.unit_phyl_order_mono_106$Mimicry.model[top5_indices]),
                 group.labels = "", ellipse.label = "",
                 levels = 0.80, center.pch = F,
                 center.cex = 1.5, draw=TRUE, plot.points = F, add = T,
                 col = c("deepskyblue","limegreen","orange", "red", "black"),
                 lwd = 2, fill = FALSE, fill.alpha = 0.3, grid = F)

arrows(x0 = rep(0, nrow(pPCA_correl)), y0 = rep(0, nrow(pPCA_correl)), x1 = pPCA_correl[,1]*4.5, y1 =  pPCA_correl[,2]*4.5, length = 0.2, code = 2, lwd = 1.5)

text(pPCA_correl[,1:2]*4.7 + cbind(c(0, -0.9, 0.0, -0.6), c(-0.3, -0.3, -0.2, -0.4)),
     labels = row.names(pPCA_correl),
     cex = 1.2, font = 2, adj = c(0, 0))

legend(legend = as.expression(bquote(bold("A"))), 
       x = "topright", inset = c(0.05, -0.03), xjust = 0.5,
       cex = 1.9, bty ="n")

### Panel B = Pseudo-F histogram

hist(x = c(pseudo_F_obs, pseudo_F_null), 
     breaks = 30,
     freq = TRUE, col = "gray",
     xlim = c(0, 6.5),
     ylim = c(0, 210),
     # main = "Distribution of pseudo-F from perMANOVA\nunder null Hypothesis",
     main = "",
     xlab = "pseudo-F",
     # axes = F,
     cex.axis = 1.7, cex.lab = 1.8, cex.main = 1.2, lwd = 2)

# axis(side = 1, lwd = 2, cex.axis = 1.7, labels = c(0, 1.5, 3, 8, 20), at = seq(-1,1,0.5))
# axis(side = 2, lwd = 2, cex.axis = 1.7)

arrows(pseudo_F_obs + 0.08, 80, pseudo_F_obs + 0.08, 5, length = 0.1, lwd = 2)  # Draw arrow above mean BC obs
abline(v = mean(c(pseudo_F_null, pseudo_F_obs)), lwd = 2, lty = 2) # Add vertical line for mean value
abline(v = quantile(c(pseudo_F_null, pseudo_F_obs), 0.95), lwd = 2, lty = 2, col = "red") # Add vertical line for 95% value

legend(legend = c(paste0(" Mean = ", format(round(mean(c(pseudo_F_null, pseudo_F_obs)),3), nsmall = 3)), 
                  paste0("CI 5% = ", round(quantile(c(pseudo_F_null, pseudo_F_obs), 0.95),3))), 
       x = "topright", inset = c(0, 0.15), 
       lty = 2 , lwd = 2, col = c("black", "red"), cex = 1.4, bty = "n")
legend(legend = c(paste0("pseudo-F obs = ", round(pseudo_F_obs, 3)),
                  paste0("                     p = 0.001")),
       x = "bottomright", inset = c(0.00, 0.40),
       cex = 1.4, bty ="n", xjust = 1)

legend(legend = as.expression(bquote(bold("B"))), 
       x = "topright", inset = c(0.03, -0.03), xjust = 0.5,
       cex = 1.9, bty ="n")

par(oma = original_ext_margins, mar = original_int_margins, mfrow = c(1,1))

dev.off()
