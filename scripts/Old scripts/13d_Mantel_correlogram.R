### Preparation ####

# Effacer l'environnement
rm(list = ls())

# Load packages
library(ape)

# Set Wd
internal.wd <- "F:/Documents/Etudes/Fac/Cours/TROPIMUNDO/Stage_S4/Projet_Papillons/R_codes"

## Phylogeny

# Load phylogeny
load(file = paste0(internal.wd, "/Phylogenies/Final_phylogeny.RData")) # N = 339
# Compute Phylogenetic/patristic distance matrix
phylo.dist.mat <- cophenetic.phylo(x = phylo.Ithomiini)

# Study distances distribution to define classes of phylogenetic distances
hist(phylo.dist.mat)
max(phylo.dist.mat) # Max = 52 My
# 10 classes from 0 to 50 My every 5 My


## Climate

# Load the env. variables reprojected with pPCA under two axis
load(file = paste0(internal.wd, "/Niche_evolution/PC.Revell.env.RData"))
# Compute the climatic distance matrix
Revell_Pairwise_climdist <- dist(x = PC.Revell.env, method = "euclidian")

## Mimicry patterns

# Load the Jaccard index for mimicry pattern distances between species
load(file = paste0(internal.wd,"/Niche_evolution/Jaccard_index.RData"))

mimicry.similarity.weights

x11()
plot(as.dist(phylo.dist.mat), Revell_Pairwise_climdist)

mod <- lm(Revell_Pairwise_climdist ~ as.dist(phylo.dist.mat))

abline(coef(mod), lwd = 3, col = "red")

plot(as.dist(phylo.dist.mat), mimicry.similarity.weights)
mod <- lm(mimicry.similarity.weights ~ as.dist(phylo.dist.mat))

abline(coef(mod), lwd = 3, col = "red")


### Multivariate Mantel's correlogram ####

library(vegan)

?mantel.correlog()

mantel.correlog(D.eco, # Matrice des distances écologiques (Jaccard de motifs, Distance niches climatiques)
                D.geo, # Matrice des distances géographiques (ou phylogénétiques)
                XY, # Coordonnées cartésiennes pour calcul automatique de distances géographiques euclidiennes
                n.class, # Nombre de classes de distances. Struges equation par défaut
                break.pts, # Vecteur pour définir manuellement les limites de classes de distances
                cutoff, # ???
                r.type = c("pearson","spearman","kendall"), # Choix de la stats de test pour les Mantel test
                nperm, # Nb de permutations pour les tests. 999 par défaut
                mult = c("holm","hochberg","bonferroni", ...), # Methode d'ajustement des p-value pour tests multiples. Voir ?p.adjust. Holm par défaut
                progressive = T/F, # Pour une correction progressive des p-value. Test of the first distance class: no correction; second distance class: correct for 2 simultaneous tests; distance class k: correct for k simultaneous tests. 
                # Valable surtout si on veut mettre en évidence l'autocorrelation dans les premières classes. Sinon, si on veut comparer toutes les classes de distances sur un pied d'égalité, ne pas utiliser.
                alpha # Niveau de significativité choisi pour discriminer les symboles sur plot
)

### Dclim ~ Dphylo

Mantel_correlog_Dclim <- mantel.correlog(D.eco = Revell_Pairwise_climdist, # Matrice des distances écologiques (Jaccard de motifs, Distance niches climatiques)
                D.geo = phylo.dist.mat, # Matrice des distances géographiques (ou phylogénétiques)
                n.class = 0, # Nombre de classes de distances. Struges equation par défaut
                break.pts = # Vecteur pour définir manuellement les limites de classes de distances
                # seq(from=0, to = 55, by = 5), 
                 seq(from=0, to = 30, by = 2),
                # NULL,
                cutoff = F, # Si T : limite les classes de distances dès que tous les points ont au moins été inclus une fois (i.e., distance maximale parmi les distances minimales des sites/espèces).
                            # En contexte de Dphylo : 2xlongueur de la branche terminale de l'espèce la plus "originiale"
                r.type = "pearson", # Choix de la stats de test pour les Mantel test
                nperm = 999, # Nb de permutations pour les tests. 999 par défaut
                mult = "holm", # Methode d'ajustement des p-value pour tests multiples. Voir ?p.adjust. Holm par défaut
                progressive = F, # Pour une correction progressive des p-value. Test of the first distance class: no correction; second distance class: correct for 2 simultaneous tests; distance class k: correct for k simultaneous tests. 
                # Valable surtout si on veut mettre en évidence l'autocorrelation dans les premières classes. Sinon, si on veut comparer toutes les classes de distances sur un pied d'égalité, ne pas utiliser.
                )
# Warning if the last break.pts is not higher than max.dist (i.e., some pairs are not evaluated, no included inthe last distance class)

# Mantel_correlog_Dclim_Struges <- Mantel_correlog_Dclim
# Mantel_correlog_Dclim_5My <- Mantel_correlog_Dclim
# Mantel_correlog_Dclim_2My <- Mantel_correlog_Dclim
# Mantel_correlog_Dclim_2My_restr <- Mantel_correlog_Dclim
# save(Mantel_correlog_Dclim_Struges, Mantel_correlog_Dclim_5My, Mantel_correlog_Dclim_2My, Mantel_correlog_Dclim_2My_restr, file = paste0(internal.wd,"/Niche_evolution/Mantel_correlogs_Dclim.RData"))

Mantel_correlog_Dclim
str(Mantel_correlog_Dclim)

plot(Mantel_correlog_Dclim, alpha = 0.05)  # Niveau de significativité choisi pour discriminer les symboles sur plot

# Export

# Chose
Mantel_correlog <- Mantel_correlog_Dclim_2My ; name <- "Mantel_correlog_Dclim_2My"
Mantel_correlog <- Mantel_correlog_Dclim_5My ; name <- "Mantel_correlog_Dclim_5My"

pdf(file = paste0(internal.wd,"/Niche_evolution/Plots/Mantel_correlogs/",name,".pdf"))
plot(Mantel_correlog, alpha = 0.05) ; title(main = paste0("Dclim ~ Dphyl\n step = ",substr(x = name, start = nchar(name) - 2, stop = nchar(name))))
dev.off()



### Dmim ~ Dphylo

Mantel_correlog_Dmim <- mantel.correlog(D.eco = mimicry.similarity.weights, # Matrice des distances écologiques (Jaccard de motifs, Distance niches climatiques)
                                         D.geo = phylo.dist.mat, # Matrice des distances géographiques (ou phylogénétiques)
                                         n.class = 0, # Nombre de classes de distances. Struges equation par défaut
                                         break.pts = # Vecteur pour définir manuellement les limites de classes de distances
                                           seq(from=0, to = 55, by = 5), 
                                         #  seq(from=0, to = 30, by = 2),
                                         # NULL,
                                         cutoff = F, # Si T : limite les classes de distances dès que tous les points ont au moins été inclus une fois (i.e., distance maximale parmi les distances minimales des sites/espèces).
                                         # En contexte de Dphylo : 2xlongueur de la branche terminale de l'espèce la plus "originiale"
                                         r.type = "spearman", # Choix de la stats de test pour les Mantel test
                                         nperm = 999, # Nb de permutations pour les tests. 999 par défaut
                                         mult = "holm", # Methode d'ajustement des p-value pour tests multiples. Voir ?p.adjust. Holm par défaut
                                         progressive = F, # Pour une correction progressive des p-value. Test of the first distance class: no correction; second distance class: correct for 2 simultaneous tests; distance class k: correct for k simultaneous tests. 
                                         # Valable surtout si on veut mettre en évidence l'autocorrelation dans les premières classes. Sinon, si on veut comparer toutes les classes de distances sur un pied d'égalité, ne pas utiliser.
)
# Warning if the last break.pts is not higher than max.dist (i.e., some pairs are not evaluated, no included inthe last distance class)

# Mantel_correlog_Dmim_Struges <- Mantel_correlog_Dmim
# Mantel_correlog_Dmim_5My <- Mantel_correlog_Dmim
# Mantel_correlog_Dmim_5My_pearson <- Mantel_correlog_Dmim
# Mantel_correlog_Dmim_2My <- Mantel_correlog_Dmim
# save(Mantel_correlog_Dmim_Struges, Mantel_correlog_Dmim_5My, Mantel_correlog_Dmim_5My_pearson, Mantel_correlog_Dmim_2My, file = paste0(internal.wd,"/Niche_evolution/Mantel_correlogs_Dmim.RData"))

Mantel_correlog_Dmim
str(Mantel_correlog_Dmim)

plot(Mantel_correlog_Dmim, alpha = 0.05) # Niveau de significativité choisi pour discriminer les symboles sur plot


# Export

# Chose
Mantel_correlog <- Mantel_correlog_Dmim_5My ; name <- "Mantel_correlog_Dmim_5My"

pdf(file = paste0(internal.wd,"/Niche_evolution/Plots/Mantel_correlogs/",name,".pdf"))
plot(Mantel_correlog, alpha = 0.05) ; title(main = paste0("Dmim ~ Dphyl\n step = ",substr(x = name, start = nchar(name) - 2, stop = nchar(name))))
dev.off()
