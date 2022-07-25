
##### Build a phylogeny including the OMUs on null distance terminal branches ####

###################################
#      Author: Maël Doré          #
#  Contact: mael.dore@gmail.com   #
###################################

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
