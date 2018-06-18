t = read.tree("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/all_genomic_data/species_tree/best_trees_tol1e-05_collapse50.tre")
t = unroot(t)
t = root(t, c("Plestiodon_laticeps", "NA_ABTC3943_Er_grac"))
t$node.label = as.numeric(t$node.label)
t$edge.length[!complete.cases(t$edge.length)] = 1

badnodes <- which(as.numeric(t$node.label) < 0.95) + length(t$tip.label)
t$edge.length[which(t$edge[,2] %in% badnodes)] = 0
tm <- di2multi(t)
tm$edge.length = NULL
tm$node.label = NULL

tm = drop.tip(tm, "NA_ABTC3943_Er_grac")
d = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv", stringsAsFactors=F)

write.tree(tm, "~/Desktop/multifurcating_tree.tre")