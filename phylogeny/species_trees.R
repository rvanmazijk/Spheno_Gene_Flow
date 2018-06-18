setwd("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/all_genomic_data/species_tree/")
t1 = read.tree("best_trees_tol1e-05_collapse30.tre")
t2 = read.tree("best_trees_tol1e-05_collapse50.tre")
t3 = read.tree("best_trees_tol1e-05_collapse80.tre")

d = read.csv('~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv', stringsAsFactors=F)

trees = list(t1, t2, t3)
for (i in 1:length(trees)) {
	trees[[i]] = unroot(trees[[i]])
	trees[[i]] = root(trees[[i]], c("Plestiodon_laticeps", "NA_ABTC3943_Er_grac"))
	trees[[i]] = drop.tip(trees[[i]], c("Plestiodon_laticeps", "NA_ABTC3943_Er_grac"))

	trees[[i]]$tip.label = d[match(trees[[i]]$tip.label, d$genomic_sample), "SPECIES"]
}
pdf("species_trees.pdf", width=14, height=8)
par(mfrow=c(1,3), mai=c(0.1, 0.1, 0.1, 0.1))
val = c(30, 50, 80)
for (i in 1:length(trees)) {
	plot(trees[[i]], cex=0.5, main=paste("collapse", val[i]))
	nodes = as.numeric(trees[[i]]$node.label)
	
	for (j in 1:length(nodes)) {
		if (!is.na(nodes[j])) {
			if (nodes[j] < 0.95) {
				nodelabels(node= j + length(trees[[i]]$tip.label), bg="red", pch=21, cex=1)
				}
			}
		}
	}
dev.off()


t1 = read.tree("RAxML_result.collapse30")
t2 = read.tree("RAxML_result.collapse50")
t3 = read.tree("RAxML_result.collapse80")

trees = list(t1, t2, t3)
for (i in 1:length(trees)) {
	trees[[i]] = unroot(trees[[i]])
	trees[[i]] = root(trees[[i]], c("Plestiodon_laticeps", "NA_ABTC3943_Er_grac"))
	trees[[i]] = drop.tip(trees[[i]], c("Plestiodon_laticeps", "NA_ABTC3943_Er_grac"))

	trees[[i]]$tip.label = d[match(trees[[i]]$tip.label, d$genomic_sample), "SPECIES"]
}
pdf("species_trees_BL_pf.pdf", width=14, height=8)
par(mfrow=c(1,3), mai=c(0.1, 0.1, 0.1, 0.1))
val = c(30, 50, 80)
for (i in 1:length(trees)) {
	plot(trees[[i]], cex=0.5, main=paste("collapse", val[i]))
	nodes = as.numeric(trees[[i]]$node.label)
	
	}
dev.off()

pdf("species_trees_BL_pf_chronopl.pdf", width=14, height=8)
par(mfrow=c(1,3), mai=c(0.1, 0.1, 0.1, 0.1))
val = c(30, 50, 80)
for (i in 1:length(trees)) {
	plot(chronopl(trees[[i]], 0.0001), cex=0.5, main=paste("collapse", val[i]))
	nodes = as.numeric(trees[[i]]$node.label)
	
	}
dev.off()