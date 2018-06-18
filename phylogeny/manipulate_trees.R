library(ggtree)

# lineage data
d = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv", stringsAsFactors=F)
s = d[complete.cases(d$sanger_sample),]

# read in tree data
td = read.beast("/Users/Sonal/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/beast_ucln_31July17/mtcds_pos1.sub.tre")
td = fortify(td)

# read in tree
t = read.nexus("/Users/Sonal/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/beast_ucln_31July17/mtcds_pos1.sub.tre")
t$node.label = td[td$node > Ntip(t), "posterior"]

# drop outgroups
outs = d[d$outgroup == TRUE, "sanger_sample"]
outs = outs[complete.cases(outs)]
t = drop.tip(t, outs)

# drop tips either matched to 2 otus or 2 species
keep = s[s$SPECIES_KEEP == TRUE, "sanger_sample"]
t1 = drop.tip(t, setdiff(t$tip.label, keep))
t1$tip.label = d[match(t1$tip.label, d$sanger_sample), "SPECIES"]
write.tree(t1, "/Users/Sonal/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/beast_ucln_31July17//species.tre")

# drop tips either matched to 2 otus or 2 species
keep = s[s$OTU_KEEP == TRUE, "sanger_sample"]
t1 = drop.tip(t, setdiff(t$tip.label, keep))
t1$tip.label = d[match(t1$tip.label, d$sanger_sample), "OTU"]
write.tree(t1, "/Users/Sonal/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/beast_ucln_31July17//otu.tre")

# get sampling probs
length(unique(d[complete.cases(d$sanger_sample), "OTU"])) / length(unique(d$OTU))
length(unique(d[complete.cases(d$sanger_sample), "SPECIES"])) / length(unique(d$SPECIES))

trees = read.nexus("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/beast_ucln_31July17/mtcds_pos1.sub.trees")
trees1 = trees[sample(1:length(trees),  100, replace=F)]
trees2 = vector('list', 100)
for (i in 1:length(trees1)) {
  t = trees1[[i]]
  t = drop.tip(t, outs)
  keep = s[s$OTU_KEEP == TRUE, "sanger_sample"]
  t1 = drop.tip(t, setdiff(t$tip.label, keep))
  t1$tip.label = d[match(t1$tip.label, d$sanger_sample), "OTU"]
  trees2[[i]] = t1
}
class(trees2) <- "multiPhylo"
write.tree(trees2, "~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/beast_ucln_31July17/otu_posterior.trees")
