library(ape)
library(phangorn)
library(treespace)

d = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv", stringsAsFactors = F)

t1 = read.tree("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/skinktree_216.tre")
t2 = read.nexus("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/RAxML_unconstrained/mtcds_pos1.tre")
t3 = read.nexus("~/Desktop/starbeast/run2/species.tre")
t3$tip.label = gsub("^[A-Z|a-z]+_[A-Z|a-z]+_?\\d?_?", "", t3$tip.label)
t4 = read.nexus("~/Desktop/mtcds_pos1.sub.tre")

outs = d[d$outgroup == TRUE, "sanger_sample"]
outs = outs[!is.na(outs)]

trees = list(t1, t2, t3, t4)
names(trees) = c("dan", "raxml", "starbeast", "beast")
tips = gsub("\\d$", "", t1$tip.label)
dan = data.frame(t1$tip.label, d[match(tips, tolower(d$SPECIES)), "SPECIES"], stringsAsFactors = F)
names(dan) = c("dan", "me")
for (i in 1:nrow(dan)) {
  if (is.na(dan[i, "me"])) {
    sp = dan[i, "dan"]
    keep = grep(sp, d$OTHER_NAMES, ignore.case = T)
    if (length(keep) > 0) {
      dan[i, "me"]  = d[keep[1], "SPECIES"]
      }
  }
}
dan = dan[complete.cases(dan), ]
dan = dan[!duplicated(dan$me), ]
keep = d[d$SPECIES_KEEP == TRUE, "sanger_sample"]
keep = keep[!is.na(keep)]
t1 = drop.tip(trees[[2]], setdiff(trees[[2]]$tip.label, keep))
t1$tip.label = d[match(t1$tip.label, d$sanger_sample), "SPECIES"]
keep = t1$tip.label[t1$tip.label %in% dan$me]
t1 = drop.tip(t1, setdiff(t1$tip.label, keep))
allkeep1 = t1$tip.label
allkeep2 = dan[dan$me %in% allkeep, "dan"]
for (i in 1:length(trees)) {
  trees[[i]] = drop.tip(trees[[i]], outs)
  if (i > 1) {
    keep = d[d$SPECIES_KEEP == TRUE, "sanger_sample"]
    keep = keep[!is.na(keep)]
    t1 = drop.tip(trees[[i]], setdiff(trees[[i]]$tip.label, keep))
    t1$tip.label = d[match(t1$tip.label, d$sanger_sample), "SPECIES"]
    t1 = drop.tip(t1, setdiff(t1$tip.label, allkeep1))
  } else {
    t1 = trees[[i]]
    t1 = drop.tip(t1, setdiff(t1$tip.label, allkeep2))
    t1$tip.label = dan[match(t1$tip.label, dan$dan), "me"]
  }
  trees[[i]] = t1
}

class(trees) = "multiPhylo"
treespace(trees, nf=2)$D
plotTreeDiff(trees$dan, trees$beast, use.edge.length=T)
