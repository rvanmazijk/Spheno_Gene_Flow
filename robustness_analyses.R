library(ape)
library(phytools)
library(phangorn)
library(phylolm)
library(nlme)

##########################
# functions
##########################

add_tips <- function(tree, tips) {
  # for depicting at least, add tip to tree next to a species in complex
  # long term, figure out a way to get these OTUs in tree
  not_in_tree = tips[!complete.cases(match(tips, tree$tip.label))]
  if (length(not_in_tree) > 0) {
    for (i in 1:length(not_in_tree)) {
      sister = grep(gsub("_\\d", "", not_in_tree[i]), tree$tip.label)
      sister = tree$tip.label[sister]
      if (length(sister) > 0) {
        sister = sister[1]
        tree <- bind.tip(tree, not_in_tree[i], where=which(tree$tip.label== sister),
                       position = 0.5 * tree$edge.length[which(tree$edge[,2]==
                       which(tree$tip.label==sister))])
      }
    }  
  }
  return(tree)
}

prune_tree <- function(tree, tips) {
  keep = intersect(tree$tip.label, tips)
  tree = drop.tip(tree, setdiff(tree$tip.label, keep))
  return(tree)
}

run_phylolm <- function(tree, data, iv, dv) {
  data2 = data[, c(iv, dv)]
  data2[, dv] = log(data2[, dv])
  data2 = data2[complete.cases(data2), ]
  keep = intersect(rownames(data2), tree$tip.label)
  tree2 = drop.tip(tree, setdiff(tree$tip.label, keep))
  data2 = data2[tree2$tip.label, ]
  m = phylolm(data2[, dv] ~ data2[, iv], data=data2, tree2, model="lambda", lower.bound=0)
  pval = summary(m)$coefficients[2, 4]
  sig = FALSE
  if (pval <= 0.05) {
    sig = TRUE
  }
  return(list(sig, pval))
}

run_diffdiv <- function(tree, data, rates) {
  data2 = data
  data2$rate = rates[rownames(data2)]
  data2$log_slope = log(data2$slope)
  data2 = data2[, c("log_slope", "rate")]
  data2 = data2[complete.cases(data2), ]
  keep = intersect(rownames(data2), tree$tip.label)
  tree2 = drop.tip(tree, setdiff(tree$tip.label, keep))
  data2 = data2[tree2$tip.label, ]
  m = gls(rate ~ log_slope, correlation=corBrownian(phy=tree2), data=data2, method="ML")
  pval = summary(m)$tTable[2, 4]
  sig = 0.00
  if (pval <= 0.05) {
    sig = 1.00
  }
  return(list(sig, pval))
} 

source("~/scripts/gene_flow/jetz_div_rates.R")

################################
# prepare data
################################

# species data
d = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv", stringsAsFactors=F)

# get in pi data
pifiles = list.files("~/Dropbox/Sphenomorphine_Gene_Flow/data/diversity/", full.names = T)
pi = do.call(rbind, lapply(pifiles, read.csv, stringsAsFactors=F))
pi = pi[pi$type == 'IND',]
pi = pi[pi$pi_denom > 50000,]
pi = aggregate(pi$pi, by=list(pi$lineage), mean, na.rm=T)
names(pi) = c("cluster", "mean_pi")

# get in morphology data
morph = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/morphology/species_OTU_means_PCAs.csv", stringsAsFactors=F)

# get in range data
r = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/envirospatial/range_data.csv", stringsAsFactors=F)

# trees
stree = read.tree("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/beast_ucln_31July17/species.tre")
otree = read.tree("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/beast_ucln_31July17/otu.tre")
srates = jetzDivRates(stree)
orates = jetzDivRates(otree)

#######################################
# NOW START
#######################################

cols = c("analysis_type", "ntips", "lambda", 
         "lambda_sig", "pi", "pi_sig", "limb", "limb_sig", "divdiff",
         "divdiff_sig")
res = data.frame(matrix(NA, nrow=13, ncol=length(cols)))
names(res) = cols

# 1 bootstrapping loci in each F_ST estimate and recalculating slopes
diff = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.csv", stringsAsFactors = F)
diff1 = diff[diff$distance_type == 'geo_dist' & diff$genetic_type == 'fst', ]
diff1$mean_pi = pi[match(diff1$cluster, pi$cluster), "mean_pi"]
diff1 = cbind(diff1, morph[match(diff1$OTU, morph$OTU), 3:ncol(morph)])
diff1 = cbind(diff1, r[match(diff1$OTU, r$OTU), 2:ncol(r)])

nrep = 100
tmpres = vector('list', length(cols))
names(tmpres) = cols
for (i in 1:length(tmpres)) {
  tmpres[[i]] = rep(NA, nrep)
}

slopes = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.bootloci.csv", stringsAsFactors = F)
for (i in 1:nrep) {
  slopes2 = slopes[slopes$rep == i, ]
  diff2 = diff1
  diff2$slope = slopes2[match(diff1$OTU, slopes2$OTU), "slope"]
  tree2 = add_tips(otree, diff2$OTU)
  tree3 = prune_tree(tree2, diff2$OTU)
  tmpres[['ntips']][i] = Ntip(tree3)
  rownames(diff2) = diff2$OTU
  lambda = phylosig(tree3, diff2[tree3$tip.label, "slope"], method="lambda", test=T)
  tmpres[['lambda']][i] = lambda$lambda
  tmpres[['lambda_sig']][i]= lambda$P
  a = run_phylolm(tree3, diff2, 'mean_pi', 'slope')
  tmpres[['pi']][i] = a[[1]]
  tmpres[['pi_sig']][i] = a[[2]]
  b = run_phylolm(tree3, diff2, 'PC2', 'slope')
  tmpres[['limb']][i] = b[[1]]
  tmpres[['limb_sig']][i] = b[[2]]
  c = run_diffdiv(tree3, diff2, orates)
  tmpres[['divdiff']][i] = c[[1]]
  tmpres[['divdiff_sig']][i] = c[[2]]
}
res[1, 'analysis_type'] = 'bootstrap loci Fst'
for (i in 2:length(tmpres)) {
  res[1, names(tmpres)[[i]]] = mean(tmpres[[i]])
}

# 2 bootstrapping all pairwise F_ST estimates and recalculating slopes
nrep = 100
tmpres = vector('list', length(cols))
names(tmpres) = cols
for (i in 1:length(tmpres)) {
  tmpres[[i]] = rep(NA, nrep)
}

slopes = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.bootinds.csv", stringsAsFactors = F)
for (i in 1:nrep) {
  slopes2 = slopes[slopes$rep == i, ]
  slopes2 = slopes2[complete.cases(slopes2), ]
  diff2 = diff1
  rownames(diff2) = diff2$OTU
  inds = diff2$OTU[diff2$OTU %in% slopes2$OTU]
  diff2[inds, "slope"] = slopes2[match(inds, slopes2$OTU), "slope"]
  tree2 = add_tips(otree, diff2$OTU)
  tree3 = prune_tree(tree2, diff2$OTU)
  tmpres[['ntips']][i] = Ntip(tree3)
  
  lambda = phylosig(tree3, diff2[tree3$tip.label, "slope"], method="lambda", test=T)
  tmpres[['lambda']][i] = lambda$lambda
  tmpres[['lambda_sig']][i]= lambda$P
  a = run_phylolm(tree3, diff2, 'mean_pi', 'slope')
  tmpres[['pi']][i] = a[[1]]
  tmpres[['pi_sig']][i] = a[[2]]
  b = run_phylolm(tree3, diff2, 'PC2', 'slope')
  tmpres[['limb']][i] = b[[1]]
  tmpres[['limb_sig']][i] = b[[2]]
  c = run_diffdiv(tree3, diff2, orates)
  tmpres[['divdiff']][i] = c[[1]]
  tmpres[['divdiff_sig']][i] = c[[2]]
}
res[2, 'analysis_type'] = 'bootstrap inds Fst'
for (i in 2:length(tmpres)) {
  res[2, names(tmpres)[[i]]] = mean(tmpres[[i]])
}

# 3 removing pairwise F_ST comparisons <20 km and recalculating slopes
diff = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.filtered_dist.csv", stringsAsFactors = F)
diff1 = diff[diff$distance_type == 'geo_dist' & diff$genetic_type == 'fst', ]
diff1$mean_pi = pi[match(diff1$cluster, pi$cluster), "mean_pi"]
diff1 = cbind(diff1, morph[match(diff1$OTU, morph$OTU), 3:ncol(morph)])
diff1 = cbind(diff1, r[match(diff1$OTU, r$OTU), 2:ncol(r)])
res[3, 'analysis_type'] = 'remove close Fst'
diff2 = diff1[diff1$ninds >= 5, ]
tree2 = add_tips(otree, diff2$OTU)
tree3 = prune_tree(tree2, diff2$OTU)
res[3, 'ntips'] = Ntip(tree3)
rownames(diff2) = diff2$OTU
lambda = phylosig(tree3, diff2[tree3$tip.label, "slope"], method="lambda", test=T)
res[3, 'lambda'] = lambda$lambda
res[3, 'lambda_sig'] = lambda$P
a = run_phylolm(tree3, diff2, 'mean_pi', 'slope')
res[3, 'pi'] = a[[1]]
res[3, 'pi_sig'] = a[[2]]
b = run_phylolm(tree3, diff2, 'PC2', 'slope')
res[3, 'limb'] = b[[1]]
res[3, 'limb_sig'] = b[[2]]
c = run_diffdiv(tree3, diff2, orates)
res[3, 'divdiff'] = c[[1]]
res[3, 'divdiff_sig'] = c[[2]]

# 4 remove pairwise F_ST comparisons where F_ST > 0.7
diff = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.filtered_fst.csv", stringsAsFactors = F)
diff1 = diff[diff$distance_type == 'geo_dist' & diff$genetic_type == 'fst', ]
diff1$mean_pi = pi[match(diff1$cluster, pi$cluster), "mean_pi"]
diff1 = cbind(diff1, morph[match(diff1$OTU, morph$OTU), 3:ncol(morph)])
diff1 = cbind(diff1, r[match(diff1$OTU, r$OTU), 2:ncol(r)])
res[4, 'analysis_type'] = 'filtered fst'
diff2 = diff1[diff1$ninds >= 5, ]
tree2 = add_tips(otree, diff2$OTU)
tree3 = prune_tree(tree2, diff2$OTU)
res[4, 'ntips'] = Ntip(tree3)
rownames(diff2) = diff2$OTU
lambda = phylosig(tree3, diff2[tree3$tip.label, "slope"], method="lambda", test=T)
res[4, 'lambda'] = lambda$lambda
res[4, 'lambda_sig'] = lambda$P
a = run_phylolm(tree3, diff2, 'mean_pi', 'slope')
res[4, 'pi'] = a[[1]]
res[4, 'pi_sig'] = a[[2]]
b = run_phylolm(tree3, diff2, 'PC2', 'slope')
res[4, 'limb'] = b[[1]]
res[4, 'limb_sig'] = b[[2]]
c = run_diffdiv(tree3, diff2, orates)
res[4, 'divdiff'] = c[[1]]
res[4, 'divdiff_sig'] = c[[2]]

# 5 restricting analyses to only pairwise F_ST estimated from the same biome (here, either desert, XXX, or XXX)
diff = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.biomes.csv", stringsAsFactors = F)
diff1 = diff[diff$biome == 'desert', ]
diff1 = diff1[complete.cases(diff1$slope),]
diff = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.csv", stringsAsFactors = F)
diff1$cluster = diff[match(diff1$OTU, diff$OTU), "cluster"]
diff1$mean_pi = pi[match(diff1$cluster, pi$cluster), "mean_pi"]
diff1 = cbind(diff1, morph[match(diff1$OTU, morph$OTU), 3:ncol(morph)])
diff1 = cbind(diff1, r[match(diff1$OTU, r$OTU), 2:ncol(r)])
res[5, 'analysis_type'] = 'only desert'
diff2 = diff1
tree2 = add_tips(otree, diff2$OTU)
tree3 = prune_tree(tree2, diff2$OTU)
res[5, 'ntips'] = Ntip(tree3)
rownames(diff2) = diff2$OTU
lambda = phylosig(tree3, diff2[tree3$tip.label, "slope"], method="lambda", test=T)
res[5, 'lambda'] = lambda$lambda
res[5, 'lambda_sig'] = lambda$P
a = run_phylolm(tree3, diff2, 'mean_pi', 'slope')
res[5, 'pi'] = a[[1]]
res[5, 'pi_sig'] = a[[2]]
b = run_phylolm(tree3, diff2, 'PC2', 'slope')
res[5, 'limb'] = b[[1]]
res[5, 'limb_sig'] = b[[2]]
c = run_diffdiv(tree3, diff2, orates)
res[5, 'divdiff'] = c[[1]]
res[5, 'divdiff_sig'] = c[[2]]

# 6 remove species with <5 individuals
diff = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.csv", stringsAsFactors = F)
diff1 = diff[diff$distance_type == 'geo_dist' & diff$genetic_type == 'fst', ]
diff1$mean_pi = pi[match(diff1$cluster, pi$cluster), "mean_pi"]
diff1 = cbind(diff1, morph[match(diff1$OTU, morph$OTU), 3:ncol(morph)])
diff1 = cbind(diff1, r[match(diff1$OTU, r$OTU), 2:ncol(r)])
res[6, 'analysis_type'] = '<5 inds'
diff2 = diff1[diff1$ninds >= 5, ]
tree2 = add_tips(otree, diff2$OTU)
tree3 = prune_tree(tree2, diff2$OTU)
res[6, 'ntips'] = Ntip(tree3)
rownames(diff2) = diff2$OTU
lambda = phylosig(tree3, diff2[tree3$tip.label, "slope"], method="lambda", test=T)
res[6, 'lambda'] = lambda$lambda
res[6, 'lambda_sig'] = lambda$P
a = run_phylolm(tree3, diff2, 'mean_pi', 'slope')
res[6, 'pi'] = a[[1]]
res[6, 'pi_sig'] = a[[2]]
b = run_phylolm(tree3, diff2, 'PC2', 'slope')
res[6, 'limb'] = b[[1]]
res[6, 'limb_sig'] = b[[2]]
c = run_diffdiv(tree3, diff2, orates)
res[6, 'divdiff'] = c[[1]]
res[6, 'divdiff_sig'] = c[[2]]

# 7 remove species with non-significant rates of differentiation
res[7, 'analysis_type'] = 'sig slope only'
diff2 = diff1[diff1$mantel_pval <= 0.05, ]
tree2 = add_tips(otree, diff2$OTU)
tree3 = prune_tree(tree2, diff2$OTU)
res[7, 'ntips'] = Ntip(tree3)
rownames(diff2) = diff2$OTU
lambda = phylosig(tree3, diff2[tree3$tip.label, "slope"], method="lambda", test=T)
res[7, 'lambda'] = lambda$lambda
res[7, 'lambda_sig'] = lambda$P
a = run_phylolm(tree3, diff2, 'mean_pi', 'slope')
res[7, 'pi'] = a[[1]]
res[7, 'pi_sig'] = a[[2]]
b = run_phylolm(tree3, diff2, 'PC2', 'slope')
res[7, 'limb'] = b[[1]]
res[7, 'limb_sig'] = b[[2]]
c = run_diffdiv(tree3, diff2, orates)
res[7, 'divdiff'] = c[[1]]
res[7, 'divdiff_sig'] = c[[2]]

# 8 only include nominal species
diff = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.csv", stringsAsFactors = F)
diff1 = diff[diff$distance_type == 'geo_dist' & diff$genetic_type == 'fst', ]
diff1$mean_pi = pi[match(diff1$cluster, pi$cluster), "mean_pi"]
diff1 = cbind(diff1, morph[match(diff1$OTU, morph$OTU), 3:ncol(morph)])
diff1 = cbind(diff1, r[match(diff1$OTU, r$OTU), 2:ncol(r)])
diff1$SPECIES = d[match(diff1$OTU, d$OTU), "SPECIES"]

res[8, 'analysis_type'] = 'nominal species'
diff2 = aggregate(diff1, by=list(diff1$SPECIES), mean, na.rm=T)
names(diff2)[1] = "SPECIES"
tree2 = add_tips(stree, diff2$SPECIES)
tree3 = prune_tree(tree2, diff2$SPECIES)
res[8, 'ntips'] = Ntip(tree3)
rownames(diff2) = diff2$SPECIES
lambda = phylosig(tree3, diff2[tree3$tip.label, "slope"], method="lambda", test=T)
res[8, 'lambda'] = lambda$lambda
res[8, 'lambda_sig'] = lambda$P
a = run_phylolm(tree3, diff2, 'mean_pi', 'slope')
res[8, 'pi'] = a[[1]]
res[8, 'pi_sig'] = a[[2]]
b = run_phylolm(tree3, diff2, 'PC2', 'slope')
res[8, 'limb'] = b[[1]]
res[8, 'limb_sig'] = b[[2]]
c = run_diffdiv(tree3, diff2, srates)
res[8, 'divdiff'] = c[[1]]
res[8, 'divdiff_sig'] = c[[2]]

# 9 remove recently-diverged species
diff = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.csv", stringsAsFactors = F)
diff1 = diff[diff$distance_type == 'geo_dist' & diff$genetic_type == 'fst', ]
diff1$mean_pi = pi[match(diff1$cluster, pi$cluster), "mean_pi"]
diff1 = cbind(diff1, morph[match(diff1$OTU, morph$OTU), 3:ncol(morph)])
diff1 = cbind(diff1, r[match(diff1$OTU, r$OTU), 2:ncol(r)])

bt = branching.times(otree)
names(bt) = seq(Ntip(otree) + 1, Nnode(otree) + Ntip(otree))
drop = bt[bt < 5]
dropTips = rep(NA, length(drop))
for (i in 1:length(drop)) {
  tips = otree$tip.label[Descendants(otree, as.numeric(names(drop)[i]), type="tips")[[1]]]
  if (length(tips) == 2) {
    droptip = sample(tips, 1)
    dropTips[i] = droptip
    if (droptip %in% diff1$OTU) {
      diff1[diff1$OTU == droptip, 'OTU'] = tips[!tips %in% droptip]
    }
  }
}
diff2 = aggregate(diff1, by=list(diff1$OTU), mean, na.rm=T)
names(diff2)[1] = "OTU"

res[9, 'analysis_type'] = 'drop young species'
otree2 = drop.tip(otree, dropTips[complete.cases(dropTips)])
tree2 = add_tips(otree2, diff2$OTU)
tree3 = prune_tree(tree2, diff2$OTU)
res[9, 'ntips'] = Ntip(tree3)
rownames(diff2) = diff2$OTU
lambda = phylosig(tree3, diff2[tree3$tip.label, "slope"], method="lambda", test=T)
res[9, 'lambda'] = lambda$lambda
res[9, 'lambda_sig'] = lambda$P
a = run_phylolm(tree3, diff2, 'mean_pi', 'slope')
res[9, 'pi'] = a[[1]]
res[9, 'pi_sig'] = a[[2]]
b = run_phylolm(tree3, diff2, 'PC2', 'slope')
res[9, 'limb'] = b[[1]]
res[9, 'limb_sig'] = b[[2]]
c = run_diffdiv(tree3, diff2, jetzDivRates(otree2))
res[9, 'divdiff'] = c[[1]]
res[9, 'divdiff_sig'] = c[[2]]

# 10 repeat analyses with a random selection of 80% of the taxa
diff = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.csv", stringsAsFactors = F)
diff1 = diff[diff$distance_type == 'geo_dist' & diff$genetic_type == 'fst', ]
diff1$mean_pi = pi[match(diff1$cluster, pi$cluster), "mean_pi"]
diff1 = cbind(diff1, morph[match(diff1$OTU, morph$OTU), 3:ncol(morph)])
diff1 = cbind(diff1, r[match(diff1$OTU, r$OTU), 2:ncol(r)])

nrep = 100
tmpres = vector('list', length(cols))
names(tmpres) = cols
for (i in 1:length(tmpres)) {
  tmpres[[i]] = rep(NA, nrep)
}

for (i in 1:nrep) {
  diff2 = diff1[sample(seq(1, nrow(diff1)), nrow(diff1) * 0.8), ]
  tmpres[['analysis_type']][i] = 'random 80%'
  tree2 = add_tips(otree, diff2$OTU)
  tree3 = prune_tree(tree2, diff2$OTU)
  tmpres[['ntips']][i] = Ntip(tree3)
  rownames(diff2) = diff2$OTU
  lambda = phylosig(tree3, diff2[tree3$tip.label, "slope"], method="lambda", test=T)
  tmpres[['lambda']][i] = lambda$lambda
  tmpres[['lambda_sig']][i]= lambda$P
  a = run_phylolm(tree3, diff2, 'mean_pi', 'slope')
  tmpres[['pi']][i] = a[[1]]
  tmpres[['pi_sig']][i] = a[[2]]
  b = run_phylolm(tree3, diff2, 'PC1', 'slope')
  tmpres[['limb']][i] = b[[1]]
  tmpres[['limb_sig']][i] = b[[2]]
  c = run_diffdiv(tree3, diff2, orates)
  tmpres[['divdiff']][i] = c[[1]]
  tmpres[['divdiff_sig']][i] = c[[2]]
}
res[10, 'analysis_type'] = 'random 80%'
for (i in 2:length(tmpres)) {
  res[10, names(tmpres)[[i]]] = mean(tmpres[[i]])
}

# 11 for phylogenetic uncertainty by repeating our analyses across 100 samples from the posterior distribution of our phylogenetic analysis
nrep = 100
tmpres = vector('list', length(cols))
names(tmpres) = cols
for (i in 1:length(tmpres)) {
  tmpres[[i]] = rep(NA, nrep)
}

otrees = read.tree("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/beast_ucln_31July17/otu_posterior.trees")
for (i in 1:nrep) {
  diff2 = diff1
  otree2 = otrees[[i]]
  tree2 = add_tips(otree2, diff2$OTU)
  tree3 = prune_tree(tree2, diff2$OTU)
  tmpres[['ntips']][i] = Ntip(tree3)
  rownames(diff2) = diff2$OTU
  lambda = phylosig(tree3, diff2[tree3$tip.label, "slope"], method="lambda", test=T)
  tmpres[['lambda']][i] = lambda$lambda
  tmpres[['lambda_sig']][i]= lambda$P
  a = run_phylolm(tree3, diff2, 'mean_pi', 'slope')
  tmpres[['pi']][i] = a[[1]]
  tmpres[['pi_sig']][i] = a[[2]]
  b = run_phylolm(tree3, diff2, 'PC2', 'slope')
  tmpres[['limb']][i] = b[[1]]
  tmpres[['limb_sig']][i] = b[[2]]
  c = run_diffdiv(tree3, diff2, jetzDivRates(otree2))
  tmpres[['divdiff']][i] = c[[1]]
  tmpres[['divdiff_sig']][i] = c[[2]]
}
res[11, 'analysis_type'] = 'tree posterior'
for (i in 2:length(tmpres)) {
  res[11, names(tmpres)[[i]]] = mean(tmpres[[i]])
}

# 12 d_xy and geodist
res[12, 'analysis_type'] = 'dxy slope across geo'
diff1 = diff[diff$distance_type == 'geo_dist' & diff$genetic_type == 'nuc_dxy', ]
diff1$mean_pi = pi[match(diff1$cluster, pi$cluster), "mean_pi"]
diff1 = cbind(diff1, morph[match(diff1$OTU, morph$OTU), 3:ncol(morph)])
diff1 = cbind(diff1, r[match(diff1$OTU, r$OTU), 2:ncol(r)])
diff2 = diff1
tree2 = add_tips(otree, diff2$OTU)
tree3 = prune_tree(tree2, diff2$OTU)
res[12, 'ntips'] = Ntip(tree3)
rownames(diff2) = diff2$OTU
lambda = phylosig(tree3, diff2[tree3$tip.label, "slope"], method="lambda", test=T)
res[12, 'lambda'] = lambda$lambda
res[12, 'lambda_sig'] = lambda$P
a = run_phylolm(tree3, diff2, 'mean_pi', 'slope')
res[12, 'pi'] = a[[1]]
res[12, 'pi_sig'] = a[[2]]
b = run_phylolm(tree3, diff2, 'PC2', 'slope')
res[12, 'limb'] = b[[1]]
res[12, 'limb_sig'] = b[[2]]
c = run_diffdiv(tree3, diff2, orates)
res[12, 'divdiff'] = c[[1]]
res[12, 'divdiff_sig'] = c[[2]]

# 13 fst and envdist
res[13, 'analysis_type'] = 'fst slope across env'
diff1 = diff[diff$distance_type == 'env_dist' & diff$genetic_type == 'fst', ]
diff1$mean_pi = pi[match(diff1$cluster, pi$cluster), "mean_pi"]
diff1 = cbind(diff1, morph[match(diff1$OTU, morph$OTU), 3:ncol(morph)])
diff1 = cbind(diff1, r[match(diff1$OTU, r$OTU), 2:ncol(r)])
diff2 = diff1
tree2 = add_tips(otree, diff2$OTU)
tree3 = prune_tree(tree2, diff2$OTU)
res[13, 'ntips'] = Ntip(tree3)
rownames(diff2) = diff2$OTU
lambda = phylosig(tree3, diff2[tree3$tip.label, "slope"], method="lambda", test=T)
res[13, 'lambda'] = lambda$lambda
res[13, 'lambda_sig'] = lambda$P
a = run_phylolm(tree3, diff2, 'mean_pi', 'slope')
res[13, 'pi'] = a[[1]]
res[13, 'pi_sig'] = a[[2]]
b = run_phylolm(tree3, diff2, 'PC2', 'slope')
res[13, 'limb'] = b[[1]]
res[13, 'limb_sig'] = b[[2]]
c = run_diffdiv(tree3, diff2, orates)
res[13, 'divdiff'] = c[[1]]
res[13, 'divdiff_sig'] = c[[2]]

write.csv(res, "~/Dropbox/Sphenomorphine_Gene_Flow/data/robustness_analyses.csv", row.names = F)