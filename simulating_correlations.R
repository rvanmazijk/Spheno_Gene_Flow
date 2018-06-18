library(mvtnorm)
library(ape)
library(phytools)
library(nlme)

mytree <- read.tree("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/beast_ucln_31July17/otu.tre")

diff = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.csv", stringsAsFactors = F)
diff = diff[diff$distance_type == 'geo_dist' & diff$genetic_type == 'fst',]
diff$log_slope = log(diff$slope)
diff = diff[complete.cases(diff$log_slope), ]
ibd = diff$log_slope
names(ibd) = diff$OTU

keep = intersect(diff$OTU, mytree$tip.label)
source("~/scripts/gene_flow/jetz_div_rates.R")
source("~/scripts/gene_flow/es-sim.R")
rates = jetzDivRates(mytree)
vv <- vcv.phylo(mytree, corr = TRUE)
# let xvals be log-transformed measures of speciation rate
#    in same order as tip labels, so
#    if mytraits is vector of rates with names
rates <- rates[mytree$tip.label]
xvals <- log(rates)

# Define a target correlation 
rcorrs = seq(0, 1, by=0.1)
nrep = c(100)
res = vector('list', length(rcorrs))
for (i in 1:length(res)) {
  res[[i]] = list(rep(NA, nrep), rep(NA, nrep), rep(NA, nrep))
  names(res[[i]]) = c("lambda", "sig", "cor")
}
names(res) = rcorrs
for (j in 1:length(rcorrs)) {
  cat(j, "\n")
  for (i in 1:nrep) {
    # Generate traits under Brownian motion:
    d1 <- as.vector(rmvnorm(1, sigma=vv))

    # rotate them by the target correlation with the speciation rates:
    newvals <- (rcorrs[j] * xvals) + sqrt(1 - rcorrs[j]^2)*d1
    # The newvals will have the approximate correlation desired 
    #   and will also have the right correlation due to phylogeny
    
    mytree2 = drop.tip(mytree, setdiff(mytree$tip.label, keep))
    newvals2 = newvals[names(newvals) %in% keep]
    xvals2 = xvals[names(xvals) %in% keep]

    # test lambda
    res[[j]]$lambda[i] = phylosig(mytree2, newvals2, method="lambda")$lambda
  
    # test significance
    d = data.frame(xvals2, newvals2)
    # g = essim(mytree2, newvals, nsim=1000, xvals, take_log = FALSE)
    g = gls(xvals2 ~ newvals2, correlation=corBrownian(phy=mytree2), data=d, method="ML")
    res[[j]]$sig[i] = summary(g)$tTable[2, 4]
    res[[j]]$cor[i] = summary(g)$tTable[2, 1]
    # res[[j]]$sig[i] = g[2]
    # res[[j]]$cor[i] = g[1]
  }
}
res1 = data.frame(unlist(lapply(res, function(x) {length(x$sig[x$sig < 0.05])})) / 100)
names(res1) = c("power")
res1$cor = rownames(res1)

write.csv(res1, "~/Dropbox/Sphenomorphine_Gene_Flow/data/power_analysis/PGLS_power.csv", row.names=F)
# write.csv(res1, "~/Dropbox/Sphenomorphine_Gene_Flow/data/power_analysis/ES-SIM_power.csv", row.names=F)
