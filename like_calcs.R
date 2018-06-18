library(ape)

v <- read.tree("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/beast_ucln_31July17/otu.tre")


# mrca node
mm <- getMRCA(v, c("Ctenotus_pantherinus", "Lerista_labialis"))

#-------------------------------
# the one-rate ML solution
#    Single overall speciation rate for tree

N <- length(v$tip.label)

# maximum likelihood estimate of lambda, using analytical
lambda_1rate <- (N - 2) / sum(v$edge.length)

# exact log-likelihood at the maximum:
like1rate <- (N-2)*log(lambda_1rate) - lambda_1rate*sum(v$edge.length)


#---------------------------
# Likelihood under 2 rate model, 
#    where Lerista and Ctenotus have separate rates:

v_ctle <- extract.clade(v, node = mm)
v_other <- drop.tip(v, v_ctle$tip.label[2:length(v_ctle$tip.label)])

time_ctle <- sum(v_ctle$edge.length)
time_other <- sum(v_other$edge.length) - max(branching.times(v_ctle))

lambda_ctle <- (length(v_ctle$tip.label)-1) / time_ctle
lambda_other <- (length(v_other$tip.label) - 2) / time_other

like_ctle <- (length(v_ctle$tip.label) - 1) * log(lambda_ctle) - lambda_ctle * time_ctle

like_other <- (length(v_other$tip.label) - 2) * log(lambda_other) - lambda_other * time_other

like2rate <- like_ctle + like_other



#---------------------------
# Likelihood under 2 rate model, again
#    but constraining to overall ML rate (to check our work)

 
v_ctle <- extract.clade(v, node = mm)
v_other <- drop.tip(v, v_ctle$tip.label[2:length(v_ctle$tip.label)])

time_ctle <- sum(v_ctle$edge.length)
time_other <- sum(v_other$edge.length) - max(branching.times(v_ctle))

lambda_ctle <- lambda_1rate
lambda_other <- lambda_1rate

like_ctle <- (length(v_ctle$tip.label) - 1) * log(lambda_ctle) - lambda_ctle * time_ctle

like_other <- (length(v_other$tip.label) - 2) * log(lambda_other) - lambda_other * time_other

like_ctle + like_other
# should give same as like1rate

#-------------------------------

# Likelihood ratio test

lr <- 2 * (like2rate - like1rate)


# p-value if we treat this as 1 parameter difference:
# 	e.g., only parameter difference is the Lerista/Ctenotus speciation rate

1 - pchisq(lr, df = 1)


# One could make the case that there are 2 parameters 
#    the extra speciation rate + location parameter, but I'd say it could 
#    go either way in this case. This is more conservative, though:

1 - pchisq(lr, df = 2)


#-------------------------------

# For the heck of it, let us test whether 
#   tree is more imbalanced than expected under a constant-rate process
#  We will use the Colless imbalance statistic:

library(geiger)

colless <- function(phy){
  
  bb <- balance(phy);
  ss <- sum(abs(bb[,1] - bb[,2]));
  n <- length(phy$tip.label);
  return((2 / ((n-1)*(n-2))) * ss);
  
}

Ntips <- length(v$tip.label)

# Generate null distribution of balance statistics
#  under a pure-birth process. Is tree more imbalanced? If so,
#   implies diversification rate heterogeneity

REPS <- 2000
null_dist <- rep(NA, REPS)

for (ii in 1:REPS){
	
	tmp <- sim.bdtree(stop = "taxa", n = Ntips)
	null_dist[ii] <- colless(tmp)
	
	
}

obs <- colless(v)


geom_hist(null_dist, breaks=50, xlim=c(0, 0.07))
d = data.frame(null_dist=null_dist)
a = ggplot(d, aes(null_dist)) + geom_histogram(bins=50, fill="gray") + 
  labs(x="null distribution") +
  geom_vline(xintercept = obs, colour="red", linetype=2)
save_plot(paste(figdir, "S8_treebalance.pdf", sep=""), a)

# Assess statistical significance: two-tailed p-value
2* sum(null_dist > obs) / (1 + length(null_dist))  
effsize = (obs - mean(null_dist)) / sd(null_dist)
# So, VERY strong evidence for rate heterogeneity 
#    from a tree-balance perspective




