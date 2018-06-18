

# phy: the tree
# trait: vector of trait data
#         with names
#    is: vector of equal-splits values (eg DR statistics).
#          if not specified, function will compute from the tree.




essim <- function(phy, trait, nsim = 1000, is, take_log = FALSE) {
	
	require(ape) 
	require(phytools)
	
	if(missing(is)) { # If inverse equal splits statistics not provided, calculate it
		rootnode <- length(phy$tip.label) + 1
		is <- numeric(length(phy$tip.label))
		for (i in 1:length(is)){
			node <- i
			index <- 1
			qx <- 0
			while (node != rootnode){
				el <- phy$edge.length[phy$edge[,2] == node]
				node <- phy$edge[,1][phy$edge[,2] == node]			
				qx <- qx + el* (1 / 2^(index-1))			
				index <- index + 1
			}
			is[i] <- 1/qx
		}		
		names(is) <- phy$tip.label
	}
	
  if (take_log ) {
	  is <- log(is[phy$tip.label]) # log transform
  } else {
    is <- is[phy$tip.label]
  }
  
	trait <- trait[phy$tip.label]
	
	# Pearson's correlation between splits statistic and trait
	res <- cor.test(is, trait, method="pearson")

	# Brownian simulations 
	sims <- fastBM(phy, nsim = nsim)
	#rownames(sims) <- rownames(vv)  # already have rownames
		
	# Pearson's correlations of simulated datasets
	sim.r <- sapply(1:nsim, function(x) cor.test(is[as.vector(rownames(sims))], sims[,x], method="pearson")$estimate)
	
	# Calculate the two-tailed p value
	corr <- res$estimate
	upper <- (length(sim.r[sim.r >= corr])+1)/(nsim+1)
	lower <- (length(sim.r[sim.r <= corr])+1)/(nsim+1)
	pval <- 2*min(c(upper,lower)) # Remove "2" for one-tailed

	result <- as.vector(c(corr, pval))
	names(result) <- c("rho", "P Value")
	return(result)

}


