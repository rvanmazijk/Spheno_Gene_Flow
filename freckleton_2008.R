
library(ape)
 

speciationRateNodeDensity <- function(phy){
	
	dmat <- vcv.phylo(phy)
	
	phy$edge.length <- rep(1, length(phy$edge.length))
	xmat <- vcv.phylo(phy)
	dg <- diag(xmat) - 1
	
	rates <- dg / diag(dmat)
	return(rates)
	
}

