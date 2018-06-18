files = list.files("~/Dropbox/Sphenomorphine_Gene_Flow/data/species_revisions/prelim_divergence/", full.names=T)

# returns the inverse of Fst
inv_fst = function(x) {
	if (x == 1) {
		return(NA)
	} else {
	return(x / (1 - x))
	}
}

# returns log distance
# note that this drops all comparisons
# between points at the same place
log_dist = function(x) {
	if (x == 0) {
		return(NA)
	} else {
		return(log(x))
	}
}

pdf("~/Desktop/prelim_Fst.pdf", height=4, width=10)
for (i in 1:length(files)) {
	d = read.csv(files[i], stringsAsFactors=F)
	
	sp = gsub(".*//", "", files[i])
	sp = gsub(".divergence_cov10.csv", "", sp)
	
	if (nrow(d) > 2) {
		d = d[complete.cases(d$fst),]
		par(mfrow=c(1,3), las=1)
		# dxy vs geo dist
		cor1 = cor.test(d$geo_dist, d$nuc_dxy, method="spearman")$estimate
		sig1 = cor.test(d$geo_dist, d$nuc_dxy, method="spearman")$p.value
		plot(d$geo_dist, d$nuc_dxy, pch=16, main=paste(sp, "; r=", round(cor1, 3), "; p=", round(sig1, 3), sep=""), xlab="distance (m)", ylab="nuc_dxy", )
	
		# fst vs geo dist
		d$inv_fst = sapply(d$fst, inv_fst)
		d$dist_log = sapply(d$geo_dist, log_dist) 
		cor2 = cor.test(d$dist_log, d$inv_fst, method="spearman")$estimate
		sig2 = cor.test(d$dist_log, d$inv_fst, method="spearman")$p.value
		plot(d$dist_log, d$inv_fst, pch=16, main=paste(sp, "; r=", round(cor2, 3), "; p=", round(sig2, 3), sep=""), xlab="log dist", ylab="inv fst", )
	
		# plot denom hist
		hist(d$fst_denom, breaks=10)
		}
}
dev.off()