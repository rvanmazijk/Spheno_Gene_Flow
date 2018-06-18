library(vegan)

# types of distance
dist = c("geo_dist", "env_dist")

# types of genetic measures
gen = c("nuc_dxy", "mt_dxy", "fst")

# things to track
## MMRR results

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

get_matrix = function(df) {
	inds = unique(c(df$ind1, df$ind2))
	m = matrix(NA, nrow=length(inds), ncol=length(inds))
	
	for (a in 1:length(inds)) {
		for (b in 1:length(inds)) {
			if (a != b) {
				tmp1 = df[(df$ind1 == inds[a]) & (df$ind2 == inds[b]),]
				tmp2 = df[(df$ind2 == inds[a]) & (df$ind1 == inds[b]),]
				if (nrow(tmp1) == 1) {
					m[a, b] = tmp1[1, 3]
					m[b, a] = tmp1[1, 3]
				} else if (nrow(tmp2) == 1) {
					m[a, b] = tmp2[1, 3]
					m[b, a] = tmp2[1, 3]
				} 
			}
		}
	}
	return(m)
}
		

envdist = readRDS("~/Dropbox/Sphenomorphine_Gene_Flow/data/envirospatial/envdist.rds")
divfiles = list.files("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/species/", full.names=T)

# set up result vector with the maximum number of comparisons
cols = c("cluster", "OTU", "ninds", "distance_type", "genetic_type", "slope", "intercept", "r2", "pval", "mantel_corr", "mantel_pval")
res = data.frame(matrix(NA, nrow=length(divfiles) * 6, ncol=length(cols)))
names(res) = cols

# remove populations that are too close?
filt_dist = FALSE
# filt_dist = TRUE

# remove populations that have too high fst
# filt_fst = TRUE
filt_fst = FALSE


row = 1
for (i in 1:length(divfiles)) {
	d = read.csv(divfiles[i], stringsAsFactors=F)
	cl = unique(d$cl)
	otu = gsub('^.*/', '', divfiles[i])
	otu = gsub('.diver.*', '', otu)

	# immediately get rid of any comparisons based on too little data
	d = d[d$fst_denom >= 1000,]
	ninds = length(unique(c(d$ind1, d$ind2)))
	
	# only keep comparisons with at least 3 inds
	if (ninds >= 3) {
		# for rousset method
		d$inv_fst = sapply(d$fst, inv_fst)
		d$dist_log = sapply(d$geo_dist, log_dist)
	
		# get env distances
		d$env_dist = rep(NA, nrow(d))
		for (x in 1:nrow(d)) {
			d[x, "env_dist"] = envdist[d[x, "ind1"], d[x, "ind2"]]
		}
		
		for (x in 1:length(dist)) {
			for (y in 1:length(gen)) {
				
				# account for differences for model
				if (gen[y] == 'fst') {
				  if (filt_fst) {
				    d = d[d$inv_fst <= 3, ]
				  }
				  
					if (dist[x] == 'geo_dist') {
						disttype = 'dist_log'
						if (filt_dist) {
							d = d[which(d$dist_log >= 8),]
						}
					} else {
						disttype = dist[x]
					}
					gentype = 'inv_fst'
				} else {
					disttype = dist[x]
					gentype = gen[y]
				}
				
				dd = d[, c("ind1", "ind2", gentype, disttype)]
				dd = dd[complete.cases(dd),]
				if (nrow(dd) >= 3) {
					# get matrix 1
					m1 = get_matrix(dd[, c("ind1", "ind2", gentype)])
					m2 = get_matrix(dd[, c("ind1", "ind2", disttype)])
					mt = try(mantel(m1, m2, na.rm=T))
					if (class(mt) == 'try-error') {
					  mt = list(NA, NA)
					  names(mt) = c('statistic', 'signif')
					}
					
					model = lm(dd[, 3] ~ dd[, 4])
				
					res[row, "cluster"] = cl
					res[row, "OTU"] = otu
					res[row, "ninds"] = ninds
					res[row, "distance_type"] = dist[x]
					res[row, "genetic_type"] = gen[y]
					res[row, "slope"] = coef(model)[2]
					res[row, "intercept"] = coef(model)[1]
					res[row, "r2"] = summary(model)$r.squared
					if (dim(summary(model)$coefficients)[1] > 1) {
						res[row, "pval"] =  summary(model)$coefficients[2, 4]
					} else {
						res[row, "pval"] = NA
					}	
					res[row, "mantel_corr"] = mt$statistic
					res[row, "mantel_pval"] = mt$signif 
					row = row + 1
					}
				}
			}
		}
	}
	
res = res[complete.cases(res$cl),]
# add in genus
res$genus = gsub("_.*", "", res$OTU)

if (filt_dist) {
	write.csv(res, "~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.filtered_dist.csv", row.names=F)
} else	{
  if (filt_fst) {
    write.csv(res, "~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.filtered_fst.csv", row.names=F)
  } else {
    write.csv(res, "~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.csv", row.names=F)
  }
}
