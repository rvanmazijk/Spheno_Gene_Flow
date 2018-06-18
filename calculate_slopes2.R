# normalize a matrix to make it more comparable
normalize <- function(mat) {
	mat = (mat - mean(mat, na.rm=T)) / sd(mat, na.rm=T)
	return(mat)
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

# types of genetic measures
gen = c("nuc_dxy", "mt_dxy", "fst")

source("/Users/sonal/scripts/eco_IBD/MMRR.R")

envdist = readRDS("~/Dropbox/Sphenomorphine_Gene_Flow/data/envirospatial/envdist.rds")
divfiles = list.files("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/species/", full.names=T)
cols = c("cluster", "OTU", "ninds", "genetic_type", "fit", "fit_pval", "geo_coeff", "env_coeff", "geo_pval", "env_pval", "env_geo_corr", "env_geo_pval")
res = data.frame(matrix(NA, nrow=length(divfiles) * 3, ncol=length(cols)))
names(res) = cols

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
		# get env distances
		d$env_dist = rep(NA, nrow(d))
		for (x in 1:nrow(d)) {
			d[x, "env_dist"] = env_dist[d[x, "ind1"], d[x, "ind2"]]
		}
		
		m_geo = get_matrix(d[, c("ind1", "ind2", "geo_dist")])
		m_geo = normalize(m_geo)
		m_env = get_matrix(d[, c("ind1", "ind2", "env_dist")])
		m_env = normalize(m_env)
		
		for (y in 1:length(gen)) {
			row = 3 * (i - 1) + y
			cat(row, "\n") 
			m_dna = get_matrix(d[, c("ind1", "ind2", gen[y])])
			m_dna = normalize(m_dna)
			
			model = try(MMRR(m_dna, list(m_geo, m_env)), silent=T)
			covar = try(MRR(m_geo, m_env), silent=T)

			res[row, "cluster"] = cl
			res[row, "OTU"] = otu
			res[row, "ninds"] = ninds
			res[row, "genetic_type"] = gen[y]

			if (class(model) != "try-error") {
				res[row, "fit"] = model$r.squared
				res[row, "fit_pval"] = model$Fpvalue
				res[row, "geo_coeff"] = model$coefficients[2]
				res[row, "env_coeff"] = model$coefficients[3]
				res[row, "geo_pval"] = model$tpvalue[2]
				res[row, "env_pval"] = model$tpvalue[3]
			}
			if (class(covar) != "try-error") {
				res[row, "env_geo_corr"] = covar$r.squared
				res[row, "env_geo_pval"] = covar$Fpvalue
			}
		}
	}		
}		

res = res[complete.cases(res$cluster),]
write.csv(res, "~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.MMRR.csv", row.names=F)
			