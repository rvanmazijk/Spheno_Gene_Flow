library(raster)
library(pcaMethods)
library(maptools)
library(rgeos)
library(ggplot2)
library(scales)
library(MASS)
library(lme4)
library(piecewiseSEM)

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

define_type <- function(x, y, per_cut) {
	# function defines the relationship between
	# two individuals w.r.t. the range
	if (x > per_cut & y > per_cut) {
		return('CEN-CEN')
	} else if (x <= per_cut & y <= per_cut) {
		return('PER-PER')
	} else {
		return('CEN-PER')
	}
}

polyEdgeDist <- function(pt, shape) {
     asPts = function(x) as(as(x, "SpatialLines"), "SpatialPoints")
     d = spDistsN1(asPts(shape), pt)
     return(min(d))
}

divdir = '~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/species/'
rangedir = '~/Dropbox/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/'
divfiles = list.files(divdir, pattern=".csv", full.names=T)
envdist = readRDS("~/Dropbox/Sphenomorphine_Gene_Flow/data/envirospatial/envdist.rds") 
  
# lat long
# get lat longs for all individuals
ll1 = read.csv("/Users/sonal/macroevolution/eco_IBD_oz/data/metadata/individual_data_nomissing22Oct15.csv", stringsAsFactors=F, na.string=c("", "NA"))
ll2 = read.csv("/Users/sonal/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/spheno_ind_data.csv", stringsAsFactors=F, na.string=c("", "NA"))
# some overlap, and ll2 is more complete
ll1 = ll1[!(ll1$sample_id %in% ll2$SAMPLE_ID), ]
ll1 = ll1[, c("sample_id", "lon", "lat")]
ll2 = ll2[, c("SAMPLE_ID", "LON", "LAT")]
names(ll2) = names(ll1)
ll = rbind(ll1, ll2)
ll = ll[complete.cases(ll$lon), ]
rownames(ll) = ll$sample_id

# get the geo & fst data
div = vector("list", length(divfiles))
for (i in 1:length(divfiles)) {
	d = read.csv(divfiles[i], stringsAsFactors=F)
	
	# get lat longs	
	d$lat1 = ll[match(d$ind1, ll$sample_id), "lat"]
	d$lon1 = ll[match(d$ind1, ll$sample_id), "lon"]
	d$lat2 = ll[match(d$ind2, ll$sample_id), "lat"]
	d$lon2 = ll[match(d$ind2, ll$sample_id), "lon"]
	
	# remove low quality, missing etc data
	d = d[!is.na(d$fst),]
	d = d[!is.na(d$geo_dist),]
	d = d[d$fst_denom > 1000,]
	d = d[!is.na(d$lat1),]
	d = d[!is.na(d$lat2),]
	
	inds = sort(unique(c(d$ind1, d$ind2)))
	ind_ll = cbind(inds, ll[match(inds, ll$sample_id), c("lon", "lat")])
	ind_ll = na.omit(ind_ll)
	
	##########################
	# range edge
	##########################
	cl = gsub("/Users/Sonal/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/species//", "", divfiles[i])
	cl = gsub(".divergence_cov10.csv", "", cl)
	d$cl = cl
	
	# get the species range and the center of its range
	range = paste(rangedir, cl, ".shp", sep="")
	if (!file.exists(range)) {
	  if (cl == 'Notoscincus_wotjulum') {
	    cl = 'Notoscincus_ornatus'
	  } else {
	    cl = gsub('_\\d', '', cl)
	  }
	  range = paste(rangedir, cl, '.shp', sep='')
	}
	range = readShapePoly(range)
	
	distances = rep(NA, 200)
	pts = spsample(range, 200, "random")
	rangeS = gSimplify(range, 0.1, topologyPreserve=T)
	for (x in 1:200) {
	 	distances[x] = polyEdgeDist(pts[x], rangeS)
		}
	per_cut = quantile(distances, 0.1)
	
	points1 = SpatialPoints(ind_ll[,c("lon", "lat")], proj4string=CRS('+proj=longlat +datum=WGS84'))
	ind_ll$distance = rep(NA, nrow(ind_ll))
	for (x in 1:nrow(ind_ll)) {
		ind_ll[x, "distance"] = polyEdgeDist(points1[x], rangeS)
	}
	
	d$distance1 = ind_ll[match(d$ind1, ind_ll$inds), "distance"]
	d$distance2 = ind_ll[match(d$ind2, ind_ll$inds), "distance"]
	d$type = mapply(define_type, d$distance1, d$distance2, per_cut)
	
	##########################
	# environmental distance
	##########################
 	d$env_dist = rep(NA)
 	for (j in 1:nrow(d)) {
 		ind1 = d[j, "ind1"]
 		ind2 = d[j, "ind2"]

		if (ind1 %in% rownames(envdist) & ind2 %in% rownames(envdist)) {
			d[j, "env_dist"] = envdist[ind1, ind2]
			}		
	 	}
	div[[i]] = d
	cat(d$cl[1], "\n")
}
keep = names(div[[1]])
div2 = vector('list', length(div))
for (i in 1:length(div)) {
  div2[[i]] = div[[i]][, keep] 
}
div2 = do.call(rbind, div2)

# get genus
div2$OTU = div2$cl
div2$cl = NULL
div2$genus = gsub("_\\S+", '', div2$OTU)

# get unique inds
inds = unique(c(div2$ind1, div2$ind2))

# get biomes
# the original biomes downloaded from the TNC
# http://maps.tnc.org/gis_data.html
biomes = readShapePoly("~/macroevolution/eco_IBD_oz/data/ecoregions/terr-ecoregions-aus_combined/terr_ecoregions-aus_combined.shp")
# simplify this so that it takes less time to process
x = gSimplify(biomes, 0.1, topologyPreserve=T)
x$biome = as.character(biomes$WWF_MHTNAM)
# rename the biomes
bnames1 = c("Tropical and Subtropical Moist Broadleaf Forests", "Montane Grasslands and Shrublands", "Temperate Broadleaf and Mixed Forests", "Tropical and Subtropical Grasslands, Savannas and Shrublands", "Deserts and Xeric Shrublands", "Mediterranean Forests, Woodlands and Scrub", "Temperate Grasslands, Savannas and Shrublands")
bnames2 = c("tropical forests", "montane", "temperate forests", "tropical grasslands", "desert", "Mediterranean forests", "temperate grasslands")
for (i in 1:length(bnames1)) {
	 x$biome = gsub(bnames1[i], bnames2[i],  x$biome)
}
rownames(ll) = ll$sample_id
ll_b = extract(x, ll[inds, c("lon", "lat")])
rownames(ll_b) = ll[inds, "sample_id"]
b1 = ll_b[div2$ind1, "biome"]
b2 = ll_b[div2$ind2, "biome"]
div2$biome = rep(NA, nrow(div2))
for (i in 1:nrow(div2)) {
	if (!is.na(b1[i]) & !is.na(b2[i])) {
		if (b1[i] == b2[i]) {
			div2[i, "biome"] = b1[i]
			} else {
			div2[i, "biome"] = "DIFF"
			}
		}	
	}

# factors
# env dist X
# geo dist X
# range location X
# biomes X
# genus X (random effect)
# species X (random nested effect)
div_all = div2
write.csv(div_all, "~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/pairwise_Fst_comparisons.csv", row.names=F)

############
# what distribution fits my data
############
fst = div2[div2$fst > 0, "fst"]
# normal
qqp(fst, "norm")
# lognormal
qqp(fst, "lnorm")
# negative binominal
# fails ... why?
nbinom <- fitdistr(fst, "Negative Binomial")
qqp(fst, "nbinom", size = nbinom$estimate[[1]], mu = nbinom$estimate[[2]])
# poisson
poisson <- fitdistr(fst, "Poisson")
qqp(fst, "pois", poisson$estimate)
# gamma
gamma <- fitdistr(fst, "gamma")
qqp(fst, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

######################################
# fst seems to best fit a normal model
######################################
div3 = div2[complete.cases(div2[,c("fst", "geo_dist", "env_dist", "type", "biome")]),]
div3$geo_dist = rescale(div3$geo_dist, to=c(0,1))
div3$env_dist = rescale(div3$env_dist, to=c(0,1))
div3$type = as.factor(div3$type)
div3$biome = as.factor(div3$biome)
lmm0 = lmer(fst ~ 1 + (1 | OTU), data=div3, REML=FALSE)
lmm1 = lmer(fst ~ geo_dist + env_dist + type + biome + (1 | OTU), data=div3, REML=FALSE)
lmm2 = lmer(fst ~ geo_dist + env_dist + type + (1 | OTU), data=div3, REML=FALSE)
lmm3 = lmer(fst ~ geo_dist + env_dist + biome + (1 | OTU), data=div3, REML=FALSE)
lmm4 = lmer(fst ~ geo_dist + type + biome + (1 | OTU), data=div3, REML=FALSE)
lmm5 = lmer(fst ~ geo_dist + env_dist + (1 | OTU), data=div3, REML=FALSE)
lmm6 = lmer(fst ~ geo_dist + biome + (1 | OTU), data=div3, REML=FALSE)
lmm7 = lmer(fst ~ geo_dist + type + (1 | OTU), data=div3, REML=FALSE)
lmm8 = lmer(fst ~ geo_dist + (1 | OTU), data=div3, REML=FALSE)
lmm9 = lmer(fst ~ env_dist + (1 | OTU), data=div3, REML=FALSE)
lmm10 = lmer(fst ~ type + (1 | OTU), data=div3, REML=FALSE)
lmm11 = lmer(fst ~ biome + (1 | OTU), data=div3, REML=FALSE)

# plot these results
# should plot these holding one constant
theme_set(theme_bw(base_size = 30))
a = ggplot(div, aes(x=env_dist, y=fst)) + geom_point(alpha=0.3, col="gray") + stat_smooth(method="lm", col="darkblue") + ylim(0,1) + theme_bw() + xlab("rescaled environmental distance") +ylab(expression('F'[ST]))
b = ggplot(div, aes(x=geo_dist, y=fst)) + geom_point(alpha=0.3, col="gray") + stat_smooth(method="lm", col="darkblue") + ylim(0,1) + theme_bw() + xlab("rescaled geographic distance") + ylab(expression('F'[ST]))


keep = intersect(t2$tip.label, unique(div3$OTU))
t3 = drop.tip(t2, setdiff(t2$tip.label, keep))
div4 = div3[div3$OTU %in% keep, ]
div4$type = as.factor(div4$type)
div4$biome = as.factor(div4$biome)

# get inverse phylogenetic matrix
inv.phylo = inverseA(t3, nodes="TIPS", scale=TRUE)

# priors
# these are fairly canonical
prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
iterations = 1e5
burnins = 1e4
m1 = MCMCglmm(fst~geo_dist,random=~OTU,
                ginverse=list(OTU=inv.phylo$Ainv),
                prior=prior,
                data=div4, nitt=iterations, burnin=burnins)
m2 = MCMCglmm(fst~env_dist, random=~OTU,
              ginverse=list(OTU=inv.phylo$Ainv),
              prior=prior,
              data=div4, nitt=iterations, burnin=burnins)
m3 = MCMCglmm(fst~type, random=~OTU,
              ginverse=list(OTU=inv.phylo$Ainv),
              prior=prior,
              data=div4, nitt=iterations, burnin=burnins)
m4 = MCMCglmm(fst~biome, random=~OTU,
              ginverse=list(OTU=inv.phylo$Ainv),
              prior=prior,
              data=div4, nitt=iterations, burnin=burnins)
m5 = MCMCglmm(fst~env_dist + geo_dist, random=~OTU,
              ginverse=list(OTU=inv.phylo$Ainv),
              prior=prior,
              data=div4, nitt=iterations, burnin=burnins)
m6 = MCMCglmm(fst~env_dist + geo_dist + biome, random=~OTU,
              ginverse=list(OTU=inv.phylo$Ainv),
              prior=prior,
              data=div4, nitt=iterations, burnin=burnins)
m7 = MCMCglmm(fst~env_dist + geo_dist + type, random=~OTU,
              ginverse=list(OTU=inv.phylo$Ainv),
              prior=prior,
              data=div4, nitt=iterations, burnin=burnins)
m8 = MCMCglmm(fst~env_dist + geo_dist + type + biome, random=~OTU,
              ginverse=list(OTU=inv.phylo$Ainv),
              prior=prior,
              data=div4, nitt=iterations, burnin=burnins)
models = list(m1, m2, m3, m4, m5, m6, m7, m8)
saveRDS(models, "~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/pairwise_Fst_models.rds")

#################################
# for best fitting model, check for convergence
##################################

m8a = MCMCglmm(fst~env_dist + geo_dist + type + biome, random=~OTU,
              ginverse=list(OTU=inv.phylo$Ainv),
              prior=prior,
              data=div4, nitt=iterations, burnin=burnins)
m8b = MCMCglmm(fst~env_dist + geo_dist + type + biome, random=~OTU,
               ginverse=list(OTU=inv.phylo$Ainv),
               prior=prior,
               data=div4, nitt=iterations, burnin=burnins)
m8c = MCMCglmm(fst~env_dist + geo_dist + type + biome, random=~OTU,
               ginverse=list(OTU=inv.phylo$Ainv),
               prior=prior,
               data=div4, nitt=iterations, burnin=burnins)

#################################
# check convergences
##################################

effectiveSize(m8$Sol)
effectiveSize(m8$VCV)
autocorr.diag(m8$Sol)
mm = mcmc.list(m8a$Sol, m8b$Sol)
gelman.plot(mm)
gelman.diag(mm)

lambda <- m8$VCV[,'OTU']/
  (m8$VCV[,'OTU']+m8$VCV[,'units'])

#################################
# summarize the results across a series of factors
##################################
div4$biome = as.character(div4$biome)
div4$type = as.character(div4$type)
div4[div4$biome == 'temperature grasslands', "biome"] = 'temperate grasslands'
div4[div4$biome == 'DIFF', "biome"] = 'different biomes'
div5 = div4[div4$biome != 'temperate grasslands', ]

means = aggregate(div5$fst, by=list(div5$biome), mean)
div5$biome <- factor(div5$biome, levels = means[order(means$x), "Group.1"], ordered = TRUE)
d = ggplot(div5, aes(biome, fst)) + geom_boxplot() + coord_flip() + ylim(0, 1) + labs(y=expression('F'[ST]))

div5[div5$type == 'CEN-CEN', "type"] = 'center-center'
div5[div5$type == 'CEN-PER', "type"] = 'center-edge'
div5[div5$type == 'PER-PER', "type"] = 'edge-edge'
means = aggregate(div5$fst, by=list(div5$type), mean)
div5$type <- factor(div5$type, levels = means[order(means$x), "Group.1"], ordered = TRUE)
c = ggplot(div5, aes(type, fst)) + geom_boxplot() + coord_flip() + ylim(0, 1) + labs(y=expression('F'[ST]))

a = ggplot(div5, aes(geo_dist, fst)) + geom_point(alpha = 0.1) + ylim(0, 1) + labs(y=expression('F'[ST]), x="scaled geographic distance")
b = ggplot(div5, aes(env_dist, fst)) + geom_point(alpha = 0.1) + ylim(0, 1) + labs(y=expression('F'[ST]), x="scaled environmental distance")
abc = plot_grid(a, b, c, d, labels = c("A", "B", "C", "D"), hjust=-1, nrow=2)
save_plot("~/Dropbox/Sphenomorphine_Gene_Flow/figures/pairwise_Fst_comparisons.pdf", abc, ncol = 2, nrow = 2, base_aspect_ratio = 1.3, base_height=3)
