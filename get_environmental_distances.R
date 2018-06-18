library(pcaMethods)
library(raster)

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

# get raster environmental data
rasters = list.files('/Users/sonal/macroevolution/eco_IBD_oz/data/AUS_5arc/', pattern=".tif", full.names=T)
rasters = rasters[-grep("slope", rasters)]
rasters = rasters[-grep("aspect", rasters)]
predictors = stack(rasters)

# get environmental divergences
# do this in one PCA so comparable across species?
vals = extract(predictors, ll[,c("lon", "lat")])
vals = cbind(ll, vals)

# because missing data, do PCA with imputation
# pca = pca(vals[,4:ncol(vals)], scale="vector", nPcs=9, method="ppca")
# only keeping the first three axes because they explain >84% and each subsequent axis only adds ~2 - 3%
# env_dist = as.matrix(dist(pca@scores[,1:3], diag=T, upper=T, method="euclidean"))
env_dist = as.matrix(dist(vals, diag=T, upper=T, method="euclidean"))
dimnames(env_dist) = list(ll$sample_id, ll$sample_id)

saveRDS(env_dist, file="/Users/sonal/Dropbox/Sphenomorphine_Gene_Flow/data/envirospatial/envdist.rds")
