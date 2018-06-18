library(sp)
library(maptools)
library(rgeos)
library(raster)
library(geosphere)
library(scales)
data(wrld_simpl)

################################
# get biomes per range
################################

rangedir = '/Users/sonal/Dropbox/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/'
d = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.csv", stringsAsFactors=F)
otus = unique(d$OTU)

# the original biomes downloaded from the TNC
# http://maps.tnc.org/gis_data.html
biomes = readShapePoly("~/macroevolution/eco_IBD_oz/data/ecoregions/terr-ecoregions-aus_combined/terr_ecoregions-aus_combined.shp")
# simplify this so that it takes less time to process
biomes1 = gSimplify(biomes, 0.1, topologyPreserve=T)
biomes1$biome = as.character(biomes$WWF_MHTNAM)

# rename the biomes
bnames1 = c("Tropical and Subtropical Moist Broadleaf Forests", "Montane Grasslands and Shrublands", "Temperate Broadleaf and Mixed Forests", "Tropical and Subtropical Grasslands, Savannas and Shrublands", "Deserts and Xeric Shrublands", "Mediterranean Forests, Woodlands and Scrub", "Temperate Grasslands, Savannas and Shrublands")
bnames2 = c("tropical forests", "montane", "temperate forests", "tropical grasslands", "desert", "Mediterranean forests", "temperature grasslands")
for (i in 1:length(bnames1)) {
   biomes1$biome = gsub(bnames1[i], bnames2[i],  biomes1$biome)
 }

biomes2 = vector('list', length(biomes1))
names(biomes2) = biomes1$biome
# helps deal with issues around illegal spatial geometries
for (i in 1:length(biomes1)) {
  a = biomes1[i, ]
  biomes2[[i]] = gBuffer(SpatialPolygons(a@polygons,proj4string=a@proj4string), width=0)
}

# identify the biomes in which
# more than (min) of the area overlaps with the biome
classify <- function(overlaps, minOverlap) {
  overlaps = overlaps / sum(overlaps)
  hab = overlaps[ overlaps > minOverlap ]
  # if nothing is significant overlapping, then call it a generalist
  if (length(hab) < 1) {
    return("generalist")
  } else {
    return(paste(sort(names(hab)), collapse=","))
  }
}

res = data.frame(OTU=otus, biome=rep(NA, length(otus)), stringsAsFactors = F)
for (i in 1:length(otus)) {
  # get the range
  range = paste(rangedir, otus[i], ".shp", sep="")
  if (!file.exists(range)) {
    if (otus[i] == 'Notoscincus_wotjulum') {
      otu = 'Notoscincus_ornatus'
    } else {
      otu = gsub('_\\d', '', otus[i])
    }
    range = paste(rangedir, otu, '.shp', sep='')
    cat(range, "\n")
  }
  
  range = readShapePoly(range)
  if (class(range) == 'SpatialPolygonsDataFrame') {
    range = SpatialPolygons(range@polygons,proj4string=range@proj4string)
  }
  
  # calculate all the overlaps in km2
  overlaps = rep(0, length(biomes2))
  names(overlaps) = names(biomes2)
  for (x in 1:length(biomes2)) {
    overlap = gIntersection(gBuffer(range, width=0), biomes2[[x]])
    if (class(overlap) != 'NULL') {
      overlaps[x] = sum(areaPolygon(overlap)) / 1e6
    }
  }	
  
  res[i, "biome"] = classify(overlaps, 0.2)
  cat(i, "\n")
}
write.csv(res, "~/Dropbox/Sphenomorphine_Gene_Flow/data/geographic_ranges/otu_biomes.csv", row.names=F)

################################
# check if biomes explain PC2 & lat midpoint patterns
################################

diff = readRDS("~/Dropbox/Sphenomorphine_Gene_Flow/data/summary.Rds")
dd = merge(diff, res, by="OTU")

bio_counts = table(dd$biome)
dd1 = dd[dd$biome %in% names(bio_counts[ bio_counts >= 3]), ]

aggregate(dd1$lat_midpoint, by=list(dd1$biome), mean, na.rm=T)
aggregate(dd1$PC2_range, by=list(dd1$biome), mean, na.rm=T)
aggregate(log(dd1$slope), by=list(dd1$biome), mean, na.rm=T)

boxplot(log(dd1$slope) ~ dd1$biome)

################################
# confirm that this is driven by the eastern species
################################

aus <- wrld_simpl[which(wrld_simpl@data$NAME == 'Australia'),]
plot(NULL, xlim=c(112, 155), ylim=c(-44, -9), xlab="", ylab="", axes=F)
plot(aus, add=T, col=alpha("gray", 0.5), border=NA)
check = c("temperate forests,tropical grasslands", "temperate forests", "tropical grasslands")
otus2 = dd1[dd1$biome %in% check, 'OTU']
for (i in 1:length(otus2)) {
  # get the range
  range = paste(rangedir, otus2[i], ".shp", sep="")
  if (!file.exists(range)) {
    if (otus[i] == 'Notoscincus_wotjulum') {
      otu = 'Notoscincus_ornatus'
    } else {
      otu = gsub('_\\d', '', otus[i])
    }
    range = paste(rangedir, otu, '.shp', sep='')
  }
  range = readShapePoly(range)
  plot(range, add=T)
}

################################
# do the analysis by species restricting to just one
# major biome at a time
################################

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

bio_ll = extract(biomes1, ll[, 2:3])
bio_ll = bio_ll[!duplicated(bio_ll$point.ID), ]
bio_ll2 = cbind(ll, bio_ll)

plot(bio_ll2$lon, bio_ll2$lat, col=as.factor(bio_ll2$biome), pch=16)

inv_fst = function(x) {
  if (x == 1) {
    return(NA)
  } else {
    return(x / (1 - x))
  }
}

log_dist = function(x) {
  if (x == 0) {
    return(NA)
  } else {
    return(log(x))
  }
}

biome_winners = c('desert', "Mediterranean forests", "tropical grasslands")
divfiles = list.files("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/species/", full.names = T)

nrow = length(biome_winners) * length(divfiles)
bslopes1 = c("OTU", "biome", "ninds", "slope")
bslopes = data.frame(matrix(NA, nrow, length(bslopes1)), stringsAsFactors = F)
names(bslopes) = bslopes1

count = 1
for (j in 1:length(biome_winners)) {
  biome_winner = biome_winners[j]
  for (i in 1:length(divfiles)) {
    div = read.csv(divfiles[i], stringsAsFactors = F)
    cl = gsub("/Users/Sonal/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/species//", "", divfiles[i])
    cl = gsub(".divergence_cov10.csv", "", cl)
    
    bslopes[count, "OTU"] = cl
    bslopes[count, "biome"] = biome_winner
    
    div = div[complete.cases(div$geo_dist),]
    div = div[div$fst_denom >= 1000, ]
    
    div$log_dist = sapply(div$geo_dist, log_dist)
    div$inv_fst = sapply(div$fst, inv_fst)
    div = div[complete.cases(div$inv_fst), ]
    div = div[complete.cases(div$log_dist), ]
    
    inds = unique(c(div$ind1, div$ind2))
    inds = bio_ll2[inds, ]
    
    keep = inds[inds$biome %in% biome_winner, 'sample_id']
    bslopes[count, "ninds"] = length(keep)
    div2 = div[div$ind1 %in% keep, ]
    div2 = div2[div2$ind2 %in% keep, ]
  
    if (nrow(div2) >= 3) {
      a = lm(div2$inv_fst ~ div2$log_dist)
      bslopes[count, "slope"] = coef(a)[2]
    }
    count = count + 1
  }
}

write.csv(bslopes, "~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.biomes.csv", row.names=F)
