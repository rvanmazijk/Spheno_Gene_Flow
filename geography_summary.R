library(maptools)
library(raster)
library(rgeos)
library(sp)
library(geosphere)

# geography
## lat_midpoint - done
## range_size - done
## elev_range - done
## PC1_range, PC2_range - done

## avg_suitability - revisit 

get_heterogeneity <- function() {
  # australia coast
  aus = readShapePoly("/Users/Sonal/scripts/eco_IBD/pascal_scripts/AusCoast.shp", proj4string=CRS("+proj=longlat +datum=WGS84"))
  pts = spsample(aus, 50000, type="regular")
  
  #  all rasters
  r_files = list.files('/Users/Sonal/macroevolution/eco_IBD_oz/data/AUS_5arc/', pattern=".tif", full.names=T)
  r_files = r_files[grep("alt", r_files, invert=T)]
  r_files = r_files[grep("slope", r_files, invert=T)]
  r_files = r_files[grep("aspect", r_files, invert=T)]
  
  rasters = stack(r_files)
  # get env data at sampled points
  vals = extract(rasters, pts)
  
  # combine with randomly sampled points
  vals = cbind(vals, data.frame(pts))
  # get rid of NA
  vals2 = vals[complete.cases(vals),]
  
  # do the pca
  end = length(r_files)
  valspca = prcomp(vals2[,1:end], scale=T, center=T)
  
  # turn the first pcas (0.41, 0.78, 0.85) into data frame
  
  pts1 = cbind(vals2[, (end+1):(end+2)], valspca$x[,1])
  names(pts1) = c("x", "y", "PC1")
  pts2 = cbind(vals2[, (end+1):(end+2)], valspca$x[,2])
  names(pts2) = c("x", "y", "PC2")
  pts3 = cbind(vals2[, (end+1):(end+2)], valspca$x[,3])
  names(pts3) = c("x", "y", "PC3")

  pts1 = SpatialPoints(pts1, proj4string=CRS("+proj=longlat +datum=WGS84"))
  pts2 = SpatialPoints(pts2, proj4string=CRS("+proj=longlat +datum=WGS84"))
  pts3 = SpatialPoints(pts3, proj4string=CRS("+proj=longlat +datum=WGS84"))
  pts = list(pts1, pts2, pts3)
  names(pts) = c("PC1", "PC2", "PC3")
  return(pts)
}

d = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.csv", stringsAsFactors=F)
rangedir = '~/Dropbox/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/'
otus = unique(d$OTU)

elev = raster("~/macroevolution/eco_IBD_oz/data/AUS_5arc/aus5min_alt.tif")
pts = get_heterogeneity()

cols = c("OTU", "range_size", "lat_midpoint", "elev_range", "PC1_range", "PC2_range", "PC3_range", "PC1_sd", "PC2_sd", "PC3_sd")
res = data.frame(matrix(NA, nrow=length(otus), ncol=length(cols)))
names(res) = cols

for (i in 1:length(otus)) {
  cat(i, "\n")
  shp = paste(rangedir, otus[i], '.shp', sep='')
  if (!file.exists(shp)) {
    if (otus[i] == 'Notoscincus_wotjulum') {
      otu = 'Notoscincus_ornatus'
    } else {
      otu = gsub('_\\d', '', otus[i])
    }
    shp = paste(rangedir, otu, '.shp', sep='')
    cat(shp, "\n")
  }
  
  res[i, "OTU"] = otus[i]
  range = readShapePoly(fn = shp, proj4string=CRS('+proj=longlat +datum=WGS84'))
  range = spTransform(range, CRS("+proj=longlat +datum=WGS84")) 
  
  # need to sum because some ranges consist of multiple ranges
  # divide by 1e6 to convert from m^2 to km^2
  res[i, "range_size"] = sum(areaPolygon(range)) / 1e6
  
  # get the latitude of the range
  res[i, "lat_midpoint"] = gCentroid(range)@coords[2]
  
  # get elevation range
  elevrange = unlist(extract(elev, range))
  res[i, "elev_range"] = max(elevrange, na.rm=T) - min(elevrange, na.rm=T)
  
  # get heterogeneity values
  for (x in 1:length(pts)) {
    hetrange =  data.frame(pts[[x]][range, ])[, 3]
    res[i, paste(names(pts)[x], "_sd", sep="")] = sd(hetrange, na.rm=T)
    res[i, paste(names(pts)[x], "_range", sep="")] = max(hetrange, na.rm=T) - min(hetrange, na.rm=T)
  }
}
write.csv(res, "~/Dropbox/Sphenomorphine_Gene_Flow/data/envirospatial/range_data.csv", row.names=F)
