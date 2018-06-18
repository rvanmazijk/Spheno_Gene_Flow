library(maptools)
data(wrld_simpl)

pdf("~/Desktop/Ctenotus_robustus_3_map.pdf")
aus <- wrld_simpl[which(wrld_simpl@data$NAME == 'Australia'),]
plot(NULL, xlim=c(112, 155), ylim=c(-44, -9), xlab="", ylab="", axes=F)
plot(aus, add=T, col=alpha("gray", 0.5), border=NA)
r = readShapePoly("~/Dropbox/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/Ctenotus_robustus_3.shp")
plot(r, add=T, col='#1b9e77', border=F)
dev.off()

ll = read.csv("/Users/sonal/macroevolution/eco_IBD_oz/data/metadata/individual_data_nomissing22Oct15.csv", stringsAsFactors=F, na.string=c("", "NA"))
d = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/species/Ctenotus_robustus_3.divergence_cov10.csv", stringsAsFactors=F)
inds = unique(c(d$ind1, d$ind2))
pts = ll[ll$sample_id %in% inds,]
pts = pts[pts$lon > 136,]

pdf("~/Desktop/Ctenotus_robustus_3_map_pts.pdf")
plot(NULL, xlim=c(112, 155), ylim=c(-44, -9), xlab="", ylab="", axes=F)
plot(aus, add=T, col=alpha("gray", 0.5), border=NA)
plot(r, add=T, col='#1b9e77', border=F)
points(pts$lon, pts$lat, pch=21, bg="white")
dev.off()

inv_fst = function(x) {
  if (x == 1) {
    return(NA)
  } else {
    return(x / (1 - x))
  }
}

d$inv_fst = mapply(inv_fst, d$fst)

pdf("~/Desktop/dist_Fst.pdf", width=4, height=4)
par(mar=c(4,4, 1, 1))
plot(NULL, axes=FALSE, xlab="", ylab="", xlim=range(log(d$geo_dist)), ylim=range(d$inv_fst), xaxs="i", yaxs="i")
points(log(d$geo_dist), d$inv_fst, pch=21, bg="#1b9e77")

xvals = c(10, 11, 12, 13, 14, 15)
xvals1 = xvals[xvals >= min(log(d$geo_dist))]
xvals1 = xvals1[xvals1 <= max(log(d$geo_dist))]
axis(1, at=xvals, labels=FALSE, tck=-0.015)
mtext(xvals1, side=1, line=0.5, at=xvals1)
mtext("log(geographic distance, m)", side=1, line=2.4)

yvals = c(0, 1, 2, 3, 4)
yvals1 = yvals[yvals >= min(d$inv_fst)]
yvals1 = yvals1[yvals1 <= max(d$inv_fst)]
axis(2, at=yvals, labels=FALSE, tck=-0.015, las=1)
mtext(yvals1, side=2, line=0.5, at=yvals1, las=1)
mtext(expression('inverse F'[ST]), side=2, line=2.4)
dev.off()


pdf("~/Desktop/dist_Fst_fit.pdf", width=4, height=4)
par(mar=c(4, 4, 1, 1))
plot(NULL, axes=FALSE, xlab="", ylab="", xlim=range(log(d$geo_dist)), ylim=range(d$inv_fst), xaxs="i", yaxs="i")
points(log(d$geo_dist), d$inv_fst, pch=16, col=alpha("#1b9e77", 0.5))

m = lm(d$inv_fst ~ log(d$geo_dist))
fst <- function(x) {
  return(m$coefficients[2] * x + m$coefficients[1])
}
lines(log(d$geo_dist), mapply(fst, log(d$geo_dist)))

xvals = c(10, 11, 12, 13, 14, 15)
xvals1 = xvals[xvals >= min(log(d$geo_dist))]
xvals1 = xvals1[xvals1 <= max(log(d$geo_dist))]
axis(1, at=xvals, labels=FALSE, tck=-0.015)
mtext(xvals1, side=1, line=0.5, at=xvals1)
mtext("log(geographic distance, m)", side=1, line=2.4)

yvals = c(0, 1, 2, 3, 4)
yvals1 = yvals[yvals >= min(d$inv_fst)]
yvals1 = yvals1[yvals1 <= max(d$inv_fst)]
axis(2, at=yvals, labels=FALSE, tck=-0.015, las=1)
mtext(yvals1, side=2, line=0.5, at=yvals1, las=1)
mtext(expression('inverse F'[ST]), side=2, line=2.4)
dev.off()
