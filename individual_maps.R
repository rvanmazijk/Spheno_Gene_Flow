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

divfiles = list.files("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/species/", full.names = T)
pdf("~/Dropbox/Sphenomorphine_Gene_Flow/figures/maps.pdf", height=3, width=10)
boots = vector('list', length(divfiles))
for (i in 1:length(boots)) {
  boots[[i]] = vector('list', 2)
}
# for (i in 25:27) {
for (i in 1:length(divfiles)) {
  cat(i, "\n")
  d = read.csv(divfiles[i], stringsAsFactors=F)
  
  otu = gsub('^.*/', '', divfiles[i])
  otu = gsub('.diver.*', '', otu)
  names(boots)[i] = otu
  inds = unique(c(d$ind1, d$ind2))
  
  if (length(inds) >= 3) {
    d = d[d$fst_denom >= 1000,]
    d$inv_fst = mapply(inv_fst, d$fst)
    d$dist_log = mapply(log_dist, d$geo_dist)
    d = d[complete.cases(d$inv_fst), ]
    d = d[complete.cases(d$dist_log), ]
    inds = unique(c(d$ind1, d$ind2))
    
    if (nrow(d) >= 3) {  
      par(mfrow=c(1, 4))
      # plot map with points
      aus <- wrld_simpl[which(wrld_simpl@data$NAME == 'Australia'),]
      plot(NULL, xlim=c(112, 155), ylim=c(-44, -9), xlab="", ylab="", axes=F, main=otu)
      plot(aus, add=T, col=alpha("gray", 0.5), border=NA)
      rangefile = paste("~/Dropbox/Sphenomorphine_Gene_Flow/data/geographic_ranges/OTU_ranges/", otu, ".shp", sep="")
      if (file.exists(rangefile)) {
        r = readShapePoly(rangefile)
        plot(r, add=T, col='#1b9e77', border=F)
      }
      pts = ll[ll$sample_id %in% inds,]
      points(pts$lon, pts$lat, pch=21, bg="white")
      
      # plot ln dist and inverse fst
      par(mar=c(4, 4, 1, 1))
      plot(d$dist_log, d$inv_fst, pch=16, col=alpha("#1b9e77", 0.5), las=1, xlab="log(geographic distance, m)", ylab=expression('inverse F'[ST]))
      
      m = lm(d$inv_fst ~ d$dist_log)
      fst <- function(x) {
        return(m$coefficients[2] * x + m$coefficients[1])
      }
      lines(d$dist_log, mapply(fst, d$dist_log))
      est_slope = m$coefficients[2]
      
      # plot histogram showing uncertainty in fst
      slopes = rep(NA, 100)
      for (x in 1:100) {
        fsts = mapply(rnorm, 1, d$fst_mean, d$fst_sd)
        inv_fsts = mapply(inv_fst, fsts)
        m = lm(inv_fsts ~ d$dist_log)
        slopes[x] = m$coefficients[2]
      }
      hist(slopes, col="gray", border=F, main="Slope error due to Fst estimation", las=1)
      abline(v=est_slope, col="red", lty=2)
      boots[[i]][[1]] = slopes
      
      # plot histogram with subsampling
      slopes = rep(NA, 1000)
      if (length(inds) > 5) {
        nsamp = 6
        trynum = 1
        while (trynum < 1000) {
          trynum = trynum + 1
          inds1 = sample(inds, nsamp)
          t = d[d$ind1 %in% inds1,]
          t = t[t$ind2 %in% inds1,]
          if (nrow(t) > 14) {
            m = lm(t$inv_fst ~ t$dist_log)
            if (nrow(summary(m)$coefficients) > 1) {
              if (summary(m)$coefficients[2, 4] <= 0.05) {
                slopes[trynum] = m$coefficients[2]
              }
            }  
          }
        }
      }  
      slopes = slopes[complete.cases(slopes)]
      if (length(slopes) > 10) {
        hist(slopes, col="gray", border=F, main="bootstrapping Fst slope", las=1)
        abline(v=est_slope, col="red", lty=2)
        boots[[i]][[2]] = slopes
      } else {
        plot(0,type='n',axes=FALSE,ann=FALSE)
      }
    }  
  } 
}

dev.off()

cols = c("OTU", "rep", "slope")
bootloci = data.frame(matrix(NA, nrow=length(boots) * 100, ncol=length(cols)))
names(bootloci) = cols
for (i in 1:length(boots)) {
  bootloci[(1 + (i - 1) * 100):(100 + (i - 1) * 100), "OTU"] = names(boots)[i]
  bootloci[(1 + (i - 1) * 100):(100 + (i - 1) * 100), "rep"] = seq(1, 100)
  slopes = boots[[i]][[1]]
  if (length(slopes) > 0) {
    bootloci[(1 + (i - 1) * 100):(100 + (i - 1) * 100), "slope"] = slopes
  } else {
    bootloci[(1 + (i - 1) * 100):(100 + (i - 1) * 100), "slope"] = NA
  }
}
write.csv(bootloci, "~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.bootloci.csv", row.names = F)

cols = c("OTU", "rep", "slope")
bootinds = data.frame(matrix(NA, nrow=length(boots) * 100, ncol=length(cols)))
names(bootinds) = cols
for (i in 1:length(boots)) {
  bootinds[(1 + (i - 1) * 100):(100 + (i - 1) * 100), "OTU"] = names(boots)[i]
  bootinds[(1 + (i - 1) * 100):(100 + (i - 1) * 100), "rep"] = seq(1, 100)
  slopes = boots[[i]][[2]]
  if (length(slopes) > 0) {
    bootinds[(1 + (i - 1) * 100):(100 + (i - 1) * 100), "slope"] = slopes[1:100]
  } else {
    bootinds[(1 + (i - 1) * 100):(100 + (i - 1) * 100), "slope"] = NA
  }
}
write.csv(bootinds, "~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.bootinds.csv", row.names = F)
