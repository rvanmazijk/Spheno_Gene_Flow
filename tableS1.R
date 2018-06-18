library(rangeBuilder)

diff = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/divergence.all_OTUs.csv", stringsAsFactors = F)
otus = unique(diff$OTU)

inds = c()
sites1 = c()
sites2 = c()
for (i in 1:length(otus)) {
  file = paste("~/Dropbox/Sphenomorphine_Gene_Flow/data/divergence/species/", otus[i], ".divergence_cov10.csv", sep="")
  d = read.csv(file, stringsAsFactors = F)
  inds = c(inds, unique(c(d$ind1, d$ind2)))
  sites1 = c(sites1, d$nuc_denom)
  sites2 = c(sites2, d$fst_denom)
}

d = read.csv("~/Desktop/activeWork/databases/new_individuals.csv", stringsAsFactors = F)
x = data.frame(inds, SAMPLE_ID=rep(NA, length(inds)), stringsAsFactors = F)
x$SAMPLE_ID = d[match(x$inds, d$SAMPLE_ID), "SAMPLE_ID"]

for (i in 1:nrow(x)) {
  if (is.na(x[i, "SAMPLE_ID"])) {
    id = x[i, "inds"]
    id = gsub("[A-Z]+", "", id)
    id = gsub("_", "", id)
    xx = d[grep(id, d$SAMPLE_ID), ]
    xx = xx[xx$REGION == "Australia", ]
    if (nrow(xx) > 1) {
      if (id == 'SAMAR46891') {
        newid = 'SAMAR_46891_Er_fasc'
      } else if (id == '60807') {
        newid = 'SAMAR_60807_Er_fasc'
      } else if (id == '21409') {
        newid = 'SAMAR_21409_Eu_tymp'
      } else if (id == '36665') {
        newid = 'SAMAR_36665_He_pero'
      }
    } else {
      newid = xx$SAMPLE_ID 
    }
  x[i, "SAMPLE_ID"] = newid
  }
}

x = merge(x, d, by="SAMPLE_ID")

sps = paste(x$GENUS, "_", x$SPECIES, sep="")
sps = table(sps)
sps1 = synonymMatch(sps, "squamates")

d1 = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv", stringsAsFactors = F)
d1 = d1[d1$outgroup == F,]
length(unique(d1$SPECIES))
length(unique(gsub("_\\S+", "", unique(d1$SPECIES))))

countries = getRepDBcountryList(unique(d1$SPECIES))
for (i in 1:length(countries)) {
  if (!"AUSTRALIA" %in% countries[[i]]) {
    cat(names(countries)[i], "\n")
    cat(countries[[i]], "\n")
    cat("***\n")
  }
}

l = x[x$GENUS == 'Lerista', c("SPECIES", "OTU")]
ll = l[grep("_\\d", l$OTU), ]
ll[order(ll$SPECIES), ]

o = x[x$GENUS != 'Lerista', ]
o = o[o$GENUS != 'Ctenotus', ]
o = o[, c("SPECIES", "OTU")]
o = o[grep("_\\d", o$OTU), ]
o[order(o$SPECIES), ]


sps2 = names(sps[sps > 2])
sps2 = sps2[grep("Ctenotus", sps2, invert=T)]
sps2 = sps2[grep("Lerista", sps2, invert=T)]


x = x[, c("inds", "SAMPLE_ID", "GENUS", "SPECIES", "LAT", "LON", "OTU")]
r = read.csv("~/Desktop/read_counts.csv", stringsAsFactors = F)
r$file = gsub("/nfs/turbo/lsa-rabosky/Lab/skink_ddrad/individual_reads/", '', r$file)
r$file = gsub(".1.noRE.fq.gz", '', r$file)
r$file = gsub("_R1.fastq.gz", '', r$file)
x$reads = r[match(x$inds, r$file), "number_seqs"]
x = x[order(x$OTU), ]
x$SAMPLE_ID = gsub("_Ct_.*$", "", x$SAMPLE_ID)
x$SAMPLE_ID = gsub("_Le_.*$", "", x$SAMPLE_ID)
x$SAMPLE_ID = gsub("_No_.*$", "", x$SAMPLE_ID)
x$SAMPLE_ID = gsub("_Sa_.*$", "", x$SAMPLE_ID)
x$SAMPLE_ID = gsub("_Eu_.*$", "", x$SAMPLE_ID)
x$SAMPLE_ID = gsub("_Er_.*$", "", x$SAMPLE_ID)
x$SAMPLE_ID = gsub("_An_.*$", "", x$SAMPLE_ID)
x$SAMPLE_ID = gsub("_Ca_.*$", "", x$SAMPLE_ID)
write.csv(x, "~/Desktop/S1_individuals.csv", row.names=F)

x = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/manuscript_drafts/tables/S1_individuals.csv", stringsAsFactors = F)
d = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv", stringsAsFactors = F)
d = d[d$outgroup == F,]
allsps = unique(d$SPECIES)

ctsps = allsps[grep("Ct", allsps)]
ctsps1 = sps1[grep("Ct", sps1)]
length(ctsps1) / length(ctsps)

lesps = allsps[grep("Le", allsps)]
lesps1 = sps1[grep("Le", sps1)]
length(lesps1) / length(lesps)

osps = allsps[grep("Ctenotus", allsps, invert=T)]
osps = osps[grep("Lerista", osps, invert=T)]
osps1 = sps1[grep("Le", sps1, invert=T)]
osps1 = osps1[grep("Ct", osps1, invert=T)]
length(osps1) / length(osps)
