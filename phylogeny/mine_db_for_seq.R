# y is skinkind database
# x is cytb

xx = merge(x, y, by="SkinkRef", all=T)
c = read.csv("~/macroevolution/eco_IBD_oz/data/clustering/Ctenotus_Lerista.clustering.revised2.csv", stringsAsFactors=F)
c$sp = gsub("C. ", "Ctenotus_", c$LatinName)
c$sp = gsub("L. ", "Lerista_", c$sp)
c$sp = gsub(" ", "_", c$sp)
c = c[complete.cases(c$sp),]
d = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv", stringsAsFactors=F)

xx$sp = paste(xx$GENUS, xx$SPECIES, sep="_")
xx$sp = gsub("\\s+", "", xx$sp)

d = d[which(d$check_genbank == FALSE),]

# check reseq
for (i in 1:nrow(d)) {
	otu = d[i, "OTU"]
	sp = d[i, "SPECIES"]
	if (otu == sp) {
		tmp = xx[xx$sp == otu,]
		tmp = tmp[complete.cases(tmp$SEQ),]
		if (nrow(tmp) > 1) {
			# tmp$seqlen = str_count(tmp$SEQ, "[A|T|C|G]")
			# tmp = tmp[which(tmp$seqlen == max(tmp$seqlen)), ]
			cat(otu, nrow(tmp), "\n", sep="\t")
			}
	} else {
		samps = c[c$sp == otu, "sample"]
		tmp = xx[xx$SAMPLE_ID.y %in% samps, ]
		tmp = tmp[complete.cases(tmp$SEQ),]
		if (nrow(tmp) > 1) {
			tmp$seqlen = str_count(tmp$SEQ, "[A|T|C|G]")
			tmp = tmp[which(tmp$seqlen == max(tmp$seqlen)), ]
			cat(otu, nrow(tmp), "\n", sep="\t")
			}
	}	
}

otu = 'Lerista_muelleri_3'
samps = c[c$sp == otu, "sample"]
tmp = xx[xx$SAMPLE_ID.y %in% samps, ]
tmp = tmp[complete.cases(tmp$SEQ),]
tmp$seqlen = str_count(tmp$SEQ, "[A|T|C|G]")
tmp = tmp[which(tmp$seqlen == max(tmp$seqlen)), ]
tmp$SAMPLE_ID.y



d = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv", stringsAsFactors=F)
dd = d[which(d$check_internal_db == TRUE),]

sink("~/Desktop/cytb_database.fasta")
for (i in 1:nrow(dd)) {
	sample = dd[i, "sanger_sample"]
	cat(paste(">", sample, sep=""), "\n")
	xx = x[x$SAMPLE_ID == sample, ]
	if (exists(xx$RESEQ)) {
		seq = xx$RESEQ
	} else {
		seq = xx$SEQ
	}
	cat(seq, "\n")
}
sink()
