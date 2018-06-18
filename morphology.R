library(ape)
library(phytools)

# read in list of species
d = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv", stringsAsFactors=F)
# read in dan's sys bio measurements
m = read.csv("~/macroevolution/eco_IBD_oz/data/morphology/skink_indmeans_june7.csv", stringsAsFactors=F)

# get phylogeny
t = read.nexus("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/beast_ucln_31July17/mtcds_pos1.sub.tre")
# run this if tree isn't ultrametric -- happens because of numerical instability
# t = nnls.tree(cophenetic(t), t,rooted=TRUE)

# get rid of outgroups
outs = d[d$outgroup == TRUE, "sanger_sample"]
outs = outs[complete.cases(outs)]
t = drop.tip(t, outs)

# make species level tree
keep = d[d$SPECIES_KEEP == TRUE, "sanger_sample"]
t1 = drop.tip(t, setdiff(t$tip.label, keep))
t1$tip.label = d[match(t1$tip.label, d$sanger_sample), "SPECIES"]

# change caps on species names so they match
d$sp1 = tolower(d$SPECIES)
d$sp2 = tolower(d$OTHER_NAMES)

# match species names to species names
# note that this ends up assigning all OTUs to the measurements for their corresponding nominal species
m$SPECIES = rep(NA, nrow(m))
for (i in 1:nrow(m)) {
	sp_m = m[i, "treename"]
	if (sp_m %in% d$sp1) {
		sp_t = d[d$sp1 == sp_m, "SPECIES"][1]
	} else {
		t = d[grep(sp_m, d$sp2),]
		sps = unique(t$SPECIES)
		if (length(sps) == 1) {
			sp_t = sps
		} else {
			sp_t = NA
		}
	}
	m[i, "SPECIES"] = sp_t
}
m = m[complete.cases(m$SPECIES),]

# take means
# does not account for sex
m1 = aggregate(m[, 1:8], by=list(m$SPECIES), mean)
names(m1)[1] = "SPECIES"

# read in toe data gathered from literature
toe = read.csv("~/macroevolution/eco_IBD_oz/data/morphology/Lerista_toes.csv", stringsAsFactors=F)
toe$species = gsub("L. ", "Lerista_", toe$species)
m1 = cbind(m1, toe[match(m1$SPECIES, toe$species), c("hand_digits", "toe_digits")])

# gather in new data on toes gathered from literature
dig = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/morphology/Sphenomorphine_digit_data.csv", stringsAsFactors=F)
for (i in 1:nrow(dig)) {
  m1[m1$SPECIES == dig[i, "SPECIES"], "hand_digits"] = dig[i, "hand_digits"]
  m1[m1$SPECIES == dig[i, "SPECIES"], "toe_digits"] = dig[i, "toe_digits"]
}

# gather in new data on SVL gathered from literature
svl = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/morphology/Sphenomorphine_SVL_data.csv", stringsAsFactors=F)
for (i in 1:nrow(svl)) {
  m1[m1$SPECIES == svl[i, "SPECIES"], "svl"] = svl[i, "svl"]
}

# do phylogenetic PCA
leg_vals = c("toe", "shank", "antebrachium", "toe_digits", "hand_digits")
legs = m1[, c("SPECIES", "svl", leg_vals)]
legs = legs[complete.cases(legs),]
rownames(legs) = legs$SPECIES

keep = intersect(t1$tip.label, legs$SPECIES)
leg_t = drop.tip(t1, setdiff(t1$tip.label, keep))
legs2 = legs[leg_t$tip.label,]
legs3 = legs2

for (i in 3:7) {
  # legs3[legs3[, i] == 0, i] = 0.5
  mod = phylolm(legs3[, i] ~ legs3[, "svl"], phy=leg_t, model="lambda")
  legs3[, i] = mod$residuals
  legs3[, i] = scale(legs3[, i], scale=T)
}

# b = phyl.pca(leg_t, legs2[, 3:7], method="lambda")
c = phyl.pca(leg_t, legs3[, 3:7], method="lambda")

legs3$PC1 = c$S[legs3$SPECIES, 1]
legs3$PC2 = c$S[legs3$SPECIES, 2]
legs3$PC3 = c$S[legs3$SPECIES, 3]

#legs2$PC1 = c$x[legs2$SPECIES, 1]
#legs2$PC2 = c$x[legs2$SPECIES, 2]
#legs2$PC3 = c$x[legs2$SPECIES, 3]

corrplot.mixed(cor(legs3[,2:10], use="complete"), upper='ellipse', lower="number", tl.col="black", tl.cex=0.5, tl.srt=45)

m1$PC1 = legs3[match(m1$SPECIES, legs3$SPECIES), "PC1"]
m1$PC2 = legs3[match(m1$SPECIES, legs3$SPECIES), "PC2"]
m1$PC3 = legs3[match(m1$SPECIES, legs3$SPECIES), "PC3"]

dm = merge(d, m1, by="SPECIES", all=TRUE)

dm$sp1 = NULL
dm$sp2 = NULL
# write it out
dm = dm[,c(1, 5, 18:30)]
write.csv(dm, "~/Dropbox/Sphenomorphine_Gene_Flow/data/morphology/species_OTU_means_PCAs.csv", row.names=F)
