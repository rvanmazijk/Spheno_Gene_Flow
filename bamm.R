library(ape)
library(coda)
library(BAMMtools)

v <- read.tree("~/Desktop/otus.tre")
is.ultrametric(v)
is.binary.tree(v)
# Now to check min branch length:
min(v$edge.length)
setBAMMpriors(v)

d = read.csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/speciation_rates/bamm_runs/run1/otu_bamm_timevarying/bamm_mcmc.txt")

burnstart <- floor(0.5 * nrow(d))
postburn <- d[burnstart:nrow(d), ]
post_probs <- table(postburn$N_shifts) / nrow(postburn)

effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
effectiveSize(postburn$logPrior)

x = getMeanBranchLengthTree(edata)
xx = data.frame(edata$meanTipLambda, tree$tip.label)

tree = read.tree("~/Dropbox/Sphenomorphine_Gene_Flow/data/speciation_rates/bamm_runs/run1/otu_bamm_timevarying/otu.tre")
edata <- getEventData(tree, nsamples=2000, eventdata = "~/Dropbox/Sphenomorphine_Gene_Flow/data/speciation_rates/bamm_runs/run2/otu_bamm_timevarying/bamm_event.txt", burnin=0.5)
sub <- subtreeBAMM(edata, tips = diff$OTU)
BAMM_rates <- getTipRates(edata)$lambda.avg
jetz_rates = jetzDivRates(tree)
rates  = data.frame(BAMM_rates, jetz_rates)
pdf("~/Desktop/test.pdf", width=6, height=4)
layout(matrix(c(1, 2, 3), nrow=1, byrow=TRUE), widths=c(2, 1, 1))
par(mar=c(4.1,1.1,1.1,0))
# plot.phylo(tree, show.tip.label=F, cex=0.6)
plot.bammdata(edata, legend=F, breaksmethod='jenks', pal="RdYlBu")
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
par(mar=c(4.1,0,1.1,2.1))
ticks = pretty(rates$BAMM_rates, 3)
plot(NULL, xlim=range(ticks), ylim=c(1, nrow(rates))-0.5, xlab="", ylab="", axes=F)
for (i in 1:nrow(rates)) {
  yloc = lastPP$yy[match(rownames(rates)[i], tree$tip.label)]
  x1 = min(ticks)
  x2 = rates[i, "BAMM_rates"]
  lines(c(x1, x2), c(yloc, yloc), lwd=0.5, col="gray60")
}
points(rates[diff$OTU, "BAMM_rates"],
       lastPP$yy[match(diff$OTU, tree$tip.label)], 
       pch=16, col="#1f78b4")
axis(1, ticks, tck=-0.02, labels=NA)
axis(1, ticks, labels=ticks, lwd = 0, line = -0.6, las = 1)
mtext("tip rates: BAMM", side=1, line=1.7)
ticks = pretty(rates$jetz_rates, 3)
plot(NULL, xlim=range(ticks), ylim=c(1, nrow(rates))-0.5, xlab="", ylab="", axes=F)
for (i in 1:nrow(rates)) {
  yloc = lastPP$yy[match(rownames(rates)[i], tree$tip.label)]
  x1 = min(ticks)
  x2 = rates[i, "jetz_rates"]
  lines(c(x1, x2), c(yloc, yloc), lwd=0.5, col="gray60")
}
points(rates[diff$OTU, "jetz_rates"],
       lastPP$yy[match(diff$OTU, tree$tip.label)], 
       pch=16, col="#1f78b4")
axis(1, ticks, tck=-0.02, labels=NA)
axis(1, ticks, labels=ticks, lwd = 0, line = -0.6, las = 1)
mtext("tip rates: equal splits", side=1, line=1.7)
dev.off()


