library(ape)
library(phangorn)

d = read.csv('~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv', stringsAsFactors=F)
a = read.dna("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/sanger_alns/all_loci.aln.phy", format="sequential")
t = read.tree("~/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/all_genomic_data/species_tree/multifurcating_tree_collapse50.tre")

setno = 1
nodes = seq(length(t$tip.label) + 1, max(t$edge))
for (i in 1:length(nodes)) {
	node = nodes[i]
	tips = Descendants(t, node, type="tips")[[1]]
	if (length(tips) > 2) {
		tips1 = t$tip.label[tips]
		line1 = paste("<distribution id=\"set", setno, ".prior\" spec=\"beast.math.distributions.MRCAPrior\" monophyletic=\"true\" tree=\"@Tree.t:ND4_pos1\">", sep="")
		cat(line1, "\n")
		line2 = paste("<taxonset id=\"set", setno, "\" spec=\"TaxonSet\">", sep="")
		cat("\t", line2, "\n")
		for (j in 1:length(tips1)) {
			line = paste("<taxon idref=\"", tips1[j], "\"/>", sep="")
			cat("\t\t", line, "\n")
		}
		line = "</taxonset>"
		cat("\t", line, "\n")
		line = "</distribution>"
		cat(line, "\n")
		setno = setno + 1
	}
}