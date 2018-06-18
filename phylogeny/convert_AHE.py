import re
import pandas as pd
import glob
import os

d = pd.read_csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv")
d = d[d.type == 'AHE']

def get_seq(locus):
	f = open(locus, 'r')
	seq = {}
	for l in f:
		if not re.search('^\d', l):
			(seqid, s) = re.split('\s+', l.rstrip())
			seq[seqid] = s
	f.close()
	return seq

replace = {
'Ctenotus_robustus_borealis': 'Ctenotus_robustus',
'Eulamprus_tenuis': 'Concinnia_tenuis',
'Ctenotus_inornatus_helenae': 'Ctenotus_inornatus',
'Ctenotus_quatt': 'Ctenotus_quattuordecimlineatus',
'Glaphyromorphus_crassicaudis': 'Glaphyromorphus_crassicaudus',
'Eulamprus_amplus': 'Concinnia_amplus',
'Notoscincus_wotjulum': 'Notoscincus_ornatus',
'Ctenotus_essingtoni': 'Ctenotus_essingtonii'
}

loci = glob.glob("/Users/sonal/Desktop/Loci/*phylip")
for locus in loci:
	out = re.sub('.phylip', '.fasta', locus)
	o = open(out, 'w')

	seq = get_seq(locus)
	for seqid, s in seq.items():
		# need to drop Lerista bipes
		# because contaminated

		if seqid in replace:
			seqid = replace[seqid]

		# dropping severus because no longer species
		if seqid not in ['Lerista_bipes', 'Ctenotus_severus']:
			name = d.ix[d.SPECIES == seqid, "genomic_sample"].tolist()[0]
			s = re.sub('-', '', s)
			o.write('>%s\n%s\n' % (name, s))
	o.close()
			