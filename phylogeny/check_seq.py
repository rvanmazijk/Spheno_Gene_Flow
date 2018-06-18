import re
import pandas as pd

def get_seq(seqfile):
	s = open(seqfile, 'r')
	id = ''
	seqs = {}

	for l in s:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			if id in seqs:
				print(id, seqfile)
			seqs[id] = ''
		else:
			seqs[id] += l.rstrip()
	s.close()

	return seqs

loci = ['12s', '16S', 'ATP', 'cytB', 'LDLR', 'nd4']
d = pd.read_csv("~/Dropbox/Sphenomorphine_Gene_Flow/data/metadata/sphenomorphine_species.csv")

for locus in loci:
	seq = '/Users/sonal/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/all_phylo_data/%s.fasta' % locus
	out = re.sub('all_phylo_data', 'sanger_alns', seq)

	seq = get_seq(seq)
	tmp = d.ix[d[locus] == True]
	inds = tmp.sanger_sample.tolist()

	for ind in inds:
		if ind not in seq:
			print('%s %s' % (locus, ind))

	'''
	for ind in seq:
		if ind not in inds:
			print('%s %s' % (locus, ind))
	'''

	o = open(out, 'w')
	for ind in seq:
		if ind in d.sanger_sample.tolist():
			o.write('>%s\n%s\n' % (ind, seq[ind]))
	o.close()
