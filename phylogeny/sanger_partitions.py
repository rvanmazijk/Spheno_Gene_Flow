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

loci = ['nd4', 'cytB', '12s', '16S', 'ATP', 'LDLR']
allseq = {}

for locus in loci:
	seq = '/Users/sonal/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/sanger_alns/%s.aln.fasta' % locus
	allseq[locus] = get_seq(seq)

inds = {}
for locus in allseq:
	for ind in allseq[locus]:
		if ind not in inds:
			inds[ind] = 1

inds = sorted(list(inds.keys()))
loclen = {}
curstart = 1
for locus in loci:
	loclen[locus] = len(list(allseq[locus].values())[0])
	print(locus, curstart, curstart + loclen[locus] - 1)
	curstart += loclen[locus]

printseq = {}
o = open("/Users/sonal/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/sanger_alns/all_loci.aln.phy", 'w')
o2 = open("/Users/sonal/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/sanger_alns/all_loci.aln.fasta", 'w')
for ind in inds:
	indseq = ''
	for locus in loci:
		if ind in allseq[locus]:
			indseq += allseq[locus][ind]
		else:
			indseq += '-' * loclen[locus]
	printseq[ind] = indseq

# get rid of plestiodon
del printseq['S8']

o.write('%s %s\n' % (len(printseq), len(indseq)))
for ind, indseq in printseq.items(): 
	o.write('%s %s\n' % (ind, indseq))
	o2.write('>%s\n%s\n' % (ind, indseq))
	# print(len(indseq))
o.close()
