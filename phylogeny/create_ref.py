import re
import pandas as pd
import glob
import os

def get_seq(locus):
	f = open(locus, 'r')
	seq = {}
	seqid = ''
	for l in f:
		if re.search('>', l):
			seqid = re.search('>(\S+)', l).group(1)
			seq[seqid] = ''
		else:
			seq[seqid] += l.rstrip()
	return seq

loci = glob.glob("/Users/sonal/Desktop/Loci/*fasta")
for locus in loci:
	seq = get_seq(locus)
	if 'Plestiodon_laticeps' in seq:
		locname = re.search('(L\d+).fasta', locus).group(1)
		print('>%s\n%s' % (locname, seq['Plestiodon_laticeps']))