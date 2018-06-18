import re
import glob
import os
import argparse

parser = argparse.ArgumentParser(description="Make gene tree files.")
parser.add_argument('--g', help="Genus for which to run, comma separated")
args = parser.parse_args()

genus = re.split(',', args.g)
genus = '_'.join(genus)


alndir = os.path.join('/scratch/drabosky_flux/sosi/gene_flow/define_clusters/', '%s_alignments' % genus)
alns = glob.glob(alndir + '/*aln')

def get_seq(aln):
	f = open(aln, 'r')
	s = {}
	seqid = ''
	for l in f:
		if re.search('>', l):
			seqid = re.search('>(\S+)', l.rstrip()).group(1)
			s[seqid] = ''
		else:
			s[seqid] += l.rstrip()
	f.close()
	return s

seq = {}
for aln in alns:
	s = get_seq(aln)
	locname = re.search('(%s_\d+).fa' % genus, aln).group(1)
	seq[locname] = s

inds = {}
seqlen = {}
for loc in seq:
	seqlen[loc] = len(list(seq[loc].values())[0])
	for ind in seq[loc]:
		if ind not in inds:
			inds[ind] = 1

loci = sorted(list(seq.keys()))
totlen = sum(list(seqlen.values()))
out = os.path.join('/scratch/drabosky_flux/sosi/gene_flow/define_clusters/', '%s.aln.phy' % genus)
o = open(out, 'w')
o.write('%s %s\n' % (len(inds), totlen))
allseq = {}
for ind in inds:
	totseq = ''
	for locus in loci:
		if ind in seq[locus]:
			totseq += seq[locus][ind]
		else:
			totseq += '-' * seqlen[locus]
	allseq[ind] =  totseq
	o.write('%s\t%s\n' % (ind, totseq))
o.close()

def getdiv(seq1, seq2):
	denom = 0
	diff = 0
	allowed = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']
	for x, y in zip(seq1, seq2):
		if x in allowed and y in allowed:
			denom += 1
			if x != y:
				diff += 1
	if denom > 0:
		dxy = float(diff) / float(denom)
	else:
		dxy = None

	return dxy, denom

out = os.path.join('/scratch/drabosky_flux/sosi/gene_flow/define_clusters/', '%s_div.csv' % genus)
o = open(out, 'w')
o.write('ind1,ind2,dxy,denom\n')
for ind1 in inds:
	for ind2 in inds:
		dxy, denom = getdiv(allseq[ind1], allseq[ind2])
		o.write('%s,%s,%s,%s\n' % (ind1, ind2, dxy, denom))
o.close()
