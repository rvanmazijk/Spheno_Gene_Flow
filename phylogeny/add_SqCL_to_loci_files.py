import re
import glob
import pandas as pd

def get_seq(seqfile):
	s = open(seqfile, 'r')
	id = ''
	seqs = {}

	for l in s:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			seqs[id] = ''
		else:
			seqs[id] += l.rstrip()
	s.close()

	return seqs

def rev_comp(seq):
	seq = seq.upper()
	seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
	return "".join([seq_dict[base] for base in reversed(seq)])

files = glob.glob("/Users/sonal/Desktop/genomic_to_add/*fasta")
for file in files:
	sample = re.sub('/Users/sonal/Desktop/genomic_to_add/', '', file)
	sample = re.sub('.fasta', '', sample)

	seq = get_seq(file)
	match = re.sub("to_add/", "to_add/matches/", file)
	match = re.sub(".fasta", "_matches.csv", match)

	d = pd.read_csv(match)
	d = d[d.status.isin(['easy_recip_match', 'complicated_recip_match'])]

	for ix, row in d.iterrows():
		s = seq[row['contig']]
		if row['orr'] == '-':
			s = rev_comp(s)
		locfile = "/Users/sonal/Desktop/Loci/%s.fasta" % row['match']
		o = open(locfile, 'a')
		o.write('>%s\n%s\n' % (sample, s))