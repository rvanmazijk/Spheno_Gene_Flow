import re

def get_seq(seqfile):
	
	f = open(seqfile, 'r')

	id = ''
	s = {}

	# first line of phylip file
	head = f.readline()

	for l in f:
		d = re.split('\s+', l.rstrip())
		s[d[0]] = d[1]
	f.close()

	return s

def print_seq(out, seq):
	numseq = len(seq)
	loci_len = len(list(seq.values())[0])

	o = open(out, 'w')
	o.write('#NEXUS\n')
	o.write('Begin data;\n')
	o.write('Dimensions ntax=%s nchar=%s\n;' % (numseq, loci_len))
	o.write('Format datatype=dna missing=N gap=-;\n')
	o.write('Matrix\n')
	for ind, s in seq.items():
		o.write('%s\t%s\n' % (ind, s))
	o.write(';\n')
	o.write('End;\n')

seqfile = '/Users/sonal/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/sanger_alns/all_loci.aln.phy'
out = '/Users/sonal/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/sanger_alns/all_loci.aln.nex'

seq = get_seq(seqfile)
print_seq(out, seq)
