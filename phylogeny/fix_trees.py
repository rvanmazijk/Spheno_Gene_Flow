import re
import glob
import os

def get_seq(aln):
	f = open(aln, 'r')
	head = f.next()

	seq = {}
	for l in f:
		d = re.split('\s+', l.rstrip())
		seq[d[0]] = d[1]
	f.close()

	return seq

dir = '/scratch/drabosky_flux/sosi/gene_flow/all_genomic_data'
resdir = os.path.join(dir, 'trimtrees')
phydir = os.path.join(dir, 'trim')
phys = glob.glob(phydir + '/*phy')

for phy in phys:
	locus = re.search('(L\d+)', phy).group(1)
	tree = os.path.join(resdir, '%s.bestTree.tre' % locus)

	if not os.path.isfile(tree):
		print(locus)

		aln = os.path.join(dir, 'trim', '%s.trim.aln.phy' % locus)
		seq = get_seq(aln)

		keep = {}
		for seqid, s in seq.items():
			comp = s.count('A') + s.count('G') + s.count('C') + s.count('T')
			comp = comp / float(len(s))
	
			if comp > 0.05:
				keep[seqid] = s

		f = open(aln, 'w')
		f.write('%s %s\n' % (len(keep), len(list(keep.values())[0])))
		for seqid, s in keep.items():
			f.write('%s\t%s\n' % (seqid, s))
		f.close()
