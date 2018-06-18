import os
import re
import pandas as pd
import subprocess
import glob
import argparse

parser = argparse.ArgumentParser(description="Make gene tree files.")
parser.add_argument('--g', help="Genus for which to run, comma separated")
parser.add_argument('--m', help="Missing level tolerated.")
args = parser.parse_args()

genus = re.split(',', args.g)
genus = '_'.join(genus)
missing = float(args.m)

outdir = '/scratch/drabosky_flux/sosi/gene_flow/define_clusters/%s_alignments/' % (genus)
seqdir = '/scratch/drabosky_flux/sosi/gene_flow/ind_assemblies/'

if not os.path.isdir(outdir):
	os.mkdir(outdir)

# homology file
hom_file = '%s%s_homology_across_species.txt' % (seqdir, genus)

# ind file
indfile = '%s%s_inds.txt' % (seqdir, genus)

def define_loci(hom_file, indfile, genus, missing):
	inds = {}
	f = open(indfile, 'r')
	for l in f:
		if re.search("KEEP", l):
			d = re.split('\t', l.rstrip())
			inds[d[0]] = 1
	num_ind = len(inds)

	d = pd.read_csv(hom_file, sep = '\t')	
	d = d[d.numMatches >= (num_ind * missing)]
	
	loci = {}
	for c, match in zip(d.contig, d.matches):
		loci[c] = {}
		matches = re.split(',', match)
		matches = [re.split(':', match) for match in matches]
		for match in matches:
			loci[c][match[0]] = match[1]

	return loci, inds


def get_sequences(seqdir, inds):
	seqs = {}

	for ind in inds:
		file = seqdir + '%s.fa' % ind
		f = open(file, 'r')
		for l in f:
			if re.search('>', l):
				id = re.search('>(\S+)', l).group(1)
				id = '%s_%s' % (ind, id)
				seq = f.next().rstrip()
				seqs[id] = seq
		f.close()

	return seqs	


def revcomp(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
	rc = "".join(complement.get(base) for base in reversed(seq))
	return rc


def make_alignments(outdir, loci, seqs, inds):
	summary_file = '%sloci_summary.csv' % outdir
	s = open(summary_file, 'w')
	s.write('locus,sequenced\n')

	for locus in loci:
		out = '%s%s.fa' % (outdir, locus)
		alnout = '%s%s.fa.aln' % (outdir, locus)

		s.write('%s,%.3f\n' % (locus, len(loci[locus]) / float(len(inds))))

		o = open(out, 'w')
		for contig, strand in loci[locus].items():
			seq = seqs[contig]
			if strand == '-':
				seq = revcomp(seq)
			species = re.search('^(\S+)_E\d+_L\d+$', contig).group(1)
			o.write('>%s\n%s\n' % (species, seq))
		o.close()

		subprocess.call("muscle -in %s -out %s -quiet" % (out, alnout), shell=True)
		# subprocess.call("mafft --adjustdirection --auto --quiet %s > %s" % (out, alnout), shell=True)
		os.remove(out)
	s.close()


loci, inds = define_loci(hom_file, indfile, genus, missing)
seqs = get_sequences(seqdir, inds)
make_alignments(outdir, loci, seqs, inds)
