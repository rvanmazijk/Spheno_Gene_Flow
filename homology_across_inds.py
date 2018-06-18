import re
import subprocess
import os
import argparse
import random
import string
import sys
import pandas as pd

parser = argparse.ArgumentParser(description='Get homology file for groups of interest.')
parser.add_argument('--g', help="genera to include, comma separated")
parser.add_argument('--c', help="cluster")
args = parser.parse_args()
g = re.split(',', args.g)
c = float(args.c)

def get_seq(seqfile):
	f = open(seqfile, 'r')
	numseq = 0
	for l in f:
		if re.search('>', l):
			numseq += 1
	f.close()
	return numseq

indfile = '/scratch/drabosky_flux/sosi/gene_flow/metadata/spheno_ind_data.csv'
d = pd.read_csv(indfile)
d = d[d.GENUS.isin(g)]
genus = '_'.join(g)

fadir = '/scratch/drabosky_flux/sosi/gene_flow/ind_assemblies/'
orgfiles = [(ind, '%s%s.fa' % (fadir, ind)) for ind in d.SAMPLE_ID.tolist()]
orgfiles = dict(orgfiles)
# cluster seqs at this value
WCLUST = c
# min number of seqs needed to include
minseq = 1000 

out = fadir + '%s_inds.txt' % genus
o = open(out, 'w')
# filter files
files = {}
for ind, file in orgfiles.items():
	if os.path.isfile(file):
		numseq = get_seq(file)
		if numseq >= minseq:
			files[ind] = file
			o.write('%s\tKEEP\n' % ind)
		else:
			o.write('%s\tdropped_too_few_seq\n' % ind)	
	else:	
		o.write('%s\tdropped_no_FASTA\n' % ind)
o.close()

def create_starter(dir, file, genus, ix, ind):
	homhash = {}

	starting = '%s%s.tmp.fa' % (dir, genus)
	f = open(file, 'r')
	o = open(starting, 'w')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			id = ind + '_' + id
			newid = '%s_%s' % (genus, ix)

			homhash[newid] = {}
			homhash[newid][id] = '+'
			
			seq = f.next().rstrip()
			o.write('>%s\n%s\n' % (newid, seq))
			ix += 1
	f.close()
	o.close()

	return starting, homhash, ix


def vsearch(dir, tmpfile, file, genus, num):
	out = '%s%s_%s_search' % (dir, genus, num)
	subprocess.call("vsearch --usearch_global %s --db %s --userout %s --id %s --userfields query+target+evalue+id+qstrand --strand both --threads 4" % (file, tmpfile, out, WCLUST), shell=True)

	return out


def create_new_tmp(dir, tmpfile, file, results, homhash, genus, ix, ind):
	matches1 = {}
	matches2 = {}

	f = open(results, 'r')
	for l in f:
		d = re.split('\s+', l.rstrip())
		# is this step necessary?
		# makes sure it is 1 to 1
		match = ind + '_' + d[0]
		if d[1] not in matches1 and match not in matches2:
			matches1[d[1]] = {'match': match, 'perc': float(d[3]), 'strand': d[4]}
			matches2[match] = {'match': d[1], 'perc': float(d[3]), 'strand': d[4]}
		elif match in matches2 and d[1] not in matches1:
			if float(d[3]) > matches2[match]['perc']:
                        	matches1[d[1]] = {'match': match, 'perc': float(d[3]), 'strand': d[4]}
                        	matches2[match] = {'match': d[1], 'perc': float(d[3]), 'strand': d[4]}
	f.close()
	os.remove(results)

	for c in matches2:
		homhash[matches2[c]['match']][c] = matches2[c]['strand']

	f = open(file, 'r')
	o = open(tmpfile, 'a')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			seq = f.next().rstrip()
			if id not in matches2:
				new_id = '%s_%s' % (genus, ix)
				ix += 1

				id = ind + '_' + id
				homhash[new_id] = {}	
				homhash[new_id][id] = '+'
		
				o.write('>%s\n%s\n' % (new_id, seq))
	f.close()
	o.close()

	return (tmpfile, homhash, ix)

inds = sorted(list(files.keys()))
ix = 0
tmpfile, homhash, ix = create_starter(fadir, files[inds[0]], genus, ix, inds[0])
for num, ind in enumerate(inds[1:]):
	file = files[inds[num + 1]]
	results = vsearch(fadir, tmpfile, file, genus, num)
	(tmpfile, homhash, ix) = create_new_tmp(fadir, tmpfile, file, results, homhash, genus, ix, ind)
os.remove(tmpfile)

o = open('%s%s_homology_across_species.txt' % (fadir, genus), 'w')
o.write('contig\tmatches\tnumMatches\n')
for c, matches in homhash.items():
	matches = ['%s:%s' % (match, homhash[c][match]) for match in matches]
	o.write('%s\t%s\t%s\n' % (c, ','.join(matches), len(matches)))
o.close()
