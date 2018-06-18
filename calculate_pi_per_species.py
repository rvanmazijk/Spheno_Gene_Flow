import re
import subprocess
import pandas as pd
import os
import numpy as np
import argparse
import gzip

parser = argparse.ArgumentParser(description="Calculate het and Fst.")
parser.add_argument('--vcf', help="VCF for which to run.")
parser.add_argument('--out', help="Name of outfile")
parser.add_argument('--lineage', help="Name of lineage file")
args = parser.parse_args()
vcf = args.vcf
out = args.out
lineage = args.lineage

def initialize_hashes(inds):
	# keep track of it all
        all = {}
        for ind in inds:
                all[ind] = {'pi': {'sum': 0, 'sites': 0}}
        all['all'] = {'pi': {'sum': 0, 'sites': 0},
                      'het': {'sum': 0, 'sites': 0}}

	return all


def get_diversity(lineage, vcf, out):
	allowed = ['0/0', '0/1', '1/1']

	f = gzip.open(vcf, 'r')
	for l in f:
		if re.search('#CHROM', l):
			inds = re.split('\t', l.rstrip())
			inds = inds[9:]
			all = initialize_hashes(inds)
		elif not re.search('#', l) and not re.search('INDEL', l):
			d = re.split('\s+', l.rstrip())
			# don't mess with multiallelics
			if len(re.split(',', d[4])) == 1:
				genos = [re.search('^(\S\/\S)', x).group(1) for x in d[9:]]
				# look at inds
				for ind, gen in zip(inds, genos):
					if gen in allowed:
						all[ind]['pi']['sites'] += 1
						if gen == '0/1':
							all[ind]['pi']['sum'] += 1
						
				genos = [x for x in genos if x in allowed]

				# only do it if there is at least one non-missing site
				if len(genos) > 0:
					# calculate proportion hets
					het_prop = genos.count('0/1') / float(len(genos))
				 	all['all']['het']['sum'] += het_prop
					all['all']['het']['sites'] += 1

					alleles = []
					for geno in genos:
						alleles += re.split('/', geno)
					alleles = dict([(x, alleles.count(x)) for x in set(alleles)])
					
					if len(alleles) > 1:
						# https://binhe.org/2011/12/29/calculate-nucleotide-diversity-per-base-pair-summation-method/
						# total alleles
						n = float(np.sum(alleles.values()))
						# minor count
						j = float(np.min(alleles.values()))
						pi_prop = (2 * j * (n - j)) / (n * (n - 1))
					else:
						pi_prop = 0
					all['all']['pi']['sum'] += pi_prop
					all['all']['pi']['sites'] += 1

	f.close()
		
	o = open(out, 'w')
	o.write('type,lineage,n_inds,ind,pi,pi_denom,het,het_denom\n')
	if all['all']['pi']['sites'] > 0:	
		pi = all['all']['pi']['sum'] / float(all['all']['pi']['sites'])
		het = all['all']['het']['sum'] / float(all['all']['het']['sites'])
	else:
		pi = np.nan
		het = np.nan	
	o.write('ALL,%s,%s,NA,%.6f,%s,%.6f,%s\n' % 
                (lineage, len(inds), pi, all['all']['pi']['sites'], het, 
		all['all']['het']['sites']))
	for ind in inds:
		if all[ind]['pi']['sites'] > 0:
			pi = all[ind]['pi']['sum'] / float(all[ind]['pi']['sites'])
			het = pi
		else:
			pi = np.nan
			het = np.nan
		o.write('IND,%s,%s,%s,%.6f,%s,%.6f,%s\n' % 
                        (lineage, len(inds), ind, pi, 
                         all[ind]['pi']['sites'], het, 
                         all[ind]['pi']['sites']))
	o.close()

get_diversity(lineage, vcf, out)
