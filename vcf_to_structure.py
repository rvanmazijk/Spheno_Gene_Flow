import argparse
import re
import pandas as pd
import gzip
import random

parser = argparse.ArgumentParser()
parser.add_argument("-v", required=True, help='vcf')
args = parser.parse_args()

vcf_dir = '/scratch/drabosky_flux/sosi/gene_flow/variants/'
out_dir = '/scratch/drabosky_flux/sosi/gene_flow/structure/'
miss = 0.25
# only take one snp per block
thin = 50

def make_file(cl, vcf):
	# get the vcf file
	f = gzip.open(vcf, 'r')
	snps = {}
	for l in f:
		if re.search('^#CHROM', l):
			inds = re.split('\t', l.rstrip())[9:]
		elif not re.search('^#', l):
			d = re.split('\t', l.rstrip())
			# only using variable sties
			if d[4] in ['A', 'T', 'C', 'G'] and d[3] in ['A', 'T', 'C', 'G']:
				genos = [re.search('^(\S\S\S)', x).group(1) for x in d[9:]]
				genos = [re.split('/', x) for x in genos]

				flat_genos = [x for inner_list in genos for x in inner_list]

				if flat_genos.count('0') > 0 and flat_genos.count('1') > 0:
					# only take it if missingness is less than 0.25
					if ( (flat_genos.count('1') + flat_genos.count('0')) / (len(inds) * 2.0) ) > (1 - miss):
						if d[0] not in snps:
							snps[d[0]] = {}
						snps[d[0]][int(d[1])] = genos
	f.close()


	def get_blocks(values, dist):
	    mi, ma = 0, 0
	    result = []
	    temp = []
	    for v in sorted(values):
		if not temp:
		    mi = ma = v
		    temp.append(v)
		else:
		    if abs(v - mi) < dist and abs(v - ma) < dist:
			temp.append(v)
			if v < mi:
			    mi = v
			elif v > ma:
			    ma = v
		    else:
			if len(temp) > 0:
			    result.append(temp)
			mi = ma = v
			temp = [v]
	    return result

	winners = {}
	for contig in snps:
		locs = sorted(snps[contig].keys())
		blocks = get_blocks(locs, thin)

		winners[contig] = {}

		for block in blocks:
			keep = random.choice(block)
			winners[contig][keep] = 1

	contigs = sorted(winners.keys())

	outfile = '%s%s.structure' % (out_dir, cl)
	o = open(outfile, 'w')
	for ix, ind in enumerate(inds):
		geno1 = ''
		geno2 = ''
		for c in contigs:
			keepers = sorted(winners[c].keys())
			for loc in keepers:
				geno1 += ' ' + ' '.join(snps[c][loc][ix][0])
				geno2 += ' ' +  ' '.join(snps[c][loc][ix][1])

		geno1 = geno1.replace('.', '-9')
		geno2 = geno2.replace('.', '-9')

		o.write('%s 1 0 1 1%s\n' % (ind, geno1))
		o.write('%s 1 0 1 1%s\n' % (ind, geno2))
	o.close()

	return outfile

cl = re.sub(vcf_dir, '', args.v)
cl = re.sub('\..*', '', cl)
make_file(cl, args.v)
