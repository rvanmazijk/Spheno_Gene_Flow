import re
import subprocess
import pandas as pd
import os
import numpy as np
import argparse
import gzip
import random
from numbers import Number
from geopy.distance import vincenty
import itertools as it

parser = argparse.ArgumentParser(description="Calculate het and Fst.")
parser.add_argument('--cl', help="Cluster for which to run.")
parser.add_argument('--cov', help="Coverage VCF.")
args = parser.parse_args()
cl = args.cl
cov = int(args.cov)
n_boot = 100

maindir = '/scratch/drabosky_flux/sosi/gene_flow/'
c_file = os.path.join(maindir, 'metadata/spheno_ind_data.csv')
vcf_dir = os.path.join(maindir, 'variants')
cov_dir = os.path.join(maindir, 'coverage')
out_dir = os.path.join(maindir, 'divergence')
dist_file = c_file
mt_file = '/Volumes/heloderma4/sonal/skink_mtDNA/cytb_alignment_30July15.fixed.aln.fa'

def get_mt(mt_file):
	mt = {}
	id = ''
	f = open(mt_file, 'r')
	for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			mt[id] = ''
		else:
			mt[id] += l.rstrip()
	f.close()
	return mt


def get_dist_hash(dist_file):
	# data file with individuals and lats / longs
	inds = pd.read_csv(dist_file)
	# get rid of undefined rows
	inds = inds[np.isfinite(inds.LAT)]
	inds = inds[np.isfinite(inds.LON)]

	# gets a dictionary where key is sample id and value is tuple of lat / long
	latlong = dict([(x, (y, z)) for x, y, z in zip(inds.SAMPLE_ID, inds.LAT, inds.LON)])
	return latlong


def get_distance(latlong, ind1, ind2):
        if ind1 in latlong and ind2 in latlong:
                # these distances are calculated based on a model
                #       in which earth is an oblate spheroid
                # more typically used is the great-circle distance
                #       which assumes spherical earth but this isn't true
                dist = vincenty(latlong[ind1], latlong[ind2]).meters
                return round(dist, 2)
        else:
                return np.nan


def get_clusters(c_file):
        d = pd.read_csv(c_file)
        d = d[d.CLUSTER.notnull()]
        d = d.groupby('CLUSTER')

        clusters = dict([(name, sorted(group['SAMPLE_ID'].tolist())) for name, group in d])

        return clusters


def get_mt_dist(mt, ind1, ind2):
	allowed = ['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']

	if ind1 in mt and ind2 in mt:
		seq1 = mt[ind1].upper()
		seq2 = mt[ind2].upper()

		diff = 0
		denom = 0

		for bp1, bp2 in zip(seq1, seq2):
			if bp1 in allowed and bp2 in allowed:
				denom += 1
				if bp1 != bp2:
					diff += 1
		diff = diff / float(denom)
		return (diff, denom)
	else:
		return (np.nan, np.nan)


def fst_estimator(counts, sample_sizes):
	'''
	modified from G. Bradford's R code in bedassle
	'calculate.pairwise.Fst'

	both inputs are arrays where each row is an individual
	and each column is a SNP
	'''

	counts = np.array(counts)
	sample_sizes = np.array(sample_sizes).astype('float')

	pop_af = counts / sample_sizes
	mean_af = np.sum(counts, axis=0) / np.sum(sample_sizes, axis = 0)

	MSP = np.sum((pop_af - mean_af) ** 2 * sample_sizes, axis=0)
	MSG = np.sum((1 - pop_af) * pop_af * sample_sizes, axis=0) \
			* (1 / np.sum(sample_sizes - 1, axis=0))
	n_c = np.sum(sample_sizes, axis = 0) - np.sum(sample_sizes ** 2, axis=0) \
			/ np.sum(sample_sizes, axis=0)

	fst = np.sum(MSP - MSG) / np.sum(MSP + (n_c - 1) * MSG)

	return fst


def fst_reich(counts, sample_sizes):
	counts1 = np.array(counts)[0]
	counts2 = np.array(counts)[1]

	sample_sizes1 = np.array(sample_sizes).astype('float')[0]
	sample_sizes2 = np.array(sample_sizes).astype('float')[1]

	h1 = counts1 * (sample_sizes1 - counts1) / (sample_sizes1 * (sample_sizes1 - 1))
	h2 = counts2 * (sample_sizes2 - counts2) / (sample_sizes2 * (sample_sizes2 - 1))
	
	N = []
	D = []

	for _a1, _a2, _n1, _n2, _h1, _h2 in zip(counts1, counts2, sample_sizes1, sample_sizes2, h1, h2):
		n = ((_a1 / _n1) - (_a2 / _n2)) ** 2 - (_h1 / _n1) - (_h2 / _n2)
		N.append(n)
		d = n + _h1 + _h2
		D.append(d)

	F = np.sum(N) / np.sum(D)

	return F

def sample_wr(population, k):         
    "used for bootstrap sampling"
    "Chooses k random elements (with replacement) from a population"
    n = len(population)
    _random, _int = random.random, int  # speed hack
    ix = [_int(_random() * n) for i in it.repeat(None, k)]
    which = [population[i] for i in ix]

    return which

def sample_wr1(population, k):
    "used for bootstrap sampling"
    "Chooses k random elements (with replacement) from a population"
    n = len(population[0])
    _random, _int = random.random, int  # speed hack
    ix = [_int(_random() * n) for i in it.repeat(None, k)]
    which1 = population[0][ix]
    which2 = population[1][ix]

    which = np.array([which1, which2])	

    return which

def format_out(header, res):
	out = []
	for key in header:
		val = res[key]
		if isinstance(val, float):
			val = '%.6f' % val
		elif isinstance(val, int):
			val = str(val)
		out.append(val)
	return out

def get_divergence(cl, inds, vcf_dir, out_dir, mt, latlong):
	diff = { '0/0': {'0/1': 0.5, '1/1': 1, '0/0': 0},
	        '0/1': {'0/1': 0, '1/1': 0.5, '0/0': 0.5},
	        '1/1': {'0/1': 0.5, '1/1': 0, '0/0': 1} }
	count = {'0/0': 0, '1/1': 2, '0/1': 1, './.': np.nan}

	file = os.path.join(vcf_dir, '%s.qual_filtered20.cov_filtered_min%s_max3x.vcf.gz' % (cl, cov))

	div = {}
	for ix, ind1 in enumerate(inds):
		div[ind1] = {}
		for ind2 in inds[(ix + 1):]:
			div[ind1][ind2] = []
		

	# for calculating fst
	counts = dict([(ind, []) for ind in inds])

	f = gzip.open(file, 'r')
	for l in f:
		if not re.search('#', l) and not re.search('INDEL', l):
			d = re.split('\s+', l.rstrip())
			# don't mess with multiallelics
			if len(re.split(',', d[4])) == 1:
				genos = d[9:]
				genos = [re.search('^(\S\/\S)', x).group(1) for x in genos]

				# variable site to be used in fst
				if d[4] in ['A', 'T', 'C', 'G']:
					for ind, geno in zip(inds, genos):
						counts[ind].append(count[geno])

				# get divergence data
				genos = dict(zip(inds, genos))
				for ind1 in div:
					for ind2 in div[ind1]:
						if genos[ind1] != './.' and genos[ind2] != './.':
							div[ind1][ind2].append(diff[genos[ind1]][genos[ind2]])
	f.close()


	out = os.path.join(out_dir, '%s.divergence_cov%s.csv' % (cl, cov))
	o = open(out, 'w')
	header = ['cl', 'ind1', 'ind2', 'geo_dist', 'nuc_dxy',
			'nuc_denom', 'dxy_mean',
			'dxy_sd', 'fst', 'fst_mean', 'fst_sd',
			'fst_denom', 'mt_dxy', 'mt_denom']

	o.write('%s\n' % ','.join(header))
	for ind1 in div:
		for ind2 in div[ind1]:
			res = {}
			for val in header:
				res[val] = np.nan
			res['ind1'] = ind1
			res['ind2'] = ind2
			res['cl'] = cl

			# have at least one comparison
			if len(div[ind1][ind2]) > 0:
				res['nuc_dxy'] = np.sum(div[ind1][ind2]) / float(len(div[ind1][ind2]))
				boots = []
				# now do the bootstraps
				# 100% of original data with replacement
				for ix in range(0, n_boot):
					n_samp = len(div[ind1][ind2])
					sample = sample_wr(div[ind1][ind2], n_samp)
					bootdxy = np.sum(sample) / float(len(sample))
					boots.append(bootdxy)
				res['dxy_mean'] = np.mean(boots)
				res['dxy_sd'] = np.std(boots)
			res['nuc_denom']  = len(div[ind1][ind2])

			alleles = np.array([counts[ind1], counts[ind2]])
			to_mask = np.any(np.isnan(alleles), axis=0)
			alleles = alleles[:, -to_mask]
			res['fst_denom'] = len(alleles[0])
			if len(alleles[0]) > 0:
				sizes = [[2] * len(alleles[0]), [2] * len(alleles[0])]
				res['fst'] = fst_reich(alleles, sizes)
				boots = []
				for ix in range(0, n_boot):
					n_samp = len(alleles[0])
					sample = sample_wr1(alleles, n_samp)
					fst = fst_reich(sample, sizes)
					boots.append(fst)
				res['fst_mean'] = np.mean(boots)
				res['fst_sd'] = np.std(boots)
					
			# get geographic distance
			res['geo_dist'] = get_distance(latlong, ind1, ind2)

			# get mito distance
			(res['mt_dxy'], res['mt_denom']) = get_mt_dist(mt, ind1, ind2)

			# pass in the result hash to get formatted output
			o.write('%s\n' % ','.join(format_out(header, res)))

	o.close()

if os.path.isfile(mt_file):
	mt = get_mt(mt_file)
else:
	mt = {}
latlong = get_dist_hash(dist_file)
clusters = get_clusters(c_file)
get_divergence(cl, clusters[cl], vcf_dir, out_dir, mt, latlong)
