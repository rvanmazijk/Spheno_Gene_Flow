import re 
import glob
import pandas as pd
import os
import subprocess
import argparse
import sys

parser = argparse.ArgumentParser(description="Run final alignment steps.")
parser.add_argument('--cl', help="Cluster for which to run script.")
args = parser.parse_args()
cl = args.cl

min_qual = 20

c_file1 = '/scratch/drabosky_flux/sosi/gene_flow/metadata/spheno_ind_data.csv'
d = pd.read_csv(c_file1)
d = d.groupby('CLUSTER')
clusters = dict([(cluster, sorted(group['SAMPLE_ID'].tolist())) for cluster, group in d])
inds = clusters[cl]

gatk = '~/bin/GenomeAnalysisTK.jar'
dir = '/scratch/drabosky_flux/sosi/gene_flow/alignments/'
seq_dir = '/scratch/drabosky_flux/sosi/gene_flow/cluster_assemblies/'

bamfiles = ['%s%s.realigned.rg.mateFixed.sorted.bam' % (dir, ind) for ind in inds]
	
raw_vcf = '%s%s.raw.vcf' % (dir, cl)
filt_vcf = '%s%s.filtered.vcf' % (dir, cl)
seq = '%s%s.fa' % (seq_dir, cl)

subprocess.call("~/bin/samtools-1.3.1/samtools mpileup -ugf %s %s | ~/bin/bcftools-1.3.1/bcftools call -vmO v -o %s" % (seq, ' '.join(bamfiles), raw_vcf), shell=True)

# filter the vcf
f = open(raw_vcf, 'r')
o = open(filt_vcf, 'w')

for l in f:
	if re.match('^#', l):
		o.write(l)
	else:
		d = re.split('\t', l.rstrip())
		if float(d[5]) >= min_qual:
			o.write(l)

f.close()
o.close()

for orig, ind in zip(bamfiles, inds):
	final = '%s%s.realigned.rg.mateFixed.sorted.final.bam' % (dir, ind)
	if not os.path.isfile(final):
		recal = '%s%s.recal.table' % (dir, ind)

		subprocess.call("~/bin/samtools-1.3.1/samtools index %s" % orig, shell=True)
		subprocess.call('java -Xmx10g -jar %s -T BaseRecalibrator -R %s -knownSites %s -I %s -o %s' % (gatk, seq, filt_vcf, orig, recal), shell=True)
		subprocess.call('java -Xmx10g -jar %s -T PrintReads -R %s -I %s --BQSR %s -o %s' % (gatk, seq, orig, recal, final), shell=True) 

		os.remove(recal)
		# os.remove(orig)
os.remove(raw_vcf)
os.remove(filt_vcf)
os.remove(filt_vcf + '.idx')
