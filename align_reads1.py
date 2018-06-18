import argparse
import os
import subprocess
import pandas as pd

parser = argparse.ArgumentParser(description='reads')
parser.add_argument('--ind', help="ind to run this on")
args = parser.parse_args()
ind = args.ind

c_file = '/scratch/drabosky_flux/sosi/gene_flow/metadata/spheno_ind_data.csv'
seq_dir = '/scratch/drabosky_flux/sosi/gene_flow/cluster_assemblies/'
read_dir = '/scratch/drabosky_flux/sosi/gene_flow/trim_reads/'
out_dir = '/scratch/drabosky_flux/sosi/gene_flow/alignments/'

def get_cluster(c_file, ind):
	d = pd.read_csv(c_file)
	cl = d[d.SAMPLE_ID == ind].CLUSTER.tolist()[0]
	return cl


def prepare_seq(cl, seq_dir):
	seq = '%s%s.fa' % (seq_dir, cl)
	if not os.path.isfile(seq + '.bwt'):
		subprocess.call("~/bin/bwa-0.7.12/bwa index %s" % seq, shell=True)
	if not os.path.isfile(seq + '.fai'):
		subprocess.call("~/bin/samtools-1.3.1/samtools faidx %s" % seq, shell=True)
	if not os.path.isfile(seq.replace('.fa', '.dict')):
		subprocess.call("java -jar ~/bin/picard-tools-2.4.1/picard.jar CreateSequenceDictionary R=%s O=%s" % (seq, seq.replace('.fa', '.dict')), shell=True)


def align_seq(ind, cl, out_dir, read_dir, seq_dir):
	r1 = '%s%s_R1.final.fq.gz' % (read_dir, ind)
	r2 = '%s%s_R2.final.fq.gz' % (read_dir, ind)
	rU = '%s%s_unpaired.final.fq.gz' % (read_dir, ind)
	seq = '%s%s.fa' % (seq_dir, cl)

	out1 = '%s%s.sam' % (out_dir, ind)
	out1b = '%s%s_u.sam' % (out_dir, ind)
	out2 = '%s%s.mateFixed.bam' % (out_dir, ind)
	out3a = '%s%s.mateFixed.sorted1.bam' % (out_dir, ind)
	out3b = '%s%s.mateFixed.sorted2.bam' % (out_dir, ind)
	out3 = '%s%s.mateFixed.sorted.bam' % (out_dir, ind)
	out4 = '%s%s.rg.mateFixed.sorted.bam' % (out_dir, ind)
	intervals = '%s%s.intervals' % (out_dir, ind)
	out5 = '%s%s.realigned.rg.mateFixed.sorted.bam' % (out_dir, ind)

	tmpdir = '%s%s/' % (out_dir, ind)
	if not os.path.isdir(tmpdir):
		os.mkdir(tmpdir)

	# align
	subprocess.call(" ~/bin/bwa-0.7.12/bwa mem -t 4 %s %s %s > %s" % (seq, r1, r2, out1), shell=True)
	subprocess.call(" ~/bin/bwa-0.7.12/bwa mem -t 4 %s %s > %s" % (seq, rU, out1b), shell=True)
	# fixmate
	subprocess.call("java -jar ~/bin/picard-tools-2.4.1/picard.jar FixMateInformation I=%s O=%s" % (out1, out2), shell=True)
	# sorted
	subprocess.call("~/bin/samtools-1.3.1/samtools sort -O bam -o %s -T %s %s" % (out3a, tmpdir, out2), shell=True)
	subprocess.call("~/bin/samtools-1.3.1/samtools sort -O bam -o %s -T %s %s" % (out3b, tmpdir, out1b), shell=True)
	# merge
	subprocess.call("~/bin/samtools-1.3.1/samtools merge %s %s %s" % (out3, out3a, out3b), shell=True)
	# readgroup
	subprocess.call("java -jar ~/bin/picard-tools-2.4.1/picard.jar AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGLB=%s RGPL=Illumina RGPU=%s RGSM=%s" % (out3, out4, ind, ind, ind), shell=True)
	subprocess.call("~/bin/samtools-1.3.1/samtools index %s" % out4, shell=True)
	# indeltarget
	subprocess.call("java -Xmx10g -jar ~/bin/GenomeAnalysisTK.jar -T RealignerTargetCreator -R %s -I %s -o %s -nt 4" % (seq, out4, intervals), shell=True)
	# indelrealigner
	subprocess.call("java -Xmx10g -jar ~/bin/GenomeAnalysisTK.jar -T IndelRealigner -R %s -I %s -targetIntervals %s -o %s" % (seq, out4, intervals, out5), shell=True)

	call = [os.remove(x) for x in [out1, out1b, out2, out3a, out3b, out3, out4, intervals, out4 + '.bai']]
	os.rmdir(tmpdir)

# get cluster
cl = get_cluster(c_file, ind)
# prepare seqfiles
prepare_seq(cl, seq_dir)
# align it all the way until time to call SNPs
align_seq(ind, cl, out_dir, read_dir, seq_dir)
