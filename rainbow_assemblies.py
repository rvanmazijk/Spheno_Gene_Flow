import re
import glob
import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='reads')
parser.add_argument('--ind', help="ind to run this on")
args = parser.parse_args()
ind = args.ind

base_dir = '/scratch/drabosky_flux/sosi/gene_flow/'
out_dir = os.path.join(base_dir, 'ind_assemblies')
file1 = os.path.join(base_dir, 'trim_reads', '%s_R1.final.fq.gz' % ind)
file2 = os.path.join(base_dir, 'trim_reads', '%s_R2.final.fq.gz' % ind)
	
out1 = '%s/%s.rainbow1' % (out_dir, ind)
out2 = '%s/%s.rainbow2' % (out_dir, ind)
out3 = '%s/%s.rainbow3' % (out_dir, ind)
out4 = '%s/%s.fa' % (out_dir, ind)	
	
subprocess.call('~/bin/rainbow_2.0.4/rainbow cluster -1 %s -2 %s > %s\n' % (file1, file2, out1), shell=True)
subprocess.call('~/bin/rainbow_2.0.4/rainbow div -i %s -o %s\n' % (out1, out2), shell=True)
subprocess.call('~/bin/rainbow_2.0.4/rainbow merge -i %s -a -o %s\n' % (out2, out3), shell=True)
subprocess.call('~/bin/rainbow_2.0.4/select_best_rbcontig_plus_read1.pl %s %s > %s\n' % (out3, out2, out4), shell=True)

os.remove(out1)
os.remove(out2)
os.remove(out3)
