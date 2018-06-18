import argparse
import glob
import os
import multiprocessing as mp
import pandas as pd
import re
import subprocess
import random

"""
Sonal Singhal
created on 23 June 2016
Written assuming:
	* mafft 7.294
	* RAxML 8.2.4
This script borrows heavily from:
https://github.com/faircloth-lab/phyluce/blob/master/bin/align/phyluce_align_seqcap_align
"""

def get_args():
    parser = argparse.ArgumentParser(
        description="Align, possibly trim, and infer gene trees for UCE loci.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--dir",
        type=str,
        default=None,
        help="The base directory if running the pipeline."
    )

    parser.add_argument(
        "--outdir",
        type=str,
	default=None,
        help="The directory with the phylogeny, if "
             "not using in context of a pipeline."
    )

    parser.add_argument(
        "--raxmltrees",
        action="store_true",
        default=False,
        help="Will infer RAxML gene trees for each alignment if "
	     "flagged."
    )

    parser.add_argument(
        "--phymltrees",
        action="store_true",
        default=False,
        help="Will infer PhyML gene trees for each alignment if "
             "flagged."
    )

    parser.add_argument(
        "--trim",
        action="store_true",
        default=False,
        help="Will trim alignments using gblocks if flagged."
    )

    parser.add_argument(
        "--b1",
        type=float,
        default=0.5,
        help="GBLOCKS -b1 proportion; min # of seqs" +
	     " required for a conserved position"
    )

    parser.add_argument(
        "--b2",
        type=float,
        default=0.85,
        help="GBLOCKS -b2 proportion; min # of seqs " +
	     " required to be at a flanking position"
    )

    parser.add_argument(
        "--b3",
        type=int,
        default=8,
        help="GBLOCKS -b3 proportion; max number of" +
             " contiguous nonconserved positions"
    )

    parser.add_argument(
        "--b4",
        type=int,
        default=10,
        help="GBLOCKS -b4 proportion;" +
             " minimum block length"
    )

    parser.add_argument(
        "--CPU",
        type=int,
        default=1,
        help="""Process alignments in parallel using --CPU for alignment. """ +
        """This is the number of PHYSICAL CPUs."""
    )

    parser.add_argument(
	"--mafft",
	type=str,
	default=None,
	help="Full path to mafft executable."
	)

    parser.add_argument(
        "--gblocks",
        type=str,
        default=None,
        help="Full path to GBLOCKS executable."
        )

    parser.add_argument(
        "--raxml",
        type=str,
        default=None,
        help="Full path to RAxML executable."
        )

    parser.add_argument(
        "--jmodel",
        type=str,
        default=None,
        help="Full path to jModelTest jar."
        )

    return parser.parse_args()


def sub_raxml(file, outdir, raxml):

	locus = re.sub('^.*/', '', file)
	locus = re.sub('\.trim.aln.*', '', locus)

	os.chdir(outdir)
	
	partitions = '/scratch/drabosky_flux/sosi/gene_flow/all_genomic_data/partitions/%s_partitions.txt' % locus
        orig_boot = 'RAxML_bootstrap.%s' % locus
        orig_tree = 'RAxML_bipartitions.%s' % locus

        new_boot = '%s.bootstrap.trees' % locus
        new_tree = '%s.bestTree.tre' % locus

	if not os.path.isfile(new_tree):
		subprocess.call('%s -x %s -# 100 -p %s -m GTRCAT -f a -n %s -s %s' % 
	                        (raxml, random.randint(0,1000), random.randint(0,1000), 
	                        locus, file), shell=True)

		os.rename(orig_boot, new_boot)
		os.rename(orig_tree, new_tree)

		subprocess.call("rm RAxML_*%s" % locus, shell=True) 

	return new_tree


result_list = []
def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    result_list.append(result)


def run_raxml(args):
	treedir = '/scratch/drabosky_flux/sosi/gene_flow/all_genomic_data/trimtrees'
	unf = open('/home/sosi/gene_flow/unfinished', 'r')
	todo = {}
	for l in unf:
		todo[l.rstrip()] = 1
	unf.close()
	phys = ['/scratch/drabosky_flux/sosi/gene_flow/all_genomic_data/trim/%s.trim.aln.phy' % locus for locus in todo]
	# phys = glob.glob('/scratch/drabosky_flux/sosi/gene_flow/all_genomic_data/trim/*phy')

	if not os.path.isdir(treedir):
		os.mkdir(treedir)

	if args.CPU > 1:
		pool = mp.Pool(args.CPU)
	
		dirs = [treedir] * len(phys)
		raxml = [args.raxml] * len(phys)
		for i in range(len(phys)):
			pool.apply_async(sub_raxml, args=(phys[i], treedir, args.raxml, ), callback=log_result)
		pool.close()
		pool.join()

def main():
	args = get_args()
	run_raxml(args)

if __name__ == "__main__":
	main()
