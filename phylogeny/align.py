import subprocess
import re
import glob

# files = glob.glob("/Users/sonal/Desktop/phylogeny/all_phylo_data/*fasta")
files = glob.glob("/Users/sonal/Dropbox/Sphenomorphine_Gene_Flow/data/phylogeny/sanger_alns/*fasta")
for file in files:
	if not re.search('aln.fasta', file):
		new = re.sub(".fasta", ".aln.fasta", file)
		subprocess.call("mafft --auto --thread 2 %s > %s" % (file, new), shell=True)
