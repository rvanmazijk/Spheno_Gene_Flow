import pandas as pd
import os
import glob
import re
import subprocess
import numpy as np

base_dir = '/Users/sonal/Dropbox/Sphenomorphine_Gene_Flow/data'
# my revised ranges based on putative OTUs
indir1 = '/Users/sonal/macroevolution/eco_IBD_oz/data/geography/new_ranges'
# Pascal's original ranges
indir2 = '/Users/sonal/macroevolution/eco_IBD_oz/data/geography/AusHerpRanges'
# outdir
outdir = os.path.join(base_dir, 'geographic_ranges')
outdir1 = os.path.join(outdir, 'species_ranges')
outdir2 = os.path.join(outdir, 'OTU_ranges')

def make_dir(outdir):
	if not os.path.isdir(outdir):
		os.mkdir(outdir)

# make output dirs
[make_dir(dir) for dir in [outdir, outdir1, outdir2]]

# file with OTUs
spfile = os.path.join(base_dir, 'metadata/sphenomorphine_species.csv')
d = pd.read_csv(spfile)

def do_species_ranges():
	for ix, row in d.iterrows():
		sp = row['SPECIES']
		stem = os.path.join(indir2, sp)
		files = glob.glob(stem + '*')
		if len(files) == 3:
			for file in files:
				new_file = re.sub(indir2, outdir1, file)
				subprocess.call("cp %s %s" % (file, new_file), shell=True)
		elif len(files) == 0:
			if pd.isnull(row['OTHER_NAMES']):
				print("%s\tno_geographic_range" % (sp))
			else:
				syns = re.split(',', row['OTHER_NAMES'])
				assigned = False
				for syn in syns:
					stem = os.path.join(indir2, syn)
					files = glob.glob(stem + '*')
					if len(files) == 3:
						for file in files:
							new_file = re.sub(indir2, outdir1, file)
							new_file = re.sub(syn, sp, new_file)
							subprocess.call("cp %s %s" % (file, new_file), shell=True)
						print("%s\tassigned_to_synonyms_range_%s" % (sp, syn))
						assigned = True
						break
				if not assigned:		
					print("%s\tno_geographic_range" % (sp))

def do_OTU_ranges():
	# the fixes because I messed up a few species ...
	revise = {'L_lineata': 'L_lineopunctulata_2', 
	          'L_micra': 'L_microtis_1',
	          'L_lineopunctulata': 'L_lineopunctulata_1',
	          'L_microtis': 'L_microtis_2'}

	for ix, row in d.iterrows():
		otu = row['OTU']
		done = False
		if re.search('Ctenotus', otu) or re.search('Lerista', otu):
			stem = re.sub('Lerista', 'L', otu)
			stem = re.sub('Ctenotus', 'C', stem)

			if stem in revise:
				stem = revise[stem]

			stem = os.path.join(indir1, stem)
			files = glob.glob(stem + '.*')
			if len(files) == 3:
				for file in files:
					new_file = re.sub(stem, os.path.join(outdir2, otu), file)
					subprocess.call("cp %s %s" % (file, new_file), shell=True)
				print("%s\tOTU_revised_geographic_range" % (otu))
				done = True

		if not done:
			stem = os.path.join(indir2, otu)
			files = glob.glob(stem + '*')
			if len(files) == 3:
				for file in files:
					new_file = re.sub(indir2, outdir2, file)
					subprocess.call("cp %s %s" % (file, new_file), shell=True)
			elif len(files) == 0:
				if pd.isnull(row['OTHER_NAMES']):
					print("%s\tno_geographic_range" % (otu))
				else:
					# only look at synonyms if I didn't revise this lineage
					if row['OTU'] == row['SPECIES']:
						syns = re.split(',', row['OTHER_NAMES'])
						assigned = False
						for syn in syns:
							stem = os.path.join(indir2, syn)
							files = glob.glob(stem + '*')
							if len(files) == 3:
								for file in files:
									new_file = re.sub(indir2, outdir2, file)
									new_file = re.sub(syn, otu, new_file)
									subprocess.call("cp %s %s" % (file, new_file), shell=True)
								print("%s\tassigned_to_synonyms_range_%s" % (otu, syn))
								assigned = True
								break
						if not assigned:		
							print("%s\tno_geographic_range" % (otu))
					else:
						print("%s\tno_geographic_range" % (otu))

# do_species_ranges()
do_OTU_ranges()