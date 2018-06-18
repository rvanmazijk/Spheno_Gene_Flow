import re
import os
import pandas as pd
import glob

name = 'demult'
nodes = 1
cpu = 1
mem = 4
hours = 10

basedir = '/scratch/drabosky_flux/sosi/gene_flow/'
indir = os.path.join(basedir, 'raw_reads')
reads = glob.glob(indir + '/*gz')

outdir = os.path.join(basedir, 'ind_reads')
if not os.path.isdir(outdir):
        os.mkdir(outdir)

d = os.path.join(basedir, 'metadata', 'spheno_ind_data.csv')
p = os.path.join(basedir, 'metadata', 'PCR_indices.csv')
a = os.path.join(basedir, 'metadata', 'inline_barcodes.csv')
d = pd.read_csv(d)
p = pd.read_csv(p)
a = pd.read_csv(a)

d = d.merge(p, left_on="PCR_INDEX", right_on="pcr_index")
d = d.merge(a, left_on="INLINE_BARCODE", right_on="inline_barcode_ID")
d = d[d.RAD_STATUS == 'gene_flow_02']

for i, g in d.groupby("PCR_INDEX"):
	index = g.index_x.tolist()[0]
	out = os.path.join(outdir, '%s_inds' % index)
	o = open(out, 'w')
	for ix, row in g.iterrows():
		o.write('%s\t%s\n' % (row['index_y'], row['SAMPLE_ID']))
	o.close()


indices = ['AAGACT', 'CCTGGA', 'GGATTC', 'TTCAAG']
cats = ['A', 'B', 'C', 'D']

for ix, (letter, index) in enumerate(zip(cats, indices)):
	subdir = os.path.join(outdir, index)
	if not os.path.isdir(subdir):
		os.mkdir(subdir)

	n = ix + 1
	r1 = [read for read in reads if re.search('%s_.*R1' % letter, read)][0]
	r2 = [read for read in reads if re.search('%s_.*R2' % letter, read)][0]
	call = "process_radtags -P -1 %s -2 %s -b %s/%s_inds -o %s -D -r -e ecoRI --renz_2 mspI -i gzfastq --disable_rad_check --barcode_dist_1 3" % (r1, r2, outdir, index, subdir)

	o = open("%s%s.sh" % (name, n), 'w')
	o.write("#PBS -N %s%s\n" % (name, n))
        o.write("#PBS -M sosi@umich.edu\n")
        o.write("#PBS -A drabosky_flux\n")
        o.write("#PBS -l qos=flux\n")
        o.write("#PBS -q flux\n")
        o.write("#PBS -l nodes=%s:ppn=%s,mem=%sgb\n" % (nodes, cpu, mem))
        o.write("#PBS -l walltime=%s:00:00\n" % hours)
        o.write("#PBS -j oe\n")
        o.write("#PBS -V\n")

	o.write("module load stacks/1.45\n")
	o.write(call + "\n")

	o.close()
