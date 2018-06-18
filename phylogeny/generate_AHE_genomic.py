import re
import glob
import os

def get_seq(aln):
	seq = {}

	f = open(aln, 'r')
	head = f.next()
	for l in f:
		d = re.split('\s+', l.rstrip())
		seq[d[0]] = d[1]
	f.close()

	return seq


alndir = '/scratch/drabosky_flux/sosi/gene_flow/all_genomic_data/trim/'
alns = glob.glob(alndir + '/*phy')

allseq = {}
loclens = {}
inds = {}

# alns = ['/scratch/drabosky_flux/sosi/gene_flow/all_genomic_data/trim/L115.trim.aln.phy', '/scratch/drabosky_flux/sosi/gene_flow/all_genomic_data/trim/L116.trim.aln.phy']

for aln in alns:
	locname = re.search('(L\d+)', aln).group(1)
	allseq[locname] = get_seq(aln)

for loc in allseq:
	loclens[loc] = len(list(allseq[loc].values())[0])
	for ind in allseq[loc]:
		if ind not in inds:
			inds[ind] = 1

loci = sorted(list(allseq.keys()))
part = {}
curpos = 1
for loc in loci:
	part[loc] = {}
	partfile = '/scratch/drabosky_flux/sosi/gene_flow/all_genomic_data/partitions/%s_partitions.txt' % loc
	if os.path.isfile(partfile):
		f = open(partfile, 'r')
		for l in f:
			info = re.search('=(\S+)', l.rstrip()).group(1)
			info = re.split(',', info)
			for bit in info:
				if not re.search('^\d+-\d+$', bit):
					pos = re.search('p(\d)', l).group(1)
					name = '%s_c%s' % (loc, pos)
				else:
					name = '%s_nc' % (loc)
				if name not in part[loc]:
					part[loc][name] = []
				start = int(re.search('^(\d+)', bit).group(1))
				end = int(re.search('-(\d+)', bit).group(1))
				part[loc][name].append([curpos + start - 1, curpos + end - 1])
		curpos += loclens[loc]
	else:
		part[loc]['%s_all' % loc] = [[curpos, curpos + loclens[loc]  - 1]]
		curpos += loclens[loc]

o = open("/scratch/drabosky_flux/sosi/gene_flow/all_genomic_data/partition_finder/all_AHE.phy", 'w')
totlen = sum([loclens[x] for x in loclens])
o.write('%s %s\n' % (len(inds), totlen))  
for ind in inds:
	totseq = ''
	for loc in loci:
		if ind in allseq[loc]:
			totseq += allseq[loc][ind]
		else:
			totseq += 'N' * loclens[loc]
	o.write('%s\t%s\n' % (ind, totseq))
o.close()

for loc in loci:
	for name in part[loc]:
		for ix, pos in enumerate(part[loc][name]):
			if re.search('_c', name):
				print('%s_%s = %s-%s\\3;' % (name, ix, pos[0], pos[1]))
			else:
				print('%s_%s = %s-%s;' % (name, ix, pos[0], pos[1]))
