import re
import pandas as pd

# see what inds are not defined in csv
# see what inds are not defined in files

d = pd.read_csv("/Users/sonal/Desktop/phylogeny/original_phylo_data/phylogeny.csv")
# remove all _
# remove all spaces
d['id'] = [re.sub('\s+', '', x) if isinstance(x, str) else x for x in d.id]
d['id'] = [re.sub('_', '', x)  if isinstance(x, str) else x for x in d.id]
d['voucher'] = [re.sub('\s+', '', x) if isinstance(x, str) else x for x in d.voucher]
d['voucher'] = [re.sub('_', '', x)  if isinstance(x, str) else x for x in d.voucher]
d['species2'] = [x.upper() for x in d.species]

loci = ['12s', '16S', 'ATP', 'cytB', 'LDLR', 'nd4']

for locus in loci:
	d[locus] = [False] * d.shape[0]

	file = '/Users/sonal/Desktop/phylogeny/original_phylo_data/nexus/x%s.nex' % locus
	read = 0
	seq = {}

	f = open(file, 'r')
	for l in f:
		if re.search('Matrix', l):
			read = 1
		elif re.search('^;', l):
			read = 0
		if read:
			vals = re.split('\s+', l.rstrip())
			if len(vals) > 1:
				# don't keep sequences with all NNNs
				if not re.search('^N+$', vals[1]):
					seq[vals[0]] = vals[1]

	new = '%s.fasta' % locus
	o = open(new, 'w')
	for ind, s in seq.items():
		# first check if matches id
		if d.ix[d.id == ind].shape[0] > 0:
			d.ix[d.id == ind, locus] = True
			newname = d.ix[d.id == ind, 'voucher'].tolist()[0]
		elif d.ix[d.species2 == ind].shape[0] > 0:
			tmp = d.ix[d.species2 == ind]
			if tmp.shape[0] > 1:
				tmp = tmp[pd.isnull(tmp.id)]
				newname =  tmp.voucher.tolist()[0]
			else:
				newname = d.ix[d.species2 == ind, 'voucher'].tolist()[0]
			d.ix[d.voucher == newname, locus] = True
		else:
			print(locus, '\t', ind)

		o.write('>%s\n%s\n' % (newname, s))
d.to_csv("phylogeny2.csv", index=False)

