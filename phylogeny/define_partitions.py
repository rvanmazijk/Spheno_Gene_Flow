import re
import glob
import subprocess
import pandas as pd
import os

match = '/Users/Sonal/Desktop/phylogeny/ref/Plestiodon_laticeps_to_Anolis_proteins.txt'
ref = 'anolis_proteins.fa'
alndir = '/Users/Sonal/Desktop/phylogeny/all_genomic_data/aligned/'

# trimming values
het_window = 5
too_much_het = 0.2

alns = glob.glob(alndir + '*aln')
f = open(match, 'r')
match = {}
for l in f:
	d = re.split('\t', l.rstrip())
	if d[0] not in match:
		match[d[0]] = d[1]
f.close()

def get_seq(seqfile):
	s = open(seqfile, 'r')
	id = ''
	seqs = {}

	for l in s:
		if re.search('>', l):
			id = re.search('>(\S+)', l.rstrip()).group(1)
			seqs[id] = ''
		else:
			seqs[id] += l.rstrip().upper()

	s.close()

	return seqs

def trim(seq):

	trimseq = {}
	for seqid in seq:
		trim = list(seq[seqid])
		s = list(seq[seqid])
		het = [0 if x in ['A', 'T', 'C', 'G', '-', 'N'] else 1 for x in s]
		
		for start in range(0, len(het), het_window):
			end = start + het_window
			if end > len(het):
				end = len(het)
			hetval = sum([het[i] for i in range(start, end)])
			hetval = hetval / float(het_window)
			
			if hetval >= too_much_het:
				for i in range(start, end):
					trim[i] = '-'
		trimseq[seqid] = ''.join(trim)
	
	return trimseq

def rev_comp(seq):
	seq = seq.upper()
	seq_dict = {'A':'T','T':'A','G':'C','C':'G', '-': '-', 
				'M': 'K', 'R': 'Y', 'W': 'W', 'S': 'S', 'Y': 'R', 'K': 'M',
				'V': 'B', 'H': 'D', 'D': 'H', 'B': 'V', 'N': 'N'}
	return "".join([seq_dict[base] for base in reversed(seq)])

for aln in alns:
	locname = re.search('(L\d+)', aln).group(1)

	seq = get_seq(aln)
	trimseq = trim(seq)
	# out = open('%s_partitions.txt' % locname, 'w')

	if locname in match:
		pass
		'''
		winner = ''
		wincount = 0
		for id, s in trimseq.items():
			count = s.count('A') + s.count('T') + s.count('C') + s.count('G')
			if count > wincount:
				winner = s
				wincount = count
		winner = re.sub('-', 'N', winner)
				
		p_out = '%s_prot.fasta' % locname
		subprocess.call("samtools faidx %s %s > %s" % (ref, match[locname], p_out), shell=True)
				
		d_out = '%s_dna.fasta' % (locname)
		o = open(d_out, 'w')
		o.write('>dna\n%s\n' % winner)
		o.close()
				
		call = subprocess.Popen("exonerate -m protein2genome -q %s -t %s --showtargetgff TRUE --showalignment NO --showvulgar 0 -n 1 | grep 'cds'" % (p_out, d_out), shell=True, stdout=subprocess.PIPE)
		lines = [line.decode('utf-8').rstrip() for line in call.stdout]
		
		orr = re.split('\s+', lines[0])[6]
		if orr == '-':
			for id, s in trimseq.items():
				trimseq[id] = rev_comp(s)
			winner = rev_comp(winner)

			d_out = '%s_dna.fasta' % (locname)
			o = open(d_out, 'w')
			o.write('>dna\n%s\n' % winner)
			o.close()

			call = subprocess.Popen("exonerate -m protein2genome -q %s -t %s --showtargetgff TRUE --showalignment NO --showvulgar 0 -n 1 | grep 'cds'" % (p_out, d_out), shell=True, stdout=subprocess.PIPE)
			lines = [line.decode('utf-8').rstrip() for line in call.stdout]

		cds = []
		if len(lines) > 0:
			for line in lines:
				d = re.split('\s+', line)
				start = int(d[3])
				end = int(d[4])

				cds.append([start, end])
		
			# print out noncoding partition
			nc = []
			if cds[0][0] != 1:
				nc.append([1, cds[0][0] - 1])
			if cds[len(cds) - 1][1] < len(winner):
				nc.append([cds[len(cds) - 1][1] + 1, len(winner)])
			if len(cds) > 1:
				for i in range(0, len(cds) - 1):
					nc.append([cds[i][1] + 1, cds[i + 1][0] - 1])

			if len(nc) > 0:
				nc = ','.join(['%s-%s' % (x[0], x[1]) for x in nc])
				out.write('DNA, pNC=%s\n' % nc)

			for i in range(0,3):
				tmpcds = ','.join(['%s-%s\\3' % (x[0] + i, x[1]) for x in cds])
				out.write('DNA, p%s=%s\n' % (i + 1, tmpcds))
			out.close()

			os.remove(d_out)
			os.remove(p_out)
		else:
			out.write('DNA, p1=1-%s\n' % len(winner))
		'''
	else:
		out = open('%s_partitions.txt' % locname, 'w')
		out.write('DNA, p1=1-%s\n' % len(list(trimseq.values())[0]))
		out.close()

	'''
	new = re.sub('.fasta.aln', '.trim.aln.phy', aln)
	o = open(new, 'w')
	o.write('%s %s\n' % (len(trimseq), len(winner)))
	for id, s in trimseq.items():
		o.write('%s\t%s\n' % (id, s))
	o.close()
	'''