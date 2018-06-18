import re
import os
import gzip
import subprocess
import argparse
import operator

parser = argparse.ArgumentParser()
parser.add_argument("-r1", "--read1", required=True, help='read 1 sample to prep.')
parser.add_argument("-r2", "--read2", required=True, help='read 2 sample to prep.')
parser.add_argument("-t1", "--trim1", help='how much to trim from read 1 start.', type=int, default=6)
parser.add_argument("-t2", "--trim2", help='how much to trim from read 2 start.', type=int, default=2)
args = parser.parse_args()

r1 = args.read1
r2 = args.read2
t1 = args.trim1
t2 = args.trim2

# https://github.com/lh3/readfq/blob/master/readfq.py
def readfq(fp): 
    # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def clean_file(read, trim, end):
        out = re.sub('.fq.gz', '.noRE.fq.gz', read)

        i = gzip.open(read, 'r')
	o = gzip.open(out, 'w')

	cuts = {}

	n, slen, qlen = 0, 0, 0
        for idname, seq, qual in readfq(i):
		idname = re.search('^(\S+)', idname).group(1)
		idname = re.sub('^.*XX:', '', idname)
		idname = re.sub(':', '_', idname)

		# two ways that reads are named
		idname = re.sub('/%s$' % end, '', idname)
		idname = re.sub('_%s$' % end, '', idname)

		cut = seq[0:trim]
		if cut not in cuts:
			cuts[cut] = 0
		cuts[cut] += 1
		seqtrim = seq[trim:]
		qualtrim = qual[trim:]
		o.write('@%s\n%s\n+\n%s\n' % (idname, seqtrim, qualtrim))
	i.close()
	o.close()
	
	out2 = re.sub('.fq.gz', '_cuts', read)
	o2 = open(out2, 'w')
	sumc = sum([cuts[x] for x in cuts])
	sorted_c = sorted(cuts.items(), key=operator.itemgetter(1))
	for x in sorted_c:
		per = x[1] / float(sumc)
		o2.write('%s\t%s\t%.3f\n' % (x[0], x[1], per))
	o2.close()

clean_file(r1, t1, '1')
clean_file(r2, t2, '2')
