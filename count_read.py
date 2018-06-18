import re
import glob
import commands

r1a = glob.glob("/nfs/turbo/lsa-rabosky/Lab/skink_ddrad/individual_reads/*R1*gz")
r1b = glob.glob("/nfs/turbo/lsa-rabosky/Lab/skink_ddrad/individual_reads/*.1.noRE.fq.gz")
r1 = r1a + r1b
print("file,number_seqs")
for r1file in r1:
	result = commands.getoutput('zcat %s | wc -l' % r1file)
	numseqs = int(result) / 4.0
	print("%s,%s" % (r1file, numseqs))
