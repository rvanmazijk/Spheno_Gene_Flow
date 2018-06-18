import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", required=True, help='in')
parser.add_argument("-o", required=True, help='out.')
args = parser.parse_args()

f = open(args.i, 'r')
tree = f.next().lstrip()

tree = re.sub(':NA', '', tree)
tree = re.sub(r"\b[0-9]+(?:\.[0-9]+)?\b", '', tree)
tree = re.sub(":", '', tree)
tree = re.sub("NA;", ';', tree)
o = open(args.o, 'w')
o.write(tree)
o.close()
