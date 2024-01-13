'''
Author: Anne de Jong 
date: Sept 2023

python3  /data/sigma54/GFF2FAA.py -s /tmpdrive/sigma54/Sigma54/86.103.33.131.gispng69rjvvk7tr6qkvvahhb4.8397 -fna query.fna -gff query.gff

'''

import argparse
import re
import pandas as pd


# ---------------------------------------------------------------- parse parameters -------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='test')
parser.add_argument('-s', dest='sessiondir', help='Full path sessiondir')
parser.add_argument('-fna', dest='fna', help='Full path and name of fna file')
parser.add_argument('-gff', dest='gff', help='Full path and name of gff file')
parser.add_argument('--version', action='version', version='Anne de Jong, version 1.0, Sept 2023')
args = parser.parse_args()



def read_fasta(filename):
	results = {}
	lines = open(filename, "r").read().split('\n')
	key=''
	for line in lines:
		if re.search("^>", line):
			if key != '': results[key] = seq
			items = line.split()
			key = items[0][1:]
			print(key)
			seq=''
		else:
			seq += line
	if key != '': results[key] = seq  # add the last record
	return results

# load fna
fna = read_fasta(args.sessiondir+'/'+args.fna)



#  load gff  
gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]
GFF = pd.read_csv(args.sessiondir+'/'+args.gff, header=None,  comment='#',sep='\t', names=gff_header)
convert_dict = { 'start': int, 'end': int }
GFF = GFF.astype(convert_dict)  # be sure that start and end are integers
CDS = GFF.loc[GFF['type'] == 'CDS']

for index, row in CDS.iterrows():
	chrom = row['chrom']
	start = row['start']-1
	end = row['end']
	seq = fna[chrom][start:end]
	print(seq)


#gff_genes.sort_values(by=['start']).to_csv(args.sessiondir+'/'+args.gff, index = False, sep ='\t', columns=gff_header, header=None)



