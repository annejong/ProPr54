'''
Author: Anne de Jong 
date: Sept 2023

python3  /data/sigma54/MergeRegulonMotifs.py -s /tmpdrive/sigma54/Sigma54/86.90.198.76.ptrlqdbml388olslisso49428b.5065 -annotation query.g2d.table -regulon query.regulon.tab -motifs query.mapped_motifs.tab -out query.merged_regulon_motifs.tab

'''

import argparse
import re
import pandas as pd


# ---------------------------------------------------------------- parse parameters -------------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(description='test')
parser.add_argument('-sessiondir', dest='sessiondir', help='Full path sessiondir')
parser.add_argument('-annotation', dest='annotation', help='default=query.g2d.table')
parser.add_argument('-regulon', dest='regulon', help='default=query.regulon.tab')
parser.add_argument('-motifs', dest='motifs', help='default=query.mapped_motifs.tab')
parser.add_argument('-out', dest='out', help='default=query.merged_regulon_motifs.tab')
parser.add_argument('--version', action='version', version='Anne de Jong, version 1.0, Sept 2023')
args = parser.parse_args()

gff_header = ["chrom","db", "type", "start", "end", "name", "strand", "score", "description"]
convert_dict = { 'start': int, 'end': int }


def read_fasta(filename):
	results = {}
	lines = open(filename, "r").read().split('\n')
	key=''
	for line in lines:
		if re.search("^>", line):
			if key != '': results[key] = seq
			items = line.split()
			key = items[0][1:]
			seq=''
		else:
			seq += line
	if key != '': results[key] = seq  # add the last record
	return results

ANNOTATION 	= pd.read_csv(args.sessiondir+'/'+args.annotation, header=0,  comment='#',sep='\t')



REGULON = pd.read_csv(args.sessiondir+'/'+args.regulon, header=0, comment='#',sep='\t')
REGULON.rename(columns={'query': 'locus-tag'}, inplace=True)
MOTIFS 	= pd.read_csv(args.sessiondir+'/'+args.motifs, header=0, comment='#',sep='\t')

merged_header = ["chrom", "locus-tag","sigma54", "motifStart", "motifEnd"]
MERGED = pd.merge(REGULON.astype(str), MOTIFS.astype(str), left_on='locus-tag', right_on='locus-tag')[merged_header]
MERGED.sort_values(by=['chrom','locus-tag']).to_csv(args.sessiondir+'/'+args.out, index = False, sep ='\t', header=True)


# add distance of motif to start of the gene
MOTIFS_ANNOTATED = pd.merge(MOTIFS, ANNOTATION, left_on='locus-tag', right_on='locus_tag')
MOTIFS_ANNOTATED['distance'] = MOTIFS_ANNOTATED.apply(lambda row: row['start'] - row['motifEnd'] if row['strand'] == '+' 
                        else row['motifStart'] - row['end'], axis=1)
MOTIFS_ANNOTATED = MOTIFS_ANNOTATED.rename(columns={'start': 'geneStart'}) 
MOTIFS_ANNOTATED = MOTIFS_ANNOTATED.rename(columns={'end': 'geneEnd'}) 
MOTIFS_ANNOTATED = MOTIFS_ANNOTATED.rename(columns={'strand': 'geneStrand'}) 

# Concatenate Motifs and Regulon to one comprehensive table
OUTER_MERGED = pd.merge(REGULON, MOTIFS_ANNOTATED, how="outer", on=["locus-tag"])
header=["locus-tag","uniprotID","bitscore","sigma54","motifStart","motifEnd","chrom_x","geneStart","geneEnd","geneStrand","distance"]
OUTER_MERGED[header].sort_values(by=['chrom_x','motifStart']).to_csv(args.sessiondir+'/query.concat_regulon_motifs.tab', float_format='%.0f', index = False, sep ='\t', header=True)


# Write MERGED motifs to Fasta file for WebLogo
fna=''
for index, row in MERGED.iterrows():
	fna+='>'+row['locus-tag']+'\n'+row['sigma54']+'\n'
f = open(args.sessiondir+'/query.merged.sigma54_motifs.fna', "w")
f.write(fna)
f.close()	










