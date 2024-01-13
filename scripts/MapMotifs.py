'''
Author: Tristan Achterberg and Anne de Jong
	
date: Sept 2023

python3 /data/sigma54/MapMotifs.py -sessiondir . -kmers query.Sig54_motifs.table -fna query.g2d.fna -out query.mapped_motifs.tab 


'''


import re
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='Add Regulon to UniProt')
parser.add_argument('-sessiondir', dest='sessiondir', help='Sessionfolder')
parser.add_argument('-kmers', dest='kmers', help='k-mers sequences from the prediction query.Sig54_motifs.table')
parser.add_argument('-fna', dest='fna', help='genome fna filename : query.regulon.tab')
parser.add_argument('-out', dest='out', help='output filename, default=query.mapped_motifs.tab')
parser.add_argument('--version', action='version', version='Anne de Jong, version 1.0, Sept 2023')
args = parser.parse_args()

fna_filename = args.sessiondir+'/'+args.fna

def readMultiFasta(fna_file):
    # all chars not A,C,T,G will be replaced by G, put sequence into dictionary
    fasta = {}
    key = ''
    with open(fna_file) as f:
        for line in f:
            if line[0] == ">":
                items = re.match(">(.*?)\s", line)
                if items:
                    key = items.group(1)
                    fasta[key] = ''
            else:
                if (key != ''):
                    fasta[key] += re.sub('[BDEFHIJKLMNOPQRSUVWXYZ]', 'G', line.strip().upper())
        return fasta


def calculate_best_subsequence(sequence, pwm1, pwm2, gap):
    max_score = -float('inf')
    best_subseq = ""
    for i in range(len(sequence) - len(pwm1['A']) - len(pwm2['A']) - gap + 1):
        subseq1 = sequence[i: i + len(pwm1['A'])]
        subseq2 = sequence[i + len(pwm1['A']) + gap: i + len(pwm1['A']) + gap + len(pwm2['A'])]
        score = sum(pwm1[base][j] for j, base in enumerate(subseq1)) + sum(
            pwm2[base][k] for k, base in enumerate(subseq2))
        if score > max_score:
            max_score = score
            best_subseq = sequence[
                          i: i + len(pwm1['A']) + gap + len(pwm2['A'])]  # includes the actual sequence during the gap
    return best_subseq


def elongate_subsequences(subsequences_file, pwm1, pwm2, gap_length, target_length=50):
    subsequences_with_locus_tags = read_sequences_from_file(subsequences_file)
    fasta = readMultiFasta(fna_filename)

    data = {
        'locus-tag': [],
        'kmer': [],
        'kmerStart': [],
        'kmerEnd': [],
        'sigma54': [],
        'motifStart': [],
        'motifEnd': [],
        'chrom': [],
        'posMotif': [],
        'posKmer': []
    }

    for locus_tag, subseq in subsequences_with_locus_tags:
        best_subseq = calculate_best_subsequence(subseq, pwm1, pwm2, gap_length)

        for seq_name, seq in fasta.items():
            rev_com_seq = reverse_complement(seq)
            # Searching in the reverse complement
            start_index_rev = rev_com_seq.find(best_subseq)

            if best_subseq in seq or start_index_rev != -1:
                if best_subseq in seq:
                    start_index = seq.find(best_subseq)
                else:
                    start_index = len(seq) - (start_index_rev + len(best_subseq))

                end_index = start_index + len(best_subseq)

                n_to_add = (target_length - len(best_subseq)) // 2
                n_to_add_right = n_to_add
                n_to_add_left = n_to_add
                if len(best_subseq) % 2 != 0:
                    n_to_add_left += 1

                # Adjustments for subsequences near sequence edges
                if start_index - n_to_add_left < 0:
                    n_to_add_right += n_to_add_left - start_index
                    n_to_add_left = start_index
                if end_index + n_to_add_right > len(seq):
                    n_to_add_left += n_to_add_right - (len(seq) - end_index)
                    n_to_add_right = len(seq) - end_index

                elongated_subseq = seq[start_index - n_to_add_left:end_index + n_to_add_right]

                data['locus-tag'].append(locus_tag)
                data['kmer'].append(elongated_subseq)
                data['sigma54'].append(best_subseq)
                data['chrom'].append(seq_name)
                data['motifStart'].append(start_index)
                data['motifEnd'].append(end_index)
                data['posMotif'].append((start_index, end_index))
                data['posKmer'].append((start_index - n_to_add_left, end_index + n_to_add_right))
                data['kmerStart'].append(start_index - n_to_add_left)
                data['kmerEnd'].append(end_index + n_to_add_right)

    df = pd.DataFrame(data)
    df = df.drop_duplicates()
    df = df.reset_index(drop=True)
    return df


def read_sequences_from_file(file_path):
    subsequences = []
    with open(file_path, 'r') as file:
        for line in file:
            locus_tag, sequence = line.strip().split()
            subsequences.append((locus_tag, sequence))
    return subsequences


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[base] for base in reversed(seq))


pwm = {
    'A': [0.000, 0.124, 0.000, 0.094, 0.279, 0.000, 0.282, 0.000, 0.000, 0.000, 0.555],
    'C': [0.155, 0.000, 0.000, 1.038, 0.070, 0.316, 0.000, 0.124, 0.000, 2.000, 0.123],
    'G': [0.078, 1.485, 1.485, 0.000, 0.070, 0.053, 0.000, 0.000, 2.000, 0.000, 0.123],
    'T': [0.776, 0.000, 0.124, 0.094, 0.035, 0.316, 0.939, 1.485, 0.000, 0.000, 0.000]
}

pwm1 = {base: scores[:6] for base, scores in pwm.items()}
pwm2 = {base: scores[6:] for base, scores in pwm.items()}

# From here need: sequences returned from prediction (sequences), original sequences from the organism, in form of .fnn or .fna (long_seqs)

gap_length = 5

#Example use:
#df = elongate_subsequences('CNN_dataset/query.Sig54_motifs (1).table', 'CNN_dataset/Borreliella_burgdorferi_B31_ASM868v2_genomic.fna', pwm1, pwm2, gap_length, target_length=50)
# Anne note: genome filename removed and set as global

df = elongate_subsequences(args.sessiondir+'/'+args.kmers, pwm1, pwm2, gap_length, target_length=50)
df.sort_values(by=['kmerStart']).to_csv(args.sessiondir+'/'+args.out, index = False, sep ='\t', header=True)

# Write motifs to Fasta file for WebLogo
fna=''
for index, row in df.iterrows():
	fna+='>'+row['locus-tag']+'\n'+row['sigma54']+'\n'
f = open(args.sessiondir+'/query.Sig54_motifs.fna', "w")
f.write(fna)
f.close()	

