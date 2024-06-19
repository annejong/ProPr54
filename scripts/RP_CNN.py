#!/usr/bin/env python
# coding: utf-8
# %%
import pandas as pd
import numpy as np
import argparse
import re
import random
import screed

import tensorflow as tf
from sklearn.preprocessing import LabelEncoder, OneHotEncoder

parser = argparse.ArgumentParser(description='Train the model')
parser.add_argument('-promlen', dest='PromLen', help='Length of Promoter', nargs='?', default=50)
parser.add_argument('-padsize', dest='padSize', help='Size of padding on both sides', nargs='?', default=12)
parser.add_argument('-sessiondir', dest='sessiondir', help='Session Dir', nargs='?', default='.')
parser.add_argument('-modelName', dest='modelName', help='Name of the model, see models.py', nargs='?',
                    default='model_lstm')
parser.add_argument('-batchsize', dest='batchSize', help='Batch size', nargs='?', default=64)
parser.add_argument('-epochs', dest='epochs', help='Number of training epochs', nargs='?', default=40)
args, unknown = parser.parse_known_args()


# Create functions which get a fasta file, one hot encode it and make a query set out of the DNA sequences.
def getCleanSeq(fna_file):
    # get DNA and replace N or n by G or g to prevent error in training; G is considered as most save replacment
    DNA = ''
    with open(fna_file) as lines:
        for line in lines:
            if line[0] != ">":
                DNA += line.strip()

    DNA = re.sub('[BDEFHIJKLMNOPQRSUVWXYZ]', 'G', DNA.upper())
    return DNA


def Anne_one_hot_encode(seq):
    # One hot encode your sequences for the CNN
    mapping = dict(zip("ACGT", range(4)))
    seq2 = [mapping[i] for i in seq]
    return np.eye(4)[seq2]


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


def makeQuerySet(DNA, window_shift):
    query_set = []
    for i in range(0, len(DNA) - args.PromLen - args.prime, window_shift):
        query_set.append(DNA[i:i + args.PromLen + args.prime])
    return query_set


def DNA(length):
    return ''.join(random.choice('CGTA') for _ in range(length))


def Anne_getPercentage(genome):
    # return a list of percentages of A,C,G,T
    percentageList = []
    for base in ['A', 'C', 'G', 'T']:    percentageList.append(round((genome.count(base) / len(genome)) * 100))
    return (percentageList)


def makeManyRandoms(count, ACGTcontent, length):
    # make 100 bases with proper GC content
    dnaset = "A" * ACGTcontent[0] + "C" * ACGTcontent[1] + "G" * ACGTcontent[2] + "T" * ACGTcontent[3]
    RandomSeq = []
    for i in range(0, count):    RandomSeq.append(randomseq(dnaset, length))
    return RandomSeq


def randomseq(dnaset, length):
    result = ''
    for i in range(length):
        r = random.randint(0, len(dnaset) - 1)
        result += dnaset[r]
    return result


def makeRandoms(count, length, sequence_file):
    sequence = getCleanSeq(sequence_file)
    percentage = Anne_getPercentage(sequence)
    randoms = makeManyRandoms(count, percentage, length)
    return randoms


def Elongate_sequence(sequence_list, ACGTcontent):
    # Elongate the your sequences to match the args.PromLen with a defined GC%
    for idx, item in enumerate(sequence_list):
        if len(item) < args.PromLen:
            def randomseq1(dnaset):
                result = ''
                for i in range(args.PromLen - len(item)):
                    r = random.randint(0, len(dnaset) - 1)
                    result += dnaset[r]
                return result

            RandomSeq = []
            dnaset = "A" * ACGTcontent[0] + "C" * ACGTcontent[1] + "G" * ACGTcontent[2] + "T" * ACGTcontent[3]
            for i in range(0, len(sequence_list)):    RandomSeq.append(randomseq1(dnaset))
            for x in RandomSeq:
                sequence_list[idx] = x + item
    return sequence_list


def read_gff(gff_file):
    # Puts the gff into a pandas dataframe
    with open(gff_file) as f:
        gff_header = ["genome", "name", "type", "start", "end", "score", "strand", "phase", "attributes"]
        gff = pd.read_csv(f, sep='\t', comment="#", lineterminator='\n', names=gff_header)
    return gff


def write_log(S):
    # write to console and logfile
    print(S)
    f = open(args.sessiondir + '/' + "log_ppp_{model_name}.log".format(model_name=args.modelName), "a")
    f.write(S + '\n')
    f.close()


def Validation_report(features, labels, name):
    # F1 = A measure that combines precision and recall is the harmonic mean of precision and recall, the traditional F-measure or balanced F-score:
    # F1 = 2 * (precision * recall) / (precision + recall)
    # precision = TP / TP + FP
    # recall = TP / TP + FN
    predicted_labels = model.predict(features)
    cm = confusion_matrix(np.argmax(labels, axis=1), np.argmax(predicted_labels, axis=1))
    TP = cm[1][1]
    FP = cm[0][1]
    FN = cm[1][0]
    precision = TP / (TP + FP)
    recall = TP / (TP + FN)
    F1 = 2 * (precision * recall) / (precision + recall)

    write_log('Confusion matrix of ' + name + '\n' + str(cm) + '| ' + args.sessiondir)
    write_log('Precision  TP / TP + FP                              ' + str(round(precision, 2)))
    write_log('Recall     TP / TP + FN                              ' + str(round(recall, 2)))
    write_log('F1 = 2 * (precision * recall) / (precision + recall) ' + str(round(F1, 2)))

    # Confusion matrix info is in the log, but one extra clean file:
    f = open(args.sessiondir + '/' + "Confusion_matrix.txt", "a")
    f.write(args.sessiondir + '\n')
    f.write('Confusion matrix of ' + name + '\n' + str(cm) + '\n')
    f.write('Precision  TP / TP + FP                              ' + str(round(precision, 2)) + '\n')
    f.write('Recall     TP / TP + FN                              ' + str(round(recall, 2)) + '\n')
    f.write('F1 = 2 * (precision * recall) / (precision + recall) ' + str(round(F1, 2)) + '\n\n')
    f.write('Precision\tRecall\tF1\n')
    f.write(str(round(precision, 2)) + '\t' + str(round(recall, 2)) + '\t' + str(round(F1, 2)) + '\n\n')
    f.close()


def prepare_sequence_file(sequence_file):
    # Makes everything uppercase, removes newlines and whitespaces
    sequence = []
    regex = r'>'
    replace_char = ['\n', ' ']
    with open(sequence_file) as f:
        for line in f:
            if regex not in line:
                sequence.append(line)
    for i in replace_char:
        for idx, item in enumerate(sequence):
            sequence[idx] = item.replace(i, '')
    for idx, item in enumerate(sequence):
        sequence[idx] = item.upper()
    sequence = list(filter(None, sequence))
    return sequence


def read_kmers_from_file(filename, ksize):
    # a library for reading in FASTA/FASTQ, ksize is length of window/kmer
    all_kmers = []
    for record in screed.open(filename):
        sequence = record.sequence

        kmers = []
        n_kmers = len(sequence) - ksize + 1

        for i in range(n_kmers):
            kmer = sequence[i:i + ksize]
            kmers.append(kmer)

        all_kmers += kmers

    return all_kmers


def makeFasta(name, regulon_list, output_file):
    # Makes a FASTA file from a list of sequences so you can input it into MEME e.g.
    fasta_key = []
    for number in range(len(regulon_list)):
        fasta_key.append(name + str(number))
    with open(output_file, 'w') as f:
        for i in range(len(fasta_key)):
            f.write('>' + fasta_key[i] + '\n' + regulon_list[i] + '\n')


def writeFasta(FASTA_dict, output_file):
    with open(output_file, 'w') as f:
        for key, value in FASTA_dict.items():
            f.write(f'>{key}\n')  # write the sequence identifier
            f.write(f'{value}\n')  # write the sequence data


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])


def Elongate_fna(FNA_file, motifs):
    # Takes an fna file containing whole genome, makes reverse of it and uses the motifs to elongate the sequences on both ends
    genome_fw = readMultiFasta(FNA_file)
    genome_fw = next(iter(genome_fw.values()))
    genome_rv = reverse_complement(genome_fw)

    elongated_sequences = []
    elongate_w_zeros = []

    for sequence in motifs:
        n_to_add = (args.PromLen + args.padSize - len(sequence)) // 2
        if len(sequence) % 2 != 0:
            n_to_add += 1
        start_index = genome_fw.find(sequence)
        end_index = start_index + len(sequence)
        elongated_sequence_fw = genome_fw[start_index - n_to_add:end_index + n_to_add]
        elongated_sequences.append(elongated_sequence_fw)

        start_index_rv = genome_rv.find(sequence)
        end_index_rv = start_index_rv + len(sequence)
        elongated_sequence_rv = genome_rv[start_index_rv - n_to_add:end_index_rv + n_to_add]
        elongated_sequences.append(elongated_sequence_rv)

        if not any(sequence in element for element in elongated_sequences):
            elongate_w_zeros.append(sequence)

    elongated_sequences = elongated_sequences + elongate_w_zeros

    return [item for item in elongated_sequences if item != '']


def Sliding_window_aug(sequence_list, step_size):
    kmers = []
    sequence_zero = []
    for sequence in sequence_list:
        if '0' not in sequence and len(sequence) == args.padSize + args.PromLen:
            n_iterations = (len(sequence) - args.PromLen) // step_size + 1

            for i in range(n_iterations):
                start_index = i * step_size
                end_index = start_index + args.PromLen
                kmer = sequence[start_index:end_index]
                kmers.append(kmer)
        else:
            sequence_zero.append(sequence)
    kmers = kmers + sequence_zero
    return kmers


def slide_window(sequence_list, window_size, step_size):
    windows = []
    sequences_short = []
    for sequence in sequence_list:
        if len(sequence) > window_size:
            seq_len = len(sequence)
            for i in range(0, seq_len - window_size + 1, step_size):
                window = sequence[i:i + window_size]
                windows.append(window)
    return windows


def mutate_sequence(sequence_list):
    sequences_final = []
    new_sequences = []
    for sequence in sequence_list:
        if len(sequence) > 52:
            # Looking at the first 3 NT
            if sequence[0] == 'T':
                new_sequences.append('C' + sequence[1:])
            elif sequence[0] == 'C':
                new_sequences.append('T' + sequence[1:])
            elif sequence[0] == 'G':
                new_sequences.append('A' + sequence[1:])
            elif sequence[0] == 'A':
                new_sequences.append('G' + sequence[1:])

            if sequence[1] == 'T':
                new_sequences.append(sequence[:1] + 'C' + sequence[2:])
            elif sequence[1] == 'C':
                new_sequences.append(sequence[:1] + 'T' + sequence[2:])
            elif sequence[1] == 'G':
                new_sequences.append(sequence[:1] + 'A' + sequence[2:])
            elif sequence[1] == 'A':
                new_sequences.append(sequence[:1] + 'G' + sequence[2:])

            if sequence[2] == 'T':
                new_sequences.append(sequence[:2] + 'C' + sequence[3:])
            elif sequence[2] == 'C':
                new_sequences.append(sequence[:2] + 'T' + sequence[3:])
            elif sequence[2] == 'G':
                new_sequences.append(sequence[:2] + 'A' + sequence[3:])
            elif sequence[2] == 'A':
                new_sequences.append(sequence[:2] + 'G' + sequence[3:])

            # Looking at the last 3 NT
            if sequence[-1] == 'T':
                new_sequences.append(sequence[:-1] + 'C')
            elif sequence[-1] == 'C':
                new_sequences.append(sequence[:-1] + 'T')
            elif sequence[-1] == 'G':
                new_sequences.append(sequence[:-1] + 'A')
            elif sequence[-1] == 'A':
                new_sequences.append(sequence[:-1] + 'G')

            if sequence[-2] == 'T':
                new_sequences.append(sequence[:-2] + 'C' + sequence[-1:])
            elif sequence[-2] == 'C':
                new_sequences.append(sequence[:-2] + 'T' + sequence[-1:])
            elif sequence[-2] == 'G':
                new_sequences.append(sequence[:-2] + 'A' + sequence[-1:])
            elif sequence[-2] == 'A':
                new_sequences.append(sequence[:-2] + 'G' + sequence[-1:])

            if sequence[-3] == 'T':
                new_sequences.append(sequence[:-3] + 'C' + sequence[-2:])
            elif sequence[-3] == 'C':
                new_sequences.append(sequence[:-3] + 'T' + sequence[-2:])
            elif sequence[-3] == 'G':
                new_sequences.append(sequence[:-3] + 'A' + sequence[-2:])
            elif sequence[-3] == 'A':
                new_sequences.append(sequence[:-3] + 'G' + sequence[-2:])

    sequence_list += new_sequences
    return sequence_list


def Anne_decode(numpy_array):
    int_to_nuc = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    decoded_sequences = []
    for one_hot_seq in numpy_array:
        # Get the index of the max value along the last axis (i.e. for each nucleotide)
        max_indices = np.argmax(one_hot_seq, axis=-1)
        # Map the integers back to nucleotides using the dictionary
        decoded_seq = ''.join([int_to_nuc[i] for i in max_indices])
        # Add the decoded sequence to the list
        decoded_sequences.append(decoded_seq)
    return decoded_sequences


# Preparing regulon from L interrogans
df_L_int = pd.read_excel('CNN_dataset/Feature_Table.xlsx', 'L_int_R')
L_int_R = df_L_int['Matching sequence'].tolist()
L_int_R = Elongate_fna('CNN_dataset/Leptospira_interrogans_serovar_Manilae_UP-MMC-NIID_HP_ASM104765v1_genomic.fna',
                       L_int_R)

# This is preperation of promoter sequences from different bacteria containing the sigma54 motif, but as I already had sequences
# from C. difficile from a different paper, the list will exclude them
df_C_mixed = pd.read_excel('CNN_dataset/C_mixed_R_1.xlsx', 'Table011 (Page 12)')
regex_C_diff = '^(?:(?!CD).)*$'
df_C_preped = df_C_mixed[df_C_mixed.First_gene_ID.str.match(regex_C_diff, na=False)]

# Makes the regulon list for L. monocytogenes
L_mon_R = prepare_sequence_file('CNN_dataset/L_mon_R.txt')
L_mon_R = [sequence.replace('(N)', '.{4}') for sequence in L_mon_R]
L_mon_G = readMultiFasta('CNN_dataset/Listeria_monocytogenes_10403S_ASM16869v2_genomic.fna')
L_mon_G = next(iter(L_mon_G.values()))
L_mon_G_rv = reverse_complement(L_mon_G)

L_mon_M = []
for sequence in L_mon_R:
    L_mon_M.extend(re.findall(sequence, L_mon_G))

for sequence in L_mon_R:
    L_mon_M.extend(re.findall(sequence, L_mon_G_rv))

L_mon_R = Elongate_fna('CNN_dataset/Listeria_monocytogenes_10403S_ASM16869v2_genomic.fna', L_mon_M)

# Makes regulon from V. cholera
df_V_chl_R = pd.read_excel('CNN_dataset/V_chl_R.xlsx')
df_V_chl_IG = df_V_chl_R[df_V_chl_R['Intergenic\nregion (IG)'] == 'IG']
V_chl_R = df_V_chl_IG['Binding motif'].tolist()
V_chl_R = Elongate_fna('CNN_dataset/GCF_008369605.1_ASM836960v1_genomic.fna', V_chl_R)

# Makes regulon from C. trachatomis

df_C_tra_R = pd.read_excel('CNN_dataset/C_tra_R.xlsx', 'C_tra_R')
C_tra_R = df_C_tra_R['Column3'].tolist()

for idx, item in enumerate(C_tra_R):
    C_tra_R[idx] = re.sub(r'[^ACTG234{}]', '', item)
for idx, item in enumerate(C_tra_R):
    C_tra_R[idx] = item.replace('{', '.{')

C_tra_G = readMultiFasta('CNN_dataset/Chlamydia_trachomatis_434_Bu_ASM6858v1_genomic.fna')
C_tra_G = next(iter(C_tra_G.values()))
C_tra_G_rv = reverse_complement(C_tra_G)
C_tra_M = []
for sequence in C_tra_R:
    C_tra_M.extend(re.findall(sequence, C_tra_G))
for sequence in C_tra_R:
    C_tra_M.extend(re.findall(sequence, C_tra_G_rv))
C_tra_R = Elongate_fna('CNN_dataset/Chlamydia_trachomatis_434_Bu_ASM6858v1_genomic.fna', C_tra_M)

# Makes regulon from V. parahaemolyticus

df_V_par_R = pd.read_excel('CNN_dataset/V_par_R.xlsx')
V_par_R = df_V_par_R['Binding motif'].dropna().tolist()
V_par_R = Elongate_fna(
    'CNN_dataset/Vibrio_parahaemolyticus_RIMD_2210633_O3_K6_substr_RIMD_2210633_ASM19609v1_genomic.fna', V_par_R)

# Makes regulon of L. plantarum

df_L_pla_R = pd.read_excel('CNN_dataset/L_pla_R.xlsx')
L_pla_R = df_L_pla_R['Sequence'].dropna().tolist()
for idx, sequence in enumerate(L_pla_R):
    L_pla_R[idx] = sequence.replace('\xa0', '')
L_pla_R = Elongate_fna('CNN_dataset/Lactobacillus_plantarum_WCFS1_NCIMB8826_ASM20385v3_genomic.fna', L_pla_R)

# Makes regulon of B. subtilis

B_sub_R = prepare_sequence_file('CNN_dataset/B_sub_R.txt')
B_sub_R = Elongate_fna('CNN_dataset/Bacillus_subtilis_subsp_subtilis_str_168_ASM904v1_genomic.fna', B_sub_R)

# Makes regulon of P. putida

P_put_R = prepare_sequence_file('CNN_dataset/P_put_R.txt')
P_put_R = Elongate_fna('CNN_dataset/Pseudomonas_putida_KT2440_ASM756v2_genomic.fna', P_put_R)

# Makes regulon of B. cereus

B_cer_R = prepare_sequence_file('CNN_dataset/B_cer_R.txt')
B_cer_R = list(set(B_cer_R))
B_cer_R = Elongate_fna('CNN_dataset/Bacillus_cereus_ATCC_14579_ASM782v1_genomic.fna', B_cer_R)

E_coli_NR = prepare_sequence_file('CNN_dataset/E_coli_NR.txt')
C_diff_NR = prepare_sequence_file('CNN_dataset/C_diff_NR.txt')
compilation_R = prepare_sequence_file('CNN_dataset/Mixed_R.txt')
S_typh_R = prepare_sequence_file('CNN_dataset/S_typh_R.txt')
S_typh_R = Elongate_fna(
    'CNN_dataset/Salmonella_enterica_subsp_enterica_serovar_Typhimurium_str_14028S_ASM2216v1_genomic.fna', S_typh_R)
Y_ptb_R = prepare_sequence_file('CNN_dataset/Y_ptb_R.txt')
Y_ptb_R = Elongate_fna('CNN_dataset/Yersinia_pseudotuberculosis_YPIII_ASM1946v1_genomic.fna', Y_ptb_R)
C_mixed_R = prepare_sequence_file('CNN_dataset/C_mixed_R.txt')

# Makes E_coli_R 62 in length
E_coli_R = prepare_sequence_file('CNN_dataset/E_coli_R.txt')
E_coli_R_2 = prepare_sequence_file('CNN_dataset/E_coli_R_2.txt')
# It is currently 81 in length with the motif being in the middle of the sequence
for idx, sequence in enumerate(E_coli_R):
    E_coli_R[idx] = sequence.replace(sequence, sequence[10:-9])

# Determine which ones to add to left, and which ones to add to right

E_coli_R_2_right = E_coli_R_2[:3]
E_coli_R_2_left = E_coli_R_2[-2:]
for sequence in E_coli_R_2_right:
    E_coli_R.append(sequence.replace(sequence, sequence[-62:]))
for sequence in E_coli_R_2_left:
    E_coli_R.append(sequence.replace(sequence, sequence[:62]))

# Removing the regulon sequences from the non regulon sequences of C. difficile
df_Cd_experiment = pd.read_excel('BLAST_MEME/C_difficile_experimental/Table_2.xlsx')
new_header = df_Cd_experiment.iloc[1]
df_Cd_experiment = df_Cd_experiment[2:]
df_Cd_experiment.columns = new_header

# Makes the C_diff_R using the excel from the paper, and gets a list of the genes controlled by the promoters under their CD name
df_Cd_Sigma54 = df_Cd_experiment.loc[df_Cd_experiment['Sigma factor'].str.contains('SigL|Sig54', na=False)]
C_diff_R = df_Cd_Sigma54['100 bases upstream of  +1'].tolist()

# Make C_diff_R only the last 62 of the sequence
for idx, item in enumerate(C_diff_R):
    C_diff_R[idx] = item.replace(item, item[- args.PromLen - args.padSize:])

# Results from BLASTp of the E. coli regulon are removed from each of the non-regulon sequences

df_CT_Blast = pd.read_csv('BLAST_MEME/Protein_Ortho/BLAST_Ecoli_Ctrac.txt', sep='\t', header=None)
df_CD_Blast = pd.read_csv('BLAST_MEME/Protein_Ortho/BLAST_Ecoli_Cdiff.txt', sep='\t', header=None)
df_ST_Blast = pd.read_csv('BLAST_MEME/Protein_Ortho/BLAST_Ecoli_Styph.txt', sep='\t', header=None)
df_PP_Blast = pd.read_csv('BLAST_MEME/Protein_Ortho/BLAST_Ecoli_Pput.txt', sep='\t', header=None)
df_PF_Blast = pd.read_csv('BLAST_MEME/Protein_Ortho/BLAST_Ecoli_Pflu.txt', sep='\t', header=None)
CT_Blast = df_CT_Blast[1].tolist()
CD_Blast = df_CD_Blast[1].tolist()
ST_Blast = df_ST_Blast[1].tolist()
PP_Blast = df_PP_Blast[1].tolist()
PF_Blast = df_PF_Blast[1].tolist()


# Removes regulon from non-regulon of P.putida

P_put_R_genes = PP_Blast
P_put_motifs = prepare_sequence_file('CNN_dataset/P_put_R.txt')
P_put_NR_dict = readMultiFasta('CNN_dataset/Pseudomonas_putida_KT2440_ASM756v2_genomic.g2d.intergenic.ffn')
P_put_NR_clean = dict()
for key, val in P_put_NR_dict.items():
    if not any(ele in key for ele in P_put_R_genes):
        P_put_NR_clean[key] = val

for key, value in list(P_put_NR_clean.items()):
    if any(substring in value for substring in P_put_motifs):
        del P_put_NR_clean[key]
P_put_NR = list(P_put_NR_clean.values())


# Removes the regulon from the non-regulon of C. difficile
C_diff_R_genes = df_Cd_Sigma54['Target Gene'].tolist()
C_diff_R_genes += CD_Blast
C_diff_NR_dict = readMultiFasta(
    'BLAST_MEME/C_difficile_experimental/Clostridioides_difficile_630_ASM920v2_genomic.g2d.intergenic.ffn')
C_diff_NR_clean = dict()
for key, val in C_diff_NR_dict.items():
    if not any(ele in key for ele in C_diff_R_genes):
        C_diff_NR_clean[key] = val
C_diff_NR = list(C_diff_NR_clean.values())

# Removes regulon from non-regulon of S_typh
S_typh_R_genes = pd.read_excel('CNN_dataset/S_typh_R.xlsx')
S_typh_R_genes = S_typh_R_genes['Associated 14028s ORFb'].dropna().tolist()
S_typh_R_genes = re.findall(r'STM14_\d{4}', ''.join(S_typh_R_genes))
S_typh_R_genes += ST_Blast
S_typh_NR_dict = readMultiFasta(
    'CNN_dataset/Salmonella_enterica_subsp_enterica_serovar_Typhimurium_str_14028S_ASM2216v1_genomic.g2d.intergenic.ffn')
S_typh_NR_clean = dict()
for key, val in S_typh_NR_dict.items():
    if not any(ele in key for ele in S_typh_R_genes):
        S_typh_NR_clean[key] = val
S_typh_NR = list(S_typh_NR_clean.values())

# Removes regulon from non-regulon of C_tra
C_tra_R_genes = pd.read_excel('CNN_dataset/C_tra_R.xlsx', 'C_tra_R_genes')
C_tra_R_genes = C_tra_R_genes['Locus tag'].dropna().tolist()
C_tra_R_genes += CT_Blast
C_tra_NR_dict = readMultiFasta('CNN_dataset/Chlamydia_trachomatis_434_Bu_ASM6858v1_genomic.g2d.intergenic.ffn')
C_tra_NR_clean = dict()
for key, val in C_tra_NR_dict.items():
    if not any(ele in key for ele in C_tra_R_genes):
        C_tra_NR_clean[key] = val
C_tra_NR = list(C_tra_NR_clean.values())

# Removes regulon from non-regulon of E_fae

E_fae_gff = read_gff('CNN_dataset/Enterococcus_faecalis_V583_ASM778v1_genomic.gff')
df_Loci = pd.read_excel('CNN_dataset/E_fae_Loci.xlsx')
df_Loci.dropna(inplace=True)
E_fae_loci = df_Loci['Former Locus Tag'].tolist()
mask = E_fae_gff['attributes'].str.contains('|'.join(E_fae_loci), case=False)
df_E_fae_reg = E_fae_gff[mask]
E_fae_R_genes = df_E_fae_reg['attributes'].tolist()
for idx, item in enumerate(E_fae_R_genes):
    E_fae_R_genes[idx] = re.findall(r'locus_tag=EF_\d+', item)
E_fae_R_genes = [item for sublist in E_fae_R_genes for item in sublist]
for idx, item in enumerate(E_fae_R_genes):
    E_fae_R_genes[idx] = item.replace('locus_tag=', '')
E_fae_NR_dict = readMultiFasta('CNN_dataset/Enterococcus_faecalis_V583_ASM778v1_genomic.g2d.intergenic.ffn')
E_fae_NR_clean = dict()
for key, val in E_fae_NR_dict.items():
    if not any(ele in key for ele in E_fae_R_genes):
        E_fae_NR_clean[key] = val
E_fae_NR = list(E_fae_NR_clean.values())

# Removes regulon from non-regulon of P_flu

df_P_flu = pd.read_excel('Third_party_predict/fmicb-12-641844-t002.xlsx')
P_flu_R_pre = df_P_flu['Locus tag\' log2 (fold change)/ Function description '].dropna().tolist()
P_flu_genes = []
for i in P_flu_R_pre:
    if 'RS' in i:
        P_flu_genes.append(i)
for idx, item in enumerate(P_flu_genes):
    P_flu_genes[idx] = item.replace('O', '0')
for idx, item in enumerate(P_flu_genes):
    P_flu_genes[idx] = item.replace('-', '')
P_flu_genes = ''.join(P_flu_genes)
P_flu_genes = re.findall(r'RS\d{5}', P_flu_genes)

P_flu_genes += PF_Blast

P_flu_NR_dict = readMultiFasta(
    'CNN_dataset/Pseudomonas_fluorescens_UK4_ASM73042v1_genomic.g2d.intergenic.ffn')
P_flu_NR_clean = dict()
for key, val in P_flu_NR_dict.items():
    if not any(ele in key for ele in P_flu_genes):
        P_flu_NR_clean[key] = val
P_flu_NR = list(P_flu_NR_clean.values())


# Make sliding window for the non-regulon sequences so that all parts of the intergenic regions are captured

C_tra_window = slide_window(C_tra_NR, args.PromLen, args.PromLen)
S_typh_window = slide_window(S_typh_NR, args.PromLen, args.PromLen)
E_coli_window = slide_window(E_coli_NR, args.PromLen, args.PromLen)
C_diff_window = slide_window(C_diff_NR, args.PromLen, args.PromLen)
E_fae_window = slide_window(E_fae_NR, args.PromLen, args.PromLen)
P_put_window = slide_window(P_put_NR, args.PromLen, args.PromLen)
P_flu_window = slide_window(P_flu_NR, args.PromLen, args.PromLen)

# Concatenates the regulons
regulon_pre = C_diff_R + V_chl_R + S_typh_R+ P_put_R + compilation_R + C_mixed_R + L_int_R + L_mon_R + V_par_R + L_pla_R + B_sub_R + C_tra_R + E_coli_R + Y_ptb_R
non_regulon = C_diff_window + S_typh_window + E_fae_window + P_flu_window + C_tra_window + E_coli_window + P_put_window


# Left out motifs:

# + B_cer_R

# Removes sequences which did not contain a motif after MEME input

regulon_clean = []
regulon = []
with open('BLAST_MEME/MEME/regulon_motifs.txt') as f:
    for line in f:
        regulon_clean = f.read().splitlines()
        regulon_clean = [element.strip() for element in regulon_clean]
for i in regulon_pre:
    for sub in regulon_clean:
        if sub in i:
            regulon.append(i)
            break
print(len(regulon))


makeFasta('regulon', regulon, 'CNN_dataset/regulon.txt')

makeFasta('non_regulon', non_regulon, 'CNN_dataset/non_regulon.txt')
#makeFasta('R', regulon, 'CNN_dataset/regulon.fasta')

# Data augmentation via making the sliding window, adding in the reverse complement, and adding base flipping

for idx, item in enumerate(regulon):
    regulon[idx] = re.sub('[BDEFHIJKLMNOPQRSUVWXYZ]', 'G', item)

# Include substitutions from C <-> T in the first 3 nucleotides of the sequence
#regulon = mutate_sequence(regulon)

# Reverse complement
#regulon_rv = []
#for i in regulon:
#    regulon_rv.append(reverse_complement(i))
#regulon = regulon_rv + regulon

# Sliding window
regulon = Sliding_window_aug(regulon, 3)

# Makes sure all the sequences are the same length and removes any non ACTG with a G
for idx, item in enumerate(regulon):
    regulon[idx] = item.replace(item, item[-args.PromLen:])
for i in non_regulon.copy():
    if len(i) < args.PromLen:
        non_regulon.remove(i)
for idx, item in enumerate(non_regulon):
    non_regulon[idx] = item.replace(item, item[-args.PromLen:])
for idx, item in enumerate(regulon):
    regulon[idx] = re.sub('[^ACTG0]', 'G', item.strip().upper())
print(len(regulon))
print(len(non_regulon))

#makeFasta('pos', regulon, 'Third_party_predict/WEKA/training_pos.fasta')
#makeFasta('neg', non_regulon, 'Third_party_predict/WEKA/training_neg.fasta')

# Sets 10% of the set to be the test set and 90% to be the training set
test_set_fraction = 0.1
test_set_percentage = 100 / (100 * test_set_fraction)

training_sequences = []
training_response = []
test_sequences = []
test_response = []

# Assignes a value to the sequences so machine can discriminate between motif and non-motif
for i in range(0, len(regulon)):
    if i % test_set_percentage == 0:
        test_sequences.append(regulon[i])
        test_response.append(1)
    else:
        training_sequences.append(regulon[i])
        training_response.append(1)

for i in range(0, len(non_regulon)):
    if i % test_set_percentage == 0:
        test_sequences.append(non_regulon[i])
        test_response.append(0)
    else:
        training_sequences.append(non_regulon[i])
        training_response.append(0)


integer_encoder = LabelEncoder()
one_hot_encoder = OneHotEncoder(categories='auto')

train_features = []
for sequence in training_sequences: train_features.append(Anne_one_hot_encode(sequence))
np.set_printoptions(threshold=40)
train_features = tf.keras.preprocessing.sequence.pad_sequences(train_features)
train_features = np.stack(train_features)


test_features = []
for sequence in test_sequences:    test_features.append(Anne_one_hot_encode(sequence))
np.set_printoptions(threshold=40)
test_features = tf.keras.preprocessing.sequence.pad_sequences(test_features)
test_features = np.stack(test_features)

train_labels = training_response
one_hot_encoder = OneHotEncoder(categories='auto')
train_labels = np.array(train_labels).reshape(-1, 1)
train_labels = one_hot_encoder.fit_transform(train_labels).toarray()

test_labels = test_response
one_hot_encoder = OneHotEncoder(categories='auto')
test_labels = np.array(test_labels).reshape(-1, 1)
test_labels = one_hot_encoder.fit_transform(test_labels).toarray()

# Save the augmented data to a NumPy array
np.save('train_features-Bc.npy', train_features)
np.save('train_labels-Bc.npy', train_labels)
np.save(args.sessiondir + '/test_features.npy', test_features)
np.save(args.sessiondir + '/test_labels.npy', test_labels)
