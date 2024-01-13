import numpy as np
import argparse
from tensorflow.keras.models import load_model
import screed
import re
import glob
import sys

parser = argparse.ArgumentParser(description='Train the model')
parser.add_argument('-promlen', dest='PromLen', help='Length of Promoter', nargs='?', default=50)
parser.add_argument('-padsize', dest='padSize', help='Size of padding on both sides', nargs='?', default=12)
parser.add_argument('-sessiondir', dest='sessiondir', help='Session Dir', nargs='?', default='.')
parser.add_argument('-modelName', dest='modelName', help='Name of the model, see models.py', nargs='?',
                    default='model_lstm-B_cer')
parser.add_argument('-batchsize', dest='batchSize', help='Batch size', nargs='?', default=64)
parser.add_argument('-epochs', dest='epochs', help='Number of training epochs', nargs='?', default=40)
args, unknown = parser.parse_known_args()

def write_log(S):
    # write to console and logfile
    print(S)
    f = open(args.sessiondir + '/' + "log_ppp_{model_name}.log".format(model_name=args.modelName), "a")
    f.write(S + '\n')
    f.close()

def read_kmers_from_file(filename, ksize, step_size):
    # a library for reading in FASTA/FASTQ, ksize is length of window/kmer
    all_kmers = []
    for record in screed.open(filename):
        sequence = record.sequence

        kmers = []
        n_kmers = len(sequence) - ksize + 1

        for i in range(0, n_kmers, step_size):
            kmer = sequence[i:i + ksize]
            kmers.append(kmer)

        all_kmers += kmers

    return all_kmers

def Anne_one_hot_encode(seq):
    # One hot encode your sequences for the CNN
    mapping = dict(zip("ACGT", range(4)))
    seq2 = [mapping[i] for i in seq]
    return np.eye(4)[seq2]

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

def readMultiFasta(fna_file):
    # all chars not A,C,T,G will be replaced by G, put sequence into dictionary
    fasta = {}
    key= ''
    with open(fna_file) as f:
        for line in f:
            if line[0] == ">":
                items = re.match(">(.*?)\s", line)
                if items:
                    key=items.group(1)
                    fasta[key]=''
            else:
                if (key != ''):
                    fasta[key] += re.sub('[BDEFHIJKLMNOPQRSUVWXYZ]','G', line.strip().upper())
        return fasta

def create_filtered_fasta(fasta_dict, sequence_list):
    # Uses a FASTA file and looks if sequence list in FASTA. If it is then appends to new dictionary
    filtered_fasta = {}
    unique_seqs = set()
    for header, seq in fasta_dict.items():
        for s in sequence_list:
            if s in seq and seq not in unique_seqs:
                filtered_fasta[header] = seq
                unique_seqs.add(seq)
                break
    return filtered_fasta

def writeFasta(FASTA_dict, output_file):
    with open(output_file, 'w') as f:
        for key, value in FASTA_dict.items():
            f.write(f'>{key}\n')  # write the sequence identifier
            f.write(f'{value}\n')  # write the sequence data

write_log(' ')


''' ==============================  LOAD MODEL ========================================================== '''
write_log(' ==>  LOAD MODEL ')
write_log('model file ='+args.modelName)

checkpoint_filepath = args.sessiondir + '/' + args.modelName + '.h5'

model = load_model(checkpoint_filepath)

write_log(' ')

# Finds files with extension .fnn, .fasta or .fa.

files_to_process = []
for ext in ['.fnn', '.fasta', '.fa']:
    files_to_process.extend(glob.glob(f'{args.query}*{ext}'))

prediction_set = read_kmers_from_file(files_to_process, args.PromLen, 3)

prediction_features = []

for sequence in prediction_set:
    prediction_features.append(Anne_one_hot_encode(sequence))
np.set_printoptions(threshold=40)
prediction_features = np.stack(prediction_features)

write_log(' ')


''' ==============================  PREDICTION ========================================================== '''
write_log(' ==>  PREDICTION RUNNING')

predicted_labels = model.predict(prediction_features)

write_log(' ==>  PREDICTION COMPLETE ')

write_log(' ')

threshold = args.pvalue
regulon_sequences = prediction_features[predicted_labels[:, 1] >= threshold]

decoded_sequences = Anne_decode(regulon_sequences)

output = create_filtered_fasta(readMultiFasta(files_to_process), decoded_sequences)

if len(decoded_sequences) == 0:
	with open(args.sessiondir+'/'+args.outPrefix+'.Sig54.txt', "w") as f:
		f.write("No promoters found")
	sys.exit()

if len(decoded_sequences) != 0:
    writeFasta(output, args.sessiondir+'/'+args.outPrefix+'.Sig54.txt')