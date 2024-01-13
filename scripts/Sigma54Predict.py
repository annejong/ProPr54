import numpy as np
import argparse
from tensorflow.keras.models import load_model
import screed
import re
import glob
import sys


from tensorflow.keras.layers import Conv1D, Dense, MaxPooling1D, Flatten, Bidirectional, LSTM, GRU, Dropout, Masking, Input, GaussianNoise, BatchNormalization
from tensorflow.keras.models import Sequential
from tensorflow.keras.callbacks import ModelCheckpoint, History, ReduceLROnPlateau, EarlyStopping



parser = argparse.ArgumentParser(description='Predict Sigma54 binding sites')
parser.add_argument('-sessiondir', dest='sessiondir', help='Session Dir', nargs='?', default='.')
parser.add_argument('-modelName', dest='modelName', help='Full path and Name of the model', nargs='?', default='model_lstm-B_cer.h5')
parser.add_argument('-pvalue', dest='pvalue', help='pvalue', nargs='?', default=0.99)
parser.add_argument('-promlen', dest='PromLen', help='Length of Promoter', nargs='?', default=50)

args, unknown = parser.parse_known_args()

def write_log(S):
    # write to console and logfile
    print(S)
    f = open(args.sessiondir + '/' + "Predict_Sigma54_binding_sites.log", "a")
    f.write(S + '\n')
    f.close()

def read_kmers_from_file(filename, ksize, step_size):
    # a library for reading in FASTA/FASTQ, ksize is length of window/kmer
    all_kmers = []
    for record in screed.open(filename):
        #sequence = record.sequence  # Anne; removed this to get clean sequences
        sequence = re.sub('[BDEFHIJKLMNOPQRSUVWXYZ]','G', record.sequence.strip().upper())
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
    # sequence_list are the identified motifs
	annotated_motifs = {} 
	filtered_fasta = {}
	unique_seqs = set()
	for header, seq in fasta_dict.items():
		for s in sequence_list:
			pos = seq.find(s)
			if pos > -1: 
				annotated_motifs[s] = header
				#print(s+' '+header+' '+str(pos))
			if s in seq and seq not in unique_seqs:
				filtered_fasta[header] = seq
				unique_seqs.add(seq)
				#break
	return filtered_fasta, annotated_motifs



write_log(' ')


''' ==============================  LOAD MODEL ========================================================== '''
write_log('\t load model file ='+args.modelName)

model = load_model(args.modelName, compile=False)
model.compile()

# Finds files with extension .fnn, .fasta or .fa.

# tristan files_to_process = []
# tristan for ext in ['.fnn', '.fasta', '.fa']:
# tristan     files_to_process.extend(glob.glob(f'{args.query}*{ext}'))
# tristan 
# tristan prediction_set = read_kmers_from_file(files_to_process, args.PromLen, 3)

# code for one intergenic fasta file called query.g2d.intergenic.ffn
prediction_set = read_kmers_from_file(args.sessiondir + '/query.g2d.intergenic.ffn', args.PromLen, 3)

prediction_features = []
for sequence in prediction_set:
    prediction_features.append(Anne_one_hot_encode(sequence))
np.set_printoptions(threshold=40)
prediction_features = np.stack(prediction_features)



''' ==============================  PREDICTION ========================================================== '''
write_log('\tPREDICTION RUNNING')
predicted_labels = model.predict(prediction_features)

write_log('\tPREDICTION COMPLETE ')
threshold = args.pvalue
regulon_sequences = prediction_features[predicted_labels[:, 1] >= threshold]
decoded_sequences = Anne_decode(regulon_sequences)


if len(decoded_sequences) == 0:
	with open(args.sessiondir+'/query.Sig54.fna', "w") as f:
		f.write("No promoters found")
	write_log(' ==>  No promoters found ')
	sys.exit()
else:
	intergenic, annotated_motifs  = create_filtered_fasta(readMultiFasta(args.sessiondir + '/query.g2d.intergenic.ffn'), decoded_sequences)
	# write intergenic regions sequences as fasta
	with open(args.sessiondir+'/query.Sig54_intergenic.fna', 'w') as f:
		for key, value in intergenic.items(): f.write(f'>{key}\n{value}\n')  
	# write motifs table
	with open(args.sessiondir+'/query.Sig54_motifs.table', 'w') as f: 
		for key, value in annotated_motifs.items():  f.write(f'{value}\t{key}\n')  
	write_log(' ==>  Results written to files ')

