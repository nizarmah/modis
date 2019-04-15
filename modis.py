import sys
import math

from threading import Thread, Lock

DEBUG		= True

mutex 		= {}
num_threads 	= 1

# reads fasta from input file
def readfasta(in_file):
	sequences = []

	with open(in_file, "r") as file:
		sequence = ""

		# go over each line in the file
		for line in file.readlines():
			# if the line starts with >
			# then it is a new fasta
			# otherwise it is a continuation
			if line.startswith(">"):
				if len(sequence) > 0:
					sequences.append(sequence.strip())
					sequence = ""
			elif len(line) > 0:
				sequence += line.strip()

	return sequences[:]

def sequential_pmatt_calc(sequences, pmatt, alphabet=[ "A", "G", "T", "C" ],
				size_vect=0, size_matt=0):
	# go through every index of each sequence and then
	# calculate the probability of each character at that index
	for i in range(size_vect):
		for sequence in sequences:
			if sequence[i] not in alphabet:
				continue

			# increase the probability of the alphabet character
			pmatt[str(sequence[i])][i] += (1 / size_matt)

def threaded_pmatt_calc(tid, sequences, pmatt, alphabet=[ "A", "G", "T", "C" ],
				size_vect=0, size_matt=0):
	# divide the segments of the sequences
	# over the different num of processors ( size )

	# each thread has a Start_Tuple and an Final_Tuple
	# Each Tuple -> ( Sequence Index, Character Index )

	# Number of Characters we have is (size_vect * size_matt), so:
	#	- len_of_segment  = (size_vect * size_matt) / num_threads
	#	- offset_segment  = len_of_segment * tid
	#	- character_index = offset_segment % size_vect
	#	- sequence_index  = (offset_segment - character_index) / size_vect

	len_of_segm	= math.ceil((size_vect * size_matt)
					/ num_threads)		# process segment length
	offset_segm	= len_of_segm * tid			# segment start offset

	# if this is last thread, then only iterate over
	# remaining length of all the segments of sequences
	if tid == (num_threads - 1):
		len_of_segm = (size_vect * size_matt) - offset_segm

	for i in range(len_of_segm):
		index_char	= (offset_segm + i) % size_vect
		index_seq	= ((offset_segm + i) - index_char) // size_vect

		if DEBUG:
			if i == 0 or i == (len_of_segm - 1):
				print(str(tid) + " : " + str(i) +
					" : " + str((index_seq, index_char)))

		# specify which sequence we are working on
		sequence = sequences[index_seq]

		# check if the character at character_index
		# is part of the alphabet or if it is not
		if sequence[index_char] not in alphabet:
			continue

		# acquire mutex for character_index
		mutex[str(sequence[index_char])][index_char].acquire()
		try:
			pmatt[str(sequence[index_char])][index_char] += (1 / size_matt)
		finally:
			mutex[str(sequence[index_char])][index_char].release()


def sequential_pmatt_intr(pmatt, p_threshold, len_motif=6,
				size_vect=0, size_matt=0):
	motif	= [ "", 0, 0 ]
	motifs	= {}

	for i in range(size_vect):
		for character in pmatt.keys():
			if motif[1] == 0 and \
				pmatt[character][i] < p_threshold:
				continue
			elif motif[1] > 0 and \
				motif[2] == len_motif and \
				pmatt[character][i] < p_threshold:
				motifs[motif[0]] = motif[1]
				motif = [ "", 0, 0 ]
			else:
				motif[0] += character
				motif[1] += pmatt[character][i]
				motif[2] += 1

	if motif[2] == len_motif and motif[1] > 0:
		motifs.append((motif[0], motif[1], motif[2]))
		motif = [ "", 0, 0 ]

	return motifs

def threaded_pmatt_intr(pmatt, p_threshold,
				size_vect=0, size_matt=0):
	motif	= [ "", 0, 0 ]
	motifs	= {}

# nucleotide motif discovery through profile matrix
def nucleotide_modis(sequences, len_motif=6, max_motifs=3, p_threshold=0.8):
	alphabet	= [ "A", "G", "T", "C" ]

	size_vect	= 0
	size_matt	= 0

	pmatt		= {}
	for character in alphabet:
		pmatt[character] = [0] * len(sequences[0])

	size_matt	= len(sequences)    # matrix size of all sequences
	size_vect	= len(sequences[0]) # vector size of all sequences

	if DEBUG:
		print("size_vect : " + str(size_vect))
		print("size_matt : " + str(size_matt))

	sequential_pmatt_calc(sequences, pmatt, alphabet=alphabet,
	 			size_vect=size_vect, size_matt=size_matt)
	motifs = sequential_pmatt_intr(pmatt, p_threshold, len_motif,
				size_vect=size_vect, size_matt=size_matt)

	for character in alphabet:
		mutex[character] = [Lock()] * size_vect

	num_threads = 0
	threads 	= [None] * num_threads
	for tid in range(num_threads):
		threads[tid] = Thread(target=threaded_pmatt_calc,
			args=(tid, sequences, pmatt),
			kwargs=dict(alphabet=alphabet,
				size_vect=size_vect, size_matt=size_matt))
		threads[tid].start()

	for thread in threads:
		thread.join()

	if DEBUG:
		print("\n> Probability Matrix : " + str(pmatt))

#	parallel_pmatt_intr(pmatt, len_motif, p_threshold,
#				size_vect=size_vect, size_matt=size_matt)

	if DEBUG:
		print("\n> Motifs : " + str(motifs))

def main(argv):
	in_file		= argv[0] 		# input file of fastas after MSA

	global num_threads			# prevents creation of local var
	num_threads	= int(argv[1]) 		# number of threads

	len_motif	= int(argv[2])		# number of characters in motif
	max_motifs	= int(argv[3])		# maximum number of motifs
	p_threshold	= float(argv[4])	# probability matt item min threshold

	# read aligned fasta from input file
	sequences = readfasta(in_file)

	nucleotide_modis(sequences, len_motif=len_motif,
				max_motifs=max_motifs, p_threshold=p_threshold)

if __name__ == "__main__":
	main(sys.argv[1:])
