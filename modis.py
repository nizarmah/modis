import sys
import math

from queue import Queue
from threading import Thread, Lock

DEBUG		= True	# DEBUG Mode - Printing

mutex 		= {}

num_threads 	= 0	# number of available threads
threads		= []	# vector of threads
results		= []	# vector of results of threads

# reads fasta from input file
def readfasta(in_file):
	'''
	read the multiple sequence alignment that are in fasta format
	'''
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
	'''
	go through every index of each sequence and then
	calculate the probability of each character at that index
	'''
	for i in range(size_vect):
		for sequence in sequences:
			if sequence[i] not in alphabet:
				continue

			# increase the probability of the alphabet character
			pmatt[str(sequence[i])][i] += (1 / size_matt)

def threaded_pmatt_calc(tid, sequences, pmatt, alphabet=[ "A", "G", "T", "C" ],
				size_vect=0, size_matt=0):
	'''
	divide the segments of the sequences
	over the different num of threads ( size )

	each thread has a Start_Tuple and an Final_Tuple
	each tuple -> ( Sequence Index, Character Index )

	Number of Characters we have is (size_vect * size_matt), so:
		- len_of_segment  = (size_vect * size_matt) / num_threads
		- offset_segment  = len_of_segment * tid
		- character_index = offset_segment % size_vect
		- sequence_index  = (offset_segment - character_index) / size_vect
	'''

	len_of_segm	= math.ceil((size_vect * size_matt)
					/ num_threads)		# thread segment length
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
	'''
	start from each character index and try all subsets until end
	keep saving the subsets as long as all the characters in it
	have a probability > p_threshold, if not save the subset and its score

	score(AGT) = p(A) * p(G) * p(T) / (1/4)**3
	'''

	motif		= [ "", 0, 0 ]	# motif [ sequence, score, length ]
	motifs		= {}		# dictionary of motifs

	size_alph	= len(pmatt)	# alphabet size

	# go over the msa profile matrix
	# get probability of each character at an index
	for i in range(size_vect):
		for character in pmatt.keys():
			if (motif[2] == 0 and				# motif is empty
				pmatt[character][i] < p_threshold):	# char prob < threshold
				continue				# ignore character

			elif (motif[2] == len_motif and			# motif is len required
				pmatt[character][i] < p_threshold):	# char prob < threshold
				motifs[motif[0]] = motif[1]		# save motif in dict
				motif = [ "", 0, 0 ]			# reset motif

			else:
				motif[0] += character			# append char to motif
				motif[1] *= (pmatt[character][i]
						/ (1 / size_alph))	# increase motif score
				motif[2] += 1				# increase motif len

	# if loop ends with an unappended motif of required len
	if motif[2] == len_motif:
		motifs[motif[0]] = motif[1]
		motif = [ "", 0, 0 ]

	return motifs

def anchored_pmatt_fltr(tid, pmatt, p_threshold, size_vect=0):
	'''
	divides the columns of the pmatt across the threads and then
	keeps only the characters with probability > p_threshold
	inside a dictionary so the possible set to enumerate is smaller
	'''

	len_of_segm	= math.ceil(size_vect
					/ num_threads)	# thread segment length
	offset_segm	= len_of_segm * tid		# segment start offset

	# if this is last thread, then only iterate over
	# remaining length of all the segments of sequences
	if tid == (num_threads - 1):
		len_of_segm = size_vect - offset_segm

	div		= {}				# div of probability vector

	# go through all the columns of the probability matrix
	for index_col in \
		range(offset_segm, offset_segm + len_of_segm):
			subset_div = [ val for val in
						# go through probability at column of index
						[ (index_row, pmatt[index_row][index_col])
							for index_row in pmatt.keys() ]

						# go through the probability of the chars
						# that are greater than the p_threshold
						if val[1] >= p_threshold ]

			if subset_div:
				div[index_col] = subset_div

	results[tid] = div

def threaded_pmatt_intr(pmatt, p_threshold, len_motif=6,
				size_vect=0, size_matt=0):
	'''
	dividing and conquering the sequence, then merging is a possible solution
	however its speedup is not as significant as expected, as the worst case
	would be having to redo the sequential code in steps

	so, we will anchor the thresholds < p_threshold in the bottom, and keep
	the characters that are above the threshold as possible permutations
	by dividing the columns across the threads and each anchoring each

	then we will cross that list once and collect all possible permutations
	'''

# nucleotide motif discovery through profile matrix
def nucleotide_modis(sequences, len_motif=6, max_motifs=3, p_threshold=0.8):
	alphabet	= [ "A", "G", "T", "C" ]

	pmatt		= {}			# probability matrix
	for character in alphabet:
		pmatt[character] = [0] * len(sequences[0])

	size_matt	= len(sequences)    	# matrix size of all sequences
	size_vect	= len(sequences[0]) 	# vector size of all sequences

	if DEBUG:
		print("size_vect : " + str(size_vect))
		print("size_matt : " + str(size_matt))

	''' sequential code
	sequential_pmatt_calc(sequences, pmatt, alphabet=alphabet,
	 			size_vect=size_vect, size_matt=size_matt)

	motifs		= sequential_pmatt_intr(pmatt, p_threshold, len_motif,
				size_vect=size_vect, size_matt=size_matt)
	'''

	for character in alphabet:
		mutex[character] = [Lock()] * size_vect

	global num_threads, threads, results	# getting global variables
#	num_threads	= 0			# uncomment to stop parallelism
	threads 	= [None] * num_threads	# vector of threads
	results		= [None] * num_threads	# vector of results

	for tid in range(num_threads):
		threads[tid] = Thread(target=threaded_pmatt_calc,
			args=(tid, sequences, pmatt),
			kwargs=dict(alphabet=alphabet,
				size_vect=size_vect, size_matt=size_matt))
		threads[tid].start()

	for thread in threads:
		thread.join()

	if DEBUG:
		print("\n> Probability Matrix : \n" + str(pmatt))

	pvect		= {}			# probability vector

	for tid in range(num_threads):
		threads[tid] = Thread(target=anchored_pmatt_fltr,
				args=(tid, pmatt, p_threshold),
				kwargs=dict(size_vect=size_vect))
		threads[tid].start()

	for thread in threads:
		thread.join()
	for result in results:
		if result: pvect.update(result)

	if DEBUG:
		print("\n> Anchored Probability Vector : \n" + str(pvect))

#	if DEBUG:
#		print("\n> Motifs : \n" + str(motifs))

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
