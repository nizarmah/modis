import sys
import math
import time

from random import shuffle

from queue import Queue
from threading import Thread, Lock

DEBUG		= True	# DEBUG Mode - Printing

mutex 		= {}	# mutex locks

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
				print("> Thread : " + str(tid) + " : " +
					("Start" if i == 0 else " End ") +
					" : " + str((index_seq, index_char)))
				if i == (len_of_segm - 1): print()

		# specify which sequence we are working on
		sequence = sequences[index_seq]

		# check if the character at character_index
		# is part of the alphabet or if it is not
		if sequence[index_char] not in alphabet:
			continue

		# while mutex is still locked, wait
		while mutex[str(sequence[index_char])][index_char].locked:
			continue

		# acquire mutex for character_index
		mutex[str(sequence[index_char])][index_char].acquire()
		try:
			pmatt[str(sequence[index_char])][index_char] += (1 / size_matt)
		finally:
			mutex[str(sequence[index_char])][index_char].release()

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
							for index_row in pmatt.keys()
							if index_col < len(pmatt[index_row]) ]

						# go through the probability of the chars
						# that are greater than the p_threshold
						if val[1] >= p_threshold ]

			if subset_div:
				div[index_col] = subset_div

	results[tid] = div

def advanced_pmatt_intr(pvect, len_motif=6, size_vect=0):
	'''
	dividing and conquering the sequence, then merging is a possible solution
	however its speedup is not as significant as expected, as the worst case
	would be having to redo the sequential code in steps

	so we will take the pvect and keep the segments that have the desired length
	through keeping track of their start and their end
	'''

	intv_motif	= [ -1, -1 ]			# interval motif [ start, end ]
	intv_motifs	= []				# list of intervals of motifs

	for index_char in pvect.keys():
		if (intv_motif[0] == -1):		# if interval starts at -1
			intv_motif[0] = index_char	# start setting interval
		elif (intv_motif[1] == -1 and index_char ==
			(intv_motif[0] + 1)):		# if interval has consecutive
			intv_motif[1] = index_char	# add first consecutive as end
		elif (intv_motif[1] > -1 and index_char ==
			(intv_motif[1] + 1)):		# end of intv has consecutive
			intv_motif[1] = index_char	# checkpoint as new end
		elif (index_char >
			(intv_motif[1] + 1)):		# if interval not consecutive
			if ((intv_motif[1] - intv_motif[0] + 1)
					>= len_motif):	# if interval > len of motif
				intv_motifs.append(	# save the interval in list
					(intv_motif[0], intv_motif[1]))

			intv_motif = [ -1, -1 ]

	if DEBUG:
		print("> Possible Indeces of Motif Sequences : \n" +
				str(intv_motifs))
		print()

	return intv_motifs

def threaded_motif_sset(tid, pvect, intv_motif, len_motif=6,
				size_vect=0, size_alph=4):
	'''
	first of all get the intervals of possible motifs from advanced_pmatt_intr

	after that we will use the threads to compute all possible subsets of motifs
	since the worst case we have if all the sequence itself is considered as one
	because it is highly conserved, and it's easier to traverse all other indeces
	'''
	motif		= [ "", 1, 0 ]			# [ sequence, score, length ]
	motifs		= []				# list of chosen motifs

	len_of_intv	= (intv_motif[1] -		# length of the interval
				intv_motif[0] + 1)
	len_of_segm	= math.ceil(len_of_intv
					/ num_threads)	# thread segment length
	offset_segm	= len_of_segm * tid		# segment start offset

	# if this is last thread, then only iterate over
	# remaining length of all the segments of sequences
	if tid == (num_threads - 1):
		len_of_segm = len_of_intv - offset_segm
		len_of_segm = len_of_segm if len_of_segm > 0 else 0

	# traverse this loop according to how the interval was
	# divided between the threads,  how big each segment is
	for i in range(len_of_segm):
		size_of_motif	= 0			# size of temp motif
		tmp_motifs	= [ motif ]		# temporary motifs
		new_motifs	= []			# new temporary motifs

		# go through the pvect len_motif times
		# hoping we will find a big enough motif
		for j in range(len_motif):
			# check if index_char still exists
			index_char = intv_motif[0] + offset_segm + j
			if (index_char > intv_motif[1]):
				break

			# for the different characters of the alphabet
			# at the index 'index_char' of pvect
			for k in pvect[index_char]:
				if DEBUG:
					print("> Thread : " + str(tid) + " : " +
						str(i) + " : " + str(j) + " : " + str(k))

				# go through each previous motif and append to it
				# the different characters, to the previous motifs
				for old_motif in tmp_motifs:
					new_motif = old_motif

					new_motif[0] += k[0]
					new_motif[1] *= (k[1] / (1 / size_alph))
					new_motif[2] += 1

					size_of_motif += 1
					new_motifs.append(new_motif)
				tmp_motifs = new_motifs
				new_motifs = []

		if size_of_motif == len_motif:
			motifs += tmp_motifs

	if DEBUG:
		print()

	results[tid] = motifs

def pairs_of_distance(tid, motifs, m_distance, e_tolerance=0.05):
	'''
	distance of two scores = abs(d1 - d2)
	'''

	size_motifs	= len(motifs)			# size of motifs
	distant_motifs	= []				# coll of distant motifs

	len_of_segm	= math.ceil(size_motifs
					/ num_threads)	# thread segment length
	offset_segm	= len_of_segm * tid		# segment start offset

	# if this is last thread, then only iterate over
	# remaining length of all the segments of sequences
	if tid == (num_threads - 1):
		len_of_segm = size_motifs - offset_segm

	# go over the segment for this thread
	for i in range(len_of_segm):
		# if index is greater than list, stop
		index_motif = offset_segm + i
		if index_motif == (size_motifs - 1):
			break

		# check every motif after this motif
		# and get the distance between them
		for other_motif in \
				range(index_motif + 1, size_motifs):
			tmp_distance = abs(motifs[index_motif][1] -
					motifs[other_motif][1])

			# if distance is as required, append
			if math.isclose(m_distance, tmp_distance,
					rel_tol=e_tolerance, abs_tol=0.0):
				distant_motifs.append((motifs[index_motif],
						motifs[other_motif]))

	results[tid] = distant_motifs

# nucleotide motif discovery through profile matrix
def nucleotide_modis(sequences, len_motif=6, max_motifs=3,
		p_threshold=0.8, m_distance=2, e_tolerance=0.05):
	alphabet	= [ "A", "G", "T", "C" ]

	pmatt		= {}			# probability matrix
	for character in alphabet:
		pmatt[character] = [0] * len(sequences[0])

	size_matt	= len(sequences)    	# matrix size of all sequences
	size_vect	= len(sequences[0]) 	# vector size of all sequences
	size_alph	= len(alphabet)		# alphabet size of all characters

	if DEBUG:
		print("size_vect : " + str(size_vect))
		print("size_matt : " + str(size_matt))
		print()

	for character in alphabet:
		mutex[character] = [Lock()] * size_vect

	global num_threads, threads, results	# getting global variables
#	num_threads	= 0			# uncomment to stop parallelism
	threads 	= [None] * num_threads	# vector of threads
	results		= [None] * num_threads	# vector of results

	start_time	= time.time()		# start time of execution time

	for tid in range(num_threads):
		threads[tid] = Thread(target=threaded_pmatt_calc,
			args=(tid, sequences, pmatt),
			kwargs=dict(alphabet=alphabet,
				size_vect=size_vect, size_matt=size_matt))
		threads[tid].start()

	for thread in threads:
		thread.join()

	if DEBUG:
		print("> Probability Matrix : \n" + str(pmatt))
		print()

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
		print("> Anchored Probability Vector : \n" + str(pvect))
		print()

	intv_motifs	= advanced_pmatt_intr(pvect,
				len_motif=len_motif, size_vect=size_vect)

	motifs		= []			# list of motifs
	for intv_motif in intv_motifs:
		for tid in range(num_threads):
			threads[tid] = Thread(target=threaded_motif_sset,
				args=(tid, pvect, intv_motif),
				kwargs=dict(len_motif=len_motif,
					size_vect=size_vect, size_alph=size_alph))
			threads[tid].start()

		for thread in threads:
			thread.join()
		for result in results:
			motifs += result

	if DEBUG:
		print("> Motifs : \n"+ str(motifs))
		print()

	for tid in range(num_threads):
		threads[tid] = Thread(target=pairs_of_distance,
				args=(tid, motifs, m_distance),
				kwargs=dict(e_tolerance=e_tolerance))
		threads[tid].start()

	for thread in threads:
		thread.join()
	distant_motifs = []
	for result in results:
		if result:
			distant_motifs += result

	if DEBUG:
		if (max_motifs > -1):
			shuffle(distant_motifs)
			distant_motifs = distant_motifs[0:max_motifs]

		print("> Distant Motifs : \n" + str(distant_motifs))
		print()

	end_time	= time.time()		# end time of execution time

	if DEBUG:
		print("> Time Taken : " + str(end_time - start_time))

def main(argv):
	in_file		= argv[0] 		# input file of fastas after MSA

	global num_threads			# prevents creation of local var
	num_threads	= int(argv[1]) 		# number of threads

	len_motif	= int(argv[2])		# number of characters in motif
	max_motifs	= int(argv[3])		# maximum number of motifs
	p_threshold	= float(argv[4])	# probability matt item min threshold
	m_distance	= float(argv[5])	# metric distance between pairs
	e_tolerance	= float(argv[6])	# error tolerance for the distance

	# read aligned fasta from input file
	sequences = readfasta(in_file)

	nucleotide_modis(sequences, len_motif=len_motif, max_motifs=max_motifs,
		p_threshold=p_threshold, m_distance=m_distance, e_tolerance=e_tolerance)

if __name__ == "__main__":
	main(sys.argv[1:])
