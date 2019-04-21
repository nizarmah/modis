import sys
import math
import time

from random import shuffle

from multiprocessing import Process, Lock, Manager

DEBUG		= True		# DEBUG Mode - Printing

num_threads 	= 0		# number of available threads
threads		= []		# vector of threads
results		= []

manager		= Manager()	# manager to comm all results

# reads fasta from input file
def readfasta(sequences, in_file):
	'''
	read the multiple sequence alignment that are in fasta format
	'''
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

def threaded_pmatt_calc(tid, sequences, alphabet=[ "A", "G", "T", "C" ],
				size_vect=0, size_matt=0):
	'''
	divide the sequence into segments
	over the different num of threads ( size )
	instead of dividing sequences into segments
	that way we don't need mutex locks

	so:
		- len_of_segment  = size_vect / num_threads
		- offset_segment  = len_of_segment * tid
	'''

	len_of_segm	= math.ceil(size_vect
					/ num_threads)	# thread segment length
	offset_segm	= len_of_segm * tid		# segment start offset

	# if this is last thread, then only iterate over
	# remaining length of all the segments of sequences
	if tid == (num_threads - 1):
		len_of_segm = size_vect - offset_segm

	pmatt_segm	= {}				# pmatt local segment
	for character in alphabet:
		pmatt_segm[character] = [0] * len_of_segm

	tmp_sid		= 0				# temp sequence id
	for sequence in sequences:
		tmp_sid += 1
		if DEBUG:
			print("# Thread : " + str(tid) + " : " +
				"Sequence : " + str(tmp_sid) + " : " +
				"Length : " + str(len_of_segm), end="\r")

		for i in range(len_of_segm):
			index_char	= offset_segm + i

			# check if the character at character_index
			# is part of the alphabet or if it is not
			if sequence[index_char] not in alphabet:
				continue

			# update pmatt_segm because it's faster
			pmatt_segm[str(sequence[
				index_char])][i] += (1 / size_matt)

	results[tid] = pmatt_segm

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

	if (intv_motif[0] > -1 and
		intv_motif[1] > -1):
		intv_motifs.append(
			(intv_motif[0], intv_motif[1]))

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
		tmp_motifs	= [ motif[:] ]		# temporary motifs
		new_motifs	= []			# new temporary motifs

		index_start	= (intv_motif[0] + offset_segm +
						i)	# start index in interval

		if (intv_motif[1] -
			index_start) < len_motif:
			break

		# go through the pvect len_motif times
		# hoping we will find a big enough motif
		for j in range(len_motif):
			index_char = index_start + i

			# for the different characters of the alphabet
			# at the index 'index_char' of pvect
			for k in pvect[index_char]:
				if DEBUG:
					print("# Thread : " + str(tid) + " : " +
						str(i) + " : " + str(j) + " : " + str(k),
						end="\r")

				# go through each previous motif and append to it
				# the different characters, to the previous motifs
				for old_motif in tmp_motifs:
					new_motif = old_motif[:]

					new_motif[0] += k[0]
					new_motif[1] *= (k[1] / (1 / size_alph))
					new_motif[2] += 1

					new_motifs.append(new_motif)

			tmp_motifs = new_motifs
			new_motifs = []

		motifs += tmp_motifs

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
	global manager
	global num_threads, threads, results	# getting global variables

	alphabet	= [ "A", "G", "T", "C" ]

	size_matt	= len(sequences)    	# matrix size of all sequences
	size_vect	= len(sequences[0]) 	# vector size of all sequences
	size_alph	= len(alphabet)		# alphabet size of all characters

	if DEBUG:
		print("# n_processors\t:\t" + str(num_threads))
		print("# size_vect\t:\t" + str(size_vect))
		print("# size_matt\t:\t" + str(size_matt))
		print()

	pmatt		= {}			# manager probability matrix
	for character in alphabet:
		pmatt[character] = []

	threads 	= [None] * num_threads	# vector of threads
	results		= manager.list(range(
				num_threads))	# manager list to store results

	start_time	= time.time()		# start time of execution time

	for tid in range(num_threads):
		threads[tid] = Process(target=threaded_pmatt_calc,
			args=(tid, sequences),
			kwargs=dict(alphabet=alphabet,
				size_vect=size_vect, size_matt=size_matt))
		threads[tid].start()

	for thread in threads:
		thread.join()
	for result in results:
		for character in alphabet:
			pmatt[character] += result[character]

	if DEBUG:
		print()
		print()

		print("> Probability Matrix : \n" + str(pmatt))
		print()

	pvect		= {}			# probability vector
	results		= manager.list(range(
				num_threads))	# manager list to store results
	for tid in range(num_threads):
		threads[tid] = Process(target=anchored_pmatt_fltr,
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
		results		= manager.list(range(
				num_threads))	# manager list to store results
		for tid in range(num_threads):
			threads[tid] = Process(target=threaded_motif_sset,
				args=(tid, pvect, intv_motif),
				kwargs=dict(len_motif=len_motif,
					size_vect=size_vect, size_alph=size_alph))
			threads[tid].start()

		for thread in threads:
			thread.join()

		for result in results:
			if isinstance(result, list):
				motifs += result

	if DEBUG:
		print()
		print()

		print("> Motifs : " + str(len(motifs)) +
				" \n" + str(motifs))
		print()

	distant_motifs	= []			# motifs of certain distance
	results		= manager.list(range(
				num_threads))	# manager list to store results
	for tid in range(num_threads):
		threads[tid] = Process(target=pairs_of_distance,
				args=(tid, motifs, m_distance),
				kwargs=dict(e_tolerance=e_tolerance))
		threads[tid].start()

	for thread in threads:
		thread.join()
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
	sequences	= manager.list()	# manager sequences list
	readfasta(sequences, in_file)

	nucleotide_modis(sequences, len_motif=len_motif, max_motifs=max_motifs,
		p_threshold=p_threshold, m_distance=m_distance, e_tolerance=e_tolerance)

if __name__ == "__main__":
	main(sys.argv[1:])
