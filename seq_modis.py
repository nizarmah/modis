import sys
import math
import time

from random import shuffle

from queue import Queue
from threading import Thread, Lock

DEBUG		= True	# DEBUG Mode - Printing

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

def sequential_pmatt_intr(pmatt, p_threshold, len_motif=6,
				size_vect=0, size_matt=0):
	'''
	start from each character index and try all subsets until end
	keep saving the subsets as long as all the characters in it
	have a probability > p_threshold, if not save the subset and its score

	example:
	score(AGT) = p(A) * p(G) * p(T) / (1/4)**3
	'''

	motif		= [ "", 1, 0 ]				# [ sequence, score, length ]
	motifs		= []					# list of motifs

	size_alph	= len(pmatt)				# alphabet size

	# go over the msa profile matrix
	# get probability of each character at an index
	for i in range(size_vect):
		for character in pmatt.keys():
			if (motif[2] == 0 and				# motif is empty
				pmatt[character][i] < p_threshold):	# char prob < threshold
				continue				# ignore character

			elif (motif[2] == len_motif and			# motif is len required
				pmatt[character][i] < p_threshold):	# char prob < threshold
				motifs.append((motif[0], motif[1]))	# save motif in list
				motif = [ "", 1, 0 ]			# reset motif

			else:
				motif[0] += character			# append char to motif
				motif[1] *= (pmatt[character][i]
						/ (1 / size_alph))	# increase motif score
				motif[2] += 1				# increase motif len

	# if loop ends with an unappended motif of required len
	if motif[2] == len_motif:
		motifs.append((motif[0], motif[1]))
		motif = [ "", 0, 0 ]

	return motifs

def sequential_distance(motifs, m_distance, e_tolerance=0.05):
	'''
	distance of two scores = abs(d1 - d2)
	'''

	size_motifs	= len(motifs)			# size of motifs
	distant_motifs	= []				# coll of distant motifs

	# go over all the motifs
	for index_motif in range(len(motifs)):
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

	return distant_motifs

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

	start_time	= time.time()		# start time of execution time

	if DEBUG:
		print("size_vect : " + str(size_vect))
		print("size_matt : " + str(size_matt))
		print()

	sequential_pmatt_calc(sequences, pmatt, alphabet=alphabet,
	 			size_vect=size_vect, size_matt=size_matt)

	if DEBUG:
		print("> Probability Matrix : \n" + str(pmatt))
		print()

	motifs		= sequential_pmatt_intr(pmatt, p_threshold, len_motif,
				size_vect=size_vect, size_matt=size_matt)

	if DEBUG:
		print("> Motifs : \n"+ str(motifs))
		print()

	distant_motifs	= sequential_distance(motifs, m_distance,
						e_tolerance=e_tolerance)

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

	len_motif	= int(argv[1])		# number of characters in motif
	max_motifs	= int(argv[2])		# maximum number of motifs
	p_threshold	= float(argv[3])	# probability matt item min threshold
	m_distance	= float(argv[4])	# metric distance between pairs
	e_tolerance	= float(argv[5])	# error tolerance for the distance

	# read aligned fasta from input file
	sequences = readfasta(in_file)

	nucleotide_modis(sequences, len_motif=len_motif, max_motifs=max_motifs,
		p_threshold=p_threshold, m_distance=m_distance, e_tolerance=e_tolerance)

if __name__ == "__main__":
	main(sys.argv[1:])
