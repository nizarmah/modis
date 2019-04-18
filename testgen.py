import sys
import math
import random

from threading import Thread

num_threads	= 0
threads		= []
results		= []

sequences	= {}

def disco_shuffle(tid, sequence, size_vect, size_matt):
	for i in range(int(math.ceil(size_matt / num_threads))):
		random.shuffle(sequence)
		str_seq = "".join(sequence)

		try:
			tmp = sequences[str_seq]
		except:
			results[tid].append(str_seq)
			sequences[str_seq] = 1

def main(argv):
	global num_threads, threads, results
	num_threads	= int(argv[0])	# number of threads
	threads 	= [None] * num_threads
	results		= [[]] * num_threads

	size_vect	= int(argv[1])	# size of a single sequence
	size_matt	= int(argv[2])	# number of sequences

	sequence	= []

	alphabet	= [ "A", "C", "G", "T" ]
	for c in range(len(alphabet)):
		sequence += [alphabet[c]] * int(size_vect *  float(argv[c + 3]))

	global sequences
	sequences	= {}
	for tid in range(num_threads):
		threads[tid] = Thread(target=disco_shuffle,
				args=(tid, sequence, size_vect, size_matt))
		threads[tid].start()

	for thread in threads:
		thread.join()
	for result in results:
		for row in result:
			print(">")
			print(row)

if __name__ == "__main__":
	main(sys.argv[1:])
