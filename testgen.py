import os
import sys
import math
import random

from multiprocessing import Process

DEBUG		= False

num_processes	= 0
processes	= []
results		= []

sequences	= {}

def disco_shuffle(tid, size_vect, size_matt,
			alphabet=["A","G","T","C"]):
	tmp_results	= []
	size_process	= int(math.ceil(size_matt / num_processes))

	if DEBUG:
		print("Processes " + str(tid) + " : Started : " +
				str(size_process))

	for i in range(size_process):
		str_seq = ""
		for j in range(size_vect):
			str_seq += random.choice(alphabet)

#		try:
#			tmp = sequences[str_seq]
#		except:
#			results[tid].append(str_seq)
#			sequences[str_seq] = 1
		tmp_results.append(str_seq)

	if DEBUG:
		print("Processes " + str(tid) + " : Ended")

	return tmp_results

def main(argv):
	global num_processes, processes, results
	num_processes	= int(argv[0])	# number of processes
	processes	= [None] * num_processes
	results		= [[]] * num_processes

	size_vect	= int(argv[1])	# size of a single sequence
	size_matt	= int(argv[2])	# number of sequences

	alphabet	= [ "A", "C", "G", "T" ]

	global sequences
	sequences	= {}
	for pid in range(num_processes):
		processes[pid] = Process(target=disco_shuffle,
				args=(pid, size_vect, size_matt),
				kwargs=dict(alphabet=alphabet))
		processes[pid].start()

	for process in processes:
		process.join()


	for result in results:
		for row in result:
			print(">")
			print(row)

if __name__ == "__main__":
	main(sys.argv[1:])
