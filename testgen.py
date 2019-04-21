import sys
import math
import random
import time

from multiprocessing import Process, Queue

DEBUG		= False

num_threads	= 0
threads		= []
queue		= None

sequences	= {}

def disco_shuffle(tid, size_vect, size_matt,
			alphabet=["A","G","T","C"]):
	tmp_results	= []
	size_thread	= int(math.ceil(size_matt / num_threads))

	start_time = time.time()
	if DEBUG:
		print("Thread : " + str(tid) + " : Started : " +
				str(size_thread))

	for i in range(size_thread):
		str_seq = ""
		for j in range(size_vect):
			str_seq += random.choice(alphabet)

#		try:
#			tmp = sequences[str_seq]
#		except:
#			results[tid].append(str_seq)
#			sequences[str_seq] = 1
		tmp_results.append(str_seq)

	queue.put(tmp_results)

	end_time = time.time()
	if DEBUG:
		print("Thread : " + str(tid) + " : Ended : " +
				str(end_time - start_time))

def main(argv):
	global num_threads, threads, queue
	num_threads	= int(argv[0])	# number of threads
	threads		= [None] * num_threads
	queue		= Queue()

	size_vect	= int(argv[1])	# size of a single sequence
	size_matt	= int(argv[2])	# number of sequences

	alphabet	= [ "A", "G", "T", "C" ]

	global sequences
	sequences	= {}

	global DEBUG
	try:
		DEBUG 	= argv[3]	# optional debug parameter
	except: pass

	for tid in range(num_threads):
		threads[tid] = Process(target=disco_shuffle,
				args=(tid, size_vect, size_matt),
				kwargs=dict(alphabet=alphabet))
		threads[tid].start()

#	for thread in threads:
#		thread.join()

	for thread in threads:
		res = queue.get()
		for row in res:
			print(">")
			print(row)

if __name__ == "__main__":
	main(sys.argv[1:])
