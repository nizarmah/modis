import os
import sys

# from Bio.Align.Applications import ClustalOmegaCommandline

def main(argv):
	# output file is the first argument
	out_file = argv[0]
	# all other arguments are input / fasta files
	in_files = " ".join(argv[1:])

	# save all fastas in one file using `cat`
	cin_file = ".cin"
	os.system("cat " + in_files + " > " + cin_file)

	# get clustal's command
# 	clustalomega_cline = ClustalOmegaCommandline(infile=cin_file, outfile=out_file,
#							verbose=True, auto=True)

	# hard-coded clustal's command in order
	# to not use the library that is required
	clustalomega_cline = "clustalo -i .cin -o " + out_file + " --auto -v"

	# run clustal's command
	os.system("./" + str(clustalomega_cline))

if __name__ == "__main__":
	main(sys.argv[1:])
