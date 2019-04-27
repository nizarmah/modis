# modis

> Project Report : https://www.overleaf.com/read/hpjdzybksrty

The reason behind this project is to create my own parallelized enumeration algorithm for all ungapped nucleotide sequence motifs in a set of sequences.
While the `testgen.py` test generator has been also parallelized, the `msa.py` multiple sequence alignment, which uses ClustalOmega, has not been parallelized, as it is not in the scope of this project.

### My Approach:
1. Divide the profile matrix over threads by columns
 
2. Distribute the profile matrix by columns over threads, and anchor/drop all the probabilities that are less than a certain threshold to the bottom.
 
3. Identify the ranges of indices that contain an ungapped sequence, possible motif, of the requested length or more.
 
4. Distribute each range across multiple threads in order to enumerate all possible sub-sets of required length.
 
5. In an embarrassingly parallel way, calculate the distance between each pair.


### How to Run:

#### Generate Random Data:
```bash
python3 testgen.py 10 1000000 1000 > fasta/test_2.py
```
##### Arguments:
argument | example | name
-------- | -------- | --------
1 | 10 | number of processors
2 | 1000000 | nucleotide sequence length
3 | 1000 | number of nucleotide sequences

#### Performing MSA:
```bash
python3 msa.py .out fasta/test_2.py
```
##### Arguments:
argument | example | name
-------- | -------- | --------
1 | .out | output file
2+ | fasta/test_2.py | input files

> Please note, this is serial. It is done by running ClustalOmega on the fastas. If your sequences are of the same length and you prefer not wasting your time running this, then don't.

#### Running Motif Enumeration:
```bash
python3 modis.py .out 64 10 20 0.3 100 0.05
```
##### Arguments:
argument | example | name
-------- | -------- | --------
1 | .out | input file
2 | 64 | number of processors
3 | 10 | motif length
4 | 20 | max number of motifs in result
5 | 0.3 | nucleotide probability threshold
6 | 100 | distance between motifs
7 | 0.05 | percentage of accepted error in distance

### Test Results:

After running the code on a dataset of size 1,000 sequences of size 1,000,000, the following speedup was achieved.

![Sequence Motif Enumeration](https://raw.githubusercontent.com/nizarmah/modis/master/sequence-motif-enumeration-testresults.png)
