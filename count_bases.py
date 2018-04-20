import sys
import numpy as np
from Bio import SeqIO

fastq = sys.argv[1]

lines = open(fastq, "r").readlines()

count = 0
line_count = 0
seq_lens = []
for record in SeqIO.parse(fastq, "fastq"):
    sl = len(record)
    count += sl
    seq_lens.append(sl)



print("Number of bases: " + str(count) + "\n")
print("Mean seq length: " + str(np.mean(seq_lens)) + "\n")
print("Number of seqs: " + str(len(seq_lens)))
