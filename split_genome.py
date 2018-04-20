####
# Split genome by sequences

# Imports
import argparse
from Bio import SeqIO

# parser = argparse.ArgumentParser()
# parser.add_argument("genome_file", help="The genome file")
# parser.add_argument("out_dir", help="The output directory")
# args = parser.parse_args()
# genome_file = args.genome_file
# out_dir = args.out_dir
out_dir = "/home/kevin/resources/genomes/GRCh38_v27_gencode/chromosomes/"
genome_file = "/home/kevin/resources/genomes/GRCh38_v27_gencode/GRCh38.primary_assembly.genome.fa"

# Load the genome
print("\n#==== Loading genome ====# \n" + genome_file)
genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
print("#==== Genome loaded ====# \n")

for i in range(1,23):
    chr = "chr" + str(i)
    seq = genome_dict[chr].seq
    outfile = out_dir + chr + ".fa"
    nf = open(outfile, "w")
    nf.write(">" + chr + "\n")
    nf.write(str(seq))
    nf.write("\n")
    nf.close()



for chr in ["chrX", "chrY"]:
    seq = genome_dict[chr].seq
    outfile = out_dir + chr + ".fa"
    nf = open(outfile, "w")
    nf.write(">" + chr + "\n")
    nf.write(str(seq))
    nf.write("\n")
    nf.close()