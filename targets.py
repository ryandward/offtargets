import os
import sys
import subprocess
from Bio import SeqIO
import gzip
import shutil

def create_working_directory():
    if not os.path.exists("working_directory"):
        os.makedirs("working_directory")
    return "working_directory"

def convert_genbank_to_fasta(genbank_file, fasta_file):
    with open(fasta_file, "w") as output_handle:
        if genbank_file.endswith(".gz"):
            with gzip.open(genbank_file, "rt") as input_handle:
                for record in SeqIO.parse(input_handle, "genbank"):
                    SeqIO.write(record, output_handle, "fasta")
        else:
            with open(genbank_file, "r") as input_handle:
                for record in SeqIO.parse(input_handle, "genbank"):
                    SeqIO.write(record, output_handle, "fasta")

def convert_fasta_to_fake_fastq(fasta_file, fastq_file):
    with open(fastq_file, "w") as output_handle:
        for record in SeqIO.parse(fasta_file, "fasta"):
            record.letter_annotations["phred_quality"] = [40] * len(record)
            SeqIO.write(record, output_handle, "fastq")

def run_bowtie(sgrna_fastq, genome_fasta, output_file):
    index_prefix = "genome_index"
    subprocess.run(["bowtie-build", genome_fasta, index_prefix])
    subprocess.run(["bowtie", "-v", "1", "-S", index_prefix, sgrna_fastq, output_file])

def main(sgrna_file, genome_file):
    working_dir = create_working_directory()
    genome_fasta = os.path.join(working_dir, os.path.splitext(os.path.basename(genome_file))[0] + ".fasta")
    sgrna_fastq = os.path.join(working_dir, os.path.splitext(os.path.basename(sgrna_file))[0] + ".fastq")
    output_folder = "results"
    output_file = os.path.join(output_folder, f"{os.path.splitext(os.path.basename(sgrna_file))[0]}_{os.path.splitext(os.path.basename(genome_file))[0]}.sam")

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    convert_genbank_to_fasta(genome_file, genome_fasta)
    convert_fasta_to_fake_fastq(sgrna_file, sgrna_fastq)
    run_bowtie(sgrna_fastq, genome_fasta, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python bowtie_wrapper.py <sgrna_fasta_file> <genome_gb_file>")
        sys.exit(1)

    sgrna_file = sys.argv[1]
    genome_file = sys.argv[2]
    main(sgrna_file, genome_file)
