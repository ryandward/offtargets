from Bio.Seq import Seq
import os
import sys
import subprocess
from Bio import SeqIO
from Bio import SeqFeature
import gzip
import shutil
import re
import pysam

# this function creates a working directory
def create_working_directory():
    if not os.path.exists("working_directory"):
        os.makedirs("working_directory")
    return "working_directory"

# this function converts genbank to fasta
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

# this function converts fasta to fake fastq
def convert_fasta_to_fake_fastq(fasta_file, fastq_file):
    with open(fastq_file, "w") as output_handle:
        for record in SeqIO.parse(fasta_file, "fasta"):
            record.letter_annotations["phred_quality"] = [40] * len(record)
            SeqIO.write(record, output_handle, "fastq")

# this function runs bowtie, which is a tool for aligning short DNA sequences
def run_bowtie(sgrna_fastq, genome_fasta, output_file):
    index_prefix = "genome_index"
    subprocess.run(["bowtie-build", genome_fasta, index_prefix])
    subprocess.run(["bowtie", "-v", "1", "-S",
                   index_prefix, sgrna_fastq, output_file])

# this function extracts the sequence, start, end, and strand of a gene
def create_locus_map(genbank_file):
    locus_map = {}
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type == "gene":
                for position in range(int(feature.location.start), int(feature.location.end) + 1):
                    locus_map[(record.id, position)] = (feature.qualifiers.get("locus_tag", [None])[
                        0], int(feature.location.start), int(feature.location.end), feature.strand)
    return locus_map

# this function reconstructs the sequence of the target, using cigar string
def reconstruct_t_seq(read):
    reference_seq = read.get_reference_sequence()
    return reference_seq

# this function parses the bowtie output, and creates a tsv file
# the tsv file contains the following columns:
# q_name: the name of the query sequence
# q_sequence: the query sequence
# q_len: the length of the query sequence
# t_locus_tag: the locus tag of the target gene
# t_chromosome: the chromosome of the target gene
# t_sequence: the sequence of the target gene
# diff_len: the number of nucleotide differences between the query sequence and the target sequence
# diff: the nucleotide differences between the query sequence and the target sequence
# coord: the start and end position of the target sequence
# offset: the offset of the query sequence from the start of the target sequence
# q_dir: the orientation of the query sequence
# t_dir: the orientation of the target sequence

def parse_sam_output(sam_file, locus_map, output_tsv):
    samfile = pysam.AlignmentFile(sam_file, "r")

    with open(output_tsv, 'w') as tsv_out:
        header = ["q_name", "q_sequence", "q_len", "t_locus_tag", "t_chromosome",
                  "t_sequence", "diff_len", "diff", "coord", "offset", "q_dir", "t_dir"]
        tsv_out.write("\t".join(header) + "\n")

        for read in samfile.fetch():
            if read.is_unmapped:
                continue

            q_name = read.query_name
            q_seq = read.query_sequence
            q_len = read.query_length
            t_chromosome = read.reference_name
            t_start = read.reference_start
            t_end = read.reference_end
            t_sequence = reconstruct_t_seq(read)

            q_dir = "F" if not read.is_reverse else "R"
            t_locus_tag, feature_start, feature_end, feature_strand = locus_map.get(
                (t_chromosome, t_start), (None, None, None, None))

            t_dir = "F" if feature_strand == 1 else "R" if feature_strand == -1 else None

            # Reverse complement the q_seq and t_sequence if orientation is reverse
            if q_dir == "R":
                q_seq = str(Seq(q_seq).reverse_complement())
                t_sequence = str(Seq(t_sequence).reverse_complement())

            # Calculate the number of nucleotide differences and the diff string
            num_nt_diff = sum(1 for a, b in zip(q_seq, t_sequence) if a != b)
            diff = ",".join(f"{a}{i+1}{b}" for i, (a, b)
                            in enumerate(zip(q_seq, t_sequence)) if a != b)

            coord = f"{t_start}-{t_start + q_len - 1}"

            # Calculate the offset and adjust for 1-based position
            if feature_start is not None:
                if t_dir == "F":
                    offset = t_start - feature_start + 1
                elif t_dir == "R":
                    offset = feature_end - (t_start + q_len + 1)
            else:
                offset = None

            row = [q_name, q_seq, str(q_len), str(t_locus_tag), t_chromosome, t_sequence, str(
                num_nt_diff), diff, coord, str(offset), q_dir, str(t_dir) if t_dir is not None else '']
            tsv_out.write("\t".join(row) + "\n")

    samfile.close()

# this function runs the entire pipeline
def main(sgrna_file, genome_file):
    working_dir = create_working_directory()
    genome_fasta = os.path.join(working_dir, os.path.splitext(
        os.path.basename(genome_file))[0] + ".fasta")
    sgrna_fastq = os.path.join(working_dir, os.path.splitext(
        os.path.basename(sgrna_file))[0] + ".fastq")
    output_folder = "results"
    output_file = os.path.join(
        output_folder, f"{os.path.splitext(os.path.basename(sgrna_file))[0]}_{os.path.splitext(os.path.basename(genome_file))[0]}.sam")

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    convert_genbank_to_fasta(genome_file, genome_fasta)
    convert_fasta_to_fake_fastq(sgrna_file, sgrna_fastq)
    run_bowtie(sgrna_fastq, genome_fasta, output_file)
    output_tsv = os.path.join(
        output_folder, f"{os.path.splitext(os.path.basename(sgrna_file))[0]}_{os.path.splitext(os.path.basename(genome_file))[0]}.tsv")
    locus_map = create_locus_map(genome_file)
    parse_sam_output(output_file, locus_map, output_tsv)

# this function is the entry point of the program
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python bowtie_wrapper.py <sgrna_fasta_file> <genome_gb_file>")
        sys.exit(1)

    sgrna_file = sys.argv[1]
    genome_file = sys.argv[2]
    main(sgrna_file, genome_file)