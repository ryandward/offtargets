from Bio.Seq import Seq
import os
import sys
import subprocess
from Bio import SeqIO
import gzip
import pysam
from rich.console import Console
from Bio.SeqFeature import CompoundLocation

console = Console()

# Utility function to open files
def open_file(file, mode):
    return gzip.open(file, mode) if file.endswith(".gz") else open(file, mode)

# Create a working directory
def create_working_directory(dir_name="working_directory"):
    os.makedirs(dir_name, exist_ok=True)
    return dir_name

def convert_genbank_to_fasta(genbank_file, fasta_file):
    # handle errors, if the file does not exist using rich.console.Console
    if not os.path.exists(genbank_file):
        console.print(f"[red]File \"{genbank_file}\" does not exist![/red]")
        sys.exit(1)

    # if it is gzipped, open with gzip.open otherwise open with open
    if genbank_file.endswith(".gz"):
        with gzip.open(genbank_file, "rt") as input_handle, open(fasta_file, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "genbank"):
                SeqIO.write(record, output_handle, "fasta")
    else:
        with open(genbank_file, "rt") as input_handle, open(fasta_file, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "genbank"):
                SeqIO.write(record, output_handle, "fasta")

# Convert fasta to fake fastq
def convert_fasta_to_fake_fastq(fasta_file, fastq_file):
    # handle errors, if the file does not exist using rich.console.Console
    if not os.path.exists(fasta_file):
        console.print(f"[red]File \"{fasta_file}\" does not exist![/red]")
        sys.exit(1)

    # if it is gzipped, open with gzip.open otherwise open with open
    if fasta_file.endswith(".gz"):
        with gzip.open(fasta_file, "rt") as input_handle, open(fastq_file, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                SeqIO.write(record, output_handle, "fastq")
    else:
        with open(fasta_file, "rt") as input_handle, open(fastq_file, "w") as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                SeqIO.write(record, output_handle, "fastq")

# This function runs Bowtie, which is a tool for aligning short DNA sequences
def run_bowtie(sgrna_fastq, genome_fasta, output_file):

    index_prefix = "genome_index"

    with open(os.devnull, "w") as devnull:
        subprocess.run(
            ["bowtie-build", genome_fasta, index_prefix],
            stdout=devnull,
            stderr=devnull,
        )
        subprocess.run(
            ["bowtie", "-v", "1", "-S", index_prefix, sgrna_fastq, output_file],
            stdout=devnull,
            stderr=devnull,
        )
        # Remove the index files, which names start with the index_prefix and end with ebwt
        for file in os.listdir("."):
            if file.startswith(index_prefix) and file.endswith(".ebwt"):
                os.remove(file)

def create_locus_map(genbank_file):
    locus_map = {}
    with open_file(genbank_file, "rt") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            for feature in record.features:
                if feature.type == "gene":
                    # Check if the feature location is a compound location
                    if isinstance(feature.location, CompoundLocation):
                        # Iterate through all parts of the compound location
                        for part_location in feature.location.parts:
                            for position in range(int(part_location.start), int(part_location.end) + 1):
                                locus_map[(record.id, position)] = (
                                    feature.qualifiers.get("locus_tag", [None])[0],
                                    int(part_location.start),
                                    int(part_location.end),
                                    feature.strand,
                                )
                    else:
                        for position in range(int(feature.location.start), int(feature.location.end) + 1):
                            locus_map[(record.id, position)] = (
                                feature.qualifiers.get("locus_tag", [None])[0],
                                int(feature.location.start),
                                int(feature.location.end),
                                feature.strand,
                            )
    return locus_map

# Reconstruct the target sequence using the cigar string
def reconstruct_t_seq(read):
    reference_seq = read.get_reference_sequence()
    return reference_seq

# Parse the SAM output and create a TSV file
def parse_sam_output(sam_file, locus_map, output_tsv):
    samfile = pysam.AlignmentFile(sam_file, "r")

    with open(output_tsv, "w") as tsv_out:
        # Write header
        header = [
            "q_name",
            "q_seq",
            "q_len",
            "t_locus_tag",
            "t_chromosome",
            "t_seq",
            "diff_len",
            "diff",
            "coord",
            "offset",
            "q_dir",
            "t_dir",
        ]
        tsv_out.write("\t".join(header) + "\n")

        for read in samfile.fetch():
            if read.is_unmapped:
                continue

            # Extract read information
            q_name = read.query_name
            q_seq = read.query_sequence
            q_len = read.query_length
            t_chromosome = read.reference_name
            t_start = read.reference_start
            t_end = read.reference_end
            t_seq = reconstruct_t_seq(read)

            q_dir = "F" if not read.is_reverse else "R"

            # Find the corresponding feature by iterating through the entire aligned region
            t_locus_tag, feature_start, feature_end, feature_strand = None, None, None, None
            for position in range(t_start, t_end + 1):
                t_locus_tag, feature_start, feature_end, feature_strand = locus_map.get(
                    (t_chromosome, position), (None, None, None, None)
                )
                if t_locus_tag is not None:
                    break

            t_dir = "F" if feature_strand == 1 else "R" if feature_strand == -1 else None

            # Reverse complement the q_seq and t_seq if orientation is reverse
            if q_dir == "R":
                q_seq = str(Seq(q_seq).reverse_complement())
                t_seq = str(Seq(t_seq).reverse_complement())

            # Calculate the number of nucleotide differences and the diff string
            num_nt_diff = sum(1 for a, b in zip(q_seq, t_seq) if a != b)
            
            # The diff string is a comma-separated list of nucleotide differences
            # Each nucleotide difference is in the format of "t_nt"+"q_pos" +"q_nt"
            diff = ",".join(
                f"{a}{i + 1}{b}"
                for i, (a, b) in enumerate(zip(t_seq, q_seq))
                if a != b
            )

            # if no diff, then diff is a dash
            if diff == "":
                diff = "-"

            # Calculate the coordinate string
            # The coordinate string is the start and end position of the aligned region, 1-based.
            coord = f"{t_start + 1}-{t_end + 1}"

            # Calculate the offset and adjust 
            if feature_start is not None:
                if t_dir == "F":
                    offset = t_start - feature_start + 1  # Add 1 to fix the off-by-1 error
                elif t_dir == "R":
                    offset = feature_end - t_end - 1  # Subtract 1 to fix the off-by-2 error
            else:
                offset = None

            row = [
                q_name,
                q_seq,
                str(q_len),
                str(t_locus_tag),
                t_chromosome,
                t_seq,
                str(num_nt_diff),
                diff,
                coord,
                str(offset),
                q_dir,
                str(t_dir) if t_dir is not None else "",
            ]
            tsv_out.write("\t".join(row) + "\n")
    
    # Close the sam file
    samfile.close()

    #then delete the sam file
    os.remove(sam_file)

# Run the entire pipeline
def main(sgrna_file, genome_file):
    working_dir = create_working_directory()
    genome_fasta = os.path.join(
        working_dir, os.path.splitext(os.path.basename(genome_file))[0] + ".fasta"
    )
    sgrna_fastq = os.path.join(
        working_dir, os.path.splitext(os.path.basename(sgrna_file))[0] + ".fastq"
    )
    output_folder = "results"
    output_file = os.path.join(
        output_folder,
        f"{os.path.splitext(os.path.basename(sgrna_file))[0]}_{os.path.splitext(os.path.basename(genome_file))[0]}.sam",
    )

    os.makedirs(output_folder, exist_ok=True)

    with console.status("[bold green][1/6] Converting GenBank to FASTA..."):
        convert_genbank_to_fasta(genome_file, genome_fasta)

    with console.status("[bold green][2/6] Converting FASTA to fake FASTQ..."):
        convert_fasta_to_fake_fastq(sgrna_file, sgrna_fastq)

    with console.status("[bold green][3/6] Running Bowtie..."):
        run_bowtie(sgrna_fastq, genome_fasta, output_file)

    with console.status("[bold green][4/6] Preparing output..."):   
        output_tsv = os.path.join(
            output_folder,
            f"{os.path.splitext(os.path.basename(sgrna_file))[0]}_{os.path.splitext(os.path.basename(genome_file))[0]}.tsv",
        )

    with console.status("[bold green][5/6] Generating map of genome..."):
        locus_map = create_locus_map(genome_file)

    with console.status("[bold green][6/6] Parsing SAM output..."):
        parse_sam_output(output_file, locus_map, output_tsv)

    console.print(f"[bold green]Process complete! Results saved in the output file: [white]{output_tsv}[bold green]")

# Entry point of the program
if __name__ == "__main__":
    if len(sys.argv) != 3:
        # print name of program using sys.argv[0] using f notation
        print(f"Usage: python {sys.argv[0]} <sgrna_fasta_file> <genome_gb_file>")
        sys.exit(1)

    sgrna_file = sys.argv[1]
    genome_file = sys.argv[2]
    main(sgrna_file, genome_file)