import sys
import gzip
import regex
from Bio import SeqIO
from Bio.Seq import Seq

# Load the list of potential off-targeting elements from a FASTA file
with open(sys.argv[1], 'r') as f:
    off_target_elements = [(record.description, str(record.seq)) for record in SeqIO.parse(f, "fasta")]

# Set a score cutoff (adjust as needed)
score_cutoff = 16  # only report alignments with a score of 16 or higher
max_mismatches = len(off_target_elements[0][1]) - score_cutoff

# Check if the position is within a gene
def within_gene(start, features):
    for feature in features:
        if feature.type == 'gene':
            gene_start = int(feature.location.start)
            gene_end = int(feature.location.end)
            if gene_start <= start < gene_end:
                gene_name = feature.qualifiers.get('gene', [''])[0]
                locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
                return True, gene_name, locus_tag
    return False, '', ''

# Get differences between query and target sequences
def get_differences(query, target):
    differences = []
    for idx, (q, t) in enumerate(zip(query, target)):
        if q != t:
            differences.append(f"{q}{idx+1}{t}")
    return ', '.join(differences)

# Print header for TSV output
header = ["query", "q_seq", "q_len", "off_target_seq", "chromosome", "position", "score", "in_gene", "gene_name", "locus_tag", "differences", "orientation"]
print('\t'.join(header))

# Iterate over the input GenBank files
for genome_file in sys.argv[2:]:
    # Read the genome file
    with gzip.open(genome_file, 'rt') as file:
        genome_records = SeqIO.parse(file, "genbank")

        # Iterate over the genome records and search for potential off-targets
        for record in genome_records:
            seq = str(record.seq)
            seq_rc = str(record.seq.reverse_complement())
            for query_description, off_target_element in off_target_elements:
                pattern = f"({off_target_element}){{s<={max_mismatches}}}"
                for orientation, sequence in [("forward", seq), ("reverse", seq_rc)]:
                    for match in regex.finditer(pattern, sequence):
                        score = len(off_target_element) - match.fuzzy_counts[0]
                        position = match.start()
                        in_gene, gene_name, locus_tag = within_gene(position, record.features)
                        differences = get_differences(off_target_element, match.group())
                        output = [query_description, off_target_element, str(len(off_target_element)), match.group(), record.id, str(position), str(score), str(in_gene), gene_name, locus_tag, differences, orientation]
                        print('\t'.join(output))
