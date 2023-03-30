import sys
import gzip
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq

# Load the list of potential off-targeting elements from a text file
with open(sys.argv[1], 'r') as f:
    off_target_elements = [line.strip() for line in f]

# Set a score cutoff (adjust as needed)
score_cutoff = 16  # only report alignments with a score of 16 or higher

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
header = ["query_seq", "off_target_seq", "record_id", "position", "score", "in_gene", "gene_name", "locus_tag", "differences", "orientation"]
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
            for i in range(len(seq) - len(off_target_elements[0]) + 1):
                subseq_fwd = seq[i:i+len(off_target_elements[0])]
                subseq_rev = seq_rc[i:i+len(off_target_elements[0])]
                for off_target_element in off_target_elements:
                    for subseq, orientation in [(subseq_fwd, "forward"), (subseq_rev, "reverse")]:
                        score = sum(1 for a, b in zip(off_target_element, subseq) if a == b)

                        if score >= score_cutoff:
                            if orientation == "reverse":
                                position = len(seq) - (i + len(off_target_element))
                            else:
                                position = i

                            in_gene, gene_name, locus_tag = within_gene(position, record.features)
                            differences = get_differences(off_target_element, subseq)
                            output = [off_target_element, subseq, record.id, str(position), str(score), str(in_gene), gene_name, locus_tag, differences, orientation]
                            print('\t'.join(output))

