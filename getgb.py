# this function downloads all nucleotide sequences, i.e., chromosomes, from NCBI, saves them as a genbank file.
# it takes one or more accession numbers as input, and returns all nucleotide sequences that are related by assembly accession.
# it also returns a list of all assembly accessions that were used to download the sequences.

def download_genbank(accession_numbers, output_file):
    assembly_accessions = []
    for accession_number in accession_numbers:
        print(f"Downloading {accession_number}...")
        handle = Entrez.efetch(db="nucleotide", id=accession_number, rettype="gbwithparts", retmode="text")
        records = list(SeqIO.parse(handle, 'genbank'))
        for record in records:
            assembly_accessions.append(record.annotations['assembly_accession'])
        SeqIO.write(records, output_file, 'genbank')
    return assembly_accessions
