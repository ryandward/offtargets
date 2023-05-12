import sys
import gzip
import zlib
import argparse
import configparser
from urllib.request import urlretrieve
from Bio import Entrez
from Bio import SeqIO
from rich.progress import Progress, BarColumn, TextColumn
from rich.table import Table
import os
from rich.progress import (
    Progress,
    BarColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)
from rich.console import Console

DESCRIPTION_WIDTH = 40

CONFIG_FILE = 'config.ini'


console = Console()

def get_email():
    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)

    if not config.has_section('Entrez'):
        config.add_section('Entrez')

    if not config.has_option('Entrez', 'email'):
        console.print("NCBI requires an email address for its API usage.", style="cyan")
        email = input("Please enter your email address: ")
        config.set('Entrez', 'email', email)

        with open(CONFIG_FILE, 'w') as configfile:
            config.write(configfile)
    else:
        email = config.get('Entrez', 'email')
        console.print(f"Using stored email address: {email}", style="cyan")

    return email


# Set your email for NCBI Entrez
Entrez.email = get_email()

from urllib.error import URLError

def download_sequences_by_accession(accession_numbers, source):
    assembly_data = []

    with Progress(
        TextColumn(f"[bold cyan]{{task.description:<{DESCRIPTION_WIDTH}}}", justify="left"),
        BarColumn(bar_width=None),
        TextColumn("{task.percentage:>3.0f}%", justify="right"),
        transient=False,
    ) as progress:
        fetch_task = progress.add_task("[1/3] Fetching assembly accessions...", total=len(accession_numbers))

        for accession_number in accession_numbers:
            db = "nuccore" if len(accession_number) == 11 else "assembly"
            handle = Entrez.esearch(db=db, term=accession_number)
            record = Entrez.read(handle)

            if db == "nuccore":
                assembly_id = record["IdList"][0]
                handle = Entrez.esummary(db="assembly", id=assembly_id)
                summary = Entrez.read(handle)
                assembly_data.append({
                    "accession": summary["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"],
                    "ftp_path": summary["DocumentSummarySet"]["DocumentSummary"][0][f"FtpPath_{source}"]
                })
            else:
                for assembly_id in record["IdList"]:
                    handle = Entrez.esummary(db="assembly", id=assembly_id)
                    summary = Entrez.read(handle)
                    assembly_data.append({
                        "accession": summary["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"],
                        "ftp_path": summary["DocumentSummarySet"]["DocumentSummary"][0][f"FtpPath_{source}"]
                    })

            progress.update(fetch_task, advance=1, refresh=True)

    assemblies_sequences = {}

    with Progress(
        TextColumn(f"[bold cyan]{{task.description:<{DESCRIPTION_WIDTH}}}", justify="left"),
        BarColumn(bar_width=None),
        TextColumn("{task.percentage:>3.0f}%", justify="right"),
        transient=False,
    ) as progress:
        download_task = progress.add_task("[2/3] Downloading sequences...", total=len(assembly_data))

        for assembly in assembly_data:
            ftp_path = assembly["ftp_path"]
            if ftp_path:
                genbank_file_name = ftp_path.split("/")[-1] + f"_genomic.gbff.gz"
                genbank_url = f"{ftp_path}/{genbank_file_name}"

                try:
                    compressed_file, _ = urlretrieve(genbank_url)
                except URLError as e:
                    console.print(
                        f"[bold red]Error while downloading sequence for accession {assembly['accession']}: {str(e)}"
                    )
                    continue

                with gzip.open(compressed_file, 'rt') as compressed_handle:
                    try:
                        sequences = list(SeqIO.parse(compressed_handle, "genbank"))
                    except ValueError as e:
                        console.print(
                            f"[bold red]Error while parsing sequence for accession {assembly['accession']}: {str(e)}"
                        )
                        continue

                    assemblies_sequences[assembly['accession']] = sequences

            progress.update(download_task, advance=1, refresh=True)

    return assemblies_sequences


def save_sequences_to_file(assemblies_sequences):
    total_sequences = sum(len(sequences) for sequences in assemblies_sequences.values())
    with Progress(
        TextColumn(f"[bold cyan]{{task.description:<{DESCRIPTION_WIDTH}}}", justify="left"),
        BarColumn(bar_width=None),
        TextColumn("{task.percentage:>3.0f}%", justify="right"),
        transient=False,
    ) as progress:
        save_task = progress.add_task("[3/3] Saving sequences...", total=total_sequences)

        header = "AssemblyAccession\tOrganism\n"
        with open("organisms.tsv", "w") as organisms_file:
            organisms_file.write(header)

        for assembly_accession, sequences in assemblies_sequences.items():
            file_name = f"{assembly_accession}.gb"

            # Check that all organisms in the assembly are the same
            organisms = {seq.annotations["organism"] for seq in sequences}
            if len(organisms) == 1:
                organism = organisms.pop()
                with open("organisms.tsv", "a") as organisms_file:
                    organisms_file.write(f"{assembly_accession}\t{organism}\n")

            # Overwrite old file
            with open(file_name, "w") as output_handle:
                for seq in sequences:
                    SeqIO.write(seq, output_handle, "genbank")
                    progress.update(save_task, advance=1)

        progress.update(save_task, completed=total_sequences)

def main():
    parser = argparse.ArgumentParser(description="Download nucleotide sequences from NCBI")
    parser.add_argument("accession_numbers", metavar="N", type=str, nargs="+",
                        help="one or more nucleotide accession numbers")
    parser.add_argument("--source", choices=["RefSeq", "GenBank"], default="GenBank",
                        help="choose between RefSeq and GenBank versions of the record (default: GenBank)")

    args = parser.parse_args()

    assemblies_sequences = download_sequences_by_accession(args.accession_numbers, args.source)
    save_sequences_to_file(assemblies_sequences)

    console.print("Please check 'organisms.tsv' for a list of the accession numbers and the organism.", style="bold yellow")

if __name__ == "__main__":
    main()