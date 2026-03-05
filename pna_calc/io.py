"""
io.py

Input/output utilities and normalization.
"""

from Bio import SeqIO
import csv


def read_fasta(path: str):
    """
    Read sequences from a FASTA file.
    """
    for record in SeqIO.parse(path, "fasta"):
        yield {
            "id": record.id,
            "sequence": str(record.seq).upper().replace("U", "T"),
        }


def write_csv(results: list[dict], path: str):
    """
    Write analysis results to CSV.
    """
    if not results:
        return

    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=results[0].keys())
        writer.writeheader()
        for row in results:
            writer.writerow(row)

