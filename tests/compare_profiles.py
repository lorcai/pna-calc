#!/usr/bin/env python3
"""
Compare Tm for test.fasta sequences under different profiles.
Run: python compare_profiles.py
"""

from pna_calc.io import read_fasta
from pna_calc.tm import dna_nn_tm, pna_tm
from pna_calc.profiles import PROFILES

FASTA = "test.fasta"

def main():
    records = list(read_fasta(FASTA))
    print("Profiles:", list(PROFILES.keys()))
    print()

    for rec in records:
        seq_id = rec["id"]
        seq = rec["sequence"]
        print(f"--- {seq_id}: {seq} ---")

        for name, profile in PROFILES.items():
            try:
                tm_dna = dna_nn_tm(seq, profile=profile)
                tm_pna = pna_tm(seq, profile=profile)
                print(f"  {name}:  dna_nn_tm={tm_dna:.2f} °C   pna_tm={tm_pna:.2f} °C")
            except Exception as e:
                print(f"  {name}:  {e}")
        print()

if __name__ == "__main__":
    main()
