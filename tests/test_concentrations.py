#!/usr/bin/env python3

"""
Script to test the melting temperature of a sequence under different concentrations of DNA strands.
Reveals the different behaviors of the melting temperature calculation under different concentrations.
Some fail!
"""

from Bio.SeqUtils import MeltingTemp as mt

# Test sequence (modify if needed)
seq = "ACGTACGTACGTACGTACGT"

# Fixed thermodynamic parameters
common_params = dict(
    nn_table=mt.DNA_NN3,
    Na=1000,
    saltcorr=0,
    selfcomp=False,
    check=True,
    strict=True,
)

# Test scenarios (nM)
scenarios = [
    ("equal_high",        200000, 200000),
    ("equal_medium",      20000,  20000),
    ("equal_low",         2000,   2000),
#    ("primer_limiting",   2000,   200000), # THIS FAILS! dnac1 must be > dnac2
    ("template_limiting", 200000, 2000),
#    ("extreme_low_1",     1,      200000), # THIS FAILS
    ("extreme_low_2",     200000, 1),
#    ("very_asymmetric",   100,    100000), # THIS FAILS
    ("reversed",          100000, 100),
]

print(f"{'scenario':20s} {'dnac1(nM)':>12s} {'dnac2(nM)':>12s} {'Tm(C)':>10s}")

for name, dnac1, dnac2 in scenarios:

    tm = mt.Tm_NN(
        seq,
        dnac1=dnac1,
        dnac2=dnac2,
        **common_params
    )

    print(f"{name:20s} {dnac1:12,.0f} {dnac2:12,.0f} {tm:10.2f}")

