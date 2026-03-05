"""
cli.py

Command-line interface for the PNA tool.
"""

import argparse
from .io import read_fasta, write_csv
from .analyze import analyze_batch
from .profiles import PROFILES


def _collect_tm_overrides(args: argparse.Namespace) -> dict:
    overrides = {}
    if args.dnac1 is not None:
        overrides["dnac1"] = args.dnac1
    if args.dnac2 is not None:
        overrides["dnac2"] = args.dnac2
    if args.na is not None:
        overrides["Na"] = args.na
    if args.k is not None:
        overrides["K"] = args.k
    if args.tris is not None:
        overrides["Tris"] = args.tris
    if args.mg is not None:
        overrides["Mg"] = args.mg
    if args.dntps is not None:
        overrides["dNTPs"] = args.dntps
    if args.saltcorr is not None:
        overrides["saltcorr"] = args.saltcorr
    if args.selfcomp:
        overrides["selfcomp"] = True
    return overrides


def _validate_tm_arguments(
    parser: argparse.ArgumentParser,
    profile_name: str,
    overrides: dict,
) -> None:
    profile = PROFILES[profile_name].with_overrides(**overrides)

    for field in ("dnac1", "dnac2", "Na", "K", "Tris", "Mg", "dNTPs"):
        value = getattr(profile, field)
        if value < 0:
            parser.error(f"--{field.lower()} must be >= 0 (got {value})")

    if profile.dnac1 < profile.dnac2:
        parser.error(
            "Invalid concentrations: dnac1 must be >= dnac2 "
            f"(got dnac1={profile.dnac1}, dnac2={profile.dnac2})."
        )

    if not 0 <= profile.saltcorr <= 7:
        parser.error(
            f"--saltcorr must be between 0 and 7 (got {profile.saltcorr})"
        )


def main():
    parser = argparse.ArgumentParser(
        description="PNA sequence analysis tool"
    )
    parser.add_argument(
        "fasta",
        help="Input FASTA file"
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output CSV file"
    )
    parser.add_argument(
        "-p",
        "--profile",
        default="pnatool",
        choices=sorted(PROFILES.keys()),
        help="Experimental profile for Tm calculation",
    )
    parser.add_argument(
        "--dnac1",
        type=float,
        help="Higher concentration strand [nM]",
    )
    parser.add_argument(
        "--dnac2",
        type=float,
        help="Lower concentration strand [nM]",
    )
    parser.add_argument(
        "--na",
        type=float,
        help="Na+ concentration [mM]",
    )
    parser.add_argument(
        "--k",
        type=float,
        help="K+ concentration [mM]",
    )
    parser.add_argument(
        "--tris",
        type=float,
        help="Tris concentration [mM]",
    )
    parser.add_argument(
        "--mg",
        type=float,
        help="Mg2+ concentration [mM]",
    )
    parser.add_argument(
        "--dntps",
        type=float,
        help="dNTP concentration [mM]",
    )
    parser.add_argument(
        "--saltcorr",
        type=int,
        help="Salt correction method index [0-7]",
    )
    parser.add_argument(
        "--selfcomp",
        action="store_true",
        help="Treat sequence as self-complementary",
    )

    args = parser.parse_args()
    tm_overrides = _collect_tm_overrides(args)
    _validate_tm_arguments(parser, args.profile, tm_overrides)

    records = read_fasta(args.fasta)
    results = analyze_batch(
        records,
        profile=args.profile,
        tm_overrides=tm_overrides,
    )
    write_csv(results, args.output)


if __name__ == "__main__":
    main()

