"""
metrics.py

Sequence-derived metrics and heuristic design checks.
Pure string logic — no thermodynamics.
"""

_VALID_BASES = {"A", "C", "G", "T"}


def _validate_sequence(seq: str) -> str:
    seq = seq.upper()
    invalid = set(seq) - _VALID_BASES
    if invalid:
        raise ValueError(f"Invalid DNA characters: {sorted(invalid)}")
    return seq


def reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp[b] for b in reversed(seq))


def longest_purine_stretch(seq: str) -> int:
    max_run = run = 0
    for b in seq:
        if b in {"A", "G"}:
            run += 1
            max_run = max(max_run, run)
        else:
            run = 0
    return max_run


def self_complementarity(seq: str, min_match: int = 2) -> int:
    """
    Detect longest internal reverse-complement match.
    """
    rc = reverse_complement(seq)
    n = len(seq)
    longest = 0

    for i in range(n):
        for j in range(n):
            k = 0
            while (
                i + k < n
                and j + k < n
                and seq[i + k] == rc[j + k]
            ):
                k += 1
            longest = max(longest, k)

    return longest


def sequence_metrics(seq: str) -> dict:
    """
    Compute all non-thermodynamic metrics for a PNA candidate.
    """
    seq = _validate_sequence(seq)
    L = len(seq)

    counts = {b: seq.count(b) for b in "ACGT"}

    return {
        "length": L,
        "base_counts": counts,
        "gc_content": (counts["G"] + counts["C"]) / L * 100,
        "purine_content": (counts["A"] + counts["G"]) / L * 100,
        "purine_stretch": longest_purine_stretch(seq),
        "reverse_complement": reverse_complement(seq),
        "self_complementarity": self_complementarity(seq),
    }

