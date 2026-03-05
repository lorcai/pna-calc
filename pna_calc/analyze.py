"""
analyze.py

Combines metrics and Tm into a unified analysis result.
"""

from .metrics import sequence_metrics
from .tm import pna_tm, dna_nn_tm

def analyze_sequence(
    seq: str,
    seq_id: str | None = None,
    profile: str | None = None,
    tm_overrides: dict | None = None,
) -> dict:
    metrics = sequence_metrics(seq)
    tm_kwargs = tm_overrides or {}
    tm_dna = dna_nn_tm(seq, profile=profile, **tm_kwargs)
    tm_pna = pna_tm(seq, profile=profile, **tm_kwargs)

    return {
        "id": seq_id,
        "sequence": seq,
        "dna_tm": round(tm_dna, 2),
        "pna_tm": round(tm_pna, 2),
        **metrics,
    }


def analyze_batch(
    records,
    profile: str | None = None,
    tm_overrides: dict | None = None,
):
    results = []
    for rec in records:
        results.append(
            analyze_sequence(
                rec["sequence"],
                rec.get("id"),
                profile=profile,
                tm_overrides=tm_overrides,
            )
        )
    return results

