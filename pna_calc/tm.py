"""
tm.py

Melting temperature calculations.

This module is responsible for:
- DNA/DNA nearest-neighbor Tm calculation (Biopython)
- PNA/DNA empirical Tm prediction (Giesen et al., 1998)

No sequence metrics are computed here beyond what is strictly
required for the PNA model.
"""

from Bio.SeqUtils import MeltingTemp as mt
from .profiles import NNProfile, ALAWI_1997, PNATOOL_LIKE, PROFILES

# SantaLucia-based empirical coefficients for PNA/DNA duplexes
_C0 = 20.79
_C1 = 0.83
_C2 = -26.13
_C3 = 0.44

_VALID_BASES = {"A", "C", "G", "T"}


class TmError(Exception):
    """Base class for Tm-related errors."""


def _validate_sequence(seq: str) -> str:
    """
    Validate and normalize a DNA sequence.

    - Uppercases input
    - Ensures only A/C/G/T are present
    """
    seq = seq.upper()
    invalid = set(seq) - _VALID_BASES
    if invalid:
        raise TmError(f"Invalid DNA characters found: {sorted(invalid)}")
    return seq


def dna_nn_tm(
    seq: str,
    profile: str | NNProfile | None = None,
    **overrides
) -> float:
    """
    Compute DNA/DNA nearest-neighbor melting temperature.

    Uses:
    - DNA_NN3: values from Allawi & SantaLucia (1997) (default)
    - Assume 5 micromolar (uM) i.e. 5000 nanomolar (nM). 
    - No salt correction
    - Assumes perfect complementarity
    
    Relevant setup info:
    From https://pnabio.com/pna-tool/:
    "The calculated Tm is for the PNA-DNA hybrid at a PNA concentration of 5 uM and the target nucleic acid in significant excess"

    Another setup option:
    - Total ssDNA concentration: 4 µM
    - dnac1 = dnac2 = 2000 nM (Primer3 convention)
    
    Relevant argument info (From Bio.SeqUtils.MeltingTemp module documentation):
    - table: Thermodynamic NN values, eight tables are implemented: For DNA/DNA hybridizations:
        DNA_NN1: values from Breslauer et al. (1986)
        DNA_NN2: values from Sugimoto et al. (1996)
        DNA_NN3: values from Allawi & SantaLucia (1997) (default)
        DNA_NN4: values from SantaLucia & Hicks (2004)
    - dnac1: Concentration of the higher concentrated strand [nM]. Typically this will be the primer (for PCR) or the probe. Default=25.
    - dnac2: Concentration of the lower concentrated strand [nM]. In PCR this is the template strand which concentration is typically very low and may be ignored (dnac2=0). In oligo/oligo hybridization experiments, dnac1 equals dnac1. Default=25. MELTING and Primer3Plus use k = [Oligo(Total)]/4 by default. To mimic this behaviour, you have to divide [Oligo(Total)] by 2 and assign this concentration to dnac1 and dnac2. E.g., Total oligo concentration of 50 nM in Primer3Plus means dnac1=25, dnac2=25.
    - selfcomp: Is the sequence self-complementary? Default=False. If ‘True’ the primer is thought binding to itself, thus dnac2 is not considered.
    - Na, K, Tris, Mg, dNTPs: See method ‘Tm_GC’ for details. Defaults: Na=50, K=0, Tris=0, Mg=0, dNTPs=0.
    - saltcorr: See method ‘Tm_GC’. Default=5. 0 means no salt correction.

    Args:
        seq: DNA sequence string
        profile: Profile name (string), NNProfile object, or None.
            If None, defaults to ALAWI_1997 profile for backward compatibility.
            If string, must be a key in PROFILES registry.
            If NNProfile, uses that profile directly.
        **overrides: Optional keyword arguments to override profile defaults.
            Valid keys: dnac1, dnac2, selfcomp, Na, K, Tris, Mg, dNTPs, saltcorr, nn_table.

    Returns:
        Tm in Celsius
    """
    seq = _validate_sequence(seq)
    
    # Resolve profile
    if profile is None:
        # Default to PNATOOL_LIKE
        active_profile = PNATOOL_LIKE
    elif isinstance(profile, str):
        if profile not in PROFILES:
            raise TmError(
                f"Unknown profile '{profile}'. Available: {list(PROFILES.keys())}"
            )
        active_profile = PROFILES[profile]
    elif isinstance(profile, NNProfile):
        active_profile = profile
    else:
        raise TmError(
            f"profile must be str, NNProfile, or None, got {type(profile)}"
        )
    
    # Apply overrides if any
    if overrides:
        active_profile = active_profile.with_overrides(**overrides)

    if active_profile.dnac1 < active_profile.dnac2:
        raise TmError(
            "Invalid concentration setup: dnac1 must be greater than or equal "
            f"to dnac2 (got dnac1={active_profile.dnac1}, "
            f"dnac2={active_profile.dnac2})."
        )
    
    return mt.Tm_NN(
        seq,                    # These parameters attempt to replicate Allawi et al 1997
        #c_seq=cseq,
        **active_profile.to_tm_nn_kwargs()
    )

def pna_tm(
    seq: str,
    profile: str | NNProfile | None = None,
    **overrides
) -> float:
    """
    Predict PNA/DNA duplex melting temperature.

    Model (Giesen et al., 1998):
        Tm = 20.79
             + 0.83 * Tm_nnDNA
             - 26.13 * f_pyr
             + 0.44 * L

    Constraints:
        Validated primarily for ~6–30 nt sequences.

    Args:
        seq: DNA sequence string
        profile: Profile name (string), NNProfile object, or None.
            If None, defaults to PNATOOL_LIKE profile (5 uM PNA concentration).
            If string, must be a key in PROFILES registry.
            If NNProfile, uses that profile directly.
        **overrides: Optional keyword arguments to override profile defaults.
            Passed through to dna_nn_tm.

    Returns:
        Predicted Tm (Celsius)
    """
    seq = _validate_sequence(seq)
    L = len(seq)

    if L < 6 or L > 30:
        raise TmError(
            f"Sequence length {L} outside validated PNA Tm model range"
        )

    # Default to PNATOOL_LIKE profile for PNA calculations if not specified
    if profile is None:
        profile = PNATOOL_LIKE
    
    tm_nn = dna_nn_tm(seq, profile=profile, **overrides)
    f_pyr = (seq.count("C") + seq.count("T")) / L

    return (
        _C0
        + _C1 * tm_nn
        + _C2 * f_pyr
        + _C3 * L
    )

