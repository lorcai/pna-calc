"""
profiles.py

Reaction profiles for Tm calculations.

This module defines experimental profiles (sets of thermodynamic parameters)
for DNA/DNA nearest-neighbor Tm calculations. Profiles encapsulate the
reaction conditions (salt concentrations, strand concentrations, etc.) that
affect melting temperature predictions.

Profiles are used by tm.py to configure Biopython's Tm_NN function.
"""

from dataclasses import dataclass, replace, asdict
from typing import Optional
from Bio.SeqUtils import MeltingTemp as mt


@dataclass(frozen=True)
class NNProfile:
    """
    Configuration profile for nearest-neighbor Tm calculations.

    This dataclass encapsulates all parameters that affect Tm_NN calculations,
    providing a clean abstraction over Biopython's Tm_NN function arguments.

    Attributes:
        nn_table: Thermodynamic NN values table (e.g., mt.DNA_NN3).
            Available tables:
            - DNA_NN1: Breslauer et al. (1986)
            - DNA_NN2: Sugimoto et al. (1996)
            - DNA_NN3: Allawi & SantaLucia (1997) (default)
            - DNA_NN4: SantaLucia & Hicks (2004)
        dnac1: Concentration of the higher concentrated strand [nM].
            Typically the primer (for PCR) or the probe. Default=25.
        dnac2: Concentration of the lower concentrated strand [nM].
            In PCR this is the template strand which concentration is typically
            very low and may be ignored (dnac2=0). In oligo/oligo hybridization
            experiments, dnac1 equals dnac2. Default=25.
            MELTING and Primer3Plus use k = [Oligo(Total)]/4 by default. To
            mimic this behaviour, divide [Oligo(Total)] by 2 and assign to both
            dnac1 and dnac2. E.g., Total oligo concentration of 50 nM means
            dnac1=25, dnac2=25.
        selfcomp: Is the sequence self-complementary? Default=False.
            If True, the primer is thought binding to itself, thus dnac2 is not
            considered.
        Na: Concentration of Na+ ions [mM]. Default=50.
        K: Concentration of K+ ions [mM]. Default=0.
        Tris: Concentration of Tris buffer [mM]. Default=0.
        Mg: Concentration of Mg2+ ions [mM]. Default=0.
        dNTPs: Concentration of dNTPs [mM]. Default=0.
        saltcorr: Salt correction method. Default=5.
            0 means no salt correction.
            Available methods (1-7) use different empirical formulas.
            See Bio.SeqUtils.MeltingTemp documentation for details.
    """

    nn_table: object  # mt.DNA_NN1, mt.DNA_NN2, etc.
    dnac1: float
    dnac2: float
    selfcomp: bool = False
    Na: float = 50.0
    K: float = 0.0
    Tris: float = 0.0
    Mg: float = 0.0
    dNTPs: float = 0.0
    saltcorr: int = 5

    def to_tm_nn_kwargs(self) -> dict:
        """
        Convert profile to keyword arguments for mt.Tm_NN.

        Returns:
            Dictionary of keyword arguments suitable for passing to mt.Tm_NN
        """
        return asdict(self)

    def with_overrides(self, **overrides) -> "NNProfile":
        """
        Create a new profile with specified overrides.

        Args:
            **overrides: Keyword arguments to override in the new profile.

        Returns:
            New NNProfile instance with overrides applied.

        Example:
            >>> profile = ALAWI_1997.with_overrides(dnac1=100000, Na=500)
        """
        return replace(self, **overrides)


# Built-in profiles

# Profile attempting to replicate Allawi et al. 1997 Tm values
# It uses Nearest neighbor parameters derived in the 1997 study AND the same chemical conditions
# Note: I think I can still get closer. Not convinced that is the actual concentration.
# Can also test the mismatched duplexes by enabling the c_seq.
ALAWI_1997 = NNProfile(
    nn_table=mt.DNA_NN3,  # The parameters derived in the 1997 study
    dnac1=200000,  # 2e-4 M -> 200,000 nM
    dnac2=200000,  # 2e-4 M
    Na=1000,  # 1M of Concentration of Na ions [mM]
    saltcorr=0,  # no salt correction
    selfcomp=False,
)

# Profile similar to PNA_tool (https://pnabio.com/pna-tool/)
# The calculated Tm is for the PNA-DNA hybrid at a PNA concentration of 5 uM
# and the target nucleic acid in significant excess
PNATOOL_LIKE = NNProfile(
    nn_table=mt.DNA_NN3,
    dnac1=5000,  # nM [PNA_tool 5uM]
    dnac2=0,  # nM (target in significant excess)
    saltcorr=0,  # no ionic correction
    # Alternative: saltcorr=3 (Schildkraut & Lifson 1965, Biopolymers 3: 195-208)
    # Alternative: Na=1000  # 1M of Concentration of Na ions [mM]
    selfcomp=False,
)

# Profile registry
# Maps string identifiers to profile objects for easy lookup
PROFILES = {
    "allawi_1997": ALAWI_1997,
    "pnatool": PNATOOL_LIKE,
}
