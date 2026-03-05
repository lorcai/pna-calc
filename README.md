# pna-calc

`pna-calc` is a lightweight sequence analysis tool for PNA design workflows. Primarily designed for screening possible PCR blockers.

It reads candidate sequences (FASTA), computes sequence metrics, predicts
DNA/DNA nearest-neighbor Tm and PNA/DNA Tm, and writes results to CSV.
One FASTA file can contain one or many sequences; each input record is
processed independently.

## What It Does

For each sequence, `pna-calc` reports:

- `dna_tm`: DNA/DNA nearest-neighbor melting temperature from Biopython
  `MeltingTemp.Tm_NN`.
- `pna_tm`: PNA/DNA melting temperature using the Giesen et al. (1998)
  empirical correction on top of `dna_tm`.
- Sequence-derived metrics:
  - length
  - base counts
  - GC content
  - purine content
  - longest purine stretch
  - reverse complement
  - internal self-complementarity score

## Installation

### Option 1: Using an existing Python environment

```bash
pip install -r requirements.txt
```

### Option 2: Conda environment (example)

```bash
conda create -n pna_calc python=3.11 -y
conda activate pna_calc
pip install -r requirements.txt
```

## CLI Usage

Run:

```bash
python -m pna_calc.cli <input.fasta> -o <output.csv> [options]
```

### Mandatory CLI parameters

- Positional `fasta`: input FASTA file path.
- `-o`, `--output`: output CSV path.

The input FASTA can contain multiple sequences; output CSV has one row per
sequence.

### Defaults

- Profile default: `pnatool`.
- Profile tuned to be as similar to [PNA tool](https://pnabio.com/pna-tool/) as possible
- No override arguments are applied unless explicitly passed.

### Profile selection

- `-p`, `--profile {allawi_1997,pnatool}`

Profiles provide preset parameter sets for DNA melting temperature calculations.
Profile definitions live in `pna_calc/profiles.py`.

### Experimental overrides (optional)

You can override profile values directly from CLI:

- `--dnac1` (float, nM)
- `--dnac2` (float, nM)
- `--na` (float, mM)
- `--k` (float, mM)
- `--tris` (float, mM)
- `--mg` (float, mM)
- `--dntps` (float, mM)
- `--saltcorr` (int, 0-7)
- `--selfcomp` (flag)

### Important constraint (`dnac1 >= dnac2`)

`dnac1` must be greater than or equal to `dnac2`. This follows Biopython
`Tm_NN` conventions, where `dnac1` is the higher-concentration strand.

If `dnac1 < dnac2`, the CLI exits with an error.

### CLI examples

Default profile (`pnatool`):

```bash
python -m pna_calc.cli test/test.fasta -o results.csv
```

Use the `allawi_1997` profile:

```bash
python -m pna_calc.cli tests/test.fasta -o results.csv --profile allawi_1997
```

Use profile plus explicit overrides:

```bash
python -m pna_calc.cli tests/test.fasta -o results.csv \
  --profile pnatool \
  --dnac1 6000 \
  --dnac2 100 \
  --saltcorr 0
```

## Python API Usage

### Import core functions

```python
from pna_calc.tm import dna_nn_tm, pna_tm
from pna_calc.analyze import analyze_sequence, analyze_batch
from pna_calc.profiles import PROFILES
```

### Mandatory API inputs

- `dna_nn_tm(seq, ...)` requires `seq`.
- `pna_tm(seq, ...)` requires `seq` (validated for lengths 6-30 nt).
- `analyze_sequence(seq, ...)` requires `seq`.
- `analyze_batch(records, ...)` requires an iterable of records with at least
  `{"sequence": ...}`.

### API examples

Simple Tm calls:

```python
dna = dna_nn_tm("CAGTCCAGTT", profile="pnatool")
pna = pna_tm("CAGTCCAGTT", profile="pnatool")
```

Tm call with overrides:

```python
dna = dna_nn_tm(
    "CAGTCCAGTT",
    profile="allawi_1997",
    dnac1=200000,
    dnac2=200000,
    Na=1000,
    saltcorr=0,
)
```

Full sequence analysis:

```python
result = analyze_sequence(
    "CAGTCCAGTT",
    seq_id="candidate_1",
    profile="pnatool",
    tm_overrides={"dnac1": 6000, "dnac2": 100, "saltcorr": 0},
)

print("DNA_tm:", result["dna_tm"])
print( "PNA_tm:", result["pna_tm"])
```

## Output Format

CSV output includes one row per input sequence with at least:

- `id`
- `sequence`
- `dna_tm`
- `pna_tm`
- metric columns from `pna_calc/metrics.py`

## Profiles

Built-in profile registry is in `pna_calc/profiles.py`:

- `allawi_1997`
- `pnatool`

Each profile is an `NNProfile` object that maps directly to Biopython
`Tm_NN` parameters. Override values can be supplied via CLI flags or API
keyword arguments.

## Notes

- This tool currently assumes perfect complementarity for Tm calculations.
- Mismatch handling via `c_seq` is not yet exposed.

## References

- Biopython `MeltingTemp` (`Tm_NN` API used for DNA Tm prediction and parameter definitions):
  [Bio.SeqUtils.MeltingTemp module (Biopython 1.79)](https://biopython.org/docs/1.79/api/Bio.SeqUtils.MeltingTemp.html)

## To do

- Fully replicate Allawi et al 1997 Tm preds
- Better approximate PNAtool
- Use [PoacV9_01](https://apsjournals.apsnet.org/doi/10.1094/PBIOMES-05-20-0040-TA) as check example
 