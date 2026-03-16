# Reproduce This Law

Reproducibility harness for the seven ED-Arch Mechanism Laws (I-VII).

Each subdirectory contains everything needed to independently verify one
law against the canonical sweep data:
- A `README.md` documenting the law, metrics, and pass/fail criteria
- A `run_law_X.py` script that loads sweep data, computes metrics,
  compares against canonical values, generates plots, and exits with
  pass (0) or fail (1) status
- An `expected_outputs/metrics.json` with canonical reference values
- A `plots/` directory populated by the verification script

## Quick Start

```bash
cd "ED Research/ED Simulations"

# Run a single law
python Reproduce_This_Law/Law_I_Binary_Partition/run_law_I.py

# Run all seven laws
for d in Reproduce_This_Law/Law_*/; do
    name=$(basename "$d")
    script=$(ls "$d"run_law_*.py 2>/dev/null | head -1)
    if [ -n "$script" ]; then
        echo "=== $name ==="
        python "$script" && echo "PASS" || echo "FAIL"
        echo
    fi
done
```

## Law Index

| Dir | Law | Script | Key Criterion |
|-----|-----|--------|---------------|
| [Law_I_Binary_Partition](Law_I_Binary_Partition/) | Binary Partition | `run_law_I.py` | DECAY=0, PBC-corner=0 for radial gates |
| [Law_II_Radial_Complementarity](Law_II_Radial_Complementarity/) | Radial Complementarity | `run_law_II.py` | IC(inw)+IC(out) = 100+/-4% for N>=4 |
| [Law_III_Monotone_Gate_Ordering](Law_III_Monotone_Gate_Ordering/) | Monotone Gate-Ordering | `run_law_III.py` | Inward IC decreasing, outward increasing |
| [Law_IV_Tangent_Fragmentation](Law_IV_Tangent_Fragmentation/) | Tangent Fragmentation | `run_law_IV.py` | Tangent IC >= 7 bands for N>=5 |
| [Law_V_N3_Total_Collapse](Law_V_N3_Total_Collapse/) | N=3 Total Collapse | `run_law_V.py` | IC=100% at (inward, N=3) |
| [Law_VI_N4_Rotational_Degeneracy](Law_VI_N4_Rotational_Degeneracy/) | N=4 Rotational Degeneracy | `run_law_VI.py` | IC=0%, DECAY=77% at (tangent, N=4) |
| [Law_VII_Tangent_Decay_Exclusivity](Law_VII_Tangent_Decay_Exclusivity/) | Tangent DECAY Exclusivity | `run_law_VII.py` | DECAY>0 only at (tangent, N=4) |

## Data Dependencies

All scripts load sweep JSON files from the parent directory
(`ED Research/ED Simulations/`). Required files:

```
n{3,4,5,6,7}_inward_sweep.json
n{3,4,5,6,7}_tangent_sweep.json
n{3,4,5,6,7}_outward_sweep.json
```

Total: 15 files, ~5,700 simulation records.

## Requirements

- Python 3.8+
- `matplotlib` (for plot generation)
- `numpy`

No modifications to the simulation engine are required or permitted.
