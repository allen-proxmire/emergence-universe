# Reproduce: Law I — Binary Partition Law

## 1. Formal Law Statement

Under the two radial gates (inward and outward), the mechanism space
partitions into exactly two dynamical outcomes — inward-collapse (IC)
and outward-PBC — for all N >= 3. No DECAY and no PBC-corner events
occur.

## 2. What This Harness Verifies

This harness loads all 10 radial sweep files (5 inward + 5 outward,
N = 3-7) and checks that:
- Every record's mechanism is one of: `inward-collapse`, `outward-PBC`,
  or `other-late`
- Zero records have mechanism `DECAY`
- Zero records have mechanism `PBC-corner`
- The `other-late` fraction does not exceed 6% for any file

## 3. Files Loaded

| File | Gate | N |
|------|------|---|
| `n3_inward_sweep.json` | inward | 3 |
| `n4_inward_sweep.json` | inward | 4 |
| `n5_inward_sweep.json` | inward | 5 |
| `n6_inward_sweep.json` | inward | 6 |
| `n7_inward_sweep.json` | inward | 7 |
| `n3_outward_sweep.json` | outward | 3 |
| `n4_outward_sweep.json` | outward | 4 |
| `n5_outward_sweep.json` | outward | 5 |
| `n6_outward_sweep.json` | outward | 6 |
| `n7_outward_sweep.json` | outward | 7 |

## 4. Metrics Computed

For each of the 10 files:
- `ic_frac`: fraction of records with mechanism = `inward-collapse`
- `opbc_frac`: fraction with mechanism = `outward-PBC`
- `decay_count`: absolute count of `DECAY` records
- `pbc_corner_count`: absolute count of `PBC-corner` records
- `other_late_frac`: fraction of `other-late` records

## 5. Pass/Fail Criteria

| Criterion | Threshold | Type |
|-----------|-----------|------|
| DECAY count across all 10 files | == 0 | Exact |
| PBC-corner count across all 10 files | == 0 | Exact |
| other-late fraction per file | <= 0.06 | Upper bound |
| IC + OPBC + OL fractions sum per file | == 1.0 | Exact |

**Overall: PASS if and only if all four criteria hold for all 10 files.**

## 6. Expected Outputs

- `expected_outputs/metrics.json`: canonical DECAY=0, PBC-corner=0 for
  each (gate, N) pair
- `plots/law_I_mechanism_heatmap.png`: 2x5 heatmap of mechanism fractions
  (rows=gates, cols=N)

## 7. How to Run

```bash
cd "ED Research/ED Simulations/Reproduce_This_Law/Law_I_Binary_Partition"
python run_law_I.py
```

Exit code 0 = PASS, exit code 1 = FAIL.

## 8. Failure Modes

- A single DECAY or PBC-corner event under any radial gate at any N
- other-late fraction exceeding 6% at any (gate, N) pair
- An unexpected mechanism label not in the known set

## 9. Related Laws

- [Law VII](../Law_VII_Tangent_Decay_Exclusivity/): Strengthens the DECAY=0
  finding by showing DECAY is exclusive to (tangent, N=4).
- [Law II](../Law_II_Radial_Complementarity/): Quantifies the IC/OPBC split
  within Law I's binary partition.
