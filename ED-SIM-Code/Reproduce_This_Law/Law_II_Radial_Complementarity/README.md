# Reproduce: Law II — Radial Complementarity Law

## 1. Formal Law Statement

For N >= 4, the IC fractions under the inward and outward gates sum to a
near-constant: IC(inward) + IC(outward) = 100 +/- 4%.

## 2. What This Harness Verifies

- Computes IC(inward) and IC(outward) for each N = 3-7
- Computes the complementarity sum for each N
- Verifies that for N >= 4, |sum - 100%| < 5%
- Verifies the mean sum across N = 4-7 falls in [98%, 106%]
- Flags N = 3 as the expected outlier (Law V interaction)

## 3. Files Loaded

| File | Gate | N |
|------|------|---|
| `n{3,4,5,6,7}_inward_sweep.json` | inward | 3-7 |
| `n{3,4,5,6,7}_outward_sweep.json` | outward | 3-7 |

## 4. Metrics Computed

- `ic_inward[N]`: IC fraction under inward gate
- `ic_outward[N]`: IC fraction under outward gate
- `comp_sum[N]`: IC(inward) + IC(outward)
- `mean_sum_N4_7`: mean of comp_sum for N = 4, 5, 6, 7
- `max_deviation_N4_7`: max |comp_sum - 100%| for N >= 4

## 5. Pass/Fail Criteria

| Criterion | Threshold | Type |
|-----------|-----------|------|
| comp_sum[N] for N >= 4 | within [95%, 105%] | Range |
| mean_sum across N = 4-7 | within [98%, 106%] | Range |
| N = 3 sum flagged as outlier | sum > 115% | Expected |

## 6. Expected Outputs

- `expected_outputs/metrics.json`: canonical sums per N
- `plots/law_II_complementarity.png`: bar chart of sums with 100% line

## 7. How to Run

```bash
cd "ED Research/ED Simulations/Reproduce_This_Law/Law_II_Radial_Complementarity"
python run_law_II.py
```

## 8. Failure Modes

- Sum deviating > 10% from 100% at any N >= 4
- Systematic drift in sum with increasing N
- N = 3 sum falling below 115% (would invalidate Law V interaction)

## 9. Related Laws

- [Law I](../Law_I_Binary_Partition/): Binary partition within which
  complementarity operates.
- [Law V](../Law_V_N3_Total_Collapse/): Explains the N=3 outlier.
- [Law III](../Law_III_Monotone_Gate_Ordering/): Individual trends whose
  sum Law II constrains.
