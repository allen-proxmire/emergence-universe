# Reproduce: Law III — Monotone Gate-Ordering Law

## 1. Formal Law Statement

The inward gate's IC fraction decreases monotonically with N; the outward
gate's IC fraction increases monotonically with N. Inward IC always
exceeds outward IC at every N tested.

## 2. What This Harness Verifies

- Computes IC(inward) and IC(outward) for N = 3-7
- Checks strict monotonic decrease of IC(inward) with N
- Checks strict monotonic increase of IC(outward) with N
- Checks IC(inward) > IC(outward) at every N
- Reports the convergence trend (spread narrowing)

## 3. Files Loaded

| File | Gate | N |
|------|------|---|
| `n{3,4,5,6,7}_inward_sweep.json` | inward | 3-7 |
| `n{3,4,5,6,7}_outward_sweep.json` | outward | 3-7 |

## 4. Metrics Computed

- `ic_inward[N]`: IC fraction per N under inward gate
- `ic_outward[N]`: IC fraction per N under outward gate
- `inward_monotone`: boolean, True if strictly decreasing
- `outward_monotone`: boolean, True if strictly increasing
- `inward_dominates`: boolean, True if IC(inward) > IC(outward) at all N
- `spread[N]`: IC(inward) - IC(outward) per N

## 5. Pass/Fail Criteria

| Criterion | Threshold | Type |
|-----------|-----------|------|
| Inward IC strictly decreasing | all step diffs < 0 | Exact |
| Outward IC strictly increasing | all step diffs > 0 | Exact |
| Inward > Outward at every N | IC(inw) > IC(out) | Exact |

## 6. Expected Outputs

- `expected_outputs/metrics.json`: canonical IC sequences
- `plots/law_III_monotone_ordering.png`: dual-axis IC trend plot

## 7. How to Run

```bash
cd "ED Research/ED Simulations/Reproduce_This_Law/Law_III_Monotone_Gate_Ordering"
python run_law_III.py
```

## 8. Failure Modes

- Monotonicity reversal at any adjacent N pair
- Outward IC exceeding inward IC at some N (curve crossing)
- Spread widening at large N (divergence instead of convergence)

## 9. Related Laws

- [Law II](../Law_II_Radial_Complementarity/): Constrains the sum of the
  two monotone sequences.
- [Law V](../Law_V_N3_Total_Collapse/): The N=3 starting point of the
  inward decrease (100%).
