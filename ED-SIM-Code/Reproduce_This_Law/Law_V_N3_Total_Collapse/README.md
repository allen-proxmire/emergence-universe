# Reproduce: Law V — N=3 Inward Total Collapse Law

## 1. Formal Law Statement

At N=3 under the inward gate, every diameter produces inward-collapse:
IC = 100% across the full 385-point sweep. This is the unique (gate, N)
pair achieving total collapse.

## 2. What This Harness Verifies

- Loads the N=3 inward sweep and checks that every single record has
  mechanism = `inward-collapse`
- Verifies IC fraction is exactly 100% (385/385)
- Checks all other (gate, N) pairs to confirm none reach 100% IC

## 3. Files Loaded

Primary: `n3_inward_sweep.json`

Uniqueness check: all 15 sweep files `n{3,4,5,6,7}_{inward,tangent,outward}_sweep.json`

## 4. Metrics Computed

- `n3_inward_ic_count`: count of IC records in N=3 inward sweep
- `n3_inward_total`: total records in N=3 inward sweep
- `n3_inward_ic_frac`: IC fraction (must be 1.0000)
- `other_100_cells`: list of any other (gate, N) with IC = 100%

## 5. Pass/Fail Criteria

| Criterion | Threshold | Type |
|-----------|-----------|------|
| N=3 inward IC fraction | == 1.0 exactly | Exact |
| No other cell at 100% IC | other_100_cells empty | Exact |

## 6. Expected Outputs

- `expected_outputs/metrics.json`: canonical IC=385/385
- `plots/law_V_total_collapse.png`: IC fraction across all 15 cells,
  highlighting the (inward, N=3) cell

## 7. How to Run

```bash
cd "ED Research/ED Simulations/Reproduce_This_Law/Law_V_N3_Total_Collapse"
python run_law_V.py
```

## 8. Failure Modes

- A single non-IC event in the N=3 inward sweep at any diameter
- Another (gate, N) pair achieving 100% IC (removes uniqueness)

## 9. Related Laws

- [Law III](../Law_III_Monotone_Gate_Ordering/): Law V is the N=3
  endpoint of the inward monotone decrease.
- [Law II](../Law_II_Radial_Complementarity/): Explains the N=3
  complementarity outlier.
- [Law VI](../Law_VI_N4_Rotational_Degeneracy/): The opposite extreme
  (0% IC at tangent N=4).
