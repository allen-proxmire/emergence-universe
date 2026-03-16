# Reproduce: Law VII — Tangent DECAY Exclusivity Law

## 1. Formal Law Statement

DECAY events occur if and only if the gate is tangent and N=4. Across all
15 (gate, N) pairs (over 5,700 simulations), DECAY appears in exactly one
cell: (tangent, N=4) at 77.4%. The 14 remaining cells show zero DECAY.

## 2. What This Harness Verifies

- Loads all 15 sweep files
- Counts DECAY events in each (gate, N) pair
- Verifies exactly one cell has DECAY > 0
- Verifies that cell is (tangent, N=4)
- Verifies DECAY fraction at (tangent, N=4) exceeds 50%

## 3. Files Loaded

All 15 sweep files: `n{3,4,5,6,7}_{inward,tangent,outward}_sweep.json`

## 4. Metrics Computed

- `decay_count[gate][N]`: absolute DECAY count per cell
- `decay_frac[gate][N]`: DECAY fraction per cell
- `nonzero_decay_cells`: list of cells with DECAY > 0
- `total_zero_decay_points`: total simulation points with DECAY = 0

## 5. Pass/Fail Criteria

| Criterion | Threshold | Type |
|-----------|-----------|------|
| Exactly one nonzero DECAY cell | count == 1 | Exact |
| That cell is (tangent, N=4) | identity match | Exact |
| DECAY fraction at (tangent, N=4) | >= 0.50 | Lower bound |
| All other 14 cells | DECAY == 0 | Exact |

## 6. Expected Outputs

- `expected_outputs/metrics.json`: canonical DECAY counts for all 15 cells
- `plots/law_VII_decay_heatmap.png`: 3x5 heatmap of DECAY fraction

## 7. How to Run

```bash
cd "ED Research/ED Simulations/Reproduce_This_Law/Law_VII_Tangent_Decay_Exclusivity"
python run_law_VII.py
```

## 8. Failure Modes

- DECAY under any radial gate at any N
- DECAY under tangent gate at N != 4
- DECAY vanishing at (tangent, N=4)

## 9. Related Laws

- [Law I](../Law_I_Binary_Partition/): Law VII strengthens the radial
  DECAY=0 into a table-wide exclusivity statement.
- [Law VI](../Law_VI_N4_Rotational_Degeneracy/): DECAY is the mechanism
  signature of the N=4 tangent anomaly.
