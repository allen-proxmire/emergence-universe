# Reproduce: Law IV — Tangent Fragmentation Law

## 1. Formal Law Statement

Under the tangent gate, the IC regime shatters into 7-10 disjoint bands
for N >= 5. No comparable fragmentation occurs under the inward or
outward gates (1-3 IC bands at all N). PBC-corner events (4-12% for
N >= 4) are exclusively tangent-gate phenomena.

## 2. What This Harness Verifies

- Counts contiguous IC bands for each (gate, N) pair
- Verifies tangent gate has >= 7 bands for N >= 5
- Verifies radial gates have <= 3 bands for all N
- Counts PBC-corner events; verifies zero for radial, nonzero for
  tangent N >= 4
- Measures largest IC band width for tangent N >= 5

## 3. Files Loaded

All 15 sweep files: `n{3,4,5,6,7}_{inward,tangent,outward}_sweep.json`

## 4. Metrics Computed

- `ic_band_count[gate][N]`: number of contiguous IC bands
- `max_band_width[gate][N]`: width of the largest IC band (d_px)
- `pbc_corner_count[gate][N]`: absolute count of PBC-corner events
- `pbc_corner_frac[gate][N]`: PBC-corner fraction

## 5. Pass/Fail Criteria

| Criterion | Threshold | Type |
|-----------|-----------|------|
| Tangent IC bands for N >= 5 | >= 7 | Lower bound |
| Radial IC bands for all N | <= 3 | Upper bound |
| PBC-corner for all radial cells | == 0 | Exact |
| PBC-corner for tangent N >= 4 | > 0 | Existence |

## 6. Expected Outputs

- `expected_outputs/metrics.json`: canonical band counts and PBC-corner
- `plots/law_IV_band_counts.png`: bar chart comparing band counts across gates

## 7. How to Run

```bash
cd "ED Research/ED Simulations/Reproduce_This_Law/Law_IV_Tangent_Fragmentation"
python run_law_IV.py
```

## 8. Failure Modes

- Tangent IC coalescing to <= 3 bands at some N >= 8
- Radial gates developing >= 5 IC bands at large N
- PBC-corner appearing under a radial gate

## 9. Related Laws

- [Law I](../Law_I_Binary_Partition/): Binary partition within which
  tangent fragmentation adds complexity.
- [Law VI](../Law_VI_N4_Rotational_Degeneracy/): N=4 tangent anomaly is
  the extreme case (zero IC bands).
