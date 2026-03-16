# Reproduce: Law VI — N=4 Rotational Degeneracy Law

## 1. Formal Law Statement

At N=4 under the tangent gate, IC vanishes entirely (0%) and DECAY
dominates (77.4%). The square's tangent velocity field produces exact
closed rotation with zero net radial drift. This degeneracy is strictly
tangent-specific: N=4 inward yields 73.5% IC and N=4 outward yields
26.8% IC, both with zero DECAY.

## 2. What This Harness Verifies

- Loads the N=4 tangent sweep and confirms zero IC events
- Confirms DECAY fraction exceeds 70%
- Verifies N=4 inward and outward show normal IC (> 20%) and zero DECAY
- Reports the full mechanism distribution for the anomalous cell

## 3. Files Loaded

| File | Gate | N |
|------|------|---|
| `n4_tangent_sweep.json` | tangent | 4 |
| `n4_inward_sweep.json` | inward | 4 |
| `n4_outward_sweep.json` | outward | 4 |

## 4. Metrics Computed

- `n4_tangent_ic_count`: IC count (expected 0)
- `n4_tangent_decay_frac`: DECAY fraction (expected ~77%)
- `n4_tangent_mechanism_dist`: full mechanism distribution
- `n4_inward_ic_frac`: IC fraction under inward gate
- `n4_outward_ic_frac`: IC fraction under outward gate
- `n4_inward_decay_count`: DECAY count under inward (expected 0)
- `n4_outward_decay_count`: DECAY count under outward (expected 0)

## 5. Pass/Fail Criteria

| Criterion | Threshold | Type |
|-----------|-----------|------|
| N=4 tangent IC count | == 0 | Exact |
| N=4 tangent DECAY fraction | >= 0.70 | Lower bound |
| N=4 inward IC fraction | >= 0.20 | Lower bound |
| N=4 outward IC fraction | >= 0.20 | Lower bound |
| N=4 inward DECAY count | == 0 | Exact |
| N=4 outward DECAY count | == 0 | Exact |

## 6. Expected Outputs

- `expected_outputs/metrics.json`: canonical values for the anomalous cell
- `plots/law_VI_n4_anomaly.png`: pie chart of N=4 tangent mechanism dist
  alongside N=4 inward/outward for contrast

## 7. How to Run

```bash
cd "ED Research/ED Simulations/Reproduce_This_Law/Law_VI_N4_Rotational_Degeneracy"
python run_law_VI.py
```

## 8. Failure Modes

- Any IC event in the N=4 tangent sweep
- DECAY dropping below 50%
- Another N showing zero-IC/high-DECAY under tangent gate
- DECAY appearing under N=4 radial gates

## 9. Related Laws

- [Law V](../Law_V_N3_Total_Collapse/): Opposite extreme (100% IC at
  inward N=3 vs. 0% IC at tangent N=4).
- [Law VII](../Law_VII_Tangent_Decay_Exclusivity/): DECAY confined to
  this single cell.
- [Law IV](../Law_IV_Tangent_Fragmentation/): Zero IC bands is the
  extreme case of tangent fragmentation.
