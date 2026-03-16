# Law VI: N=4 Rotational Degeneracy Law

## Formal Statement

At N=4 under the tangent gate, IC vanishes entirely (0%) and DECAY
dominates (77.4%). The square's tangent velocity field produces exact
closed rotation with zero net radial drift, trapping particles in
perpetual orbit. This degeneracy is strictly tangent-specific: the same
N=4 under the inward gate yields 73.5% IC and under the outward gate
yields 26.8% IC, both with zero DECAY.

## Intuitive Explanation

A square has four vertices whose tangent directions form a closed
rotational cycle: each vertex's tangent vector points toward the next
edge, not toward or away from any other vertex. The result is that all
four particles orbit the center in perfect lockstep, never approaching
or separating. This is a geometric degeneracy specific to the square —
the only regular polygon whose tangent directions form exact 90-degree
rotations. No net radial component exists, so no collapse occurs, and
the simulation reaches the DECAY horizon.

## Supporting Evidence

### The Anomalous Cell

| Gate | N | IC % | OPBC % | DECAY % | PBC-corner % | OL % |
|------|---|------|--------|---------|--------------|------|
| inward | 4 | 73.5 | 26.5 | 0.0 | 0.0 | 0.0 |
| **tangent** | **4** | **0.0** | **10.6** | **77.4** | **4.5** | **7.4** |
| outward | 4 | 26.8 | 67.5 | 0.0 | 0.0 | 5.7 |

### Gate Specificity

N=4 under radial gates shows normal IC fractions:
- Inward: 73.5% IC (second highest after N=3's 100%)
- Outward: 26.8% IC (consistent with outward trend)

The anomaly is purely a property of the tangent velocity geometry
interacting with the square's four-fold rotational symmetry.

### DECAY Concentration

DECAY accounts for 240 of 310 points (77.4%) in the N=4 tangent sweep.
The remaining non-DECAY events are:
- outward-PBC: 33 (10.6%)
- other-late: 23 (7.4%)
- PBC-corner: 14 (4.5%)
- inward-collapse: 0 (0.0%)

### Uniqueness of N=4 DECAY

Across all 15 (gate, N) pairs, DECAY > 0 occurs only at (tangent, N=4).
See Law VII for the formal exclusivity statement.

## Reproducibility

### Sweep File

Load: `n4_tangent_sweep.json` (note: contains 310 records in
`{results: [...], intervals: [...]}` format)

### Metrics to Compute

Extract the `results` list. Count mechanism labels.

### Success Criterion

Law VI is reproduced if:
1. Zero records have `mechanism == "inward-collapse"` in the N=4 tangent
   sweep
2. DECAY fraction exceeds 70%
3. N=4 inward and outward sweeps both show IC > 20% and DECAY = 0%

## Failure Modes

### What Would Falsify Law VI

- Any IC event in the N=4 tangent sweep (at finer diameter resolution
  or different merge threshold)
- DECAY dropping below 50% under parameter variation, suggesting the
  closed orbit is not robust
- Another N (e.g., N=8) showing the same zero-IC, high-DECAY pattern
  under the tangent gate

### Scientific Meaning of Failure

IC appearing at N=4 tangent would mean the closed rotational orbit is
not perfectly degenerate — some diameters produce a slight radial drift
that accumulates into collapse. This would downgrade Law VI from an
exact structural zero to an approximate suppression. Another N showing
the same pattern would generalize the degeneracy, suggesting that
regular polygons whose tangent cycle length divides 2*pi may all exhibit
rotational trapping (e.g., N=4 with 90-degree steps, hypothetically
N=8 with 45-degree steps).

## Related Laws

- [Law V](../Law_V_N3_Total_Collapse/): The complementary singular cell
  (100% IC at inward N=3 vs. 0% IC at tangent N=4). Together they define
  the extremes of the mechanism table.
- [Law VII](../Law_VII_Tangent_Decay_Exclusivity/): Formalizes that Law
  VI's DECAY phenomenon is confined to this single cell.
- [Law IV](../Law_IV_Tangent_Fragmentation/): Law VI is the most extreme
  manifestation of tangent-gate anomalous behavior — zero IC bands is
  the limiting case of fragmentation.
