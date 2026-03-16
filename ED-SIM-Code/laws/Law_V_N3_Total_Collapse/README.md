# Law V: N=3 Inward Total Collapse Law

## Formal Statement

At N=3 under the inward gate, every diameter produces inward-collapse:
IC = 100% across the full 385-point diameter sweep. This is the unique
(gate, N) pair achieving total collapse. No other N achieves total
collapse under any gate, and no other gate achieves it at N=3 (tangent:
43.1%, outward: 25.7%).

## Intuitive Explanation

The equilateral triangle's three inward-pointing velocities converge to
the centroid by construction. At small diameters, the particles are
already close and collapse immediately. At large diameters, PBC wrapping
shifts the positions but the three-fold symmetric inward velocity field
still drives convergence. Unlike higher-N polygons, where inward
velocities can create competing approach channels that some particles
escape, the triangle has too few particles for any to avoid the
converging flow. The result is total collapse at every ring size.

## Supporting Evidence

### The Singular Cell

| Gate | N | IC % | OPBC % | DECAY % | Total Points |
|------|---|------|--------|---------|--------------|
| **inward** | **3** | **100.0** | **0.0** | **0.0** | **385** |
| tangent | 3 | 43.1 | 26.8 | 0.0 | 385 |
| outward | 3 | 25.7 | 68.8 | 0.0 | 385 |

### Uniqueness Verification

IC fraction for all 15 (gate, N) pairs:

| N | Inward | Tangent | Outward |
|---|--------|---------|---------|
| 3 | **100.0%** | 43.1% | 25.7% |
| 4 | 73.5% | 0.0% | 26.8% |
| 5 | 67.5% | 36.6% | 36.4% |
| 6 | 60.0% | 36.9% | 40.3% |
| 7 | 59.5% | 20.8% | 44.7% |

No other cell reaches 100%. The closest is inward N=4 at 73.5%.

### Chi Profile at N=3 Inward

The chi values span the full range [0.010, 2.050], with 21 distinct
chi levels and 103 plateaus. The total-collapse regime is not trivial
uniform collapse — it exhibits rich internal chi structure while
maintaining 100% IC mechanism assignment.

## Reproducibility

### Sweep File

Load: `n3_inward_sweep.json`

### Metrics to Compute

Count mechanism labels. Compute IC fraction.

### Success Criterion

Law V is reproduced if and only if:
1. Every record in `n3_inward_sweep.json` has `mechanism == "inward-collapse"`
2. IC fraction = 100.0% exactly (385/385)

## Failure Modes

### What Would Falsify Law V

- A single non-IC event in the N=3 inward sweep at any diameter
- Achieving 100% IC at a different (gate, N) pair, which would remove
  the uniqueness claim

### Scientific Meaning of Failure

If a non-IC event appeared at N=3 inward, it would imply that PBC
geometry can deflect the triangular convergence at some critical
diameter — revealing a previously unknown resonance between the
triangle's symmetry and the square PBC domain. If another (gate, N)
pair achieved 100% IC, it would generalize the total-collapse phenomenon
beyond the triangle, suggesting a broader class of geometries with
inescapable convergence.

## Related Laws

- [Law III](../Law_III_Monotone_Gate_Ordering/): Law V is the N=3
  endpoint of the inward gate's monotone IC decrease (100% is the
  maximum).
- [Law II](../Law_II_Radial_Complementarity/): Law V explains the N=3
  complementarity outlier (125.7%) — the inward gate's total collapse
  inflates the sum above 100%.
- [Law VI](../Law_VI_N4_Rotational_Degeneracy/): The other singular cell
  in the mechanism table, at the opposite extreme (0% IC vs. 100% IC).
