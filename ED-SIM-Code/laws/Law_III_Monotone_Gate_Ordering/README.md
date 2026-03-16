# Law III: Monotone Gate-Ordering Law

## Formal Statement

The inward gate is IC-dominant and its IC fraction decreases monotonically
with N; the outward gate is OPBC-dominant and its IC fraction increases
monotonically with N. Specifically: IC(inward) = 100%, 74%, 68%, 60%, 59%
and IC(outward) = 26%, 27%, 36%, 40%, 45% for N = 3, 4, 5, 6, 7. The
inward gate always ranks first in IC fraction at every N tested. The two
radial trends converge: their spread narrows from 74 percentage points at
N=3 to 15 at N=7, consistent with an asymptotic meeting point near 50-55%.

## Intuitive Explanation

Under the inward gate, all velocities point toward the center, maximizing
direct collapse probability. As N grows, the inward momentum distributes
across more directions, reducing the probability of simultaneous pairwise
convergence — hence IC decreases. Under the outward gate, particles escape
radially and must wrap through PBC boundaries to encounter neighbors; more
particles create more pairwise wrap-around approach channels, so IC
increases with N. The inward gate always dominates because direct approach
is geometrically more efficient than wrap-around approach at every N.

## Supporting Evidence

### IC Fraction by Gate and N

| N | IC(inward) | IC(outward) | Spread | Ranking |
|---|-----------|------------|--------|---------|
| 3 | 100.0% | 25.7% | 74.3 pp | inward > tangent > outward |
| 4 | 73.5% | 26.8% | 46.7 pp | inward > outward > tangent |
| 5 | 67.5% | 36.4% | 31.2 pp | inward > tangent ~ outward |
| 6 | 60.0% | 40.3% | 19.7 pp | inward > outward > tangent |
| 7 | 59.5% | 44.7% | 14.8 pp | inward > outward > tangent |

### Monotonicity Verification

Inward IC step differences (N to N+1):
- N=3 to 4: -26.5 pp (decrease)
- N=4 to 5: -6.0 pp (decrease)
- N=5 to 6: -7.5 pp (decrease)
- N=6 to 7: -0.5 pp (decrease)

All steps negative. **Strictly monotonically decreasing: confirmed.**

Outward IC step differences (N to N+1):
- N=3 to 4: +1.0 pp (increase)
- N=4 to 5: +9.6 pp (increase)
- N=5 to 6: +3.9 pp (increase)
- N=6 to 7: +4.4 pp (increase)

All steps positive. **Strictly monotonically increasing: confirmed.**

### Convergence Trend

The inward-outward spread narrows consistently:
74.3 -> 46.7 -> 31.2 -> 19.7 -> 14.8 pp

Linear extrapolation suggests convergence near N = 10-12.

## Reproducibility

### Sweep Files

Load all 10 radial sweep files:
- `n{3,4,5,6,7}_inward_sweep.json`
- `n{3,4,5,6,7}_outward_sweep.json`

### Metrics to Compute

For each file, compute IC fraction. Verify:
1. Inward IC sequence is strictly decreasing
2. Outward IC sequence is strictly increasing
3. IC(inward) > IC(outward) at every N

### Success Criterion

Law III is reproduced if all three conditions hold.

## Failure Modes

### What Would Falsify Law III

- A reversal in either monotone trend at any N (e.g., IC(inward) at N=8
  exceeding IC(inward) at N=7)
- Outward IC exceeding inward IC at some N (crossing of the two curves)
- The convergence trend reversing (spread widening at large N)

### Scientific Meaning of Failure

A monotonicity reversal would imply a non-trivial resonance at specific N
where polygon geometry creates an anomalous enhancement or suppression of
collapse probability — analogous to the N=4 tangent anomaly (Law VI) but
under radial conditions. A crossing point where outward IC exceeds inward
IC would overturn the intuition that direct approach is always more
efficient than wrap-around approach.

## Related Laws

- [Law II](../Law_II_Radial_Complementarity/): Constrains the sum of the
  two monotone trends to ~100%.
- [Law V](../Law_V_N3_Total_Collapse/): The N=3 endpoint of the inward
  monotone decrease (100% IC = total collapse).
- [Law I](../Law_I_Binary_Partition/): Guarantees that IC and OPBC are the
  only two mechanisms whose fractions Law III tracks.
