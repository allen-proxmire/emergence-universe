# Law II: Radial Complementarity Law

## Formal Statement

For N >= 4, the IC fractions under the inward and outward gates sum to a
near-constant: IC(inward) + IC(outward) = 100 +/- 4%. The measured sums
are 100.3% (N=4), 103.9% (N=5), 100.3% (N=6), 104.2% (N=7). The N=3
outlier (125.7%) arises from the inward gate's total-collapse singularity
(Law V).

## Intuitive Explanation

At any given diameter, a particle pair that collapses via direct radial
approach under the inward gate will instead collapse via PBC wrap-around
under the outward gate, and vice versa. The two radial gates probe the
same underlying distance-closing process from opposite initial conditions.
Their IC fractions are therefore approximately complementary: what one
gate counts as IC, the other counts as outward-PBC, and the sum is
conserved near 100%.

## Supporting Evidence

### Complementarity Sums

| N | IC(inward) | IC(outward) | Sum | Deviation from 100% |
|---|-----------|------------|-----|---------------------|
| 3 | 100.0% | 25.7% | 125.7% | +25.7% (Law V outlier) |
| 4 | 73.5% | 26.8% | 100.3% | +0.3% |
| 5 | 67.5% | 36.4% | 103.9% | +3.9% |
| 6 | 60.0% | 40.3% | 100.3% | +0.3% |
| 7 | 59.5% | 44.7% | 104.2% | +4.2% |

### Statistics (N >= 4)

- Mean sum: 102.2%
- Standard deviation: 2.1%
- Maximum deviation: 4.2% (N=7)

### Geometric Interpretation

The complementarity reflects PBC mirror symmetry. A particle at position
x moving with velocity +v reaches the same PBC-wrapped configuration as
a particle at x moving with -v, but via the complementary route (direct
vs. wrap-around). The ~2-4% excess above 100% represents configurations
where both inward and outward initial conditions lead to IC via different
particle pairs.

## Reproducibility

### Sweep Files

Load paired files for each N:
- `n{4,5,6,7}_inward_sweep.json`
- `n{4,5,6,7}_outward_sweep.json`

### Metrics to Compute

For each file, compute IC fraction = (count of `inward-collapse`) / total.
Sum IC(inward) + IC(outward) for each N.

### Success Criterion

Law II is reproduced if:
1. For all N >= 4: |IC(inward) + IC(outward) - 100%| < 5%
2. The mean sum across N=4-7 falls in [98%, 106%]

## Failure Modes

### What Would Falsify Law II

- A sum deviating more than 10% from 100% at any N >= 4
- A systematic drift in the sum with increasing N (e.g., sums growing
  toward 120% or shrinking toward 80%)
- A second total-collapse singularity at N > 3 inflating the sum

### Scientific Meaning of Failure

If the complementarity sum diverged significantly from 100%, it would
imply that the inward and outward gates access qualitatively different
regions of the collapse landscape — not merely opposite views of the
same PBC geometry. This would break the mirror-symmetry interpretation
and suggest that radial direction carries dynamical content beyond
simple sign inversion.

## Related Laws

- [Law I](../Law_I_Binary_Partition/): Establishes the binary IC/OPBC
  partition within which complementarity operates.
- [Law III](../Law_III_Monotone_Gate_Ordering/): Specifies the individual
  monotone trends whose sum Law II constrains.
- [Law V](../Law_V_N3_Total_Collapse/): Explains the N=3 outlier where
  inward total collapse inflates the sum to 126%.
