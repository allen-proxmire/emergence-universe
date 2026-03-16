# Law I: Binary Partition Law

## Formal Statement

Under the two radial gates (inward and outward), the mechanism space
partitions into exactly two dynamical outcomes — inward-collapse (IC)
and outward-PBC — for all N >= 3. No DECAY and no PBC-corner events
occur. A thin other-late strip (at most 6% of the diameter range)
appears at sub-threshold diameters but carries no dynamical content.

## Intuitive Explanation

When particles move purely radially (either toward or away from the ring
center), they can only do two things: collapse onto each other via direct
approach or via PBC wrap-around. The rotational complexity needed to
produce DECAY (perpetual orbiting) or PBC-corner events (diagonal
wrap-around collapse) is geometrically inaccessible to radial initial
conditions. The tangent gate, by contrast, injects rotational momentum
that opens these additional mechanism channels.

## Supporting Evidence

### Gates and N Values

Law I is verified across all 10 radial (gate, N) pairs:

| Gate | N | IC % | OPBC % | DECAY % | PBC-corner % | OL % |
|------|---|------|--------|---------|--------------|------|
| inward | 3 | 100.0 | 0.0 | **0.0** | **0.0** | 0.0 |
| inward | 4 | 73.5 | 26.5 | **0.0** | **0.0** | 0.0 |
| inward | 5 | 67.5 | 32.5 | **0.0** | **0.0** | 0.0 |
| inward | 6 | 60.0 | 40.0 | **0.0** | **0.0** | 0.0 |
| inward | 7 | 59.5 | 36.9 | **0.0** | **0.0** | 3.6 |
| outward | 3 | 25.7 | 68.8 | **0.0** | **0.0** | 5.5 |
| outward | 4 | 26.8 | 67.5 | **0.0** | **0.0** | 5.7 |
| outward | 5 | 36.4 | 57.9 | **0.0** | **0.0** | 5.7 |
| outward | 6 | 40.3 | 54.0 | **0.0** | **0.0** | 5.7 |
| outward | 7 | 44.7 | 49.6 | **0.0** | **0.0** | 5.7 |

### Structural Zeros

- DECAY = 0 across all 10 radial cells (3,850 simulation points)
- PBC-corner = 0 across all 10 radial cells

### Contrast with Tangent Gate

The tangent gate violates binary partition:
- N=4 tangent: 77.4% DECAY, 4.5% PBC-corner
- N=5 tangent: 11.2% PBC-corner
- N=6 tangent: 11.9% PBC-corner
- N=7 tangent: 7.3% PBC-corner

## Reproducibility

### Sweep Files

Load any of the 10 radial sweep files:
- `n{3,4,5,6,7}_inward_sweep.json`
- `n{3,4,5,6,7}_outward_sweep.json`

### Metrics to Compute

For each file, count occurrences of each mechanism label.

### Success Criterion

Law I is reproduced if and only if:
1. Every record has `mechanism` in `{inward-collapse, outward-PBC, other-late}`
2. Zero records have `mechanism == "DECAY"`
3. Zero records have `mechanism == "PBC-corner"`

## Failure Modes

### What Would Falsify Law I

A single DECAY or PBC-corner event under the inward or outward gate at
any N would falsify the law. Potential sources:

- **N >= 8:** Higher polygon orders might create near-tangent velocity
  components from PBC-reflected trajectories, opening rotational channels.
- **Finer diameter resolution:** Sub-pixel diameter steps could reveal
  narrow mechanism pockets invisible at d\_px step = 1.
- **Different merge threshold:** Changing MERGE\_THR\_PX could alter the
  collapse detection boundary, potentially reclassifying edge cases.

### Scientific Meaning of Failure

If DECAY appeared under a radial gate, it would imply that PBC geometry
alone (without rotational initial conditions) can trap particles in
non-collapsing orbits — a fundamentally different dynamical regime than
the current understanding where DECAY requires tangent momentum.

## Related Laws

- [Law VII](../Law_VII_Tangent_Decay_Exclusivity/): Strengthens Law I by
  localizing DECAY to the single cell (tangent, N=4).
- [Law II](../Law_II_Radial_Complementarity/): Quantifies the IC/OPBC
  partition within the binary mechanism space established by Law I.
- [Law IV](../Law_IV_Tangent_Fragmentation/): Describes the richer
  mechanism space that arises when Law I's binary constraint is relaxed
  under the tangent gate.
