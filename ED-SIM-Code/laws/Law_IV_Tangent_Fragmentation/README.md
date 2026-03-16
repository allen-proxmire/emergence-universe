# Law IV: Tangent Fragmentation Law

## Formal Statement

Under the tangent gate, the IC regime is contiguous for N <= 4 but
shatters into 7-10 disjoint bands for N >= 5, with the largest band
shrinking from 166 px (N=3) to 23 px (N=7). No comparable fragmentation
occurs under the inward or outward gates, which maintain 1-3 contiguous
IC regions at all N. Bowl-shaped chi sub-plateaus within IC bands and
PBC-corner events (4-12% for N >= 4) are exclusively tangent-gate
phenomena.

## Intuitive Explanation

The tangent gate injects rotational momentum that interacts with the PBC
lattice in complex, N-dependent ways. At small N (3-4), the rotational
geometry is simple enough to produce a single contiguous IC band. At
N >= 5, the increasing number of tangent velocity directions creates
interference patterns with the PBC boundaries, producing alternating
zones of collapse and escape. The radial gates avoid this complexity
because their velocities project purely along the radial direction,
bypassing the rotational-PBC coupling entirely.

## Supporting Evidence

### IC Band Count by Gate and N

| N | Inward | Tangent | Outward |
|---|--------|---------|---------|
| 3 | 1 | 1 | 1 |
| 4 | 1 | 0 | 1 |
| 5 | 3 | **7** | 2 |
| 6 | 1 | **7** | 1 |
| 7 | 3 | **10** | 2 |

### Largest IC Band Width (d_px)

| N | Inward | Tangent | Outward |
|---|--------|---------|---------|
| 3 | 385 | 166 | 99 |
| 4 | 283 | 0 | 103 |
| 5 | 235 | 70 | 138 |
| 6 | 231 | 69 | 155 |
| 7 | 173 | **23** | 171 |

### Tangent-Only Complexity Features

**PBC-corner events:**
- N=4: 4.5%
- N=5: 11.2%
- N=6: 11.9%
- N=7: 7.3%
- Inward/outward: 0% at all N

**Bowl-shaped sub-plateaus:**
- Detected at N=3, N=6, N=7 under tangent gate
- Never detected under inward or outward gates

### Mechanism Transitions (total band count - 1)

| N | Inward | Tangent | Outward |
|---|--------|---------|---------|
| 3 | 0 | 11 | 3 |
| 4 | 1 | 5 | 2 |
| 5 | 5 | 18 | 4 |
| 6 | 1 | 19 | 2 |
| 7 | 6 | **24** | 5 |

The tangent gate's transition count grows roughly linearly with N,
reaching 24 at N=7.

## Reproducibility

### Sweep Files

Load all tangent sweep files:
- `n{3,4,5,6,7}_tangent_sweep.json`

For comparison, also load radial sweeps for the same N values.

### Metrics to Compute

1. Count contiguous IC bands (runs of consecutive d\_px with
   mechanism = "inward-collapse")
2. Measure width of largest IC band
3. Count PBC-corner events
4. Check for bowl-shaped sub-plateaus (chi descends then ascends within
   a single plateau of width >= 10)

### Success Criterion

Law IV is reproduced if:
1. Tangent gate shows >= 7 IC bands for all N >= 5
2. Radial gates show <= 3 IC bands for all N
3. PBC-corner count is zero for all radial sweeps and positive for
   tangent N >= 4

## Failure Modes

### What Would Falsify Law IV

- Tangent IC coalescing back into <= 3 bands at some N >= 8
  (defragmentation)
- Radial gates developing >= 5 IC bands at large N (radial fragmentation)
- PBC-corner events appearing under a radial gate
- Bowl sub-plateaus appearing under a radial gate

### Scientific Meaning of Failure

Tangent defragmentation at large N would imply that the rotational-PBC
interference pattern simplifies at high polygon order, possibly due to
the polygon approaching a smooth circle where tangent directions become
indistinguishable. Radial fragmentation would imply that PBC geometry
alone can create interference effects without rotational input, which
would challenge the interpretation that fragmentation requires
tangent-PBC coupling.

## Related Laws

- [Law I](../Law_I_Binary_Partition/): Establishes the binary mechanism
  space within which tangent fragmentation adds complexity.
- [Law VI](../Law_VI_N4_Rotational_Degeneracy/): The N=4 tangent anomaly
  (zero IC bands) is the extreme case of tangent-specific behavior.
- [Law VII](../Law_VII_Tangent_Decay_Exclusivity/): DECAY is another
  tangent-only phenomenon, concentrated at the same N=4 that shows zero
  IC bands.
