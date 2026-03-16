# Law VII: Tangent DECAY Exclusivity Law

## Formal Statement

DECAY events occur if and only if the gate is tangent and N=4. Across all
15 (gate, N) pairs tested — three gates times five particle counts,
totalling over 5,700 simulations — DECAY appears in exactly one cell:
(tangent, N=4) at 77.4%. The 14 remaining cells show identically zero
DECAY.

## Intuitive Explanation

DECAY means the simulation reached its time horizon without any pair of
particles approaching within the merge threshold. This requires that all
pairwise distances remain bounded away from zero for the entire evolution.
Under radial gates, particles are always either approaching each other
(inward) or wrapping through PBC to approach from the far side (outward),
so collapse is inevitable. Only the tangent gate at N=4 produces a
velocity field where particles orbit without any net approach — the
square's perfect rotational symmetry. At all other N under the tangent
gate, slight asymmetries in the tangent-PBC interaction break the
rotational degeneracy and restore collapse.

## Supporting Evidence

### Complete DECAY Census

| Gate | N=3 | N=4 | N=5 | N=6 | N=7 |
|------|-----|-----|-----|-----|-----|
| inward | 0% | 0% | 0% | 0% | 0% |
| tangent | 0% | **77.4%** | 0% | 0% | 0% |
| outward | 0% | 0% | 0% | 0% | 0% |

### Point Counts

Total simulation points with DECAY = 0:
- All inward sweeps: 1,925 points, 0 DECAY
- Tangent N=3: 385 points, 0 DECAY
- Tangent N=5: 385 points, 0 DECAY
- Tangent N=6: 385 points, 0 DECAY
- Tangent N=7: 385 points, 0 DECAY
- All outward sweeps: 1,925 points, 0 DECAY

Total zero-DECAY points: 5,390

The sole nonzero cell:
- Tangent N=4: 310 points, 240 DECAY (77.4%)

### If-and-Only-If Structure

This law has a biconditional form:
- **If** gate=tangent AND N=4, **then** DECAY > 0 (specifically 77.4%)
- **If** DECAY > 0, **then** gate=tangent AND N=4

Both directions are verified by the data.

## Reproducibility

### Sweep Files

Load all 15 sweep files:
- `n{3,4,5,6,7}_{inward,tangent,outward}_sweep.json`

### Metrics to Compute

For each file, count records with `mechanism == "DECAY"`.

### Success Criterion

Law VII is reproduced if and only if:
1. Exactly one (gate, N) pair has nonzero DECAY count
2. That pair is (tangent, N=4)
3. DECAY fraction at (tangent, N=4) exceeds 50%

## Failure Modes

### What Would Falsify Law VII

- Any DECAY event under the inward or outward gate at any N
- Any DECAY event under the tangent gate at N != 4
- DECAY vanishing at (tangent, N=4) under parameter variation

### Scientific Meaning of Failure

DECAY under a radial gate would falsify both Law VII and Law I
simultaneously, implying that PBC geometry alone can produce non-collapsing
trapped orbits. DECAY at tangent N != 4 would generalize the rotational
degeneracy beyond the square, suggesting a family of DECAY-supporting
polygon orders (possibly N=4k for integer k). DECAY vanishing at tangent
N=4 would downgrade Law VI's rotational degeneracy to an approximate
resonance rather than an exact symmetry.

## Related Laws

- [Law I](../Law_I_Binary_Partition/): Law VII strengthens Law I's
  "no DECAY under radial gates" by localizing DECAY to a single cell
  in the entire table.
- [Law VI](../Law_VI_N4_Rotational_Degeneracy/): Law VII is the DECAY
  component of the N=4 tangent anomaly; Law VI describes the full
  mechanism distribution including the 0% IC.
- [Law IV](../Law_IV_Tangent_Fragmentation/): DECAY exclusivity is part
  of the broader pattern where exotic mechanisms (DECAY, PBC-corner,
  bowls) are confined to the tangent gate.
