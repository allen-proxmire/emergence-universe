# ED-Arch Mechanism Laws

## Overview

The ED-Arch Mechanism Laws are seven formal statements that encode the
complete structure of the micro-event mechanism space across the three
canonical velocity gates (inward, tangent, outward) and particle counts
N = 3 through 7. They were derived from exhaustive diameter sweeps
(d\_px = 1..385, totalling over 5,700 simulations) and represent the
first systematic classification of how ring-particle collapse mechanisms
depend on initial velocity geometry and polygon order. The laws are
organized in a hierarchy: universal invariants that hold across all gates
and N, gate-specific behaviors unique to each velocity geometry, and
N-specific anomalies tied to particular polygon symmetries.

## Law Index

| Law | Name | Type | One-Line Summary |
|-----|------|------|------------------|
| [I](Law_I_Binary_Partition/) | Binary Partition | Universal | Radial gates produce only IC and outward-PBC; no DECAY, no PBC-corner. |
| [II](Law_II_Radial_Complementarity/) | Radial Complementarity | Universal | IC(inward) + IC(outward) ≈ 100% for N ≥ 4. |
| [III](Law_III_Monotone_Gate_Ordering/) | Monotone Gate-Ordering | Universal | Inward IC decreases with N; outward IC increases; inward always ranks first. |
| [IV](Law_IV_Tangent_Fragmentation/) | Tangent Fragmentation | Gate-specific | Tangent IC shatters into 7-10 bands for N ≥ 5; bowls and PBC-corner are tangent-only. |
| [V](Law_V_N3_Total_Collapse/) | N=3 Total Collapse | N-specific | N=3 inward achieves 100% IC across all diameters. |
| [VI](Law_VI_N4_Rotational_Degeneracy/) | N=4 Rotational Degeneracy | N-specific | N=4 tangent produces 0% IC and 77% DECAY. |
| [VII](Law_VII_Tangent_Decay_Exclusivity/) | Tangent DECAY Exclusivity | Universal | DECAY occurs if and only if gate=tangent and N=4. |

## Law Hierarchy

```
UNIVERSAL INVARIANTS
├── Law I:   Binary Partition Law
│             (constrains which mechanisms can appear)
├── Law II:  Radial Complementarity Law
│             (constrains numerical IC values across radial gates)
├── Law III: Monotone Gate-Ordering Law
│             (constrains IC trends with N)
└── Law VII: Tangent DECAY Exclusivity Law
              (constrains DECAY to a single cell)

GATE-SPECIFIC LAWS
└── Law IV:  Tangent Fragmentation Law
              (describes tangent-only complexity class)

N-SPECIFIC ANOMALIES
├── Law V:   N=3 Inward Total Collapse Law
│             (singular cell: inward × N=3)
└── Law VI:  N=4 Rotational Degeneracy Law
              (singular cell: tangent × N=4)
```

## Relationship Diagram

The seven laws collectively encode the full Gate x N mechanism table:

```
              INWARD         TANGENT        OUTWARD
         ┌──────────────┬──────────────┬──────────────┐
   N=3   │ 100% IC [V]  │  43% IC      │  26% IC      │  Law I: only IC + OPBC
         ├──────────────┼──────────────┼──────────────┤    in radial columns
   N=4   │  74% IC      │  0% IC [VI]  │  27% IC      │
         ├──────────────┼──────────────┼──────────────┤  Law II: each row's
   N=5   │  68% IC      │  37% IC [IV] │  36% IC      │    inward + outward
         ├──────────────┼──────────────┼──────────────┤    ≈ 100%
   N=6   │  60% IC      │  37% IC [IV] │  40% IC      │
         ├──────────────┼──────────────┼──────────────┤  Law III: inward col
   N=7   │  60% IC      │  21% IC [IV] │  45% IC      │    decreases ↓
         └──────────────┴──────────────┴──────────────┘    outward col
              ↓ decreasing    fragmented    ↑ increasing     increases ↑
            [Law III]         [Law IV]      [Law III]

         Law VII: DECAY appears only at (tangent, N=4)
```

## Data Foundation

All laws are derived from the following sweep files located in
`ED Research/ED Simulations/`:

| File | Gate | N | Records |
|------|------|---|---------|
| `n3_inward_sweep.json` | inward | 3 | 385 |
| `n3_tangent_sweep.json` | tangent | 3 | 385 |
| `n3_outward_sweep.json` | outward | 3 | 385 |
| `n4_inward_sweep.json` | inward | 4 | 385 |
| `n4_tangent_sweep.json` | tangent | 4 | 310 |
| `n4_outward_sweep.json` | outward | 4 | 385 |
| `n5_inward_sweep.json` | inward | 5 | 385 |
| `n5_tangent_sweep.json` | tangent | 5 | 385 |
| `n5_outward_sweep.json` | outward | 5 | 385 |
| `n6_inward_sweep.json` | inward | 6 | 385 |
| `n6_tangent_sweep.json` | tangent | 6 | 385 |
| `n6_outward_sweep.json` | outward | 6 | 385 |
| `n7_inward_sweep.json` | inward | 7 | 385 |
| `n7_tangent_sweep.json` | tangent | 7 | 385 |
| `n7_outward_sweep.json` | outward | 7 | 385 |

## Exact vs. Empirical Status

- **Exact (structural zeros):** Laws I, V, VI, VII. These describe
  absolute presences or absences (100% IC, 0% IC, 0% DECAY) that hold
  without exception across every simulation point. Falsifiable by a
  single counterexample at untested N or finer diameter resolution.

- **Empirical (quantitative trends):** Laws II, III, IV. These describe
  numerical trends (complementarity sums, monotonicity, fragmentation
  thresholds) established over N = 3-7. Extension to N = 8-12 would
  strengthen their status.
