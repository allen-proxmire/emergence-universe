ED-ARCH CANON — LAW PIPELINE
═══════════════════════════════════════════════════════════════════════════════════

  ┌─────────────────────────────────┐
  │         ED-ARCH-20              │
  │      GLOBAL LAW SURFACES        │
  │                                 │
  │  Three γ-sheets over (d, N)     │
  │  Analytic χ and mechanism maps  │
  └───────────────┬─────────────────┘
                  │
                  ▼
  ┌─────────────────────────────────┐
  │         ED-ARCH-21              │
  │   COMPOSITIONALITY & γ-SWITCHING│
  │                                 │
  │  Independent clusters;          │
  │  memoryless sheet transitions   │
  └───────────────┬─────────────────┘
                  │
    ┌─────────────┼─────────────┐
    │             │             │
    ▼             ▼             ▼
  ┌───────┐   ┌────────┐   ┌────────┐
  │ γ = +1│   │ γ =  0 │   │ γ = −1 │
  │INWARD │   │TANGENT │   │OUTWARD │
  └───┬───┘   └───┬────┘   └───┬────┘
      │             │             │
      └─────────────┼─────────────┘
                    │
                    ▼
  ┌─────────────────────────────────┐
  │            LAW I                │
  │   SHEETS · MIRRORS · LADDERS    │
  │                                 │
  │  Collapse geometry of each      │
  │  γ-sheet; bilateral structure   │
  └───────────────┬─────────────────┘
                  │
                  ▼
  ┌─────────────────────────────────┐
  │            LAW II               │
  │   ISLANDS · RIDGES · RESONANCE  │
  │                                 │
  │  Resonant timing; ladder        │
  │  interruptions; DECAY islands   │
  └───────────────┬─────────────────┘
                  │
                  ▼
  ┌─────────────────────────────────┐
  │            LAW III              │
  │ WINDOWS · CYCLES · RECURRENCE   │
  │                                 │
  │  PBC quantization; periodic     │
  │  DECAY windows; cycle timing    │
  └───────────────┬─────────────────┘
                  │
                  ▼
  ┌─────────────────────────────────┐
  │            LAW IV               │
  │ SCALING · SENSITIVITY · STRUCTURE│
  │                                 │
  │  N-scaling of α; window width;  │
  │  K = 2 sin(π/N) closing rate    │
  └───────────────┬─────────────────┘
                  │
                  ▼
  ┌─────────────────────────────────┐
  │            LAW V                │
  │    COMMENSURABILITY CLASSES     │
  │                                 │
  │  Ring-to-box ratio partitions   │
  │  constructive/destructive bands │
  └─────────────────────────────────┘


VALIDATION HARNESS
───────────────────────────────────────────────────────────────────────────────

  run_arch_harness.py
      │
      ├── --suite law-surfaces  →  arch20_law_surfaces.py
      │       • sample_law_surface()
      │       • closing_rate()
      │       • decay_window()
      │       • engine χ comparison via micro_event_operator
      │
      ├── --suite invariants     →  arch21_invariants.py
      │       • INV-21-1  Compositionality
      │       • INV-21-2  Memoryless switching
      │       • INV-21-3  Pythagorean drift
      │       • INV-21-4  Programmed collapse
      │       • INV-21-5  Perturbation hardness
      │       • INV-21-6  DECAY subregime
      │       • INV-21-7  Ghost compositionality
      │
      └── outputs/
              arch_harness_report.txt
              law_surface_results.json
              invariant_results.json
              law_surface_engine_comparison.png
              invariant_summary.png
