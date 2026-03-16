# Reproduce This Paper — ED‑Arch Laws I–V and Precursor Studies (20 & 21)

This directory contains the **complete reproducibility suite** for the Event Density Architecture (ED‑Arch) canon.  
It provides **one‑command, deterministic reproduction** of every figure, table, mechanism map, and law‑surface used in:

- **ED‑Arch‑20** — Global Law Surfaces  
- **ED‑Arch‑21** — Compositionality & y‑Switching  
- **ED‑Arch‑22** — Law I: Sheets, Mirrors, Ladders  
- **ED‑Arch‑23** — Law II: Islands, Ridges, Resonance  
- **ED‑Arch‑24** — Law III: Windows, Cycles, Recurrence  
- **ED‑Arch‑25** — Law IV: Scaling, Sensitivity, Structure  
- **ED‑Arch‑26** — Law V: Commensurability Classes  

The goal is simple:

> **Anyone should be able to reproduce the entire ED‑Arch canon with zero guesswork.**

This directory provides the structure, wrappers, and expected outputs to make that possible.

---

## 📦 What’s in this directory?

```
Reproduce_This_Paper/
│
├── 20_Manifold_Global_Law_Surfaces/
├── 21_Compositionality_y_Switching/
│
├── Law_I_Sheets_Mirrors_Ladders/
├── Law_II_Islands_Ridges_Resonance/
├── Law_III_Windows_Cycles_Recurrence/
├── Law_IV_Scaling_Sensitivity_Structure/
├── Law_V_Commensurability_Classes/
│
└── README.md   ← (this file)
```

Each law folder contains:

- a **run_law_XX.py** wrapper  
- an **expected_outputs/** directory  
- a **README** describing the law  
- the harness configuration for that law  

All wrappers call the same unified harness:

```
python run_arch_harness.py --suite <suite-name>
```

---

## 🚀 Quickstart — Reproduce Everything

To reproduce **all seven laws and both precursor studies**, run:

```
python run_arch_harness.py --suite all
```

This will:

- generate all law surfaces  
- compute all invariants  
- produce all mechanism maps  
- regenerate all figures and JSON reports  
- compare them to the canonical expected outputs  

If anything diverges, the harness will tell you exactly where and why.

---

## ▶️ Run a Single Law

Each law has a tiny wrapper script:

```
python 20_Manifold_Global_Law_Surfaces/run_law_20.py
python 21_Compositionality_y_Switching/run_law_21.py

python Law_I_Sheets_Mirrors_Ladders/run_law_I.py
python Law_II_Islands_Ridges_Resonance/run_law_II.py
python Law_III_Windows_Cycles_Recurrence/run_law_III.py
python Law_IV_Scaling_Sensitivity_Structure/run_law_IV.py
python Law_V_Commensurability_Classes/run_law_V.py
```

Each wrapper:

- locates the repo root  
- calls the unified harness  
- passes the correct suite name  
- writes outputs into that law’s folder  

---

## 📁 Expected Outputs

Every law folder contains:

```
expected_outputs/
    json/
    figures/
```

These hold the **canonical ground‑truth outputs** used for:

- regression testing  
- correctness verification  
- scientific reproducibility  
- student comparison  
- paper‑figure regeneration  

The harness automatically compares new outputs to these canonical ones.

---

## 🧭 Orientation Map — How the Laws Fit Together

| Component | Purpose | Output |
|----------|---------|--------|
| **ED‑Arch‑20** | Global law surfaces | 3‑sheet manifold, mechanism boundaries |
| **ED‑Arch‑21** | Invariants & y‑switching | 7 invariants, compositionality, robustness |
| **Law I** | Sheets, mirrors, ladders | Geometry of collapse surfaces |
| **Law II** | Islands, ridges, resonance | Resonance‑driven structure |
| **Law III** | Windows, cycles, recurrence | Timing laws & recurrence windows |
| **Law IV** | Scaling laws | N‑scaling of all structures |
| **Law V** | Commensurability classes | Constructive/destructive regimes |

Together, these form the **ED‑Arch canon**:  
a complete architectural description of collapse dynamics.

---

## 🔬 Reproducibility Philosophy

This suite is designed so that:

- every figure in every paper can be regenerated  
- every mechanism map is deterministic  
- every JSON report is stable  
- every law surface is reproducible  
- every scaling law is testable  

The harness enforces:

- fixed seeds  
- fixed sampling grids  
- fixed numerical tolerances  
- fixed figure styles  
- fixed output paths  

This ensures **scientific invariance** across machines, OSes, and Python versions.

---

## 🧱 Folder‑Level READMEs

Each law folder contains its own README describing:

- the purpose of the law  
- the structure of the expected outputs  
- the figures and JSON files generated  
- how to interpret the results  

These READMEs are intentionally short and practical.

---

## 🧩 Extending the Suite

To add a new law:

1. Create a new folder  
2. Add a wrapper script  
3. Add an `expected_outputs/` scaffold  
4. Add a law‑specific config to the harness  
5. Add a short README  

The harness will automatically pick it up.

---

## 🏁 Final Notes

This directory is the **operational core** of ED‑Arch reproducibility.  
It is designed to be:

- minimal  
- deterministic  
- transparent  
- educational  
- architecturally honest  

If you run this suite, you reproduce the entire canon.

