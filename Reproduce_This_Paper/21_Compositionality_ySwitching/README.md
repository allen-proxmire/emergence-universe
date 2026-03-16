# 📘 ED‑Arch‑21 — Compositionality, y‑Switching, and Architectural Robustness  
### *Reproducibility Folder*

This folder reproduces the results of:

**ED‑Arch‑21 — Compositionality, y‑Switching, and Architectural Robustness**  
(see: `../../docs/ED‑Arch‑21_NOFIGS_Compositionality.pdf`)

ED‑Arch‑21 extends the ED‑Arch‑20 law surfaces into the **multi‑cluster**, **time‑varying**, and **perturbed** regimes.  
It identifies **seven architectural invariants** that show the ED manifold is compositional, memoryless, geometry‑aware, and robust.

---

# 1. 🔍 What ED‑Arch‑21 Establishes

### **Seven Architectural Invariants**

#### **INV‑21‑1 — 3×3 Sheet Compositionality**
Clusters evolve independently on their own y‑sheets, even under deep physical overlap.

#### **INV‑21‑2 — Memoryless Law‑Surface Jumps**
A y→y′ switch produces a trajectory identical to a fresh start at the same geometric state.

#### **INV‑21‑3 — Pythagorean Tangential Drift**
Under y = 0, the effective diameter grows as:


\[
d_{\text{sw}} = \sqrt{d_{\text{init}}^2 + (t_{\text{sw}}\cdot DT)^2}
\]



#### **INV‑21‑4 — Geometry‑Aware Programmed Collapse**
Multi‑phase y‑programs (+1→0→+1, etc.) produce exact collapse times when geometry is accounted for.

#### **INV‑21‑5 — Perturbation Hardness Hierarchy**
Robustness ordering:


\[
y = -1 \;\ge\; y = +1 \;>\; y = 0
\]



#### **INV‑21‑6 — DECAY Angular Sub‑Regime**
Inside the DECAY window, angular perturbations can trigger late PBC‑corner collapse.

#### **INV‑21‑7 — Ghost Compositionality**
Multi‑cluster y‑programming remains fully compositional even under deep pass‑through.

---

# 2. 🧪 What This Folder Tests

Running the harness reproduces:

- all seven invariants  
- all 21A sheet‑interaction tests  
- all 21B y‑switching tests  
- all 21C perturbation‑robustness tests  
- all 21D multi‑cluster y‑programming tests  

Expected outcomes:

- **all invariants pass**  
- deviations match canonical values (0.0 or the documented tiny tolerances)  
- mechanisms match the ED‑Arch‑21 paper exactly  

---

# 3. ▶️ How to Run the Reproduction

From the **repo root**:

```bash
python run_arch_harness.py --suite invariants
```

This command:

- runs all ED‑Arch‑21 invariant tests  
- validates compositionality  
- validates y‑switching  
- validates geometry‑aware programming  
- validates perturbation robustness  
- validates multi‑cluster independence  
- writes JSON + TXT reports  
- generates canonical figures  

Outputs appear in:

```
arch_harness_output/
```

If you prefer a local wrapper:

```bash
python run_law_21.py
```

---

# 4. 📂 Expected Outputs

### **JSON Reports**
- `invariant_results.json`  
- `arch_harness_report.txt`  

### **Figures**
- compositionality matrix  
- y‑switching residuals  
- programmed‑collapse diagrams  
- perturbation‑robustness heatmaps  
- DECAY sub‑regime scatterplots  
- ghost‑compositionality traces  

Canonical copies live in:

```
expected_outputs/
    json/
    figures/
```

---

# 5. 🧭 Mapping to the Paper

This reproducibility folder corresponds to the following sections of ED‑Arch‑21:

| Paper Section | Reproduced Here |
|---------------|-----------------|
| 21A — Sheet Compositionality | ✔ INV‑21‑1 |
| 21B — y‑Switching | ✔ INV‑21‑2, INV‑21‑3, INV‑21‑4 |
| 21C — Perturbation Robustness | ✔ INV‑21‑5, INV‑21‑6 |
| 21D — Multi‑Cluster y‑Programming | ✔ INV‑21‑7 |
| Architectural Synthesis | ✔ all invariants validated |

This folder is the **computational companion** to the paper.

---

# 6. 📌 Notes for Researchers

- All results are deterministic.  
- No random seeds are used.  
- Deviations should match the canonical JSON exactly.  
- Ghost compositionality is exact even under deep physical overlap.  
- Angular perturbations are the only source of mechanism variation inside DECAY.  

If deviations appear, check:

- environment consistency  
- local code modifications  
- harness version  

---

# 7. 📜 Citation

If you use this reproduction in research:

**Proxmire, A. (2026).  
*ED‑Arch‑21 — Compositionality, y‑Switching, and Architectural Robustness*.**
