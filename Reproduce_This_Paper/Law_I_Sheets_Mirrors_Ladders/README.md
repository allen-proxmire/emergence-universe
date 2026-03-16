# 📘 Law I — Sheets, Mirrors, and Ladders  
### *Reproducibility Folder*

This folder reproduces the results of:

**ED‑Arch‑22 — Law I: Sheets, Mirrors, and Ladders**  
(see: `../../docs/ED‑Arch‑22_Law I_Sheets Mirrors and Ladders.pdf`)

Law I reveals the **three fundamental geometric law structures** of ED‑Arch collapse:  
a continuous sheet, a reflected mirror, and a quantized ladder.

---

# 1. 🔍 What Law I Establishes

### **1. The Sheet (Inward Gate)**
A smooth, continuous, metric‑driven law:
- collapse time increases linearly with diameter  
- no discontinuities  
- no resonances  
- most stable and predictable regime  

### **2. The Mirror (Outward Gate)**
A reflected linear law:
- slope magnitude identical to the inward sheet  
- collapse driven by PBC wrap‑around  
- bilateral symmetry between inward and outward  

### **3. The Ladder (Tangent Gate)**
A quantized timing law:


\[
\chi = \chi_0 + k\Delta
\]


with:
- discrete rungs  
- recurrence‑driven timing  
- commensurability‑dependent rung selection  

### **4. The Broken Ridgeline**
The ladder is interrupted by:
- DECAY islands  
- outward‑PBC pockets  
- quasi‑periodic gaps  

This produces the characteristic serrated structure of the tangent regime.

---

# 2. 🧪 What This Folder Tests

Running the harness reproduces:

### ✔ Inward Sheet  
- linear collapse law  
- smooth dependence on \( d \)  

### ✔ Outward Mirror  
- reflected linear law  
- matching slope magnitude  

### ✔ Tangent Ladder  
- discrete rungs  
- correct rung spacing  
- correct rung index assignment  

### ✔ DECAY Islands  
- correct spacing  
- correct mechanism classification  

### ✔ Mechanism Map  
- inward‑collapse  
- outward‑PBC  
- DECAY  
- other‑late  

These are the core empirical signatures of Law I.

---

# 3. ▶️ How to Run the Reproduction

From the **repo root**:

```bash
python run_arch_harness.py --suite law-surfaces
```

This command:

- samples the law surfaces for all three gates  
- identifies ladder rungs  
- detects DECAY islands  
- classifies mechanisms  
- writes JSON + TXT reports  
- generates canonical figures  

Outputs appear in:

```
arch_harness_output/
```

If you prefer a local wrapper:

```bash
python run_law_I.py
```

---

# 4. 📂 Expected Outputs

### **JSON Reports**
- `law_surface_results.json`  
- `decay_window_results.json`  
- `angular_rate_results.json`  

### **Figures**
- inward sheet plot  
- outward mirror plot  
- tangent ladder plot  
- broken‑ridgeline mechanism map  

Canonical copies live in:

```
expected_outputs/
    json/
    figures/
```

---

# 5. 🧭 Mapping to the Paper

This reproducibility folder corresponds to the following sections of ED‑Arch‑22:

| Paper Section | Reproduced Here |
|---------------|-----------------|
| Inward sheet law | ✔ sheet validation |
| Outward mirror law | ✔ mirror validation |
| Tangent ladder | ✔ rung detection |
| DECAY islands | ✔ DECAY identification |
| Broken ridgeline | ✔ mechanism map |
| Three‑geometry structure | ✔ all gates reproduced |

This folder is the **computational companion** to the paper.

---

# 6. 📌 Notes for Researchers

- Ladder rungs must be **discrete**, not continuous.  
- Mirror slope must match inward slope in magnitude.  
- DECAY islands appear at predictable quasi‑periodic spacing.  
- Mechanism classification should match the canonical map exactly.  
- All results are deterministic.  

If deviations appear, check:
- environment consistency  
- local code modifications  
- harness version  

---

# 7. 📜 Citation

If you use this reproduction in research:

**Proxmire, A. (2026).  
*ED‑Arch‑22 — Law I: Sheets, Mirrors, and Ladders*.**
