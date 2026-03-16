# 📘 ED‑Arch‑20 — Global Law Surfaces  
### *Reproducibility Folder*

This folder reproduces the results of:

**ED‑Arch‑20 — Global Law Surfaces of the Event Density Manifold**  
(see: `../../docs/ED‑Arch‑20_NOFIGS_Global Law Surfaces_Manifold.pdf`)

ED‑Arch‑20 is the chapter where the ED manifold becomes a **closed architectural object**:  
all three sheets, all boundaries, all mechanisms, and the collapse‑rate law are mapped, bounded, and derived.

---

# 1. 🔍 What ED‑Arch‑20 Establishes

### **The Three‑Sheet ED Manifold**
- **y = +1** → inward‑collapse sheet  
- **y = 0** → tangential saddle (inward‑spiral, DECAY, PBC‑corner, satellite‑PBC)  
- **y = −1** → outward‑PBC sheet  

### **Thirteen Mechanism Boundaries**
All resolved to ±0.5 px, including:
- MERGE_THR  
- topology‑dependent `d_lower`, `d_upper`  
- `chi_zero_d`  
- saddle‑sheet transitions  

### **Derived Angular‑Rate Law**


\[
K = 2\sin(\pi / n_{cl})
\]


and the invariant:


\[
\omega_N(d)\, d = K
\]



### **Linearity on y = ±1**
Collapse time is strictly linear in \( d \), with slopes determined solely by \( n_{cl} \).

### **Full Manifold Closure**
The ED manifold is:
- mapped  
- bounded  
- derived  
- self‑consistent  

This is the chapter where ED transitions from empirical mapping to **predictive theory**.

---

# 2. 🧪 What This Folder Tests

Running the harness from this folder reproduces:

### ✔ Analytic Law‑Surface Validation
- 56 sampled points  
- `max_deviation = 0.0`  
- mechanism predictions match exactly  

### ✔ DECAY Window Boundaries
- correct `d_lower`, `d_upper` for each topology  
- `boundary_accuracy = 1.0`  

### ✔ Angular‑Rate Law Validation
- predicted vs empirical \( K \)  
- `max_deviation = 0.0`  

### ✔ Mechanism Classification
- inward‑collapse  
- outward‑PBC  
- DECAY  
- PBC‑corner / satellite‑PBC  

These are the core empirical signatures of ED‑Arch‑20.

---

# 3. ▶️ How to Run the Reproduction

From the **repo root**:

```bash
python run_arch_harness.py --suite law-surfaces
```

This command:

- samples the analytic law surfaces  
- computes DECAY windows  
- validates the angular‑rate law  
- classifies mechanisms  
- writes JSON + TXT reports  
- generates canonical figures  

Outputs appear in:

```
arch_harness_output/
```

If you prefer a local wrapper:

```bash
python run_law_20.py
```

---

# 4. 📂 Expected Outputs

After running the harness, you should see:

### **JSON Reports**
- `law_surface_results.json`  
- `decay_window_results.json`  
- `angular_rate_results.json`  

### **Figures**
- `fig_law_surface_mechanisms.png`  
- `fig_decay_windows.png`  
- `fig_angular_rate_law.png`  

### **Summary**
- `arch_harness_report.txt`  

Canonical copies of these outputs are stored in:

```
expected_outputs/
    json/
    figures/
```

Use these to verify correctness.

---

# 5. 🧭 Mapping to the Paper

This reproducibility folder corresponds to the following sections of ED‑Arch‑20:

| Paper Section | Reproduced Here |
|---------------|-----------------|
| Three‑sheet mechanism map | ✔ mechanism classification |
| Boundary catalog (13 boundaries) | ✔ DECAY windows |
| Derived angular‑rate law | ✔ angular‑rate validation |
| Linearity on y = ±1 | ✔ law‑surface sampling |
| Full manifold closure | ✔ all tests pass |

This folder is the **computational companion** to the paper.

---

# 6. 📌 Notes for Researchers

- All results are deterministic.  
- No random seeds are used.  
- The analytic oracle is exact for the tested domain.  
- Deviations should be **0.0** for all analytic tests.  
- Boundary accuracy should be **1.0**.  
- Mechanism classifications should match the paper exactly.  

If any deviation appears, it indicates:
- a local code modification  
- a mismatch in environment  
- or a regression in the harness  

---

# 7. 📜 Citation

If you use this reproduction in research:

**Proxmire, A. (2026).  
*ED‑Arch‑20 — Global Law Surfaces of the Event Density Manifold*.**
