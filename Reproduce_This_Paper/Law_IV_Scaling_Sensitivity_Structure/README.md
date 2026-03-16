# 📘 Law IV — Scaling, Sensitivity, and Structure  
### *Reproducibility Folder*

This folder reproduces the results of:

**ED‑Arch‑25 — Law IV: Scaling, Sensitivity, and Structure**  
(see: `../../docs/ED‑Arch‑25_Law IV_Scaling Sensitivity and Structure.pdf`)

Law IV derives the **N‑scaling laws** that govern how collapse surfaces transform as the number of particles increases.  
It shows that timing, resonance, recurrence, and sensitivity all scale predictably with \( N \), revealing universal geometric patterns.

---

# 1. 🔍 What Law IV Establishes

### **1. Circumference Scaling**
For a ring of \( N \) particles:


\[
C(d, N) \propto N
\]


This linear scaling drives all downstream effects.

### **2. Drift‑Phase Velocity Scaling**


\[
\dot{\phi}(N) \propto N
\]


Larger rings drift faster → phase accumulates faster → windows narrow.

### **3. Recurrence Time Scaling**


\[
T_r(N) \propto \frac{1}{N}
\]


This explains why ladder spacing shrinks as \( N \) increases.

### **4. Window‑Width Scaling**


\[
W(N) \propto \frac{1}{N}
\]


Narrower windows → greater sensitivity → more DECAY islands.

### **5. DECAY Spacing Scaling**


\[
\Delta d(N) \propto \frac{1}{N}
\]


More particles → more resonance → denser DECAY structure.

### **6. Ladder Density Scaling**


\[
\text{rung density} \propto N
\]


Larger rings produce more rungs in the same time interval.

### **7. Universal Collapse Surfaces**
When collapse surfaces are rescaled by \( N \), they converge to **universal curves**.  
ED‑Arch collapse is a **self‑similar geometric family** across particle count.

---

# 2. 🧪 What This Folder Tests

Running the harness reproduces:

### ✔ N‑Scaling of Recurrence  
- correct \( T_r(N) \propto 1/N \)  
- correct ladder spacing shrinkage  

### ✔ N‑Scaling of Window Width  
- windows narrow as \( 1/N \)  
- sensitivity increases with \( N \)

### ✔ N‑Scaling of DECAY  
- DECAY islands become more frequent  
- spacing shrinks as \( 1/N \)

### ✔ Collapse‑Surface Transformation  
- surfaces sharpen  
- boundaries shift predictably  
- ridgeline becomes more intricate  

### ✔ Universal Curves  
- rescaled collapse surfaces align  
- self‑similar structure emerges  

These are the empirical signatures of Law IV.

---

# 3. ▶️ How to Run the Reproduction

From the **repo root**:

```bash
python run_arch_harness.py --suite law-surfaces
```

This command:

- samples collapse surfaces across multiple \( N \)  
- measures recurrence times  
- detects ladder rungs  
- identifies DECAY islands  
- writes JSON + TXT reports  
- generates canonical figures  

Outputs appear in:

```
arch_harness_output/
```

If you prefer a local wrapper:

```bash
python run_law_IV.py
```

---

# 4. 📂 Expected Outputs

### **JSON Reports**
- `law_surface_results.json`  
- `decay_window_results.json`  

### **Figures**
- recurrence‑time scaling  
- window‑width scaling  
- DECAY‑density scaling  
- ladder‑density scaling  
- universal‑curve overlays  
- mechanism maps across \( N \)

Canonical copies live in:

```
expected_outputs/
    json/
    figures/
```

---

# 5. 🧭 Mapping to the Paper

This reproducibility folder corresponds to the following sections of ED‑Arch‑25:

| Paper Section | Reproduced Here |
|---------------|-----------------|
| Circumference scaling | ✔ drift‑phase scaling |
| Recurrence scaling | ✔ ladder spacing |
| Window‑width scaling | ✔ sensitivity validation |
| DECAY spacing scaling | ✔ resonance‑density validation |
| Collapse‑surface sharpening | ✔ mechanism maps |
| Universal curves | ✔ rescaled‑surface overlays |

This folder is the **computational companion** to the paper.

---

# 6. 📌 Notes for Researchers

- Ladder spacing must shrink as \( 1/N \).  
- DECAY islands must become denser with increasing \( N \).  
- Recurrence windows must narrow as \( 1/N \).  
- Collapse surfaces must sharpen and converge under rescaling.  
- All results are deterministic.  

If deviations appear, check:
- environment consistency  
- local code modifications  
- harness version  

---

# 7. 📜 Citation

If you use this reproduction in research:

**Proxmire, A. (2026).  
*ED‑Arch‑25 — Law IV: Scaling, Sensitivity, and Structure*.**
