# 📘 Law II — Islands, Ridges, and Resonance  
### *Reproducibility Folder*

This folder reproduces the results of:

**ED‑Arch‑23 — Law II: Islands, Ridges, and Resonance**  
(see: `../../docs/ED‑Arch‑23_Law II_Islands Ridges and Resonance.pdf`)

Law II explains the **resonance structure** that interrupts the tangent‑gate ladder:  
the periodic DECAY islands that give the tangent regime its serrated, broken‑ridgeline geometry.

---

# 1. 🔍 What Law II Establishes

### **1. DECAY as a Resonance Phenomenon**
DECAY occurs when the circumference–box ratio satisfies a rational relationship:


\[
\frac{C(d)}{L} \approx \frac{p}{q}
\]


Small‑denominator rationals produce strong destructive alignment → DECAY.

### **2. DECAY Island Spacing**
For N = 3:


\[
\Delta d \approx 0.04
\]


This spacing is:
- stable  
- quasi‑periodic  
- resonance‑driven  
- predictable from geometry  

### **3. Ladder + DECAY = Broken Ridgeline**
The tangent gate is not a surface — it is:
- quantized ladder rungs  
- interrupted by DECAY islands  
- interleaved with narrow outward‑PBC pockets  

This produces the characteristic serrated structure.

### **4. Resonance as a Secondary Law Structure**
Law I gives the ladder.  
Law II explains why the ladder is *not contiguous*.

---

# 2. 🧪 What This Folder Tests

Running the harness reproduces:

### ✔ DECAY Island Structure
- correct spacing  
- correct positions  
- correct mechanism classification  

### ✔ Resonance Geometry
- rational‑ratio alignment  
- destructive vs constructive alignment  

### ✔ Ladder Interruptions
- DECAY islands appear between rungs  
- outward‑PBC pockets appear at edges  

### ✔ Mechanism Map
- inward‑collapse  
- DECAY  
- outward‑PBC  
- other‑late  

These are the empirical signatures of Law II.

---

# 3. ▶️ How to Run the Reproduction

From the **repo root**:

```bash
python run_arch_harness.py --suite law-surfaces
```

This command:

- samples tangent‑gate collapse  
- detects DECAY islands  
- identifies ladder rungs  
- classifies mechanisms  
- writes JSON + TXT reports  
- generates canonical figures  

Outputs appear in:

```
arch_harness_output/
```

If you prefer a local wrapper:

```bash
python run_law_II.py
```

---

# 4. 📂 Expected Outputs

### **JSON Reports**
- `law_surface_results.json`  
- `decay_window_results.json`  

### **Figures**
- DECAY island map  
- broken‑ridgeline plot  
- resonance‑spacing visualization  
- mechanism classification map  

Canonical copies live in:

```
expected_outputs/
    json/
    figures/
```

---

# 5. 🧭 Mapping to the Paper

This reproducibility folder corresponds to the following sections of ED‑Arch‑23:

| Paper Section | Reproduced Here |
|---------------|-----------------|
| DECAY island identification | ✔ DECAY detection |
| DECAY spacing (~0.04) | ✔ spacing validation |
| Resonance condition \( C(d)/L \approx p/q \) | ✔ resonance geometry |
| Ladder interruptions | ✔ broken‑ridgeline reproduction |
| Tangent‑gate mechanism map | ✔ mechanism classification |

This folder is the **computational companion** to the paper.

---

# 6. 📌 Notes for Researchers

- DECAY islands must appear at the canonical quasi‑periodic spacing.  
- Ladder rungs must be discrete and correctly interrupted.  
- Outward‑PBC pockets appear at predictable resonance edges.  
- All results are deterministic.  

If deviations appear, check:
- environment consistency  
- local code modifications  
- harness version  

---

# 7. 📜 Citation

If you use this reproduction in research:

**Proxmire, A. (2026).  
*ED‑Arch‑23 — Law II: Islands, Ridges, and Resonance*.**
