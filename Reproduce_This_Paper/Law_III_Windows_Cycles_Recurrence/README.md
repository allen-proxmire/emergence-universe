# 📘 Law III — Windows, Cycles, and Recurrence  
### *Reproducibility Folder*

This folder reproduces the results of:

**ED‑Arch‑24 — Law III: Windows, Cycles, and Recurrence**  
(see: `../../docs/ED‑Arch‑24_Law III_Windows Cycles and Recurrence.pdf`)

Law III derives the **recurrence geometry** that creates the tangent‑gate ladder and explains why collapse times are quantized.

---

# 1. 🔍 What Law III Establishes

### **1. Recurrence Windows**
Collapse can occur **only** when the drifting ring re‑enters a PBC‑equivalent configuration:


\[
t = k T_r
\]


These discrete windows are the *only* times collapse is geometrically possible.

### **2. Ladder Step Size**
Law I observed:


\[
\chi = \chi_0 + k\Delta
\]


Law III proves:


\[
\Delta = T_r
\]


The ladder spacing *is* the recurrence time.

### **3. Window Width**
Recurrence windows have width:


\[
W \propto \frac{1}{\dot{\phi}}
\]


and since drift‑phase velocity grows with \( N \):


\[
W(N) \propto \frac{1}{N}
\]


Larger rings → narrower windows → greater sensitivity.

### **4. Interaction with Resonance (Law II)**
- Constructive alignment → ladder rung  
- Destructive alignment → DECAY island  

Recurrence determines **when** collapse can occur.  
Resonance determines **whether** collapse occurs at that time.

### **5. Temporal Scaffold of the Tangent Gate**
The tangent gate is governed by:
- recurrence (Law III)  
- quantized timing (Law I)  
- rational resonance (Law II)  

Together they form the full architecture of tangent‑gate collapse.

---

# 2. 🧪 What This Folder Tests

Running the harness reproduces:

### ✔ Recurrence Timing
- correct recurrence period \( T_r \)  
- correct ladder spacing \( \Delta = T_r \)

### ✔ Window Geometry
- correct window width scaling  
- correct sensitivity to \( d \) and \( N \)

### ✔ Ladder Structure
- discrete rungs  
- correct rung index assignment  

### ✔ Recurrence–Resonance Interaction
- DECAY when windows are missed  
- collapse when windows are hit  

### ✔ Mechanism Map
- inward‑collapse  
- DECAY  
- outward‑PBC  
- other‑late  

These are the empirical signatures of Law III.

---

# 3. ▶️ How to Run the Reproduction

From the **repo root**:

```bash
python run_arch_harness.py --suite law-surfaces
```

This command:

- samples tangent‑gate collapse  
- identifies recurrence windows  
- detects ladder rungs  
- classifies mechanisms  
- writes JSON + TXT reports  
- generates canonical figures  

Outputs appear in:

```
arch_harness_output/
```

If you prefer a local wrapper:

```bash
python run_law_III.py
```

---

# 4. 📂 Expected Outputs

### **JSON Reports**
- `law_surface_results.json`  
- `decay_window_results.json`  

### **Figures**
- recurrence‑window diagram  
- ladder spacing plot  
- window‑width scaling plot  
- mechanism classification map  

Canonical copies live in:

```
expected_outputs/
    json/
    figures/
```

---

# 5. 🧭 Mapping to the Paper

This reproducibility folder corresponds to the following sections of ED‑Arch‑24:

| Paper Section | Reproduced Here |
|---------------|-----------------|
| Recurrence time \( T_r \) | ✔ ladder spacing |
| Window width | ✔ window‑geometry validation |
| Constructive vs destructive alignment | ✔ recurrence–resonance interaction |
| Ladder quantization | ✔ rung detection |
| Tangent‑gate mechanism map | ✔ mechanism classification |

This folder is the **computational companion** to the paper.

---

# 6. 📌 Notes for Researchers

- Ladder spacing must equal the recurrence time.  
- Recurrence windows must be narrow and \( \propto 1/N \).  
- DECAY occurs when windows are repeatedly missed.  
- All results are deterministic.  

If deviations appear, check:
- environment consistency  
- local code modifications  
- harness version  

---

# 7. 📜 Citation

If you use this reproduction in research:

**Proxmire, A. (2026).  
*ED‑Arch‑24 — Law III: Windows, Cycles, and Recurrence*.**
