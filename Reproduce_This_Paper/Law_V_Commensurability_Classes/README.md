# 📘 Law V — Commensurability Classes  
### *Reproducibility Folder*

This folder reproduces the results of:

**ED‑Arch‑26 — Law V: Commensurability Classes**  
(see: `../../docs/ED‑Arch‑26_Law V_Commensurability Classes.pdf`)

Law V is the **classification law** of ED‑Arch.  
It identifies the rational and irrational structures that determine which particle diameters collapse early, which drift indefinitely, which produce DECAY, and which lock into periodic orbits.

This is the capstone of the ED‑Arch law series.

---

# 1. 🔍 What Law V Establishes

### **1. Commensurability as the Organizing Principle**
Collapse behavior clusters according to the rational structure of:


\[
\frac{C(d, N)}{L}
\]


where \( C(d, N) \) is the ring circumference and \( L \) is the box length.

This ratio determines whether drift phase:
- returns to alignment  
- nearly returns  
- or never returns  

### **2. Four Commensurability Regimes**

#### **Constructive Commensurability**
- low‑denominator rationals  
- strong periodicity  
- early collapse  
- robust alignment  

#### **Destructive Commensurability**
- high‑denominator rationals  
- fragile periodicity  
- repeated window misses  
- DECAY  

#### **Drift‑Locked Classes**
- near‑rational ratios  
- stable periodic orbits  
- collapse only on specific rungs  

#### **Quasi‑Periodic Classes**
- irrational ratios  
- no exact return  
- long drift  
- collapse only via near‑miss  

### **3. The Commensurability Atlas**
Diameter space partitions into a structured atlas of:
- constructive bands  
- destructive bands  
- drift‑locked islands  
- quasi‑periodic seas  

### **4. Scaling with N**
As \( N \) increases:
- class density grows like \( N \)  
- class width shrinks like \( 1/N \)  
- the atlas becomes more intricate  
- universal patterns emerge under rescaling  

### **5. Law V as the Capstone**
Law V unifies Laws I–IV by explaining **which diameters activate which mechanisms**.

---

# 2. 🧪 What This Folder Tests

Running the harness reproduces:

### ✔ Commensurability Classification
- constructive vs destructive  
- drift‑locked vs quasi‑periodic  

### ✔ Collapse‑Behavior Clustering
- early collapse bands  
- DECAY bands  
- periodic‑orbit bands  

### ✔ Resonance & Recurrence Interaction
- rational alignment → constructive  
- rational misalignment → destructive  
- irrational drift → quasi‑periodic  

### ✔ N‑Scaling of Classes
- class density ∝ \( N \)  
- class width ∝ \( 1/N \)  

### ✔ Mechanism Map
- inward‑collapse  
- DECAY  
- drift‑locked periodicity  
- quasi‑periodic drift  

These are the empirical signatures of Law V.

---

# 3. ▶️ How to Run the Reproduction

From the **repo root**:

```bash
python run_arch_harness.py --suite law-surfaces
```

This command:

- samples collapse behavior across diameter space  
- identifies commensurability classes  
- detects constructive/destructive regimes  
- identifies drift‑locked and quasi‑periodic classes  
- writes JSON + TXT reports  
- generates canonical figures  

Outputs appear in:

```
arch_harness_output/
```

If you prefer a local wrapper:

```bash
python run_law_V.py
```

---

# 4. 📂 Expected Outputs

### **JSON Reports**
- `law_surface_results.json`  
- `decay_window_results.json`  

### **Figures**
- commensurability atlas  
- constructive/destructive band map  
- drift‑locked periodicity map  
- quasi‑periodic drift map  
- N‑scaling overlays  

Canonical copies live in:

```
expected_outputs/
    json/
    figures/
```

---

# 5. 🧭 Mapping to the Paper

This reproducibility folder corresponds to the following sections of ED‑Arch‑26:

| Paper Section | Reproduced Here |
|---------------|-----------------|
| Constructive classes | ✔ early‑collapse bands |
| Destructive classes | ✔ DECAY bands |
| Drift‑locked classes | ✔ periodic‑orbit detection |
| Quasi‑periodic classes | ✔ long‑drift detection |
| Rational vs irrational structure | ✔ resonance/recurrence interaction |
| N‑scaling of class density | ✔ scaling validation |
| Commensurability atlas | ✔ full mechanism map |

This folder is the **computational companion** to the paper.

---

# 6. 📌 Notes for Researchers

- Constructive classes correspond to low‑denominator rationals.  
- Destructive classes correspond to high‑denominator rationals.  
- Drift‑locked classes appear at near‑rational ratios.  
- Quasi‑periodic classes correspond to irrational ratios.  
- Class density must grow with \( N \).  
- Class width must shrink with \( 1/N \).  
- All results are deterministic.  

If deviations appear, check:
- environment consistency  
- local code modifications  
- harness version  

---

# 7. 📜 Citation

If you use this reproduction in research:

**Proxmire, A. (2026).  
*ED‑Arch‑26 — Law V: Commensurability Classes*.**
