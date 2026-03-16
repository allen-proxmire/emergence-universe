#!/usr/bin/env python
"""
check_consolidation_sweep.py
==============================
Post-Law VIII Consolidation Analysis: determines whether the unified manifold
remains coherent beyond the canonical parameter grid.

Loads ALL data (canonical + extended), then runs 11 sections:
  1.  Inventory & coverage
  2.  Complexity flow at higher N (N=24, 28, 32)
  3.  Monotonicity test: does complexity decrease for every metric?
  4.  N=8 resonance persistence check
  5.  DECAY confinement (still N=4 only?)
  6.  Fine-angle structure near tangent (new folds or cusps?)
  7.  Fine-drift structure (new instabilities?)
  8.  Temporal fold persistence at N>=8
  9.  Unified manifold coherence at extended N
 10.  Taxonomy collapse: do sub-mechanisms simplify further?
 11.  Summary: deviations from Law VIII
"""

import json
import os
import glob
import math
import statistics
from collections import defaultdict, Counter

SWEEP_DIR = os.path.dirname(__file__)

# Full parameter lists (canonical + extended)
NS_CANONICAL = [4, 8, 12, 20]
NS_HIGH      = [24, 28, 32]
NS_ALL       = NS_CANONICAL + NS_HIGH

ANGLES_CANONICAL = [-10, -5, -2, -1, 0, 1, 2, 5, 10]
ANGLES_FINE      = [-0.75, -0.5, -0.25, 0.25, 0.5, 0.75]
ANGLES_ALL       = sorted(set(ANGLES_CANONICAL + ANGLES_FINE))

RADII = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 0.90, 0.95]

DRIFTS_CANONICAL = [0.00, 0.01, 0.02, 0.05, 0.10]
DRIFTS_FINE      = [0.005, 0.015, 0.025]
DRIFTS_ALL       = sorted(set(DRIFTS_CANONICAL + DRIFTS_FINE))

MECHS = ["inward-collapse", "outward-PBC", "DECAY", "PBC-corner", "other-late"]
SHORT = {"inward-collapse": "IC", "outward-PBC": "OP", "DECAY": "DE",
         "PBC-corner": "CR", "other-late": "OL"}


# ============================================================
# Data loading
# ============================================================

def load_all_sweeps():
    """Load ALL sweep files (canonical + extended)."""
    cube = {}
    pattern = os.path.join(SWEEP_DIR, "n*_angle_*_r*_d*_sweep.json")
    for path in sorted(glob.glob(pattern)):
        with open(path) as f:
            data = json.load(f)
        key = (data["N"], data["angle"], data["radius"], data["drift"])
        cube[key] = data
    return cube


def load_all_atlas():
    """Load ALL atlas files (canonical + extended)."""
    cube = {}
    pattern = os.path.join(SWEEP_DIR, "n*_angle_*_r*_d*_atlas.json")
    for path in sorted(glob.glob(pattern)):
        with open(path) as f:
            data = json.load(f)
        if "atlas" not in data:
            continue
        key = (data["N"], data["angle"], data["radius"], data["drift"])
        cube[key] = data
    return cube


# ============================================================
# Sub-mechanism classification (identical to unification script)
# ============================================================

def classify_sub(rec):
    mech = rec["mechanism"]
    a = rec["atlas"]
    pbc = a["n_pbc_crossings"]
    steps = a["n_steps"]
    seq = a["corner_sequence"]
    period = a["periodicity_score"]
    multi = a["is_multi_stage"]
    tcross = a["tangent_crossings"]

    if mech == "inward-collapse":
        if pbc > 0:           return "IC-bounced"
        elif steps > 50:      return "IC-delayed"
        else:                 return "IC-direct"
    elif mech == "outward-PBC":
        if pbc == 0:          return "OPBC-single"
        elif pbc <= 5:        return "OPBC-glancing"
        else:                 return "OPBC-winding"
    elif mech == "PBC-corner":
        distinct = len(set(seq))
        diag_pairs = ({0, 3}, {1, 2})
        seq_set = set(seq)
        if distinct <= 1:     return "CR-single"
        elif any(dp <= seq_set for dp in diag_pairs):
                              return "CR-diagonal"
        else:                 return "CR-multi"
    elif mech == "other-late":
        if multi:             return "OL-multi"
        elif period > 0.3:    return "OL-periodic"
        elif steps < 5:       return "OL-stalled"
        else:                 return "OL-wandering"
    elif mech == "DECAY":
        if tcross > 0:        return "DE-tangent"
        else:                 return "DE-frozen"
    return mech


# ============================================================
# Derived quantities
# ============================================================

def compute_ic_onset(sweep_cube, ns, angles):
    """For each (N, radius, drift), find onset angle of IC."""
    sorted_angles = sorted(angles)
    onset = {}
    for n in ns:
        for r in RADII:
            for d in DRIFTS_ALL:
                ic_angs = [ang for ang in sorted_angles
                           if (n, ang, r, d) in sweep_cube
                           and sweep_cube[(n, ang, r, d)]["mechanism"] == "inward-collapse"]
                if not ic_angs:
                    onset[(n, r, d)] = None
                else:
                    onset[(n, r, d)] = max(ic_angs)
    return onset


def signed_boundary_dist(angle, onset_angle):
    if onset_angle is None:
        return float('nan')
    return angle - onset_angle


def drift_sensitivity_class(onset_d0, onset_d):
    if onset_d0 is None or onset_d is None:
        return "no-boundary"
    shift = abs(onset_d - onset_d0)
    if shift == 0:    return "rigid"
    elif shift <= 2:  return "stable"
    elif shift <= 5:  return "moderate"
    else:             return "fragile"


def temporal_gradient_class(chi, chi_neighbor_mean):
    if chi_neighbor_mean is None or chi_neighbor_mean < 0.001:
        return "isolated"
    ratio = chi / chi_neighbor_mean
    if ratio > 5.0:   return "cliff"
    elif ratio > 2.0:  return "steep"
    elif ratio > 1.2:  return "ramp"
    elif ratio > 0.8:  return "flat"
    elif ratio > 0.5:  return "ramp-down"
    else:              return "cliff-down"


def local_curvature(onset, n, d, r_idx):
    if r_idx < 1 or r_idx >= len(RADII) - 1:
        return 0.0
    r0, r1, r2 = RADII[r_idx - 1], RADII[r_idx], RADII[r_idx + 1]
    a0 = onset.get((n, r0, d))
    a1 = onset.get((n, r1, d))
    a2 = onset.get((n, r2, d))
    if a0 is None or a1 is None or a2 is None:
        return 0.0
    max_ang = max(ANGLES_ALL)
    if a0 > max_ang or a1 > max_ang or a2 > max_ang:
        return 0.0
    dx1, dy1 = r1 - r0, float(a1 - a0)
    dx2, dy2 = r2 - r1, float(a2 - a1)
    ddx, ddy = dx2 - dx1, dy2 - dy1
    dxa, dya = (dx1 + dx2) / 2, (dy1 + dy2) / 2
    denom = (dxa**2 + dya**2)**1.5
    if denom < 1e-15:
        return 0.0
    return abs(dxa * ddy - dya * ddx) / denom


def local_torsion(onset, n, r_idx):
    """Torsion: how much the boundary curve twists through the drift dimension."""
    if r_idx < 1 or r_idx >= len(RADII) - 1:
        return 0.0
    r1 = RADII[r_idx]
    # Need at least 3 drift values
    dvals = sorted(d for d in DRIFTS_ALL if onset.get((n, r1, d)) is not None)
    if len(dvals) < 3:
        return 0.0
    torsion_sum = 0.0
    count = 0
    for i in range(1, len(dvals) - 1):
        d0, d1, d2 = dvals[i-1], dvals[i], dvals[i+1]
        a0 = onset.get((n, r1, d0))
        a1 = onset.get((n, r1, d1))
        a2 = onset.get((n, r1, d2))
        max_ang = max(ANGLES_ALL)
        if a0 is None or a1 is None or a2 is None:
            continue
        if a0 > max_ang or a1 > max_ang or a2 > max_ang:
            continue
        dd1 = float(a1 - a0) / (d1 - d0) if d1 != d0 else 0
        dd2 = float(a2 - a1) / (d2 - d1) if d2 != d1 else 0
        torsion_sum += abs(dd2 - dd1)
        count += 1
    return torsion_sum / count if count > 0 else 0.0


# ============================================================
# Assemble unified tensor (all data)
# ============================================================

def assemble_tensor(sweep_cube, atlas_cube):
    """Build unified record for every cell (canonical + extended)."""
    # Determine the set of angles and N values present
    all_angles_present = sorted(set(k[1] for k in sweep_cube))
    all_ns_present = sorted(set(k[0] for k in sweep_cube))

    onset = compute_ic_onset(sweep_cube, all_ns_present, all_angles_present)

    # Precompute chi neighbor means
    chi_neighbors = {}
    sorted_angs = sorted(all_angles_present)
    for key, rec in sweep_cube.items():
        n, ang, r, d = key
        neighbors = []
        if ang in sorted_angs:
            ang_idx = sorted_angs.index(ang)
            for da in [-1, +1]:
                ni = ang_idx + da
                if 0 <= ni < len(sorted_angs):
                    nkey = (n, sorted_angs[ni], r, d)
                    if nkey in sweep_cube:
                        neighbors.append(sweep_cube[nkey]["chi_emp"])
        r_idx = RADII.index(r) if r in RADII else -1
        for dr in [-1, +1]:
            ri = r_idx + dr
            if 0 <= ri < len(RADII):
                nkey = (n, ang, RADII[ri], d)
                if nkey in sweep_cube:
                    neighbors.append(sweep_cube[nkey]["chi_emp"])
        chi_neighbors[key] = statistics.mean(neighbors) if neighbors else None

    tensor = {}
    for key in sweep_cube:
        n, ang, r, d = key
        rec = sweep_cube[key]
        mech = rec["mechanism"]
        chi = rec["chi_emp"]

        subtype = "unknown"
        if key in atlas_cube:
            subtype = classify_sub(atlas_cube[key])

        o = onset.get((n, r, d))
        bdist = signed_boundary_dist(ang, o)

        o_base = onset.get((n, r, 0.0))
        o_this = onset.get((n, r, d))
        dc = drift_sensitivity_class(o_base, o_this)

        tc = temporal_gradient_class(chi, chi_neighbors[key])

        r_idx = RADII.index(r) if r in RADII else -1
        kappa = local_curvature(onset, n, d, r_idx)
        tau = local_torsion(onset, n, r_idx)

        lifetime = 0
        if key in atlas_cube:
            lifetime = atlas_cube[key]["atlas"]["n_steps"]

        tensor[key] = {
            "N": n, "angle": ang, "radius": r, "drift": d,
            "mechanism": mech,
            "subtype": subtype,
            "chi": chi,
            "lifetime": lifetime,
            "boundary_dist": bdist,
            "drift_class": dc,
            "temporal_class": tc,
            "curvature": kappa,
            "torsion": tau,
        }

    return tensor, onset


# ============================================================
# 1. Inventory & coverage
# ============================================================

def sec01_inventory(tensor):
    print("=" * 72)
    print("SECTION 1: INVENTORY & COVERAGE")
    print("=" * 72)

    all_ns = sorted(set(v["N"] for v in tensor.values()))
    print(f"  Total cells: {len(tensor)}")
    print(f"  N values present: {all_ns}")

    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        angles = sorted(set(v["angle"] for v in cells))
        drifts = sorted(set(v["drift"] for v in cells))
        mechs = Counter(v["mechanism"] for v in cells)
        print(f"\n  N={n}: {len(cells)} cells, "
              f"{len(angles)} angles, {len(drifts)} drifts")
        print(f"    Angles: {angles}")
        print(f"    Drifts: {drifts}")
        print(f"    Mechanisms: {dict(mechs.most_common())}")

    # How many are canonical vs extended?
    canonical = sum(1 for v in tensor.values()
                    if v["N"] in NS_CANONICAL
                    and v["angle"] in ANGLES_CANONICAL
                    and v["drift"] in DRIFTS_CANONICAL)
    extended = len(tensor) - canonical
    print(f"\n  Canonical grid cells: {canonical}")
    print(f"  Extended cells:       {extended}")


# ============================================================
# 2. Complexity flow at higher N
# ============================================================

def sec02_complexity_flow(tensor):
    print("\n" + "=" * 72)
    print("SECTION 2: COMPLEXITY FLOW AT HIGHER N")
    print("=" * 72)

    all_ns = sorted(set(v["N"] for v in tensor.values()))

    print(f"\n  {'N':>4s}  {'cells':>5s}  {'mechs':>5s}  {'subs':>5s}  "
          f"{'mean_chi':>8s}  {'med_chi':>8s}  {'max_chi':>8s}  "
          f"{'IC%':>5s}  {'OL%':>5s}")

    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        nc = len(cells)
        n_mechs = len(set(v["mechanism"] for v in cells))
        n_subs = len(set(v["subtype"] for v in cells if v["subtype"] != "unknown"))
        chis = [v["chi"] for v in cells]
        ic_pct = 100 * sum(1 for v in cells if v["mechanism"] == "inward-collapse") / nc
        ol_pct = 100 * sum(1 for v in cells if v["mechanism"] == "other-late") / nc
        print(f"  {n:4d}  {nc:5d}  {n_mechs:5d}  {n_subs:5d}  "
              f"{statistics.mean(chis):8.2f}  {statistics.median(chis):8.2f}  "
              f"{max(chis):8.2f}  {ic_pct:5.1f}  {ol_pct:5.1f}")

    # Compare high-N mechanism distributions
    print(f"\n  Mechanism distribution at high N (canonical angles, drift=0):")
    for n in all_ns:
        cells = [v for v in tensor.values()
                 if v["N"] == n and v["angle"] in ANGLES_CANONICAL
                 and v["drift"] == 0.0]
        if not cells:
            continue
        mechs = Counter(v["mechanism"] for v in cells)
        total = len(cells)
        parts = [f"{SHORT.get(k,k)}={100*ct/total:.0f}%" for k, ct in mechs.most_common()]
        print(f"    N={n}: {', '.join(parts)}")


# ============================================================
# 3. Monotonicity test
# ============================================================

def sec03_monotonicity(tensor):
    print("\n" + "=" * 72)
    print("SECTION 3: MONOTONICITY TEST — COMPLEXITY FLOW")
    print("=" * 72)

    all_ns = sorted(set(v["N"] for v in tensor.values()))
    if len(all_ns) < 2:
        print("  Not enough N values for monotonicity test.")
        return

    # Metric 1: Mean chi
    print(f"\n  Metric: Mean chi (collapse time)")
    prev_val = None
    monotone_chi = True
    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        val = statistics.mean(v["chi"] for v in cells)
        direction = ""
        if prev_val is not None:
            if val > prev_val:
                direction = "  INCREASE"
                if n > 4:  # N=4→8 increase is known (N=8 anomaly)
                    pass
            else:
                direction = ""
        print(f"    N={n:2d}: {val:8.2f}{direction}")
        prev_val = val

    # Metric 2: Number of distinct mechanisms
    print(f"\n  Metric: Distinct mechanisms")
    prev_nm = None
    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        nm = len(set(v["mechanism"] for v in cells))
        direction = ""
        if prev_nm is not None:
            if nm > prev_nm:
                direction = "  INCREASE"
        print(f"    N={n:2d}: {nm}{direction}")
        prev_nm = nm

    # Metric 3: Number of distinct sub-mechanisms
    print(f"\n  Metric: Distinct sub-mechanisms")
    prev_ns = None
    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        ns = len(set(v["subtype"] for v in cells if v["subtype"] != "unknown"))
        direction = ""
        if prev_ns is not None:
            if ns > prev_ns:
                direction = "  INCREASE"
        print(f"    N={n:2d}: {ns}{direction}")
        prev_ns = ns

    # Metric 4: Fraction of drift-rigid cells
    print(f"\n  Metric: Drift-rigid fraction")
    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        rigid = sum(1 for v in cells if v["drift_class"] == "rigid")
        print(f"    N={n:2d}: {100*rigid/len(cells):.0f}%")

    # Metric 5: Mean boundary curvature
    print(f"\n  Metric: Mean boundary curvature")
    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        kappas = [v["curvature"] for v in cells if v["curvature"] > 0]
        if kappas:
            print(f"    N={n:2d}: mean={statistics.mean(kappas):.1f}, "
                  f"max={max(kappas):.1f}, n={len(kappas)}")
        else:
            print(f"    N={n:2d}: no curvature (flat boundary)")

    # Metric 6: Mean torsion
    print(f"\n  Metric: Mean boundary torsion")
    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        taus = [v["torsion"] for v in cells if v["torsion"] > 0]
        if taus:
            print(f"    N={n:2d}: mean={statistics.mean(taus):.1f}, "
                  f"max={max(taus):.1f}")
        else:
            print(f"    N={n:2d}: no torsion")

    # Overall verdict
    print(f"\n  MONOTONICITY VERDICT:")
    # Check: does mean chi decrease monotonically from N=8 onward?
    chi_by_n = {}
    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        chi_by_n[n] = statistics.mean(v["chi"] for v in cells)

    ns_from_8 = [n for n in all_ns if n >= 8]
    chi_mono = all(chi_by_n[ns_from_8[i]] >= chi_by_n[ns_from_8[i+1]]
                   for i in range(len(ns_from_8) - 1)) if len(ns_from_8) > 1 else True
    print(f"    Chi monotone-decreasing for N>=8: {'YES' if chi_mono else 'NO'}")

    sub_by_n = {}
    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        sub_by_n[n] = len(set(v["subtype"] for v in cells if v["subtype"] != "unknown"))
    ns_from_4 = [n for n in all_ns if n >= 4]
    sub_mono = all(sub_by_n[ns_from_4[i]] >= sub_by_n[ns_from_4[i+1]]
                   for i in range(len(ns_from_4) - 1)) if len(ns_from_4) > 1 else True
    print(f"    Sub-mechanism count monotone-decreasing: {'YES' if sub_mono else 'NO'}")


# ============================================================
# 4. N=8 resonance persistence
# ============================================================

def sec04_n8_resonance(tensor):
    print("\n" + "=" * 72)
    print("SECTION 4: N=8 RESONANCE PERSISTENCE CHECK")
    print("=" * 72)

    all_ns = sorted(set(v["N"] for v in tensor.values()))

    # OL-wandering fraction by N
    print(f"\n  OL-wandering fraction by N:")
    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        olw = sum(1 for v in cells if v["subtype"] == "OL-wandering")
        pct = 100 * olw / len(cells) if cells else 0
        marker = "  <-- RESONANCE" if pct > 10 else ""
        print(f"    N={n:2d}: {olw:4d} / {len(cells):4d} = {pct:5.1f}%{marker}")

    # Chi statistics: is N=8 still the slow plateau?
    print(f"\n  Mean chi by N (canonical angles, drift=0):")
    for n in all_ns:
        cells = [v for v in tensor.values()
                 if v["N"] == n and v["angle"] in ANGLES_CANONICAL
                 and v["drift"] == 0.0]
        if cells:
            chis = [v["chi"] for v in cells]
            print(f"    N={n:2d}: mean={statistics.mean(chis):8.2f}, "
                  f"max={max(chis):8.2f}")

    # Check: any new resonances at N=24, 28, 32?
    print(f"\n  New resonance check (N=24, 28, 32):")
    for n in NS_HIGH:
        cells = [v for v in tensor.values() if v["N"] == n]
        if not cells:
            print(f"    N={n}: no data")
            continue
        olw = sum(1 for v in cells if v["subtype"] == "OL-wandering")
        mean_chi = statistics.mean(v["chi"] for v in cells)
        n_mechs = len(set(v["mechanism"] for v in cells))
        pct_olw = 100 * olw / len(cells)
        anomalous = pct_olw > 5 or mean_chi > 1.0
        status = "ANOMALOUS" if anomalous else "FLAT"
        print(f"    N={n}: OLw={pct_olw:.1f}%, chi={mean_chi:.3f}, "
              f"mechs={n_mechs} -> {status}")


# ============================================================
# 5. DECAY confinement
# ============================================================

def sec05_decay_confinement(tensor):
    print("\n" + "=" * 72)
    print("SECTION 5: DECAY CONFINEMENT CHECK")
    print("=" * 72)

    decay_by_n = defaultdict(int)
    total_by_n = defaultdict(int)
    for v in tensor.values():
        total_by_n[v["N"]] += 1
        if v["mechanism"] == "DECAY":
            decay_by_n[v["N"]] += 1

    all_ns = sorted(total_by_n.keys())
    print(f"\n  DECAY cells by N:")
    for n in all_ns:
        dc = decay_by_n[n]
        tc = total_by_n[n]
        print(f"    N={n:2d}: {dc:4d} / {tc:4d} = "
              f"{100*dc/tc:.1f}%")

    decay_ns = set(n for n, c in decay_by_n.items() if c > 0)
    if decay_ns == {4}:
        print(f"\n  VERDICT: DECAY remains CONFINED to N=4. Law VII holds.")
    elif len(decay_ns) == 0:
        print(f"\n  VERDICT: No DECAY found anywhere (unexpected if N=4 included).")
    else:
        print(f"\n  VERDICT: DECAY found at N={sorted(decay_ns)}. "
              f"Law VII VIOLATED if N != {{4}}.")

    # DECAY at fine angles / fine drifts
    if decay_by_n[4] > 0:
        decay_cells = [v for v in tensor.values()
                       if v["mechanism"] == "DECAY" and v["N"] == 4]
        fine_angle_decay = [v for v in decay_cells if v["angle"] in ANGLES_FINE]
        fine_drift_decay = [v for v in decay_cells if v["drift"] in DRIFTS_FINE]
        print(f"\n  N=4 DECAY at fine angles: {len(fine_angle_decay)} cells")
        if fine_angle_decay:
            angs = sorted(set(v["angle"] for v in fine_angle_decay))
            print(f"    Angles: {angs}")
        print(f"  N=4 DECAY at fine drifts: {len(fine_drift_decay)} cells")
        if fine_drift_decay:
            drifts = sorted(set(v["drift"] for v in fine_drift_decay))
            print(f"    Drifts: {drifts}")


# ============================================================
# 6. Fine-angle structure near tangent
# ============================================================

def sec06_fine_angle(tensor, onset):
    print("\n" + "=" * 72)
    print("SECTION 6: FINE-ANGLE STRUCTURE NEAR TANGENT")
    print("=" * 72)

    all_ns = sorted(set(v["N"] for v in tensor.values()))
    tangent_angles = sorted(a for a in ANGLES_ALL if -1.5 <= a <= 1.5)

    for n in all_ns:
        print(f"\n  N={n}, drift=0, tangent-region angles {tangent_angles}:")
        for r in RADII:
            row = []
            for a in tangent_angles:
                key = (n, a, r, 0.0)
                if key in tensor:
                    v = tensor[key]
                    row.append(f"{SHORT.get(v['mechanism'], v['mechanism'][:2]):>2s}")
                else:
                    row.append(" .")
            print(f"    r={r:.2f}: {' '.join(row)}")

    # Check for new folds: onset non-monotonicity in the fine-angle region
    print(f"\n  Fine-angle fold detection (onset reversal in [-1, +1] range):")
    for n in all_ns:
        fold_count = 0
        for d in DRIFTS_ALL:
            for r_idx in range(1, len(RADII) - 1):
                r = RADII[r_idx]
                # Check if onset is non-monotonic vs radius
                o_prev = onset.get((n, RADII[r_idx - 1], d))
                o_curr = onset.get((n, r, d))
                o_next = onset.get((n, RADII[r_idx + 1], d))
                if o_prev is not None and o_curr is not None and o_next is not None:
                    d1 = o_curr - o_prev
                    d2 = o_next - o_curr
                    if d1 * d2 < 0:
                        fold_count += 1
        print(f"    N={n:2d}: {fold_count} folds detected")

    # New feature: does the IC/non-IC boundary now resolve to sub-degree precision?
    print(f"\n  IC boundary location at fine angles (drift=0, r=0.50):")
    for n in all_ns:
        ic_angs = []
        non_ic_angs = []
        for a in tangent_angles:
            key = (n, a, 0.50, 0.0)
            if key in tensor:
                if tensor[key]["mechanism"] == "inward-collapse":
                    ic_angs.append(a)
                else:
                    non_ic_angs.append(a)
        if ic_angs and non_ic_angs:
            boundary_lo = max(ic_angs)
            boundary_hi = min(a for a in non_ic_angs if a > boundary_lo) if [a for a in non_ic_angs if a > boundary_lo] else "?"
            print(f"    N={n:2d}: IC up to {boundary_lo:+.2f}°, "
                  f"non-IC from {boundary_hi}°")
        elif ic_angs:
            print(f"    N={n:2d}: all IC in range")
        elif non_ic_angs:
            print(f"    N={n:2d}: no IC in range")
        else:
            print(f"    N={n:2d}: no data at r=0.50")


# ============================================================
# 7. Fine-drift structure
# ============================================================

def sec07_fine_drift(tensor, onset):
    print("\n" + "=" * 72)
    print("SECTION 7: FINE-DRIFT STRUCTURE")
    print("=" * 72)

    all_ns = sorted(set(v["N"] for v in tensor.values()))

    # Onset angle as a function of drift (all drift values)
    print(f"\n  IC onset angle vs drift at r=0.50:")
    for n in all_ns:
        row = []
        for d in DRIFTS_ALL:
            o = onset.get((n, 0.50, d))
            if o is not None:
                row.append(f"{o:+6.2f}")
            else:
                row.append("  none")
        dstr = "  ".join(f"{d:.3f}" for d in DRIFTS_ALL)
        print(f"    N={n:2d} drift: {dstr}")
        print(f"    {'':4s} onset: {'  '.join(row)}")

    # Drift sensitivity at fine drifts
    print(f"\n  Drift class distribution (fine drifts only):")
    for n in all_ns:
        cells = [v for v in tensor.values()
                 if v["N"] == n and v["drift"] in DRIFTS_FINE]
        if not cells:
            continue
        dc = Counter(v["drift_class"] for v in cells)
        total = len(cells)
        parts = [f"{k}={100*ct/total:.0f}%" for k, ct in dc.most_common()]
        print(f"    N={n:2d}: {', '.join(parts)}")

    # Check: any new instabilities at fine drift?
    print(f"\n  Mechanism stability across fine drift steps (angle=+1, r=0.50):")
    for n in all_ns:
        mechs_at_drift = []
        for d in DRIFTS_ALL:
            key = (n, 1, 0.50, d)
            if key in tensor:
                mechs_at_drift.append(
                    f"d={d:.3f}:{SHORT.get(tensor[key]['mechanism'], '??')}")
        if mechs_at_drift:
            print(f"    N={n:2d}: {', '.join(mechs_at_drift)}")


# ============================================================
# 8. Temporal fold persistence at N>=8
# ============================================================

def sec08_temporal_folds(tensor):
    print("\n" + "=" * 72)
    print("SECTION 8: TEMPORAL FOLD PERSISTENCE AT N>=8")
    print("=" * 72)

    all_ns = sorted(set(v["N"] for v in tensor.values()))

    for n in all_ns:
        # Temporal folds: non-monotonicity of chi vs radius at positive angles, d=0
        chi_seq = []
        for r in RADII:
            chis_at_r = [v["chi"] for v in tensor.values()
                         if v["N"] == n and v["radius"] == r
                         and v["drift"] == 0.0 and v["angle"] > 0
                         and v["angle"] in ANGLES_CANONICAL]
            if chis_at_r:
                chi_seq.append(statistics.mean(chis_at_r))

        chi_folds = 0
        for i in range(1, len(chi_seq) - 1):
            d1 = chi_seq[i] - chi_seq[i - 1]
            d2 = chi_seq[i + 1] - chi_seq[i]
            if d1 * d2 < 0:
                chi_folds += 1

        # Temporal gradient distribution
        cells = [v for v in tensor.values() if v["N"] == n]
        tc = Counter(v["temporal_class"] for v in cells)
        total = len(cells)
        cliff_pct = 100 * tc.get("cliff", 0) / total
        steep_pct = 100 * tc.get("steep", 0) / total

        print(f"  N={n:2d}: chi_folds={chi_folds}, "
              f"cliff={cliff_pct:.1f}%, steep={steep_pct:.1f}%")

    # Fine-angle temporal structure: does sub-degree resolution reveal new structure?
    print(f"\n  Temporal structure at fine angles (drift=0, r=0.50):")
    for n in all_ns:
        tangent_angles = sorted(a for a in ANGLES_ALL if -1.5 <= a <= 1.5)
        chi_at_angle = []
        for a in tangent_angles:
            key = (n, a, 0.50, 0.0)
            if key in tensor:
                chi_at_angle.append((a, tensor[key]["chi"]))
        if chi_at_angle:
            parts = [f"{a:+.2f}°:{chi:.2f}" for a, chi in chi_at_angle]
            print(f"    N={n:2d}: {', '.join(parts)}")


# ============================================================
# 9. Unified manifold coherence at extended N
# ============================================================

def sec09_coherence(tensor, onset):
    print("\n" + "=" * 72)
    print("SECTION 9: UNIFIED MANIFOLD COHERENCE (EXTENDED N)")
    print("=" * 72)

    all_ns = sorted(set(v["N"] for v in tensor.values()))

    # Same 4 coherence tests as the unification script
    print(f"\n  Test 1: IC cells have boundary_dist <= 0")
    for n in all_ns:
        ic_outside = sum(1 for v in tensor.values()
                         if v["N"] == n and v["mechanism"] == "inward-collapse"
                         and not math.isnan(v["boundary_dist"])
                         and v["boundary_dist"] > 0)
        ic_total = sum(1 for v in tensor.values()
                       if v["N"] == n and v["mechanism"] == "inward-collapse"
                       and not math.isnan(v["boundary_dist"]))
        status = "HOLDS" if ic_outside == 0 else f"VIOLATED ({ic_outside} exceptions)"
        print(f"    N={n:2d}: {status}  (IC with boundary: {ic_total})")

    # Test 2: DECAY only at N=4
    decay_ns = set(v["N"] for v in tensor.values() if v["mechanism"] == "DECAY")
    print(f"\n  Test 2: DECAY exclusive to N=4: "
          f"{'HOLDS' if decay_ns <= {4} else 'VIOLATED'} "
          f"(found at N={sorted(decay_ns) if decay_ns else 'none'})")

    # Test 3: Cross-layer agreement per N
    print(f"\n  Test 3: Cross-layer agreement")
    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        consistent = 0
        for v in cells:
            is_ok = True
            if v["mechanism"] == "inward-collapse" and not math.isnan(v["boundary_dist"]):
                if v["boundary_dist"] > 0:
                    is_ok = False
            if v["chi"] < 0.1 and v["mechanism"] == "DECAY":
                is_ok = False
            if v["chi"] > 50 and v["mechanism"] == "inward-collapse":
                is_ok = False
            if is_ok:
                consistent += 1
        pct = 100 * consistent / len(cells) if cells else 0
        print(f"    N={n:2d}: {consistent}/{len(cells)} = {pct:.1f}%")

    # Test 4: Co-transition rate
    print(f"\n  Test 4: Mechanism/subtype co-transition rate")
    for n in all_ns:
        mech_trans = 0
        co_trans = 0
        for d in DRIFTS_ALL:
            for r in RADII:
                prev_m = None
                prev_s = None
                for a in sorted(set(v["angle"] for v in tensor.values() if v["N"] == n)):
                    key = (n, a, r, d)
                    if key not in tensor:
                        continue
                    cm = tensor[key]["mechanism"]
                    cs = tensor[key]["subtype"]
                    if prev_m is not None and cm != prev_m:
                        mech_trans += 1
                        if cs != prev_s:
                            co_trans += 1
                    prev_m = cm
                    prev_s = cs
        rate = 100 * co_trans / mech_trans if mech_trans > 0 else 0
        print(f"    N={n:2d}: {co_trans}/{mech_trans} = {rate:.0f}%")


# ============================================================
# 10. Taxonomy collapse
# ============================================================

def sec10_taxonomy(tensor):
    print("\n" + "=" * 72)
    print("SECTION 10: TAXONOMY COLLAPSE AT HIGH N")
    print("=" * 72)

    all_ns = sorted(set(v["N"] for v in tensor.values()))

    print(f"\n  Sub-mechanism inventory by N:")
    for n in all_ns:
        cells = [v for v in tensor.values()
                 if v["N"] == n and v["subtype"] != "unknown"]
        subs = Counter(v["subtype"] for v in cells)
        total = len(cells)
        print(f"\n  N={n:2d} ({total} cells, {len(subs)} sub-mechanisms):")
        for sub, ct in subs.most_common():
            print(f"    {sub:<16s}: {ct:4d} ({100*ct/total:5.1f}%)")

    # Does the taxonomy flatten further at N=24+?
    print(f"\n  Distinct sub-mechanism count progression:")
    for n in all_ns:
        cells = [v for v in tensor.values()
                 if v["N"] == n and v["subtype"] != "unknown"]
        n_subs = len(set(v["subtype"] for v in cells))
        bar = "#" * n_subs
        print(f"    N={n:2d}: {n_subs:2d} {bar}")

    # Check: do any NEW sub-mechanisms appear?
    canonical_subs = set()
    for v in tensor.values():
        if v["N"] in NS_CANONICAL and v["subtype"] != "unknown":
            canonical_subs.add(v["subtype"])

    high_n_subs = set()
    for v in tensor.values():
        if v["N"] in NS_HIGH and v["subtype"] != "unknown":
            high_n_subs.add(v["subtype"])

    new_subs = high_n_subs - canonical_subs
    lost_subs = canonical_subs - high_n_subs

    print(f"\n  Canonical sub-mechanisms: {sorted(canonical_subs)}")
    print(f"  High-N sub-mechanisms:   {sorted(high_n_subs)}")
    print(f"  New at high N:           {sorted(new_subs) if new_subs else 'none'}")
    print(f"  Lost at high N:          {sorted(lost_subs) if lost_subs else 'none'}")


# ============================================================
# 11. Summary: deviations from Law VIII
# ============================================================

def sec11_summary(tensor, onset):
    print("\n" + "=" * 72)
    print("SECTION 11: SUMMARY — DEVIATIONS FROM LAW VIII")
    print("=" * 72)

    all_ns = sorted(set(v["N"] for v in tensor.values()))
    deviations = []

    # Check 1: Single coherent manifold?
    for n in all_ns:
        ic_outside = sum(1 for v in tensor.values()
                         if v["N"] == n and v["mechanism"] == "inward-collapse"
                         and not math.isnan(v["boundary_dist"])
                         and v["boundary_dist"] > 0)
        if ic_outside > 0:
            deviations.append(f"IC outside boundary at N={n} ({ic_outside} cells)")

    # Check 2: DECAY confinement
    decay_ns = set(v["N"] for v in tensor.values() if v["mechanism"] == "DECAY")
    if not decay_ns <= {4}:
        deviations.append(f"DECAY found at N={sorted(decay_ns - {4})}")

    # Check 3: Monotone complexity flow (chi) for N>=8
    chi_by_n = {}
    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        chi_by_n[n] = statistics.mean(v["chi"] for v in cells)
    ns_post_8 = [n for n in all_ns if n >= 8]
    for i in range(len(ns_post_8) - 1):
        if chi_by_n[ns_post_8[i+1]] > chi_by_n[ns_post_8[i]]:
            deviations.append(
                f"Chi INCREASES from N={ns_post_8[i]} to N={ns_post_8[i+1]}")

    # Check 4: Sub-mechanism count monotone decreasing
    sub_by_n = {}
    for n in all_ns:
        cells = [v for v in tensor.values()
                 if v["N"] == n and v["subtype"] != "unknown"]
        sub_by_n[n] = len(set(v["subtype"] for v in cells))
    for i in range(len(all_ns) - 1):
        # Allow N=4→N=8 increase (known N=8 anomaly)
        if all_ns[i] >= 8 and sub_by_n.get(all_ns[i+1], 0) > sub_by_n.get(all_ns[i], 0):
            deviations.append(
                f"Sub-mech count INCREASES from N={all_ns[i]} to N={all_ns[i+1]}")

    # Check 5: N=8 resonance persists
    for n in NS_HIGH:
        cells = [v for v in tensor.values() if v["N"] == n]
        if cells:
            olw = sum(1 for v in cells if v["subtype"] == "OL-wandering")
            if 100 * olw / len(cells) > 10:
                deviations.append(f"N={n} shows N=8-like resonance ({100*olw/len(cells):.0f}% OLw)")

    # Check 6: New resonance at N=24, 28, or 32?
    for n in NS_HIGH:
        cells = [v for v in tensor.values() if v["N"] == n]
        if cells:
            mean_chi = statistics.mean(v["chi"] for v in cells)
            # If chi is anomalously high (> 2x the trend)
            prev_n = max(nn for nn in all_ns if nn < n)
            prev_cells = [v for v in tensor.values() if v["N"] == prev_n]
            if prev_cells:
                prev_chi = statistics.mean(v["chi"] for v in prev_cells)
                if mean_chi > 2 * prev_chi:
                    deviations.append(f"N={n} chi anomaly: {mean_chi:.3f} vs prev {prev_chi:.3f}")

    # Report
    print(f"\n  LAW VIII DEVIATION REPORT:")
    if not deviations:
        print(f"    NO DEVIATIONS DETECTED.")
        print(f"    The unified manifold remains single-sheeted, the complexity")
        print(f"    flow continues monotonically, DECAY remains confined to N=4,")
        print(f"    and no new resonances appear at N=24-32.")
        print(f"\n    Law VIII is CONFIRMED by the consolidation sweep.")
    else:
        print(f"    {len(deviations)} deviation(s) detected:")
        for i, d in enumerate(deviations, 1):
            print(f"      {i}. {d}")
        print(f"\n    Law VIII requires AMENDMENT to account for these deviations.")

    # Compact summary table
    print(f"\n  CONSOLIDATION SUMMARY TABLE:")
    print(f"  {'N':>4s}  {'cells':>5s}  {'mechs':>5s}  {'subs':>5s}  "
          f"{'chi':>8s}  {'OLw%':>5s}  {'rigid%':>6s}  {'kappa':>6s}  {'agree%':>7s}")
    for n in all_ns:
        cells = [v for v in tensor.values() if v["N"] == n]
        nc = len(cells)
        nm = len(set(v["mechanism"] for v in cells))
        ns = len(set(v["subtype"] for v in cells if v["subtype"] != "unknown"))
        chi = statistics.mean(v["chi"] for v in cells)
        olw = 100 * sum(1 for v in cells if v["subtype"] == "OL-wandering") / nc
        rigid = 100 * sum(1 for v in cells if v["drift_class"] == "rigid") / nc
        kappas = [v["curvature"] for v in cells if v["curvature"] > 0]
        kappa_mean = statistics.mean(kappas) if kappas else 0

        consistent = 0
        for v in cells:
            is_ok = True
            if v["mechanism"] == "inward-collapse" and not math.isnan(v["boundary_dist"]):
                if v["boundary_dist"] > 0:
                    is_ok = False
            if v["chi"] > 50 and v["mechanism"] == "inward-collapse":
                is_ok = False
            if is_ok:
                consistent += 1
        agree = 100 * consistent / nc

        print(f"  {n:4d}  {nc:5d}  {nm:5d}  {ns:5d}  "
              f"{chi:8.3f}  {olw:5.1f}  {rigid:6.0f}  {kappa_mean:6.1f}  {agree:7.1f}")


# ============================================================
# Main
# ============================================================

def main():
    print("Post-Law VIII Consolidation Analysis")
    print("=" * 72)

    print("Loading ALL sweep data (canonical + extended)...")
    sweep_cube = load_all_sweeps()
    print(f"  {len(sweep_cube)} sweep records")

    print("Loading ALL atlas data (canonical + extended)...")
    atlas_cube = load_all_atlas()
    print(f"  {len(atlas_cube)} atlas records")

    print("Assembling unified tensor...")
    tensor, onset = assemble_tensor(sweep_cube, atlas_cube)
    print(f"  {len(tensor)} unified cells")

    all_ns = sorted(set(v["N"] for v in tensor.values()))
    print(f"  N values: {all_ns}")
    print()

    sec01_inventory(tensor)
    sec02_complexity_flow(tensor)
    sec03_monotonicity(tensor)
    sec04_n8_resonance(tensor)
    sec05_decay_confinement(tensor)
    sec06_fine_angle(tensor, onset)
    sec07_fine_drift(tensor, onset)
    sec08_temporal_folds(tensor)
    sec09_coherence(tensor, onset)
    sec10_taxonomy(tensor)
    sec11_summary(tensor, onset)


if __name__ == "__main__":
    main()
