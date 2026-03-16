#!/usr/bin/env python
"""
check_high_n_asymptotic.py
============================
High-N Asymptotic Probe: determines whether ED-Arch exhibits a true
asymptotic law as N -> infinity.

Loads ALL data (canonical N=4..32 + high-N N=40..128), then runs 11 sections:
  1.  Inventory & coverage
  2.  Phase boundary convergence (does onset angle approach a limit?)
  3.  Collapse time scaling (chi ~ 1/N^p or chi ~ exp(-kN)?)
  4.  Taxonomy stabilization (does the sub-mechanism set freeze?)
  5.  Drift sensitivity convergence (does sensitivity -> 0?)
  6.  Curvature and torsion convergence (do they -> 0?)
  7.  Planarity test (is the manifold exactly flat at high N?)
  8.  Resonance scan (any new anomalies at N=40..128?)
  9.  IC-direct dominance curve
 10.  Scaling law fits (power law vs exponential)
 11.  Summary: asymptotic law statement
"""

import json
import os
import glob
import math
import statistics
from collections import defaultdict, Counter

SWEEP_DIR = os.path.dirname(__file__)

# All N values we want to analyze (canonical + consolidation + high-N)
NS_CANONICAL     = [4, 8, 12, 20]
NS_CONSOLIDATION = [24, 28, 32]
NS_HIGH          = [40, 48, 56, 64, 80, 96, 128]
NS_ALL           = NS_CANONICAL + NS_CONSOLIDATION + NS_HIGH

ANGLES = [-10, -5, -2, -1, 0, 1, 2, 5, 10]
RADII  = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 0.90, 0.95]

# For high-N we only generated 3 drifts; canonical had 5
DRIFTS_CANONICAL = [0.00, 0.01, 0.02, 0.05, 0.10]
DRIFTS_HIGH      = [0.00, 0.01, 0.02]
# For cross-N comparison, use the intersection
DRIFTS_COMMON    = [0.00, 0.01, 0.02]

MECHS = ["inward-collapse", "outward-PBC", "DECAY", "PBC-corner", "other-late"]
SHORT = {"inward-collapse": "IC", "outward-PBC": "OP", "DECAY": "DE",
         "PBC-corner": "CR", "other-late": "OL"}

SORTED_ANGLES = sorted(ANGLES)


# ============================================================
# Data loading
# ============================================================

def load_all_sweeps():
    cube = {}
    pattern = os.path.join(SWEEP_DIR, "n*_angle_*_r*_d*_sweep.json")
    for path in sorted(glob.glob(pattern)):
        with open(path) as f:
            data = json.load(f)
        key = (data["N"], data["angle"], data["radius"], data["drift"])
        cube[key] = data
    return cube


def load_all_atlas():
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
# Sub-mechanism classification
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

def compute_ic_onset(sweep_cube, ns_list):
    onset = {}
    for n in ns_list:
        for r in RADII:
            for d in DRIFTS_COMMON:
                ic_angs = [ang for ang in SORTED_ANGLES
                           if (n, ang, r, d) in sweep_cube
                           and sweep_cube[(n, ang, r, d)]["mechanism"] == "inward-collapse"]
                onset[(n, r, d)] = max(ic_angs) if ic_angs else None
    return onset


def local_curvature(onset, n, d, r_idx):
    if r_idx < 1 or r_idx >= len(RADII) - 1:
        return 0.0
    r0, r1, r2 = RADII[r_idx - 1], RADII[r_idx], RADII[r_idx + 1]
    a0 = onset.get((n, r0, d))
    a1 = onset.get((n, r1, d))
    a2 = onset.get((n, r2, d))
    if a0 is None or a1 is None or a2 is None:
        return 0.0
    if a0 > max(ANGLES) or a1 > max(ANGLES) or a2 > max(ANGLES):
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
    if r_idx < 1 or r_idx >= len(RADII) - 1:
        return 0.0
    r1 = RADII[r_idx]
    dvals = sorted(d for d in DRIFTS_COMMON if onset.get((n, r1, d)) is not None)
    if len(dvals) < 3:
        return 0.0
    torsion_sum = 0.0
    count = 0
    for i in range(1, len(dvals) - 1):
        d0, d1, d2 = dvals[i-1], dvals[i], dvals[i+1]
        a0 = onset.get((n, r1, d0))
        a1 = onset.get((n, r1, d1))
        a2 = onset.get((n, r1, d2))
        if a0 is None or a1 is None or a2 is None:
            continue
        if a0 > max(ANGLES) or a1 > max(ANGLES) or a2 > max(ANGLES):
            continue
        dd1 = float(a1 - a0) / (d1 - d0) if d1 != d0 else 0
        dd2 = float(a2 - a1) / (d2 - d1) if d2 != d1 else 0
        torsion_sum += abs(dd2 - dd1)
        count += 1
    return torsion_sum / count if count > 0 else 0.0


# ============================================================
# 1. Inventory
# ============================================================

def sec01_inventory(sweep_cube, atlas_cube):
    print("=" * 72)
    print("SECTION 1: INVENTORY & COVERAGE")
    print("=" * 72)

    all_ns = sorted(set(k[0] for k in sweep_cube))
    print(f"  Total sweep records: {len(sweep_cube)}")
    print(f"  Total atlas records: {len(atlas_cube)}")
    print(f"  N values present: {all_ns}")

    for n in all_ns:
        cells = {k: v for k, v in sweep_cube.items() if k[0] == n}
        drifts = sorted(set(k[3] for k in cells))
        mechs = Counter(v["mechanism"] for v in cells.values())
        print(f"\n  N={n:3d}: {len(cells)} cells, drifts={drifts}")
        parts = [f"{SHORT.get(m,m)}={c}" for m, c in mechs.most_common()]
        print(f"    Mechanisms: {', '.join(parts)}")

    # How many are high-N?
    high_n_count = sum(1 for k in sweep_cube if k[0] in NS_HIGH)
    print(f"\n  High-N (N>=40) cells: {high_n_count}")


# ============================================================
# 2. Phase boundary convergence
# ============================================================

def sec02_boundary_convergence(sweep_cube, onset, ns_present):
    print("\n" + "=" * 72)
    print("SECTION 2: PHASE BOUNDARY CONVERGENCE")
    print("=" * 72)

    # For each N, compute mean onset angle across all (radius, drift=0)
    print(f"\n  Mean IC onset angle by N (drift=0):")
    print(f"  {'N':>5s}  {'mean_onset':>10s}  {'min':>6s}  {'max':>6s}  "
          f"{'distinct':>8s}  {'IC_frac':>8s}")

    onset_by_n = {}
    for n in ns_present:
        onsets = []
        for r in RADII:
            o = onset.get((n, r, 0.0))
            if o is not None:
                onsets.append(o)
        if not onsets:
            print(f"  {n:5d}  {'no IC':>10s}")
            continue
        onset_by_n[n] = onsets
        distinct = len(set(onsets))
        # IC fraction at this N
        total_cells = sum(1 for k, v in sweep_cube.items()
                          if k[0] == n and k[3] == 0.0)
        ic_cells = sum(1 for k, v in sweep_cube.items()
                       if k[0] == n and k[3] == 0.0
                       and v["mechanism"] == "inward-collapse")
        ic_frac = ic_cells / total_cells if total_cells > 0 else 0
        print(f"  {n:5d}  {statistics.mean(onsets):10.2f}  {min(onsets):6.1f}  "
              f"{max(onsets):6.1f}  {distinct:8d}  {ic_frac:8.2f}")

    # Check: does onset converge to a single value?
    print(f"\n  Onset angle convergence check:")
    for n in ns_present:
        onsets = onset_by_n.get(n, [])
        if not onsets:
            continue
        distinct = len(set(onsets))
        is_flat = distinct == 1
        marker = "  FLAT" if is_flat else ""
        print(f"    N={n:3d}: {distinct} distinct onset angle(s): "
              f"{sorted(set(onsets))}{marker}")

    # Does the limiting onset angle stabilize?
    print(f"\n  Onset angle at r=0.50, drift=0 (boundary-tracking radius):")
    for n in ns_present:
        o = onset.get((n, 0.50, 0.0))
        val = f"{o:+.1f}" if o is not None else "no-IC"
        print(f"    N={n:3d}: {val}")


# ============================================================
# 3. Collapse time scaling
# ============================================================

def sec03_chi_scaling(sweep_cube, ns_present):
    print("\n" + "=" * 72)
    print("SECTION 3: COLLAPSE TIME SCALING (chi vs N)")
    print("=" * 72)

    # Mean chi by N (drift=0, all angles, all radii)
    chi_by_n = {}
    print(f"\n  Mean chi by N (drift=0, canonical angles):")
    print(f"  {'N':>5s}  {'mean_chi':>10s}  {'med_chi':>10s}  {'max_chi':>10s}  "
          f"{'min_chi':>10s}  {'frac_floor':>10s}")

    for n in ns_present:
        cells = [v for k, v in sweep_cube.items()
                 if k[0] == n and k[3] == 0.0 and k[1] in ANGLES]
        if not cells:
            continue
        chis = [v["chi_emp"] for v in cells]
        chi_by_n[n] = statistics.mean(chis)
        frac_floor = sum(1 for c in chis if c <= 0.01) / len(chis)
        print(f"  {n:5d}  {statistics.mean(chis):10.4f}  "
              f"{statistics.median(chis):10.4f}  {max(chis):10.4f}  "
              f"{min(chis):10.4f}  {frac_floor:10.2f}")

    # Scaling law: log-log (power law) and semi-log (exponential)
    # Use N >= 8 to avoid the N=4 outlier
    ns_fit = [n for n in ns_present if n >= 8 and n in chi_by_n and chi_by_n[n] > 0.01]
    if len(ns_fit) >= 3:
        import numpy as np

        log_n = [math.log(n) for n in ns_fit]
        log_chi = [math.log(chi_by_n[n]) for n in ns_fit]
        lin_n = [float(n) for n in ns_fit]

        # Power law fit: log(chi) = a - p*log(N) => chi ~ N^(-p)
        n_pts = len(ns_fit)
        sum_x = sum(log_n)
        sum_y = sum(log_chi)
        sum_xy = sum(x * y for x, y in zip(log_n, log_chi))
        sum_x2 = sum(x**2 for x in log_n)
        denom = n_pts * sum_x2 - sum_x**2
        if abs(denom) > 1e-15:
            slope_pow = (n_pts * sum_xy - sum_x * sum_y) / denom
            intercept_pow = (sum_y - slope_pow * sum_x) / n_pts
            # Residuals
            ss_res = sum((y - (intercept_pow + slope_pow * x))**2
                         for x, y in zip(log_n, log_chi))
            ss_tot = sum((y - sum_y / n_pts)**2 for y in log_chi)
            r2_pow = 1 - ss_res / ss_tot if ss_tot > 0 else 0
        else:
            slope_pow, intercept_pow, r2_pow = 0, 0, 0

        # Exponential fit: log(chi) = a - k*N => chi ~ exp(-kN)
        sum_x_e = sum(lin_n)
        sum_y_e = sum(log_chi)
        sum_xy_e = sum(x * y for x, y in zip(lin_n, log_chi))
        sum_x2_e = sum(x**2 for x in lin_n)
        denom_e = n_pts * sum_x2_e - sum_x_e**2
        if abs(denom_e) > 1e-15:
            slope_exp = (n_pts * sum_xy_e - sum_x_e * sum_y_e) / denom_e
            intercept_exp = (sum_y_e - slope_exp * sum_x_e) / n_pts
            ss_res_e = sum((y - (intercept_exp + slope_exp * x))**2
                           for x, y in zip(lin_n, log_chi))
            r2_exp = 1 - ss_res_e / ss_tot if ss_tot > 0 else 0
        else:
            slope_exp, intercept_exp, r2_exp = 0, 0, 0

        print(f"\n  Scaling law fits (N >= 8, excluding chi = floor):")
        print(f"    Power law:  chi ~ N^({slope_pow:.2f}),  R² = {r2_pow:.4f}")
        print(f"      => exponent p = {-slope_pow:.2f}")
        print(f"    Exponential: chi ~ exp({slope_exp:.4f} * N),  R² = {r2_exp:.4f}")
        print(f"      => decay rate k = {-slope_exp:.4f}")

        if r2_pow > r2_exp:
            print(f"\n    BEST FIT: Power law (R² = {r2_pow:.4f} > {r2_exp:.4f})")
            print(f"    chi ~ N^({slope_pow:.2f})")
        else:
            print(f"\n    BEST FIT: Exponential (R² = {r2_exp:.4f} > {r2_pow:.4f})")
            print(f"    chi ~ exp({slope_exp:.4f} * N)")

        # Check for floor: does chi hit a hard lower bound?
        print(f"\n  Floor detection:")
        for n in ns_present:
            if n not in chi_by_n:
                continue
            chi = chi_by_n[n]
            at_floor = chi <= 0.011
            print(f"    N={n:3d}: chi={chi:.4f} "
                  f"{'<-- AT FLOOR' if at_floor else ''}")
    else:
        print(f"\n  Insufficient data for scaling fit (need N>=8 with chi > floor)")


# ============================================================
# 4. Taxonomy stabilization
# ============================================================

def sec04_taxonomy(atlas_cube, ns_present):
    print("\n" + "=" * 72)
    print("SECTION 4: TAXONOMY STABILIZATION")
    print("=" * 72)

    print(f"\n  Sub-mechanism inventory by N:")
    sub_sets = {}
    for n in ns_present:
        cells = {k: v for k, v in atlas_cube.items()
                 if k[0] == n and k[3] in DRIFTS_COMMON}
        if not cells:
            continue
        subs = Counter(classify_sub(v) for v in cells.values())
        total = len(cells)
        sub_sets[n] = set(subs.keys())
        print(f"\n  N={n:3d} ({total} cells, {len(subs)} sub-mechanisms):")
        for sub, ct in subs.most_common():
            print(f"    {sub:<16s}: {ct:4d} ({100*ct/total:5.1f}%)")

    # Track the set progression
    print(f"\n  Sub-mechanism count progression:")
    for n in ns_present:
        if n not in sub_sets:
            continue
        count = len(sub_sets[n])
        bar = "#" * count
        print(f"    N={n:3d}: {count:2d} {bar}")

    # When does the set stabilize?
    print(f"\n  Set stability check (is set[N] == set[N-1]?):")
    ns_sorted = sorted(sub_sets.keys())
    stable_from = None
    for i in range(1, len(ns_sorted)):
        n_prev = ns_sorted[i - 1]
        n_curr = ns_sorted[i]
        same = sub_sets[n_curr] == sub_sets[n_prev]
        gained = sub_sets[n_curr] - sub_sets[n_prev]
        lost = sub_sets[n_prev] - sub_sets[n_curr]
        status = "SAME" if same else f"CHANGED (+{gained or '{}'} -{lost or '{}'})"
        print(f"    N={n_prev:3d} -> N={n_curr:3d}: {status}")
        if same and stable_from is None:
            stable_from = n_prev

    # What is the limiting set?
    if ns_sorted:
        last_n = ns_sorted[-1]
        print(f"\n  Taxonomy at N={last_n}: {sorted(sub_sets.get(last_n, set()))}")

    # Does the IC-direct fraction approach 100%?
    print(f"\n  IC-direct dominance curve:")
    for n in ns_present:
        cells = {k: v for k, v in atlas_cube.items()
                 if k[0] == n and k[3] in DRIFTS_COMMON}
        if not cells:
            continue
        ic_direct = sum(1 for v in cells.values()
                        if classify_sub(v) == "IC-direct")
        print(f"    N={n:3d}: IC-direct = {ic_direct}/{len(cells)} "
              f"= {100*ic_direct/len(cells):.1f}%")


# ============================================================
# 5. Drift sensitivity convergence
# ============================================================

def sec05_drift_sensitivity(sweep_cube, onset, ns_present):
    print("\n" + "=" * 72)
    print("SECTION 5: DRIFT SENSITIVITY CONVERGENCE")
    print("=" * 72)

    print(f"\n  Onset angle shift: drift=0 vs drift=0.02")
    print(f"  {'N':>5s}  {'rigid':>6s}  {'stable':>6s}  {'moderate':>8s}  "
          f"{'fragile':>7s}  {'no-bdy':>7s}")

    for n in ns_present:
        rigid = stable = moderate = fragile = no_bdy = 0
        for r in RADII:
            o0 = onset.get((n, r, 0.0))
            o2 = onset.get((n, r, 0.02))
            if o0 is None or o2 is None:
                no_bdy += 1
                continue
            shift = abs(o2 - o0)
            if shift == 0:    rigid += 1
            elif shift <= 2:  stable += 1
            elif shift <= 5:  moderate += 1
            else:             fragile += 1
        total = rigid + stable + moderate + fragile + no_bdy
        print(f"  {n:5d}  {rigid:6d}  {stable:6d}  {moderate:8d}  "
              f"{fragile:7d}  {no_bdy:7d}")

    # Chi sensitivity to drift
    print(f"\n  Mean chi at drift=0 vs drift=0.02:")
    for n in ns_present:
        chi_d0 = [v["chi_emp"] for k, v in sweep_cube.items()
                  if k[0] == n and k[3] == 0.0 and k[1] in ANGLES]
        chi_d2 = [v["chi_emp"] for k, v in sweep_cube.items()
                  if k[0] == n and k[3] == 0.02 and k[1] in ANGLES]
        if chi_d0 and chi_d2:
            m0 = statistics.mean(chi_d0)
            m2 = statistics.mean(chi_d2)
            ratio = m2 / m0 if m0 > 0 else float('inf')
            print(f"    N={n:3d}: d=0: {m0:.4f}, d=0.02: {m2:.4f}, "
                  f"ratio={ratio:.3f}")


# ============================================================
# 6. Curvature and torsion convergence
# ============================================================

def sec06_curvature_torsion(onset, ns_present):
    print("\n" + "=" * 72)
    print("SECTION 6: CURVATURE AND TORSION CONVERGENCE")
    print("=" * 72)

    print(f"\n  {'N':>5s}  {'mean_kappa':>10s}  {'max_kappa':>10s}  "
          f"{'n_curved':>8s}  {'mean_tau':>10s}  {'max_tau':>10s}")

    for n in ns_present:
        kappas = []
        taus = []
        for d in DRIFTS_COMMON:
            for r_idx in range(len(RADII)):
                k = local_curvature(onset, n, d, r_idx)
                t = local_torsion(onset, n, r_idx)
                if k > 0:
                    kappas.append(k)
                if t > 0:
                    taus.append(t)

        mk = statistics.mean(kappas) if kappas else 0
        xk = max(kappas) if kappas else 0
        mt = statistics.mean(taus) if taus else 0
        xt = max(taus) if taus else 0
        print(f"  {n:5d}  {mk:10.1f}  {xk:10.1f}  {len(kappas):8d}  "
              f"{mt:10.1f}  {xt:10.1f}")

    # Does curvature -> 0?
    print(f"\n  Curvature convergence check:")
    for n in ns_present:
        kappas = []
        for d in DRIFTS_COMMON:
            for r_idx in range(len(RADII)):
                k = local_curvature(onset, n, d, r_idx)
                if k > 0:
                    kappas.append(k)
        is_zero = len(kappas) == 0
        print(f"    N={n:3d}: {'ZERO (flat boundary)' if is_zero else f'{len(kappas)} curved segments, mean={statistics.mean(kappas):.1f}'}")


# ============================================================
# 7. Planarity test
# ============================================================

def sec07_planarity(onset, ns_present):
    print("\n" + "=" * 72)
    print("SECTION 7: PLANARITY TEST")
    print("=" * 72)

    print(f"\n  A manifold slice is 'planar' if the IC onset angle is the same")
    print(f"  at every radius for a given (N, drift).")

    for n in ns_present:
        planar_count = 0
        total_slices = 0
        for d in DRIFTS_COMMON:
            onsets_at = set()
            has_data = False
            for r in RADII:
                o = onset.get((n, r, d))
                if o is not None:
                    onsets_at.add(o)
                    has_data = True
            if has_data:
                total_slices += 1
                if len(onsets_at) <= 1:
                    planar_count += 1

        is_planar = planar_count == total_slices if total_slices > 0 else False
        print(f"    N={n:3d}: {planar_count}/{total_slices} drift slices planar"
              f"  {'=> EXACTLY PLANAR' if is_planar else ''}")

    # Fold count per N
    print(f"\n  Boundary folds (onset non-monotonicity vs radius, drift=0):")
    for n in ns_present:
        onset_seq = []
        for r in RADII:
            o = onset.get((n, r, 0.0))
            if o is not None and o <= max(ANGLES):
                onset_seq.append(o)
        folds = 0
        for i in range(1, len(onset_seq) - 1):
            d1 = onset_seq[i] - onset_seq[i - 1]
            d2 = onset_seq[i + 1] - onset_seq[i]
            if d1 * d2 < 0:
                folds += 1
        print(f"    N={n:3d}: {folds} folds")


# ============================================================
# 8. Resonance scan
# ============================================================

def sec08_resonance_scan(atlas_cube, sweep_cube, ns_present):
    print("\n" + "=" * 72)
    print("SECTION 8: RESONANCE SCAN (NEW ANOMALIES AT HIGH N?)")
    print("=" * 72)

    # OL-wandering fraction
    print(f"\n  OL-wandering fraction by N:")
    for n in ns_present:
        cells = {k: v for k, v in atlas_cube.items()
                 if k[0] == n and k[3] in DRIFTS_COMMON}
        if not cells:
            continue
        olw = sum(1 for v in cells.values()
                  if classify_sub(v) == "OL-wandering")
        pct = 100 * olw / len(cells)
        marker = "  <-- RESONANCE" if pct > 10 else ""
        marker = "  <-- TRACE" if 2 < pct <= 10 else marker
        print(f"    N={n:3d}: {olw:4d}/{len(cells):4d} = {pct:5.1f}%{marker}")

    # Mean chi anomaly: compare each N to the monotone trend
    print(f"\n  Mean chi by N (looking for non-monotone bumps):")
    prev_chi = None
    for n in ns_present:
        cells = [v for k, v in sweep_cube.items()
                 if k[0] == n and k[3] == 0.0 and k[1] in ANGLES]
        if not cells:
            continue
        chi = statistics.mean(v["chi_emp"] for v in cells)
        bump = ""
        if prev_chi is not None and chi > prev_chi and n > 8:
            bump = "  <-- BUMP"
        print(f"    N={n:3d}: chi = {chi:.4f}{bump}")
        prev_chi = chi

    # Mechanism diversity at high N
    print(f"\n  Distinct mechanisms at high N:")
    for n in NS_HIGH:
        cells = {k: v for k, v in sweep_cube.items()
                 if k[0] == n and k[3] in DRIFTS_COMMON}
        if not cells:
            continue
        mechs = set(v["mechanism"] for v in cells.values())
        print(f"    N={n:3d}: {len(mechs)} mechanisms: "
              f"{sorted(SHORT.get(m, m) for m in mechs)}")


# ============================================================
# 9. IC-direct dominance curve
# ============================================================

def sec09_ic_direct(atlas_cube, ns_present):
    print("\n" + "=" * 72)
    print("SECTION 9: IC-DIRECT DOMINANCE CURVE")
    print("=" * 72)

    print(f"\n  IC-direct fraction by N (does it approach 100%?):")
    print(f"  {'N':>5s}  {'IC-dir':>8s}  {'IC-bnc':>8s}  {'OP-s':>8s}  "
          f"{'OL-s':>8s}  {'other':>8s}")

    for n in ns_present:
        cells = {k: v for k, v in atlas_cube.items()
                 if k[0] == n and k[3] in DRIFTS_COMMON}
        if not cells:
            continue
        subs = Counter(classify_sub(v) for v in cells.values())
        total = len(cells)
        ic_dir = 100 * subs.get("IC-direct", 0) / total
        ic_bnc = 100 * subs.get("IC-bounced", 0) / total
        op_s = 100 * subs.get("OPBC-single", 0) / total
        ol_s = 100 * subs.get("OL-stalled", 0) / total
        other = 100 - ic_dir - ic_bnc - op_s - ol_s
        print(f"  {n:5d}  {ic_dir:7.1f}%  {ic_bnc:7.1f}%  {op_s:7.1f}%  "
              f"{ol_s:7.1f}%  {other:7.1f}%")

    # IC-direct as fraction of all IC
    print(f"\n  IC-direct as fraction of all IC:")
    for n in ns_present:
        cells = {k: v for k, v in atlas_cube.items()
                 if k[0] == n and k[3] in DRIFTS_COMMON}
        if not cells:
            continue
        ic_all = [v for v in cells.values()
                  if v["mechanism"] == "inward-collapse"]
        if not ic_all:
            continue
        ic_dir = sum(1 for v in ic_all if classify_sub(v) == "IC-direct")
        print(f"    N={n:3d}: {ic_dir}/{len(ic_all)} = "
              f"{100*ic_dir/len(ic_all):.1f}%")


# ============================================================
# 10. Scaling law fits (detailed)
# ============================================================

def sec10_scaling_fits(sweep_cube, atlas_cube, ns_present):
    print("\n" + "=" * 72)
    print("SECTION 10: DETAILED SCALING LAW FITS")
    print("=" * 72)

    # Fit multiple quantities vs N
    metrics = {}
    for n in ns_present:
        cells_s = [v for k, v in sweep_cube.items()
                   if k[0] == n and k[3] == 0.0 and k[1] in ANGLES]
        cells_a = {k: v for k, v in atlas_cube.items()
                   if k[0] == n and k[3] in DRIFTS_COMMON}
        if not cells_s or not cells_a:
            continue

        mean_chi = statistics.mean(v["chi_emp"] for v in cells_s)
        ic_frac = sum(1 for v in cells_s
                      if v["mechanism"] == "inward-collapse") / len(cells_s)
        n_subs = len(set(classify_sub(v) for v in cells_a.values()))
        ic_dir_frac = (sum(1 for v in cells_a.values()
                           if classify_sub(v) == "IC-direct") / len(cells_a))

        metrics[n] = {
            "chi": mean_chi,
            "ic_frac": ic_frac,
            "n_subs": n_subs,
            "ic_dir_frac": ic_dir_frac,
        }

    # Print raw data
    print(f"\n  Raw metrics by N:")
    print(f"  {'N':>5s}  {'chi':>10s}  {'IC%':>8s}  {'#subs':>6s}  {'ICdir%':>8s}")
    for n in sorted(metrics):
        m = metrics[n]
        print(f"  {n:5d}  {m['chi']:10.4f}  {100*m['ic_frac']:7.1f}%  "
              f"{m['n_subs']:6d}  {100*m['ic_dir_frac']:7.1f}%")

    # Fit IC fraction vs N (does it approach a limit?)
    print(f"\n  IC fraction limit:")
    ns_big = [n for n in sorted(metrics) if n >= 20]
    if ns_big:
        ic_fracs_big = [metrics[n]["ic_frac"] for n in ns_big]
        print(f"    IC% at N>=20: {[f'{100*f:.1f}%' for f in ic_fracs_big]}")
        if len(ic_fracs_big) >= 2:
            mean_tail = statistics.mean(ic_fracs_big[-3:])
            print(f"    Tail mean (last 3): {100*mean_tail:.1f}%")

    # Fit IC-direct fraction vs N
    print(f"\n  IC-direct fraction limit:")
    if ns_big:
        icd_fracs = [metrics[n]["ic_dir_frac"] for n in ns_big]
        print(f"    ICdir% at N>=20: {[f'{100*f:.1f}%' for f in icd_fracs]}")
        if len(icd_fracs) >= 2:
            mean_tail = statistics.mean(icd_fracs[-3:])
            print(f"    Tail mean (last 3): {100*mean_tail:.1f}%")

    # Check: at what N does chi first hit the floor (0.01)?
    print(f"\n  Chi floor onset:")
    for n in sorted(metrics):
        if metrics[n]["chi"] <= 0.011:
            print(f"    Chi hits floor at N={n}")
            break
    else:
        print(f"    Chi has not hit floor yet")


# ============================================================
# 11. Summary: asymptotic law statement
# ============================================================

def sec11_summary(sweep_cube, atlas_cube, onset, ns_present):
    print("\n" + "=" * 72)
    print("SECTION 11: ASYMPTOTIC LAW SUMMARY")
    print("=" * 72)

    findings = []

    # 1. Boundary convergence
    boundary_flat_from = None
    for n in sorted(ns_present):
        all_flat = True
        for d in DRIFTS_COMMON:
            onsets_at = set()
            for r in RADII:
                o = onset.get((n, r, d))
                if o is not None:
                    onsets_at.add(o)
            if len(onsets_at) > 1:
                all_flat = False
                break
        if all_flat and boundary_flat_from is None:
            boundary_flat_from = n
        elif not all_flat:
            boundary_flat_from = None

    if boundary_flat_from:
        findings.append(f"Boundary is exactly planar from N={boundary_flat_from} onward")
    else:
        findings.append("Boundary has not fully converged to planar")

    # 2. Chi floor
    chi_floor_n = None
    for n in sorted(ns_present):
        cells = [v for k, v in sweep_cube.items()
                 if k[0] == n and k[3] == 0.0 and k[1] in ANGLES]
        if cells:
            mean_chi = statistics.mean(v["chi_emp"] for v in cells)
            if mean_chi <= 0.011 and chi_floor_n is None:
                chi_floor_n = n

    if chi_floor_n:
        findings.append(f"Chi hits floor (0.01) at N={chi_floor_n}")
    else:
        findings.append("Chi has not reached the floor")

    # 3. DECAY confinement
    decay_ns = set(k[0] for k, v in sweep_cube.items()
                   if v["mechanism"] == "DECAY")
    if decay_ns <= {4}:
        findings.append("DECAY remains exclusively N=4")
    else:
        findings.append(f"DECAY found at N={sorted(decay_ns)}")

    # 4. Taxonomy freeze
    last_three = sorted(ns_present)[-3:]
    sub_sets_tail = []
    for n in last_three:
        cells_a = {k: v for k, v in atlas_cube.items()
                   if k[0] == n and k[3] in DRIFTS_COMMON}
        if cells_a:
            sub_sets_tail.append(set(classify_sub(v) for v in cells_a.values()))

    if len(sub_sets_tail) >= 2 and all(s == sub_sets_tail[0] for s in sub_sets_tail):
        findings.append(f"Taxonomy frozen at {sorted(sub_sets_tail[0])} "
                        f"from N={last_three[0]}")
    else:
        findings.append(f"Taxonomy still fluctuating at N={last_three}")

    # 5. Drift irrelevance
    drift_rigid_from = None
    for n in sorted(ns_present):
        all_rigid = True
        for r in RADII:
            o0 = onset.get((n, r, 0.0))
            o2 = onset.get((n, r, 0.02))
            if o0 is not None and o2 is not None and o0 != o2:
                all_rigid = False
                break
        if all_rigid and drift_rigid_from is None:
            drift_rigid_from = n
        elif not all_rigid:
            drift_rigid_from = None

    if drift_rigid_from:
        findings.append(f"Drift becomes irrelevant from N={drift_rigid_from}")
    else:
        findings.append("Drift still has residual effect")

    # 6. Curvature/torsion convergence
    curv_zero_from = None
    for n in sorted(ns_present):
        has_curv = False
        for d in DRIFTS_COMMON:
            for r_idx in range(len(RADII)):
                if local_curvature(onset, n, d, r_idx) > 0:
                    has_curv = True
                    break
            if has_curv:
                break
        if not has_curv and curv_zero_from is None:
            curv_zero_from = n
        elif has_curv:
            curv_zero_from = None

    if curv_zero_from:
        findings.append(f"Curvature/torsion identically zero from N={curv_zero_from}")
    else:
        findings.append("Some residual curvature persists")

    # 7. New resonances
    new_resonances = []
    for n in NS_HIGH:
        cells_a = {k: v for k, v in atlas_cube.items()
                   if k[0] == n and k[3] in DRIFTS_COMMON}
        if cells_a:
            olw = sum(1 for v in cells_a.values()
                      if classify_sub(v) == "OL-wandering")
            if 100 * olw / len(cells_a) > 5:
                new_resonances.append(n)

    if new_resonances:
        findings.append(f"New resonance-like features at N={new_resonances}")
    else:
        findings.append("No new resonances at N=40..128")

    # Print findings
    print(f"\n  ASYMPTOTIC FINDINGS:")
    for i, f in enumerate(findings, 1):
        print(f"    {i}. {f}")

    # Compact summary table
    print(f"\n  ASYMPTOTIC SUMMARY TABLE:")
    print(f"  {'N':>5s}  {'chi':>10s}  {'mechs':>5s}  {'subs':>5s}  "
          f"{'IC%':>6s}  {'ICdir%':>7s}  {'rigid%':>7s}  "
          f"{'kappa':>7s}  {'planar':>7s}")

    for n in ns_present:
        cells_s = [v for k, v in sweep_cube.items()
                   if k[0] == n and k[3] in DRIFTS_COMMON and k[1] in ANGLES]
        cells_a = {k: v for k, v in atlas_cube.items()
                   if k[0] == n and k[3] in DRIFTS_COMMON}
        if not cells_s:
            continue

        chi = statistics.mean(v["chi_emp"] for v in cells_s)
        mechs = len(set(v["mechanism"] for v in cells_s))
        subs = len(set(classify_sub(v) for v in cells_a.values())) if cells_a else 0
        ic_pct = 100 * sum(1 for v in cells_s
                           if v["mechanism"] == "inward-collapse") / len(cells_s)
        ic_dir_pct = (100 * sum(1 for v in cells_a.values()
                                if classify_sub(v) == "IC-direct") / len(cells_a)
                      if cells_a else 0)

        # Drift rigidity
        rigid = 0
        total_r = 0
        for r in RADII:
            o0 = onset.get((n, r, 0.0))
            o2 = onset.get((n, r, 0.02))
            if o0 is not None and o2 is not None:
                total_r += 1
                if o0 == o2:
                    rigid += 1
        rigid_pct = 100 * rigid / total_r if total_r > 0 else 100

        # Curvature
        kappas = []
        for d in DRIFTS_COMMON:
            for r_idx in range(len(RADII)):
                k = local_curvature(onset, n, d, r_idx)
                if k > 0:
                    kappas.append(k)
        mean_k = statistics.mean(kappas) if kappas else 0

        # Planarity
        all_planar = True
        for d in DRIFTS_COMMON:
            onsets_at = set()
            for r in RADII:
                o = onset.get((n, r, d))
                if o is not None:
                    onsets_at.add(o)
            if len(onsets_at) > 1:
                all_planar = False
                break

        print(f"  {n:5d}  {chi:10.4f}  {mechs:5d}  {subs:5d}  "
              f"{ic_pct:5.1f}%  {ic_dir_pct:6.1f}%  {rigid_pct:6.0f}%  "
              f"{mean_k:7.1f}  {'YES' if all_planar else 'no':>7s}")

    # Final verdict
    print(f"\n  ASYMPTOTIC LAW VERDICT:")
    all_good = (chi_floor_n is not None and
                boundary_flat_from is not None and
                not new_resonances)
    if all_good:
        print(f"    ED-Arch exhibits a TRUE ASYMPTOTIC LAW.")
        print(f"    As N -> infinity:")
        print(f"      - The phase boundary converges to a fixed planar surface")
        print(f"      - Collapse time chi converges to the measurement floor")
        print(f"      - The taxonomy freezes to a finite set of sub-mechanisms")
        print(f"      - Drift becomes irrelevant (100% rigid)")
        print(f"      - Curvature and torsion vanish identically")
        print(f"      - No new resonances appear")
        print(f"      - The unified manifold becomes exactly flat")
    else:
        print(f"    ED-Arch shows STRONG CONVERGENCE but details vary:")
        if chi_floor_n:
            print(f"      - Chi hits floor at N={chi_floor_n}")
        else:
            print(f"      - Chi still decreasing (has not hit floor)")
        if boundary_flat_from:
            print(f"      - Boundary planar from N={boundary_flat_from}")
        else:
            print(f"      - Boundary not yet fully planar")


# ============================================================
# Main
# ============================================================

def main():
    print("High-N Asymptotic Probe Analysis")
    print("=" * 72)

    print("Loading ALL sweep data...")
    sweep_cube = load_all_sweeps()
    print(f"  {len(sweep_cube)} sweep records")

    print("Loading ALL atlas data...")
    atlas_cube = load_all_atlas()
    print(f"  {len(atlas_cube)} atlas records")

    # Determine which N values are actually present
    ns_present = sorted(set(k[0] for k in sweep_cube if k[0] in NS_ALL))
    print(f"  N values for analysis: {ns_present}")

    # Compute onset angles
    print("Computing IC onset angles...")
    onset = compute_ic_onset(sweep_cube, ns_present)
    print()

    sec01_inventory(sweep_cube, atlas_cube)
    sec02_boundary_convergence(sweep_cube, onset, ns_present)
    sec03_chi_scaling(sweep_cube, ns_present)
    sec04_taxonomy(atlas_cube, ns_present)
    sec05_drift_sensitivity(sweep_cube, onset, ns_present)
    sec06_curvature_torsion(onset, ns_present)
    sec07_planarity(onset, ns_present)
    sec08_resonance_scan(atlas_cube, sweep_cube, ns_present)
    sec09_ic_direct(atlas_cube, ns_present)
    sec10_scaling_fits(sweep_cube, atlas_cube, ns_present)
    sec11_summary(sweep_cube, atlas_cube, onset, ns_present)


if __name__ == "__main__":
    main()
