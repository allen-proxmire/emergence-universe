#!/usr/bin/env python
"""
check_law_ix_consolidation.py
================================
Law IX Consolidation Sweep: tests whether lambda(N) = d_min_init / merge_thr
predicts all five descriptive layers of the ED-Arch unified manifold at
non-canonical parameter values.

Lambda is the ratio of the initial minimum pairwise distance (after PBC
wrapping) to the merge threshold.  This quantity decreases with N because
larger rings wrap around the PBC box, bringing distant ring-neighbors into
physical proximity.

Sections:
  1.  Inventory & lambda computation
  2.  Lambda-predicted mechanism vs actual mechanism
  3.  Lambda-predicted chi vs actual chi (temporal layer)
  4.  Lambda-predicted drift class vs actual drift class
  5.  Lambda-predicted sub-mechanism vs actual sub-mechanism
  6.  Lambda-predicted boundary distance vs actual boundary
  7.  Regime verification: N=4 folded, N=8 resonance, N>=48 asymptotic
  8.  Temporal scaling exponent test (-7/3)
  9.  Drift-rigidity transition test
 10.  Taxonomy freeze-out test
 11.  Asymptotic hyperplane test
 12.  Cross-layer agreement summary
"""

import json
import os
import glob
import math
import statistics
from collections import defaultdict, Counter


SWEEP_DIR = os.path.dirname(__file__)

NS_ALL = [4, 8, 12, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 128]

# Merge threshold in engine units
BOX_PX = 400
MERGE_THR_PX = 23.5
SCALE = 1.0 / BOX_PX
MERGE_THR_ENGINE = MERGE_THR_PX * SCALE   # 0.05875
DT_ENGINE = 1.0 * SCALE                   # 0.0025
T_DECAY = 100


# ============================================================
# Data loading
# ============================================================

def load_all_atlas():
    """Load all atlas JSON files into a dict keyed by (N, angle, radius, drift)."""
    cube = {}
    for path in sorted(glob.glob(os.path.join(SWEEP_DIR, "n*_angle_*_r*_d*_atlas.json"))):
        try:
            with open(path) as f:
                data = json.load(f)
        except Exception:
            continue
        if "atlas" in data:
            key = (data["N"], data["angle"], data["radius"], data["drift"])
            cube[key] = data
    return cube


def load_all_sweep():
    """Load all sweep JSON files into a dict keyed by (N, angle, radius, drift)."""
    cube = {}
    for path in sorted(glob.glob(os.path.join(SWEEP_DIR, "n*_angle_*_r*_d*_sweep.json"))):
        try:
            with open(path) as f:
                data = json.load(f)
        except Exception:
            continue
        key = (data["N"], data["angle"], data["radius"], data["drift"])
        cube[key] = data
    return cube


# ============================================================
# Lambda computation -- PBC-aware initial d_min
# ============================================================

def compute_lambda(n, r_frac):
    """
    Compute lambda(N, r) = d_min_initial / merge_threshold.

    Places N particles equi-spaced on a ring of circumradius R_engine
    inside a unit PBC box, applies PBC wrapping, and computes the
    minimum pairwise distance.  This captures the key N-dependence:
    at large N, the ring exceeds the box and PBC wrapping brings
    distant ring-neighbors close together, reducing d_min and lambda.
    """
    d_px = max(1, int(round(r_frac * BOX_PX)))
    sin_pi_n = math.sin(math.pi / n)
    diameter_engine = (d_px / sin_pi_n) * SCALE
    R_engine = diameter_engine / 2.0
    box = 1.0

    # Place particles on ring, apply PBC
    theta = [2.0 * math.pi * k / n for k in range(n)]
    cx, cy = 0.5, 0.5
    positions = []
    for t in theta:
        x = (cx + R_engine * math.cos(t)) % box
        y = (cy + R_engine * math.sin(t)) % box
        positions.append((x, y))

    # Compute minimum PBC-aware pairwise distance
    d_min_val = float('inf')
    for i in range(n):
        for j in range(i + 1, n):
            dx = abs(positions[i][0] - positions[j][0])
            dy = abs(positions[i][1] - positions[j][1])
            if dx > box / 2:
                dx = box - dx
            if dy > box / 2:
                dy = box - dy
            d = math.sqrt(dx * dx + dy * dy)
            if d < d_min_val:
                d_min_val = d

    return d_min_val / MERGE_THR_ENGINE


# Lambda cache for performance
_lambda_cache = {}

def get_lambda(n, r_frac):
    """Cached lambda computation."""
    key = (n, round(r_frac, 4))
    if key not in _lambda_cache:
        _lambda_cache[key] = compute_lambda(n, r_frac)
    return _lambda_cache[key]


# ============================================================
# Sub-mechanism classification (unchanged from prior scripts)
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
# Lambda-based predictions
# ============================================================

def predict_mechanism_from_lambda(lam, angle, drift):
    """
    Predict mechanism from lambda, angle, and drift.

    Lambda > 1:    particles start far apart relative to merge threshold.
                   Complex trajectories possible; mechanism depends on angle.
                   At angle=0, drift=0 -> DECAY (N=4 only, lambda >> 1).
    Lambda ~ 1:    crossover regime, similar to lambda < 1 but resonance possible.
    Lambda < 1:    collapse is instantaneous.  Mechanism is a simple partition:
                   angle > 0 -> IC, angle = 0 -> OL, angle < 0 -> OPBC.
    """
    if lam > 5.0:
        # Lambda >> 1 regime (N=4 style): complex
        if abs(angle) < 0.01 and abs(drift) < 0.001:
            return "DECAY"
        elif angle > 0:
            return "inward-collapse"
        elif angle < 0:
            return "outward-PBC"
        else:
            # angle = 0 with nonzero drift
            if drift > 0:
                return "inward-collapse"
            else:
                return "other-late"
    elif lam > 1.0:
        # Lambda ~ 1-5 regime: transitional
        if angle > 0:
            return "inward-collapse"
        elif angle < 0:
            return "outward-PBC"
        else:
            return "other-late"
    else:
        # Lambda < 1: instantaneous collapse
        if angle > 0:
            return "inward-collapse"
        elif angle < 0:
            return "outward-PBC"
        else:
            return "other-late"


def predict_chi_from_lambda(lam, n):
    """
    Predict chi from lambda regime.
    Lambda < 1: chi at floor (0.01)
    Lambda ~ 1-3: chi low (0.01-1)
    Lambda >> 1: chi follows power-law ~ 318.49 * N^(-2.328)
    """
    if lam < 1.0:
        return 0.01
    else:
        chi_pred = 318.49 * (n ** (-2.328))
        return max(chi_pred, 0.01)


def predict_drift_class(lam):
    """
    Predict drift class:
      lambda >> 1 -> drift-fragile
      lambda <= 1 -> drift-rigid
    """
    if lam > 5.0:
        return "fragile"
    else:
        return "rigid"


def predict_subtype_from_lambda(lam, angle, drift, mech_pred):
    """
    Predict sub-mechanism from lambda regime.
    Lambda < 1: IC-direct, OPBC-single, OL-stalled
    Lambda ~ 1-5: IC-direct or IC-bounced, OPBC-single/glancing
    Lambda >> 5: IC-bounced, OPBC-winding, DE-tangent, etc.
    """
    if mech_pred == "DECAY":
        return "DE-tangent"
    elif mech_pred == "inward-collapse":
        if lam < 1.0:
            return "IC-direct"
        elif lam < 5.0:
            return "IC-direct"
        else:
            return "IC-bounced"
    elif mech_pred == "outward-PBC":
        if lam < 1.0:
            return "OPBC-single"
        elif lam < 5.0:
            return "OPBC-single"
        else:
            return "OPBC-winding"
    elif mech_pred == "other-late":
        return "OL-stalled"
    return mech_pred


def predict_boundary_sign(angle, lam):
    """
    Predict whether a cell is inside IC domain or outside.
    Simple rule: angle > 0 -> IC, angle <= 0 -> non-IC.
    At lambda >> 1 (N=4), the boundary is complex and this breaks down.
    """
    if angle > 0:
        return "inside-IC"
    else:
        return "outside-IC"


# ============================================================
# Section 1: Inventory & lambda computation
# ============================================================

def section_1(atlas_cube, sweep_cube):
    print("=" * 70)
    print("SECTION 1: Inventory & Lambda Computation")
    print("=" * 70)

    total_atlas = len(atlas_cube)
    total_sweep = len(sweep_cube)
    print(f"  Atlas records loaded:  {total_atlas}")
    print(f"  Sweep records loaded:  {total_sweep}")

    # Count per N
    n_counts = Counter()
    for key in atlas_cube:
        n_counts[key[0]] += 1
    print(f"\n  Records per N:")
    for n in sorted(n_counts):
        print(f"    N={n:3d}: {n_counts[n]:5d} atlas cells")

    # Lambda statistics per N
    print(f"\n  Lambda statistics per N (PBC-aware d_min / merge_thr):")
    print(f"    {'N':>4s}  {'lam_min':>8s}  {'lam_mean':>9s}  {'lam_max':>8s}  {'regime':>12s}")
    for n in NS_ALL:
        lambdas = []
        for key, rec in atlas_cube.items():
            if key[0] == n:
                lam = get_lambda(n, rec["radius"])
                lambdas.append(lam)
        if lambdas:
            lmin = min(lambdas)
            lmax = max(lambdas)
            lmean = statistics.mean(lambdas)
            if lmin > 3.0:
                regime = "folded"
            elif lmax < 1.0:
                regime = "planar"
            elif lmean > 2.0:
                regime = "mixed"
            elif lmean > 0.5:
                regime = "transitional"
            else:
                regime = "planar"
            print(f"    {n:4d}  {lmin:8.3f}  {lmean:9.3f}  {lmax:8.3f}  {regime:>12s}")

    # Lambda at representative radius (r=0.50) per N -- shows the N-dependence
    print(f"\n  Lambda at r=0.50 per N (shows N-dependence via PBC wrapping):")
    for n in NS_ALL:
        lam = get_lambda(n, 0.50)
        print(f"    N={n:3d}: lambda={lam:.4f}")

    print()


# ============================================================
# Section 2: Mechanism prediction
# ============================================================

def section_2(atlas_cube):
    print("=" * 70)
    print("SECTION 2: Lambda-Predicted Mechanism vs Actual")
    print("=" * 70)

    correct = 0
    total = 0
    errors_by_n = defaultdict(int)
    totals_by_n = defaultdict(int)
    confusion = defaultdict(int)

    for key, rec in atlas_cube.items():
        n, angle, radius, drift = key
        lam = get_lambda(n, radius)
        pred = predict_mechanism_from_lambda(lam, angle, drift)
        actual = rec["mechanism"]
        total += 1
        totals_by_n[n] += 1
        if pred == actual:
            correct += 1
        else:
            errors_by_n[n] += 1
            confusion[(pred, actual)] += 1

    pct = 100.0 * correct / total if total > 0 else 0
    print(f"  Overall mechanism accuracy: {correct}/{total} ({pct:.1f}%)")

    print(f"\n  Accuracy by N:")
    for n in sorted(set(k[0] for k in atlas_cube)):
        tot_n = totals_by_n[n]
        err_n = errors_by_n.get(n, 0)
        cor_n = tot_n - err_n
        pct_n = 100.0 * cor_n / tot_n if tot_n > 0 else 0
        lam_rep = get_lambda(n, 0.50)
        print(f"    N={n:3d} (lam@r0.5={lam_rep:.3f}): {cor_n}/{tot_n} ({pct_n:.1f}%)")

    if confusion:
        print(f"\n  Top confusion pairs (predicted -> actual):")
        for (pred, actual), count in sorted(confusion.items(), key=lambda x: -x[1])[:10]:
            print(f"    {pred:20s} -> {actual:20s}: {count}")

    # Accuracy by lambda bin
    print(f"\n  Mechanism accuracy by lambda bin:")
    bins = [(0, 0.3), (0.3, 1.0), (1.0, 3.0), (3.0, 8.0), (8.0, 200.0)]
    for lo, hi in bins:
        cor_bin = 0
        tot_bin = 0
        for key, rec in atlas_cube.items():
            lam = get_lambda(key[0], key[2])
            if lo <= lam < hi:
                pred = predict_mechanism_from_lambda(lam, key[1], key[3])
                tot_bin += 1
                if pred == rec["mechanism"]:
                    cor_bin += 1
        if tot_bin > 0:
            pct_bin = 100.0 * cor_bin / tot_bin
            print(f"    lambda in [{lo:.1f}, {hi:.1f}): {cor_bin}/{tot_bin} ({pct_bin:.1f}%)")

    print()
    return correct, total


# ============================================================
# Section 3: Temporal (chi) prediction
# ============================================================

def section_3(atlas_cube):
    print("=" * 70)
    print("SECTION 3: Lambda-Predicted Chi vs Actual Chi")
    print("=" * 70)

    residuals_by_n = defaultdict(list)
    all_residuals = []

    for key, rec in atlas_cube.items():
        n, angle, radius, drift = key
        lam = get_lambda(n, radius)
        chi_actual = rec["chi_emp"]
        chi_pred = predict_chi_from_lambda(lam, n)
        resid = abs(chi_actual - chi_pred)
        residuals_by_n[n].append(resid)
        all_residuals.append(resid)

    mean_resid = statistics.mean(all_residuals)
    median_resid = statistics.median(all_residuals)
    print(f"  Overall |chi_pred - chi_actual|: mean={mean_resid:.3f}, median={median_resid:.3f}")

    print(f"\n  Residuals by N:")
    print(f"    {'N':>4s}  {'lam@r0.5':>9s}  {'chi_pred':>9s}  {'chi_mean':>9s}  {'|resid|':>9s}  {'max_resid':>10s}")
    for n in sorted(residuals_by_n):
        chi_actuals = [rec["chi_emp"] for key, rec in atlas_cube.items() if key[0] == n]
        lam_rep = get_lambda(n, 0.50)
        chi_pred = predict_chi_from_lambda(lam_rep, n)
        chi_mean = statistics.mean(chi_actuals)
        res_mean = statistics.mean(residuals_by_n[n])
        res_max = max(residuals_by_n[n])
        print(f"    {n:4d}  {lam_rep:9.3f}  {chi_pred:9.3f}  {chi_mean:9.3f}  {res_mean:9.3f}  {res_max:10.3f}")

    # At-floor accuracy (lambda < 1)
    floor_cells = [(key, rec) for key, rec in atlas_cube.items()
                   if get_lambda(key[0], key[2]) < 1.0]
    if floor_cells:
        at_floor = sum(1 for _, rec in floor_cells if rec["chi_emp"] <= 0.02)
        pct_floor = 100.0 * at_floor / len(floor_cells)
        print(f"\n  Lambda < 1 cells at floor (chi<=0.02): {at_floor}/{len(floor_cells)} ({pct_floor:.1f}%)")

    # N >= 32 at floor
    n32_cells = [(key, rec) for key, rec in atlas_cube.items() if key[0] >= 32]
    at_floor_32 = sum(1 for _, rec in n32_cells if rec["chi_emp"] <= 0.02)
    pct_floor_32 = 100.0 * at_floor_32 / len(n32_cells) if n32_cells else 0
    print(f"  N>=32 cells at floor (chi<=0.02): {at_floor_32}/{len(n32_cells)} ({pct_floor_32:.1f}%)")

    print()


# ============================================================
# Section 4: Drift class prediction
# ============================================================

def section_4(atlas_cube):
    print("=" * 70)
    print("SECTION 4: Lambda-Predicted Drift Class vs Actual")
    print("=" * 70)

    # Group cells by (N, angle, radius), compare across drifts
    by_nar = defaultdict(dict)
    for key, rec in atlas_cube.items():
        n, angle, radius, drift = key
        by_nar[(n, angle, radius)][drift] = rec["mechanism"]

    rigid_correct = 0
    rigid_total = 0
    fragile_correct = 0
    fragile_total = 0

    for (n, angle, radius), drift_mechs in by_nar.items():
        if len(drift_mechs) < 2:
            continue
        lam = get_lambda(n, radius)
        pred_class = predict_drift_class(lam)
        mechs = list(drift_mechs.values())
        actual_rigid = len(set(mechs)) == 1

        if pred_class == "rigid":
            rigid_total += 1
            if actual_rigid:
                rigid_correct += 1
        else:
            fragile_total += 1
            if not actual_rigid:
                fragile_correct += 1

    total = rigid_total + fragile_total
    correct = rigid_correct + fragile_correct
    pct = 100.0 * correct / total if total > 0 else 0

    print(f"  Overall drift-class accuracy: {correct}/{total} ({pct:.1f}%)")
    if rigid_total:
        print(f"    Predicted rigid,  actually rigid:  {rigid_correct}/{rigid_total} "
              f"({100*rigid_correct/rigid_total:.1f}%)")
    if fragile_total:
        print(f"    Predicted fragile, actually fragile: {fragile_correct}/{fragile_total} "
              f"({100*fragile_correct/fragile_total:.1f}%)")

    # Per-N drift rigidity
    print(f"\n  Drift rigidity per N:")
    for n in NS_ALL:
        groups = [(k, v) for k, v in by_nar.items() if k[0] == n and len(v) >= 2]
        if not groups:
            continue
        n_rigid = sum(1 for _, dm in groups if len(set(dm.values())) == 1)
        pct_r = 100.0 * n_rigid / len(groups)
        lam_rep = get_lambda(n, 0.50)
        print(f"    N={n:3d} (lam@r0.5={lam_rep:.3f}): {n_rigid}/{len(groups)} rigid ({pct_r:.1f}%)")

    # By lambda bin
    print(f"\n  Drift rigidity by lambda bin:")
    bins = [(0, 0.3), (0.3, 1.0), (1.0, 3.0), (3.0, 8.0), (8.0, 200.0)]
    for lo, hi in bins:
        groups_in_bin = []
        for (n, angle, radius), dm in by_nar.items():
            if len(dm) < 2:
                continue
            lam = get_lambda(n, radius)
            if lo <= lam < hi:
                groups_in_bin.append(len(set(dm.values())) == 1)
        if groups_in_bin:
            n_rigid = sum(groups_in_bin)
            pct = 100.0 * n_rigid / len(groups_in_bin)
            print(f"    lambda in [{lo:.1f}, {hi:.1f}): {n_rigid}/{len(groups_in_bin)} ({pct:.1f}%) rigid")

    print()
    return correct, total


# ============================================================
# Section 5: Sub-mechanism prediction
# ============================================================

def section_5(atlas_cube):
    print("=" * 70)
    print("SECTION 5: Lambda-Predicted Sub-mechanism vs Actual")
    print("=" * 70)

    correct = 0
    total = 0
    errors_by_n = defaultdict(int)
    totals_by_n = defaultdict(int)

    for key, rec in atlas_cube.items():
        n, angle, radius, drift = key
        lam = get_lambda(n, radius)
        mech_pred = predict_mechanism_from_lambda(lam, angle, drift)
        sub_pred = predict_subtype_from_lambda(lam, angle, drift, mech_pred)
        sub_actual = classify_sub(rec)
        total += 1
        totals_by_n[n] += 1
        if sub_pred == sub_actual:
            correct += 1
        else:
            errors_by_n[n] += 1

    pct = 100.0 * correct / total if total > 0 else 0
    print(f"  Overall subtype accuracy: {correct}/{total} ({pct:.1f}%)")

    print(f"\n  Accuracy by N:")
    for n in sorted(totals_by_n):
        tot_n = totals_by_n[n]
        err_n = errors_by_n.get(n, 0)
        cor_n = tot_n - err_n
        pct_n = 100.0 * cor_n / tot_n if tot_n > 0 else 0
        print(f"    N={n:3d}: {cor_n}/{tot_n} ({pct_n:.1f}%)")

    # Sub-mech distribution for N>=48
    print(f"\n  Sub-mechanism distribution for N>=48:")
    sub_counts = Counter()
    for key, rec in atlas_cube.items():
        if key[0] >= 48:
            sub_counts[classify_sub(rec)] += 1
    total_48 = sum(sub_counts.values())
    for sub, cnt in sub_counts.most_common():
        print(f"    {sub:20s}: {cnt:5d} ({100*cnt/total_48:.1f}%)")

    print()
    return correct, total


# ============================================================
# Section 6: Boundary distance prediction
# ============================================================

def section_6(atlas_cube):
    print("=" * 70)
    print("SECTION 6: Lambda-Predicted Boundary Sign vs Actual")
    print("=" * 70)

    correct = 0
    total = 0
    false_ic = 0
    missed_ic = 0

    for key, rec in atlas_cube.items():
        n, angle, radius, drift = key
        lam = get_lambda(n, radius)
        pred_sign = predict_boundary_sign(angle, lam)
        actual_is_ic = rec["mechanism"] == "inward-collapse"
        total += 1
        if pred_sign == "inside-IC" and actual_is_ic:
            correct += 1
        elif pred_sign == "outside-IC" and not actual_is_ic:
            correct += 1
        elif pred_sign == "inside-IC" and not actual_is_ic:
            false_ic += 1
        else:
            missed_ic += 1

    pct = 100.0 * correct / total if total > 0 else 0
    print(f"  Overall boundary-sign accuracy: {correct}/{total} ({pct:.1f}%)")
    print(f"  False IC (predicted IC, actual not): {false_ic}")
    print(f"  Missed IC (predicted not IC, actual IC): {missed_ic}")

    # By lambda bin
    print(f"\n  Boundary accuracy by lambda bin:")
    bins = [(0, 0.3), (0.3, 1.0), (1.0, 3.0), (3.0, 8.0), (8.0, 200.0)]
    for lo, hi in bins:
        cor_bin = 0
        tot_bin = 0
        for key, rec in atlas_cube.items():
            lam = get_lambda(key[0], key[2])
            if lo <= lam < hi:
                pred = predict_boundary_sign(key[1], lam)
                actual_ic = rec["mechanism"] == "inward-collapse"
                tot_bin += 1
                if (pred == "inside-IC") == actual_ic:
                    cor_bin += 1
        if tot_bin > 0:
            pct_bin = 100.0 * cor_bin / tot_bin
            print(f"    lambda in [{lo:.1f}, {hi:.1f}): {cor_bin}/{tot_bin} ({pct_bin:.1f}%)")

    # By N
    print(f"\n  Boundary accuracy by N:")
    for n in NS_ALL:
        cells = [(key, rec) for key, rec in atlas_cube.items() if key[0] == n]
        if not cells:
            continue
        cor_n = sum(1 for key, rec in cells
                    if (predict_boundary_sign(key[1], get_lambda(n, key[2])) == "inside-IC")
                    == (rec["mechanism"] == "inward-collapse"))
        pct_n = 100.0 * cor_n / len(cells)
        print(f"    N={n:3d}: {cor_n}/{len(cells)} ({pct_n:.1f}%)")

    print()
    return correct, total


# ============================================================
# Section 7: Regime verification
# ============================================================

def section_7(atlas_cube):
    print("=" * 70)
    print("SECTION 7: Regime Verification (Folded / Resonant / Asymptotic)")
    print("=" * 70)

    # N=4 folded regime
    print("  --- N=4 Folded Regime ---")
    n4 = {k: v for k, v in atlas_cube.items() if k[0] == 4}
    if n4:
        mechs_4 = Counter(v["mechanism"] for v in n4.values())
        subs_4 = Counter(classify_sub(v) for v in n4.values())
        chi_4 = [v["chi_emp"] for v in n4.values()]
        lam_4 = [get_lambda(4, k[2]) for k in n4]
        print(f"    Cells: {len(n4)}")
        print(f"    Lambda range: [{min(lam_4):.3f}, {max(lam_4):.3f}]")
        print(f"    Mechanisms: {dict(mechs_4)}")
        print(f"    Sub-types: {len(subs_4)} distinct")
        print(f"    DECAY present: {'DECAY' in mechs_4}")
        print(f"    Chi range: [{min(chi_4):.3f}, {max(chi_4):.3f}], mean={statistics.mean(chi_4):.3f}")
        print(f"    FOLDED regime confirmed: DECAY={'DECAY' in mechs_4}, "
              f"diversity>7={len(subs_4) >= 7}, max_chi>10={max(chi_4) > 10}")

    # N=8 resonance
    print(f"\n  --- N=8 Resonant Regime ---")
    n8 = {k: v for k, v in atlas_cube.items() if k[0] == 8}
    if n8:
        subs_8 = Counter(classify_sub(v) for v in n8.values())
        total_8 = sum(subs_8.values())
        ol_wander_pct = 100.0 * subs_8.get("OL-wandering", 0) / total_8
        chi_8_pos = [v["chi_emp"] for k, v in n8.items() if k[1] > 0]
        lam_8 = [get_lambda(8, k[2]) for k in n8]
        print(f"    Cells: {len(n8)}")
        print(f"    Lambda range: [{min(lam_8):.3f}, {max(lam_8):.3f}]")
        print(f"    OL-wandering: {subs_8.get('OL-wandering', 0)}/{total_8} ({ol_wander_pct:.1f}%)")
        if chi_8_pos:
            print(f"    Chi mean (positive angles): {statistics.mean(chi_8_pos):.3f}")
        print(f"    RESONANCE confirmed (OL-wandering > 10%): {ol_wander_pct > 10}")

    # N>=48 asymptotic
    print(f"\n  --- N>=48 Asymptotic Regime ---")
    n48plus = {k: v for k, v in atlas_cube.items() if k[0] >= 48}
    if n48plus:
        chi_48 = [v["chi_emp"] for v in n48plus.values()]
        at_floor = sum(1 for c in chi_48 if c <= 0.02)
        pct_floor = 100.0 * at_floor / len(chi_48)
        subs_48 = Counter(classify_sub(v) for v in n48plus.values())
        total_48 = sum(subs_48.values())
        ic_direct_pct = 100.0 * subs_48.get("IC-direct", 0) / total_48
        lam_48 = [get_lambda(k[0], k[2]) for k in n48plus]
        print(f"    Cells: {len(n48plus)}")
        print(f"    Lambda range: [{min(lam_48):.3f}, {max(lam_48):.3f}]")
        print(f"    Chi at floor (<=0.02): {at_floor}/{len(chi_48)} ({pct_floor:.1f}%)")
        print(f"    IC-direct: {subs_48.get('IC-direct',0)}/{total_48} ({ic_direct_pct:.1f}%)")
        print(f"    Sub-types: {len(subs_48)}")
        for sub, cnt in subs_48.most_common():
            print(f"      {sub:20s}: {cnt:5d} ({100*cnt/total_48:.1f}%)")
        print(f"    ASYMPTOTIC confirmed (chi at floor > 95%): {pct_floor > 95}")

    print()


# ============================================================
# Section 8: Temporal scaling exponent
# ============================================================

def section_8(atlas_cube):
    print("=" * 70)
    print("SECTION 8: Temporal Scaling Exponent Test (-7/3)")
    print("=" * 70)

    chi_by_n = defaultdict(list)
    for key, rec in atlas_cube.items():
        chi_by_n[key[0]].append(rec["chi_emp"])

    fit_ns = []
    fit_chis = []
    print(f"  Mean chi per N:")
    for n in sorted(chi_by_n):
        mean_chi = statistics.mean(chi_by_n[n])
        print(f"    N={n:3d}: mean_chi={mean_chi:.4f} (n_cells={len(chi_by_n[n])})")
        if mean_chi > 0.02:
            fit_ns.append(n)
            fit_chis.append(mean_chi)

    # Log-log least squares
    if len(fit_ns) >= 3:
        log_ns = [math.log(n) for n in fit_ns]
        log_chis = [math.log(c) for c in fit_chis]
        n_pts = len(log_ns)
        sum_x = sum(log_ns)
        sum_y = sum(log_chis)
        sum_xx = sum(x*x for x in log_ns)
        sum_xy = sum(x*y for x, y in zip(log_ns, log_chis))
        denom = n_pts * sum_xx - sum_x * sum_x
        if abs(denom) > 1e-15:
            alpha = (n_pts * sum_xy - sum_x * sum_y) / denom
            ln_c = (sum_y - alpha * sum_x) / n_pts
            c = math.exp(ln_c)
            y_mean = sum_y / n_pts
            ss_tot = sum((y - y_mean)**2 for y in log_chis)
            ss_res = sum((y - (ln_c + alpha * x))**2 for x, y in zip(log_ns, log_chis))
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

            print(f"\n  Power-law fit (N where mean_chi > 0.02):")
            print(f"    chi ~ {c:.2f} * N^({alpha:.3f})")
            print(f"    R^2 = {r2:.4f}")
            print(f"    Expected: -2.333 (-7/3)")
            print(f"    Measured:  {alpha:.3f}")
            print(f"    Deviation: {abs(alpha - (-7/3)):.3f}")
            within = abs(alpha - (-7/3)) < 0.5
            print(f"    Exponent within 0.5 of -7/3: {within}")

    print()


# ============================================================
# Section 9: Drift-rigidity transition
# ============================================================

def section_9(atlas_cube):
    print("=" * 70)
    print("SECTION 9: Drift-Rigidity Transition Test")
    print("=" * 70)

    by_nar = defaultdict(dict)
    for key, rec in atlas_cube.items():
        n, angle, radius, drift = key
        by_nar[(n, angle, radius)][drift] = rec["mechanism"]

    print(f"  Drift rigidity per N:")
    for n in NS_ALL:
        groups = [(k, v) for k, v in by_nar.items() if k[0] == n and len(v) >= 2]
        if not groups:
            continue
        n_rigid = sum(1 for _, dm in groups if len(set(dm.values())) == 1)
        pct = 100.0 * n_rigid / len(groups)
        lam_rep = get_lambda(n, 0.50)
        print(f"    N={n:3d} (lam@r0.5={lam_rep:.3f}): {n_rigid}/{len(groups)} ({pct:.1f}%) drift-rigid")

    # By lambda bin
    print(f"\n  Drift rigidity by lambda bin:")
    bins = [(0, 0.3), (0.3, 1.0), (1.0, 3.0), (3.0, 8.0), (8.0, 200.0)]
    for lo, hi in bins:
        groups_in_bin = []
        for (n, angle, radius), dm in by_nar.items():
            if len(dm) < 2:
                continue
            lam = get_lambda(n, radius)
            if lo <= lam < hi:
                groups_in_bin.append(len(set(dm.values())) == 1)
        if groups_in_bin:
            n_rigid = sum(groups_in_bin)
            pct = 100.0 * n_rigid / len(groups_in_bin)
            print(f"    lambda in [{lo:.1f}, {hi:.1f}): {n_rigid}/{len(groups_in_bin)} ({pct:.1f}%) rigid")

    print()


# ============================================================
# Section 10: Taxonomy freeze-out
# ============================================================

def section_10(atlas_cube):
    print("=" * 70)
    print("SECTION 10: Taxonomy Freeze-Out Test")
    print("=" * 70)

    print(f"  Sub-mechanism taxonomy per N:")
    print(f"    {'N':>4s}  {'n_types':>7s}  {'H_sub':>7s}  {'dominant':>20s}  {'dom_pct':>8s}  {'lam@r0.5':>9s}")

    entropies = {}
    for n in NS_ALL:
        cells = [rec for key, rec in atlas_cube.items() if key[0] == n]
        if not cells:
            continue
        subs = Counter(classify_sub(rec) for rec in cells)
        total_n = sum(subs.values())
        n_types = len(subs)
        h = 0.0
        for cnt in subs.values():
            p = cnt / total_n
            if p > 0:
                h -= p * math.log(p)
        entropies[n] = h
        dominant = subs.most_common(1)[0]
        dom_pct = 100.0 * dominant[1] / total_n
        lam_rep = get_lambda(n, 0.50)
        print(f"    {n:4d}  {n_types:7d}  {h:7.3f}  {dominant[0]:>20s}  {dom_pct:7.1f}%  {lam_rep:9.3f}")

    if len([n for n in entropies if n >= 20]) >= 3:
        tail = [entropies[n] for n in sorted(entropies) if n >= 20]
        tail_mean = statistics.mean(tail)
        tail_std = statistics.stdev(tail) if len(tail) > 1 else 0
        cv = tail_std / tail_mean if tail_mean > 0 else 0
        print(f"\n  Entropy convergence (N>=20):")
        print(f"    Mean H_sub: {tail_mean:.3f}")
        print(f"    Std H_sub:  {tail_std:.3f}")
        print(f"    CV:         {100*cv:.1f}%")
        print(f"    Taxonomy frozen (CV < 20%): {cv < 0.20}")

    types_48 = set()
    for key, rec in atlas_cube.items():
        if key[0] >= 48:
            types_48.add(classify_sub(rec))
    print(f"\n  Sub-types at N>=48: {sorted(types_48)} ({len(types_48)} types)")

    print()


# ============================================================
# Section 11: Asymptotic hyperplane test
# ============================================================

def section_11(atlas_cube):
    print("=" * 70)
    print("SECTION 11: Asymptotic Hyperplane Test")
    print("=" * 70)

    # Non-canonical angles: does IC dominate at positive angles for N>=20?
    print("  --- Boundary flatness at non-canonical angles ---")
    for test_angle in [7, 3, 0.5, -0.5, -3, -7]:
        for n in [20, 32, 48, 64, 128]:
            cells = [(key, rec) for key, rec in atlas_cube.items()
                     if key[0] == n and abs(key[1] - test_angle) < 0.01]
            if not cells:
                continue
            ic_count = sum(1 for _, rec in cells if rec["mechanism"] == "inward-collapse")
            pct = 100.0 * ic_count / len(cells)
            print(f"    N={n:3d}, angle={test_angle:+5.1f}: {ic_count}/{len(cells)} IC ({pct:.0f}%)")

    # Non-canonical radii at positive angles
    print(f"\n  --- IC at non-canonical radii (positive angles, N>=20) ---")
    for r in [0.025, 0.125, 0.275, 0.525, 0.675, 0.825]:
        for n in [20, 48, 128]:
            cells = [(key, rec) for key, rec in atlas_cube.items()
                     if key[0] == n and abs(key[2] - r) < 0.005 and key[1] == 10]
            if not cells:
                continue
            ic_count = sum(1 for _, rec in cells if rec["mechanism"] == "inward-collapse")
            chi_vals = [rec["chi_emp"] for _, rec in cells]
            lam = get_lambda(n, r)
            print(f"    N={n:3d}, r={r:.3f} (lam={lam:.3f}): {ic_count}/{len(cells)} IC, "
                  f"chi_mean={statistics.mean(chi_vals):.3f}")

    # Non-canonical drifts
    print(f"\n  --- Drift invariance at non-canonical drifts (N>=20) ---")
    for d in [0.003, 0.007, 0.013, 0.017, 0.030]:
        for n in [20, 48, 128]:
            cells = [(key, rec) for key, rec in atlas_cube.items()
                     if key[0] == n and abs(key[3] - d) < 0.0005 and key[1] == 10]
            if not cells:
                continue
            ic_count = sum(1 for _, rec in cells if rec["mechanism"] == "inward-collapse")
            chi_vals = [rec["chi_emp"] for _, rec in cells]
            print(f"    N={n:3d}, drift={d:.3f}: {ic_count}/{len(cells)} IC, "
                  f"chi_mean={statistics.mean(chi_vals):.3f}")

    # Global planarity: positive angle -> IC for N>=48
    print(f"\n  --- Global planarity: positive angles -> IC for N>=48 ---")
    for n in [48, 64, 80, 96, 128]:
        pos_cells = [(key, rec) for key, rec in atlas_cube.items()
                     if key[0] == n and key[1] > 0]
        if not pos_cells:
            continue
        ic = sum(1 for _, rec in pos_cells if rec["mechanism"] == "inward-collapse")
        pct = 100.0 * ic / len(pos_cells)
        print(f"    N={n:3d}: {ic}/{len(pos_cells)} ({pct:.0f}%)")

    print()


# ============================================================
# Section 12: Cross-layer agreement summary
# ============================================================

def section_12(mech_correct, mech_total, drift_correct, drift_total,
               sub_correct, sub_total, bound_correct, bound_total, atlas_cube):
    print("=" * 70)
    print("SECTION 12: Cross-Layer Agreement Summary")
    print("=" * 70)

    mech_pct = 100.0 * mech_correct / mech_total if mech_total else 0
    drift_pct = 100.0 * drift_correct / drift_total if drift_total else 0
    sub_pct = 100.0 * sub_correct / sub_total if sub_total else 0
    bound_pct = 100.0 * bound_correct / bound_total if bound_total else 0

    print(f"  Layer                    Correct / Total     Accuracy")
    print(f"  {'-'*58}")
    print(f"  Mechanism (lam->mech)     {mech_correct:6d} / {mech_total:<6d}    {mech_pct:.1f}%")
    print(f"  Drift class (lam->DC)     {drift_correct:6d} / {drift_total:<6d}    {drift_pct:.1f}%")
    print(f"  Sub-mechanism (lam->sub)  {sub_correct:6d} / {sub_total:<6d}    {sub_pct:.1f}%")
    print(f"  Boundary sign (lam->bnd)  {bound_correct:6d} / {bound_total:<6d}    {bound_pct:.1f}%")

    # Per-regime breakdown
    print(f"\n  Per-regime mechanism accuracy:")
    for label, n_range in [("Folded (N=4)", [4]),
                           ("Resonant (N=8)", [8]),
                           ("Transitional (N=12-32)", [12, 20, 24, 28, 32]),
                           ("Asymptotic (N>=48)", [n for n in NS_ALL if n >= 48])]:
        correct_r = 0
        total_r = 0
        for key, rec in atlas_cube.items():
            if key[0] in n_range:
                n, angle, radius, drift = key
                lam = get_lambda(n, radius)
                pred = predict_mechanism_from_lambda(lam, angle, drift)
                total_r += 1
                if pred == rec["mechanism"]:
                    correct_r += 1
        if total_r > 0:
            pct_r = 100.0 * correct_r / total_r
            print(f"    {label:30s}: {correct_r}/{total_r} ({pct_r:.1f}%)")

    # Accuracy excluding N=4
    n4_cells = sum(1 for k in atlas_cube if k[0] == 4)
    n4_errors = sum(1 for k, v in atlas_cube.items()
                    if k[0] == 4 and predict_mechanism_from_lambda(
                        get_lambda(4, k[2]), k[1], k[3]) != v["mechanism"])
    non_n4_errors = (mech_total - mech_correct) - n4_errors
    non_n4_total = mech_total - n4_cells
    non_n4_pct = 100.0 * (non_n4_total - non_n4_errors) / non_n4_total if non_n4_total else 0

    print(f"\n  Excluding N=4:")
    print(f"    N>4 mechanism accuracy: {non_n4_total - non_n4_errors}/{non_n4_total} ({non_n4_pct:.1f}%)")

    # Verdict
    print(f"\n  ===== LAW IX CONSOLIDATION VERDICT =====")

    # The real test: does the N-regime structure hold at non-canonical parameters?
    # Check: positive angle -> IC at N>=20 (non-canonical angles, radii, drifts)
    pos_n20 = [(k, v) for k, v in atlas_cube.items() if k[0] >= 20 and k[1] > 0]
    ic_pos = sum(1 for _, v in pos_n20 if v["mechanism"] == "inward-collapse")
    pct_pos = 100.0 * ic_pos / len(pos_n20) if pos_n20 else 0

    neg_n20 = [(k, v) for k, v in atlas_cube.items() if k[0] >= 20 and k[1] < 0]
    not_ic_neg = sum(1 for _, v in neg_n20 if v["mechanism"] != "inward-collapse")
    pct_neg = 100.0 * not_ic_neg / len(neg_n20) if neg_n20 else 0

    zero_n20 = [(k, v) for k, v in atlas_cube.items() if k[0] >= 20 and abs(k[1]) < 0.01]
    ol_zero = sum(1 for _, v in zero_n20 if v["mechanism"] == "other-late")
    pct_zero = 100.0 * ol_zero / len(zero_n20) if zero_n20 else 0

    n48_floor = [(k, v) for k, v in atlas_cube.items() if k[0] >= 48]
    at_floor = sum(1 for _, v in n48_floor if v["chi_emp"] <= 0.02)
    pct_floor = 100.0 * at_floor / len(n48_floor) if n48_floor else 0

    print(f"  Regime structure holds at non-canonical params:")
    print(f"    N>=20 positive angle -> IC:       {ic_pos}/{len(pos_n20)} ({pct_pos:.1f}%)")
    print(f"    N>=20 negative angle -> non-IC:   {not_ic_neg}/{len(neg_n20)} ({pct_neg:.1f}%)")
    print(f"    N>=20 angle=0 -> other-late:      {ol_zero}/{len(zero_n20)} ({pct_zero:.1f}%)")
    print(f"    N>=48 chi at floor:               {at_floor}/{len(n48_floor)} ({pct_floor:.1f}%)")

    n8_cells = {k: v for k, v in atlas_cube.items() if k[0] == 8}
    ol_wander = sum(1 for v in n8_cells.values() if classify_sub(v) == "OL-wandering")
    pct_wander = 100.0 * ol_wander / len(n8_cells) if n8_cells else 0
    print(f"    N=8 OL-wandering resonance:       {ol_wander}/{len(n8_cells)} ({pct_wander:.1f}%)")

    print(f"\n  Lambda discriminates regimes:")
    # Check if lambda < 1 predicts chi at floor
    lam_lt1 = [(k, v) for k, v in atlas_cube.items() if get_lambda(k[0], k[2]) < 1.0]
    lam_lt1_floor = sum(1 for _, v in lam_lt1 if v["chi_emp"] <= 0.02)
    pct_lt1 = 100.0 * lam_lt1_floor / len(lam_lt1) if lam_lt1 else 0
    print(f"    lambda < 1 -> chi at floor:  {lam_lt1_floor}/{len(lam_lt1)} ({pct_lt1:.1f}%)")

    lam_gt5 = [(k, v) for k, v in atlas_cube.items() if get_lambda(k[0], k[2]) > 5.0]
    lam_gt5_slow = sum(1 for _, v in lam_gt5 if v["chi_emp"] > 1.0)
    pct_gt5 = 100.0 * lam_gt5_slow / len(lam_gt5) if lam_gt5 else 0
    print(f"    lambda > 5 -> chi > 1.0:     {lam_gt5_slow}/{len(lam_gt5)} ({pct_gt5:.1f}%)")

    print()


# ============================================================
# Main
# ============================================================

def main():
    print("Law IX Consolidation Sweep Analysis")
    print("=" * 70)
    print()

    atlas_cube = load_all_atlas()
    sweep_cube = load_all_sweep()

    section_1(atlas_cube, sweep_cube)
    mech_correct, mech_total = section_2(atlas_cube)
    section_3(atlas_cube)
    drift_correct, drift_total = section_4(atlas_cube)
    sub_correct, sub_total = section_5(atlas_cube)
    bound_correct, bound_total = section_6(atlas_cube)
    section_7(atlas_cube)
    section_8(atlas_cube)
    section_9(atlas_cube)
    section_10(atlas_cube)
    section_11(atlas_cube)
    section_12(mech_correct, mech_total, drift_correct, drift_total,
               sub_correct, sub_total, bound_correct, bound_total, atlas_cube)


if __name__ == "__main__":
    main()
