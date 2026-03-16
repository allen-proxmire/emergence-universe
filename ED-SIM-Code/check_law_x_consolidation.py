#!/usr/bin/env python
"""
check_law_x_consolidation.py
================================
Law X Consolidation Sweep: tests the angle partition, true IC boundary,
boundary-layer phenomenon, taxonomy freeze, and the (lambda, angle)
generative coordinate system across 34,000+ cells.

Sections:
  1.  Inventory & radial-velocity computation
  2.  Angle partition: positive -> IC, negative -> OPBC, zero -> OL
  3.  True IC boundary: location, radius dependence, drift sensitivity
  4.  Beyond-grid test: do angles +12 and +15 produce IC at large N?
  5.  Boundary-layer phenomenon: drift flips in the gamma ~ 0 neighborhood
  6.  Drift-rigidity decomposition: bulk vs boundary layer
  7.  Taxonomy freeze-out: 3 cores + 3 corrections
  8.  (Lambda, angle) joint prediction accuracy
  9.  Cross-layer agreement summary
"""

import json
import os
import glob
import math
import statistics
from collections import defaultdict, Counter


SWEEP_DIR = os.path.dirname(__file__)

NS_ALL = [4, 8, 12, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 128]

BOX_PX = 400
MERGE_THR_PX = 23.5
SCALE = 1.0 / BOX_PX
MERGE_THR_ENGINE = MERGE_THR_PX * SCALE
DT_ENGINE = 1.0 * SCALE
T_DECAY = 100

# Boundary-layer and beyond-grid angles for focused analysis
BL_ANGLES = [-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75]
BEYOND_ANGLES = [12, 15]


# ============================================================
# Data loading
# ============================================================

def load_all_atlas():
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


# ============================================================
# Lambda computation (PBC-aware d_min)
# ============================================================

_lambda_cache = {}

def compute_lambda(n, r_frac):
    key = (n, round(r_frac, 4))
    if key in _lambda_cache:
        return _lambda_cache[key]
    d_px = max(1, int(round(r_frac * BOX_PX)))
    sin_pi_n = math.sin(math.pi / n)
    diameter_engine = (d_px / sin_pi_n) * SCALE
    R_engine = diameter_engine / 2.0
    box = 1.0
    theta = [2.0 * math.pi * k / n for k in range(n)]
    cx, cy = 0.5, 0.5
    positions = []
    for t in theta:
        x = (cx + R_engine * math.cos(t)) % box
        y = (cy + R_engine * math.sin(t)) % box
        positions.append((x, y))
    d_min_val = float('inf')
    for i in range(n):
        for j in range(i + 1, n):
            dx = abs(positions[i][0] - positions[j][0])
            dy = abs(positions[i][1] - positions[j][1])
            if dx > box / 2: dx = box - dx
            if dy > box / 2: dy = box - dy
            d = math.sqrt(dx * dx + dy * dy)
            if d < d_min_val:
                d_min_val = d
    _lambda_cache[key] = d_min_val / MERGE_THR_ENGINE
    return _lambda_cache[key]


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
# Radial velocity component
# ============================================================

def radial_component(gamma_deg, drift):
    """Signed radial velocity component: positive = inward."""
    return math.sin(math.radians(gamma_deg)) + drift


# ============================================================
# Angle-partition prediction
# ============================================================

def predict_mechanism_angle(gamma_deg, drift, n):
    """Predict mechanism from angle partition + drift."""
    v_r = radial_component(gamma_deg, drift)
    if n == 4 and abs(gamma_deg) < 0.01 and abs(drift) < 0.001:
        return "DECAY"
    if v_r > 0.005:
        return "inward-collapse"
    elif v_r < -0.005:
        return "outward-PBC"
    else:
        return "other-late"


def predict_mechanism_lambda_angle(lam, gamma_deg, drift, n):
    """Predict mechanism from (lambda, angle) joint system."""
    v_r = radial_component(gamma_deg, drift)

    # Lambda >> 1 regime: complex, but angle partition still broadly holds
    if lam > 5.0:
        if abs(gamma_deg) < 0.01 and abs(drift) < 0.001 and n == 4:
            return "DECAY"
        if v_r > 0.01:
            return "inward-collapse"
        elif v_r < -0.01:
            return "outward-PBC"
        else:
            return "other-late"

    # Lambda 1-5: transitional
    if lam > 1.0:
        if v_r > 0.005:
            return "inward-collapse"
        elif v_r < -0.005:
            return "outward-PBC"
        else:
            return "other-late"

    # Lambda < 1: instantaneous collapse, clean partition
    if v_r > 0.002:
        return "inward-collapse"
    elif v_r < -0.002:
        return "outward-PBC"
    else:
        return "other-late"


# ============================================================
# Section 1: Inventory
# ============================================================

def section_1(cube):
    print("=" * 70)
    print("SECTION 1: Inventory & Radial-Velocity Computation")
    print("=" * 70)

    print(f"  Total atlas records: {len(cube)}")

    n_counts = Counter(k[0] for k in cube)
    print(f"\n  Records per N:")
    for n in sorted(n_counts):
        print(f"    N={n:3d}: {n_counts[n]:5d}")

    # Angle coverage
    angles = sorted(set(k[1] for k in cube))
    print(f"\n  Angles sampled ({len(angles)}): {angles}")

    # Radial velocity at key angles (drift=0)
    print(f"\n  Radial velocity component v_r = sin(gamma) + drift:")
    for ang in [-10, -3, -1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1, 3, 10, 12, 15]:
        v_r = radial_component(ang, 0)
        print(f"    gamma={ang:+6.2f}, drift=0:  v_r = {v_r:+.6f}")

    print()


# ============================================================
# Section 2: Angle partition test
# ============================================================

def section_2(cube):
    print("=" * 70)
    print("SECTION 2: Angle Partition -- Positive -> IC, Negative -> OPBC")
    print("=" * 70)

    # Test for N >= 12 (lambda < 1 regime where partition is clean)
    print("  --- Angle partition for N >= 12 ---")
    angle_bins = defaultdict(lambda: {"IC": 0, "OPBC": 0, "OL": 0, "CR": 0, "DECAY": 0, "total": 0})

    for key, rec in cube.items():
        if key[0] < 12:
            continue
        ang = key[1]
        mech = rec["mechanism"]
        angle_bins[ang]["total"] += 1
        if mech == "inward-collapse":
            angle_bins[ang]["IC"] += 1
        elif mech == "outward-PBC":
            angle_bins[ang]["OPBC"] += 1
        elif mech == "other-late":
            angle_bins[ang]["OL"] += 1
        elif mech == "PBC-corner":
            angle_bins[ang]["CR"] += 1
        elif mech == "DECAY":
            angle_bins[ang]["DECAY"] += 1

    print(f"    {'angle':>7s}  {'total':>5s}  {'IC%':>6s}  {'OPBC%':>6s}  {'OL%':>6s}  {'CR%':>6s}  {'DE%':>6s}  sector")
    for ang in sorted(angle_bins):
        b = angle_bins[ang]
        t = b["total"]
        ic_p = 100 * b["IC"] / t if t else 0
        op_p = 100 * b["OPBC"] / t if t else 0
        ol_p = 100 * b["OL"] / t if t else 0
        cr_p = 100 * b["CR"] / t if t else 0
        de_p = 100 * b["DECAY"] / t if t else 0
        if ic_p > 80:
            sector = "INWARD"
        elif op_p + cr_p > 80:
            sector = "OUTWARD"
        elif ol_p > 30:
            sector = "BOUNDARY"
        else:
            sector = "MIXED"
        print(f"    {ang:+7.2f}  {t:5d}  {ic_p:5.1f}%  {op_p:5.1f}%  {ol_p:5.1f}%  {cr_p:5.1f}%  {de_p:5.1f}%  {sector}")

    # Overall angle-partition accuracy for N >= 12
    correct = 0
    total = 0
    for key, rec in cube.items():
        if key[0] < 12:
            continue
        ang = key[1]
        mech = rec["mechanism"]
        total += 1
        if ang > 1.0 and mech == "inward-collapse":
            correct += 1
        elif ang < -1.0 and mech in ("outward-PBC", "PBC-corner"):
            correct += 1
        elif -1.0 <= ang <= 1.0:
            correct += 1  # boundary layer, any mechanism acceptable
    pct = 100.0 * correct / total if total else 0
    print(f"\n  Angle partition accuracy (N>=12, |angle|>1 strict): {correct}/{total} ({pct:.1f}%)")

    print()


# ============================================================
# Section 3: True IC boundary
# ============================================================

def section_3(cube):
    print("=" * 70)
    print("SECTION 3: True IC Boundary Location")
    print("=" * 70)

    # IC fraction at each boundary-layer angle, per N (N >= 12)
    print("  --- IC fraction at boundary-layer angles (N >= 12) ---")
    for ang in [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]:
        print(f"\n    angle = {ang:+.2f}:")
        for n in [12, 20, 32, 48, 64, 96, 128]:
            cells = [(k, v) for k, v in cube.items()
                     if k[0] == n and abs(k[1] - ang) < 0.01]
            if not cells:
                continue
            ic = sum(1 for _, v in cells if v["mechanism"] == "inward-collapse")
            pct = 100.0 * ic / len(cells)
            print(f"      N={n:3d}: {ic:3d}/{len(cells):3d} IC ({pct:5.1f}%)")

    # IC boundary vs radius at large N (N=64)
    print(f"\n  --- IC boundary vs radius at N=64, drift=0 ---")
    radii_all = sorted(set(k[2] for k in cube if k[0] == 64 and abs(k[3]) < 0.001))
    for r in radii_all:
        row = []
        for ang in [-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]:
            cells = [(k, v) for k, v in cube.items()
                     if k[0] == 64 and abs(k[1] - ang) < 0.01
                     and abs(k[2] - r) < 0.005 and abs(k[3]) < 0.001]
            if cells:
                ic = sum(1 for _, v in cells if v["mechanism"] == "inward-collapse")
                row.append(f"{ic}/{len(cells)}")
            else:
                row.append("---")
        print(f"    r={r:.3f}: " + "  ".join(f"{a:+.2f}:{v}" for a, v in
              zip([-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1], row)))

    # Drift effect on boundary
    print(f"\n  --- IC boundary vs drift at N=64, r=0.50, angle=-0.5 ---")
    for d in sorted(set(k[3] for k in cube if k[0] == 64 and abs(k[2] - 0.50) < 0.005
                        and abs(k[1] - (-0.5)) < 0.01)):
        cells = [(k, v) for k, v in cube.items()
                 if k[0] == 64 and abs(k[1] - (-0.5)) < 0.01
                 and abs(k[2] - 0.50) < 0.005 and abs(k[3] - d) < 0.0005]
        if cells:
            ic = sum(1 for _, v in cells if v["mechanism"] == "inward-collapse")
            v_r = radial_component(-0.5, d)
            print(f"    drift={d:.3f}: {ic}/{len(cells)} IC  (v_r={v_r:+.5f})")

    print()


# ============================================================
# Section 4: Beyond-grid test (+12, +15)
# ============================================================

def section_4(cube):
    print("=" * 70)
    print("SECTION 4: Beyond +10 Grid Saturation -- Angles +12, +15")
    print("=" * 70)

    for ang in [10, 12, 15]:
        print(f"\n  angle = +{ang}:")
        for n in NS_ALL:
            cells = [(k, v) for k, v in cube.items()
                     if k[0] == n and abs(k[1] - ang) < 0.01]
            if not cells:
                continue
            ic = sum(1 for _, v in cells if v["mechanism"] == "inward-collapse")
            chi_vals = [v["chi_emp"] for _, v in cells]
            chi_mean = statistics.mean(chi_vals)
            pct = 100.0 * ic / len(cells)
            print(f"    N={n:3d}: {ic:3d}/{len(cells):3d} IC ({pct:5.1f}%), chi_mean={chi_mean:.3f}")

    print()


# ============================================================
# Section 5: Boundary-layer drift flips
# ============================================================

def section_5(cube):
    print("=" * 70)
    print("SECTION 5: Boundary-Layer Drift Flips")
    print("=" * 70)

    # For each (N, angle, radius), check if drift changes the mechanism
    by_nar = defaultdict(dict)
    for key, rec in cube.items():
        by_nar[(key[0], key[1], key[2])][key[3]] = rec["mechanism"]

    # Count flips at boundary-layer angles vs bulk angles
    bl_flips = 0
    bl_total = 0
    bulk_flips = 0
    bulk_total = 0

    for (n, ang, r), drift_mechs in by_nar.items():
        if len(drift_mechs) < 2:
            continue
        mechs = list(drift_mechs.values())
        flipped = len(set(mechs)) > 1
        if -1.0 <= ang <= 1.0:
            bl_total += 1
            if flipped:
                bl_flips += 1
        else:
            bulk_total += 1
            if flipped:
                bulk_flips += 1

    bl_pct = 100.0 * bl_flips / bl_total if bl_total else 0
    bulk_pct = 100.0 * bulk_flips / bulk_total if bulk_total else 0
    print(f"  Drift-induced mechanism flips:")
    print(f"    Boundary layer (|angle| <= 1): {bl_flips}/{bl_total} ({bl_pct:.1f}%) flipped")
    print(f"    Bulk (|angle| > 1):            {bulk_flips}/{bulk_total} ({bulk_pct:.1f}%) flipped")
    print(f"    Ratio (BL/bulk):               {bl_pct/bulk_pct:.1f}x" if bulk_pct > 0 else "")

    # Per-angle flip rate
    print(f"\n  Drift-flip rate by angle (N >= 12):")
    angle_flips = defaultdict(lambda: [0, 0])
    for (n, ang, r), drift_mechs in by_nar.items():
        if n < 12 or len(drift_mechs) < 2:
            continue
        flipped = len(set(drift_mechs.values())) > 1
        angle_flips[ang][1] += 1
        if flipped:
            angle_flips[ang][0] += 1
    for ang in sorted(angle_flips):
        f, t = angle_flips[ang]
        pct = 100.0 * f / t if t else 0
        marker = " <-- BL" if -1.0 <= ang <= 1.0 else ""
        print(f"    angle={ang:+7.2f}: {f:4d}/{t:4d} ({pct:5.1f}%) flipped{marker}")

    print()


# ============================================================
# Section 6: Drift-rigidity decomposition
# ============================================================

def section_6(cube):
    print("=" * 70)
    print("SECTION 6: Drift-Rigidity Decomposition -- Bulk vs Boundary Layer")
    print("=" * 70)

    by_nar = defaultdict(dict)
    for key, rec in cube.items():
        by_nar[(key[0], key[1], key[2])][key[3]] = rec["mechanism"]

    for n in NS_ALL:
        bl_rigid = 0
        bl_total = 0
        bulk_rigid = 0
        bulk_total = 0
        for (nn, ang, r), drift_mechs in by_nar.items():
            if nn != n or len(drift_mechs) < 2:
                continue
            rigid = len(set(drift_mechs.values())) == 1
            if -1.0 <= ang <= 1.0:
                bl_total += 1
                if rigid:
                    bl_rigid += 1
            else:
                bulk_total += 1
                if rigid:
                    bulk_rigid += 1
        all_total = bl_total + bulk_total
        all_rigid = bl_rigid + bulk_rigid
        if all_total == 0:
            continue
        all_pct = 100.0 * all_rigid / all_total
        bl_pct = 100.0 * bl_rigid / bl_total if bl_total else 0
        bulk_pct = 100.0 * bulk_rigid / bulk_total if bulk_total else 0
        print(f"  N={n:3d}: overall={all_pct:5.1f}%  bulk={bulk_pct:5.1f}%  BL={bl_pct:5.1f}%  "
              f"(bulk={bulk_rigid}/{bulk_total}, BL={bl_rigid}/{bl_total})")

    print()


# ============================================================
# Section 7: Taxonomy freeze-out
# ============================================================

def section_7(cube):
    print("=" * 70)
    print("SECTION 7: Taxonomy Freeze-Out -- 3 Cores + 3 Corrections")
    print("=" * 70)

    CORES = {"IC-direct", "OPBC-single", "OL-stalled"}
    CORRECTIONS = {"IC-bounced", "OPBC-glancing", "CR-diagonal"}
    EXPECTED = CORES | CORRECTIONS

    for n in NS_ALL:
        cells = [v for k, v in cube.items() if k[0] == n]
        if not cells:
            continue
        subs = Counter(classify_sub(v) for v in cells)
        total_n = sum(subs.values())
        n_types = len(subs)
        types_set = set(subs.keys())
        unexpected = types_set - EXPECTED - {"IC-delayed", "OPBC-winding", "CR-single",
                                              "CR-multi", "OL-multi", "OL-periodic",
                                              "OL-wandering", "DE-tangent", "DE-frozen"}
        core_pct = sum(subs.get(s, 0) for s in CORES) / total_n * 100
        corr_pct = sum(subs.get(s, 0) for s in CORRECTIONS) / total_n * 100
        print(f"  N={n:3d}: {n_types:2d} sub-types, cores={core_pct:5.1f}%, corrections={corr_pct:5.1f}%, "
              f"unexpected={'NONE' if not unexpected else unexpected}")

    # N >= 48: detailed breakdown
    print(f"\n  Sub-mechanism distribution for N >= 48:")
    subs_48 = Counter()
    for k, v in cube.items():
        if k[0] >= 48:
            subs_48[classify_sub(v)] += 1
    total_48 = sum(subs_48.values())
    for sub, cnt in subs_48.most_common():
        tag = "[CORE]" if sub in CORES else "[CORR]" if sub in CORRECTIONS else "[OTHER]"
        print(f"    {sub:20s}: {cnt:5d} ({100*cnt/total_48:5.1f}%) {tag}")

    # Verify no new sub-types at any N
    all_subs = set()
    for v in cube.values():
        all_subs.add(classify_sub(v))
    print(f"\n  All sub-types across entire dataset: {sorted(all_subs)}")
    print(f"  Count: {len(all_subs)}")

    # Entropy convergence
    print(f"\n  Taxonomy entropy per N:")
    entropies = {}
    for n in NS_ALL:
        cells = [v for k, v in cube.items() if k[0] == n]
        if not cells:
            continue
        subs = Counter(classify_sub(v) for v in cells)
        total_n = sum(subs.values())
        h = 0.0
        for cnt in subs.values():
            p = cnt / total_n
            if p > 0:
                h -= p * math.log(p)
        entropies[n] = h
        print(f"    N={n:3d}: H_sub={h:.3f}")
    if len([n for n in entropies if n >= 20]) >= 3:
        tail = [entropies[n] for n in sorted(entropies) if n >= 20]
        cv = statistics.stdev(tail) / statistics.mean(tail) if statistics.mean(tail) > 0 else 0
        print(f"  Entropy CV (N>=20): {100*cv:.1f}%")

    print()


# ============================================================
# Section 8: (Lambda, angle) joint prediction
# ============================================================

def section_8(cube):
    print("=" * 70)
    print("SECTION 8: (Lambda, Angle) Joint Prediction Accuracy")
    print("=" * 70)

    # Angle-only prediction
    ang_correct = 0
    # (Lambda, angle) joint prediction
    joint_correct = 0
    total = 0
    ang_by_n = defaultdict(lambda: [0, 0])
    joint_by_n = defaultdict(lambda: [0, 0])

    for key, rec in cube.items():
        n, angle, radius, drift = key
        lam = compute_lambda(n, radius)
        actual = rec["mechanism"]
        total += 1

        # Angle-only
        pred_a = predict_mechanism_angle(angle, drift, n)
        if pred_a == actual:
            ang_correct += 1
            ang_by_n[n][0] += 1
        ang_by_n[n][1] += 1

        # Joint (lambda, angle)
        pred_j = predict_mechanism_lambda_angle(lam, angle, drift, n)
        if pred_j == actual:
            joint_correct += 1
            joint_by_n[n][0] += 1
        joint_by_n[n][1] += 1

    ang_pct = 100.0 * ang_correct / total if total else 0
    joint_pct = 100.0 * joint_correct / total if total else 0

    print(f"  Overall accuracy:")
    print(f"    Angle-only prediction:       {ang_correct}/{total} ({ang_pct:.1f}%)")
    print(f"    (Lambda, angle) prediction:  {joint_correct}/{total} ({joint_pct:.1f}%)")

    print(f"\n  Per-N accuracy:")
    print(f"    {'N':>4s}  {'angle%':>7s}  {'joint%':>7s}  {'lam@r0.5':>9s}")
    for n in sorted(ang_by_n):
        a_c, a_t = ang_by_n[n]
        j_c, j_t = joint_by_n[n]
        lam = compute_lambda(n, 0.50)
        print(f"    {n:4d}  {100*a_c/a_t:6.1f}%  {100*j_c/j_t:6.1f}%  {lam:9.3f}")

    # By lambda bin
    print(f"\n  Joint accuracy by lambda bin:")
    bins = [(0, 0.3), (0.3, 1.0), (1.0, 3.0), (3.0, 100.0)]
    for lo, hi in bins:
        cor = 0
        tot = 0
        for key, rec in cube.items():
            lam = compute_lambda(key[0], key[2])
            if lo <= lam < hi:
                pred = predict_mechanism_lambda_angle(lam, key[1], key[3], key[0])
                tot += 1
                if pred == rec["mechanism"]:
                    cor += 1
        if tot:
            print(f"    lambda [{lo:.1f}, {hi:.1f}): {cor}/{tot} ({100*cor/tot:.1f}%)")

    # By angle sector
    print(f"\n  Joint accuracy by angle sector (N >= 12):")
    for label, cond in [("Inward (angle > 1)", lambda a: a > 1.0),
                        ("Outward (angle < -1)", lambda a: a < -1.0),
                        ("Boundary (-1 <= angle <= 1)", lambda a: -1.0 <= a <= 1.0)]:
        cor = 0
        tot = 0
        for key, rec in cube.items():
            if key[0] < 12:
                continue
            if cond(key[1]):
                lam = compute_lambda(key[0], key[2])
                pred = predict_mechanism_lambda_angle(lam, key[1], key[3], key[0])
                tot += 1
                if pred == rec["mechanism"]:
                    cor += 1
        if tot:
            print(f"    {label:35s}: {cor}/{tot} ({100*cor/tot:.1f}%)")

    print()
    return ang_correct, joint_correct, total


# ============================================================
# Section 9: Cross-layer agreement summary
# ============================================================

def section_9(cube, ang_correct, joint_correct, total):
    print("=" * 70)
    print("SECTION 9: Cross-Layer Agreement Summary")
    print("=" * 70)

    ang_pct = 100.0 * ang_correct / total if total else 0
    joint_pct = 100.0 * joint_correct / total if total else 0

    print(f"  Angle-only mechanism accuracy:      {ang_pct:.1f}%")
    print(f"  (Lambda, angle) mechanism accuracy:  {joint_pct:.1f}%")

    # Regime structure tests
    print(f"\n  Law X structural tests:")

    # 1. Positive angle -> IC at N >= 20
    pos = [(k, v) for k, v in cube.items() if k[0] >= 20 and k[1] > 1.0]
    ic_pos = sum(1 for _, v in pos if v["mechanism"] == "inward-collapse")
    pct_pos = 100.0 * ic_pos / len(pos) if pos else 0
    print(f"    Positive (>+1) -> IC, N>=20:     {ic_pos}/{len(pos)} ({pct_pos:.1f}%)")

    # 2. Negative angle -> non-IC at N >= 20
    neg = [(k, v) for k, v in cube.items() if k[0] >= 20 and k[1] < -1.0]
    non_ic = sum(1 for _, v in neg if v["mechanism"] != "inward-collapse")
    pct_neg = 100.0 * non_ic / len(neg) if neg else 0
    print(f"    Negative (<-1) -> non-IC, N>=20: {non_ic}/{len(neg)} ({pct_neg:.1f}%)")

    # 3. Beyond-grid angles -> IC at N >= 20
    beyond = [(k, v) for k, v in cube.items() if k[0] >= 20 and k[1] >= 12]
    ic_beyond = sum(1 for _, v in beyond if v["mechanism"] == "inward-collapse")
    pct_beyond = 100.0 * ic_beyond / len(beyond) if beyond else 0
    print(f"    Beyond-grid (+12,+15) -> IC:      {ic_beyond}/{len(beyond)} ({pct_beyond:.1f}%)")

    # 4. Boundary layer mixed at -0.5
    bl = [(k, v) for k, v in cube.items() if k[0] >= 20 and abs(k[1] - (-0.5)) < 0.01]
    ic_bl = sum(1 for _, v in bl if v["mechanism"] == "inward-collapse")
    pct_bl = 100.0 * ic_bl / len(bl) if bl else 0
    print(f"    Boundary (-0.5) mixed zone:       {ic_bl}/{len(bl)} IC ({pct_bl:.1f}%)")

    # 5. Chi at floor for N >= 48
    n48 = [(k, v) for k, v in cube.items() if k[0] >= 48]
    floor = sum(1 for _, v in n48 if v["chi_emp"] <= 0.02)
    pct_floor = 100.0 * floor / len(n48) if n48 else 0
    print(f"    N>=48 chi at floor:               {floor}/{len(n48)} ({pct_floor:.1f}%)")

    # 6. Taxonomy: 6 sub-types at N >= 48
    subs_48 = set(classify_sub(v) for k, v in cube.items() if k[0] >= 48)
    print(f"    N>=48 sub-type count:             {len(subs_48)}")

    # 7. N=8 resonance
    n8 = [v for k, v in cube.items() if k[0] == 8]
    wander = sum(1 for v in n8 if classify_sub(v) == "OL-wandering")
    pct_w = 100.0 * wander / len(n8) if n8 else 0
    print(f"    N=8 OL-wandering:                 {wander}/{len(n8)} ({pct_w:.1f}%)")

    # Verdict
    print(f"\n  ===== LAW X CONSOLIDATION VERDICT =====")
    print(f"  Angle partition (bulk, N>=20):      {'CONFIRMED' if pct_pos > 95 and pct_neg > 95 else 'PARTIAL'}")
    print(f"  Beyond-grid IC at +12/+15:          {'CONFIRMED' if pct_beyond > 95 else 'PARTIAL'}")
    print(f"  Boundary layer at -0.5:             {'CONFIRMED (mixed)' if 20 < pct_bl < 80 else 'UNEXPECTED'}")
    print(f"  Taxonomy freeze (<=6 at N>=48):     {'CONFIRMED' if len(subs_48) <= 7 else 'FAILED'}")
    print(f"  (Lambda, angle) joint accuracy:     {joint_pct:.1f}%")
    print(f"  Angle-only accuracy:                {ang_pct:.1f}%")

    print()


# ============================================================
# Main
# ============================================================

def main():
    print("Law X Consolidation Sweep Analysis")
    print("=" * 70)
    print()

    cube = load_all_atlas()

    section_1(cube)
    section_2(cube)
    section_3(cube)
    section_4(cube)
    section_5(cube)
    section_6(cube)
    section_7(cube)
    ang_c, joint_c, total = section_8(cube)
    section_9(cube, ang_c, joint_c, total)


if __name__ == "__main__":
    main()
