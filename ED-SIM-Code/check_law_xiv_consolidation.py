#!/usr/bin/env python
"""
check_law_xiv_consolidation.py
================================
Law XIV Consolidation Sweep: tests the collapse-time scaling law
chi ~ lambda^p, the temporal hierarchy, temporal invisibility of the
boundary band, and radius-modulation of the temporal exponent.

Sections:
  1.  Inventory & lambda computation
  2.  Canonical-grid exponent (9 canonical radii)
  3.  Extended-grid exponent (all radii)
  4.  Radius-dependent exponent decomposition
  5.  Universal chi-floor verification
  6.  Floor-onset by radius
  7.  Temporal invisibility of the boundary band
  8.  Strict temporal hierarchy verification
  9.  Temporal hierarchy invariance (angle, radius, drift)
 10.  Cross-layer summary & verdict
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

RADII_CANONICAL = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 0.90, 0.95]


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
# Power-law fit: log(chi) = p * log(lambda) + c
# ============================================================

def fit_power_law(lambdas, chis):
    """Fit chi ~ lambda^p using log-log OLS on cells with chi > 0.01."""
    xs = []
    ys = []
    for lam, chi in zip(lambdas, chis):
        if lam > 0.05 and chi > 0.011:  # above floor, above noise
            xs.append(math.log(lam))
            ys.append(math.log(chi))
    if len(xs) < 3:
        return None, None, 0
    n = len(xs)
    sx = sum(xs)
    sy = sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-15:
        return None, None, 0
    p = (n * sxy - sx * sy) / denom
    c = (sy - p * sx) / n
    # R^2
    y_mean = sy / n
    ss_tot = sum((y - y_mean) ** 2 for y in ys)
    ss_res = sum((y - (p * x + c)) ** 2 for x, y in zip(xs, ys))
    r2 = 1 - ss_res / ss_tot if ss_tot > 1e-15 else 0
    return p, c, r2


# ============================================================
# Section 1: Inventory
# ============================================================

def section_1(cube):
    print("=" * 70)
    print("SECTION 1: Inventory & Lambda Computation")
    print("=" * 70)

    print(f"  Total atlas records: {len(cube)}")

    n_counts = Counter(k[0] for k in cube)
    print(f"\n  Records per N:")
    for n in sorted(n_counts):
        print(f"    N={n:3d}: {n_counts[n]:5d}")

    radii = sorted(set(k[2] for k in cube))
    print(f"\n  Radii sampled ({len(radii)}): {[round(r, 3) for r in radii]}")

    angles = sorted(set(k[1] for k in cube))
    print(f"  Angles sampled ({len(angles)}): {angles}")

    drifts = sorted(set(k[3] for k in cube))
    print(f"  Drifts sampled ({len(drifts)}): {drifts}")

    # Lambda values at r=0.50
    print(f"\n  Lambda at r=0.50:")
    for n in NS_ALL:
        lam = compute_lambda(n, 0.50)
        regime = "folded" if lam > 3 else "resonant" if lam > 0.8 else "planar"
        print(f"    N={n:3d}: lambda={lam:.3f}  ({regime})")

    print()


# ============================================================
# Section 2: Canonical-grid exponent
# ============================================================

def section_2(cube):
    print("=" * 70)
    print("SECTION 2: Canonical-Grid Exponent (9 canonical radii)")
    print("=" * 70)

    canonical_set = set(round(r, 3) for r in RADII_CANONICAL)

    # Per-cell lambda using actual radius (not reference r=0.50)
    all_lambdas = []
    all_chis = []
    chi_by_n = defaultdict(list)
    lam_by_n = defaultdict(list)
    for key, rec in cube.items():
        r_rounded = round(key[2], 3)
        if r_rounded not in canonical_set:
            continue
        n = key[0]
        lam = compute_lambda(n, key[2])
        chi = rec["chi_emp"]
        all_lambdas.append(lam)
        all_chis.append(chi)
        chi_by_n[n].append(chi)
        lam_by_n[n].append(lam)

    # Print per-N summary (mean lambda and mean chi)
    print(f"  Per-N summary (canonical radii, per-cell lambda):")
    ns_sorted = sorted(chi_by_n.keys())
    for n in ns_sorted:
        mc = statistics.mean(chi_by_n[n])
        ml = statistics.mean(lam_by_n[n])
        print(f"    N={n:3d}: mean_lambda={ml:.3f}, mean_chi={mc:.4f}, n_cells={len(chi_by_n[n])}")

    # Fit on per-cell (lambda, chi) scatter
    p, c, r2 = fit_power_law(all_lambdas, all_chis)
    print(f"\n  Per-cell power-law fit (chi ~ lambda^p, above floor):")
    print(f"    Data points above floor: {sum(1 for l, ch in zip(all_lambdas, all_chis) if l > 0.05 and ch > 0.011)}")
    if p is not None:
        print(f"    p_lambda = {p:.2f}")
        print(f"    R^2 = {r2:.3f}")
    else:
        print(f"    Insufficient data for lambda fit")

    # N-scaling fit (chi ~ N^alpha) for comparison with Law IX convention
    ns_fit = []
    chis_fit = []
    for n in ns_sorted:
        mc = statistics.mean(chi_by_n[n])
        if mc > 0.011:
            ns_fit.append(n)
            chis_fit.append(mc)
    alpha, c_n, r2_n = fit_power_law(ns_fit, chis_fit)
    print(f"\n  N-scaling fit (chi ~ N^alpha, mean chi > floor):")
    print(f"    N values used: {ns_fit}")
    if alpha is not None:
        print(f"    alpha = {alpha:.2f}")
        print(f"    R^2 = {r2_n:.3f}")
        print(f"    Target: alpha ~ -2.33, R^2 >= 0.97")
        print(f"    Match: |alpha - (-2.33)| = {abs(alpha - (-2.33)):.2f}")
    else:
        print(f"    Insufficient data for N-scaling fit")

    print()
    return alpha, r2_n


# ============================================================
# Section 3: Extended-grid exponent
# ============================================================

def section_3(cube):
    print("=" * 70)
    print("SECTION 3: Extended-Grid Exponent (all radii)")
    print("=" * 70)

    # Per-cell lambda using actual radius
    all_lambdas = []
    all_chis = []
    chi_by_n = defaultdict(list)
    lam_by_n = defaultdict(list)
    for key, rec in cube.items():
        n = key[0]
        lam = compute_lambda(n, key[2])
        chi = rec["chi_emp"]
        all_lambdas.append(lam)
        all_chis.append(chi)
        chi_by_n[n].append(chi)
        lam_by_n[n].append(lam)

    print(f"  Per-N summary (all radii, per-cell lambda):")
    ns_sorted = sorted(chi_by_n.keys())
    for n in ns_sorted:
        mc = statistics.mean(chi_by_n[n])
        ml = statistics.mean(lam_by_n[n])
        print(f"    N={n:3d}: mean_lambda={ml:.3f}, mean_chi={mc:.4f}, n_cells={len(chi_by_n[n])}")

    # Fit on per-cell (lambda, chi) scatter
    p, c, r2 = fit_power_law(all_lambdas, all_chis)
    print(f"\n  Per-cell power-law fit (chi ~ lambda^p, above floor):")
    print(f"    Data points above floor: {sum(1 for l, ch in zip(all_lambdas, all_chis) if l > 0.05 and ch > 0.011)}")
    if p is not None:
        print(f"    p_lambda = {p:.2f}")
        print(f"    R^2 = {r2:.3f}")
    else:
        print(f"    Insufficient data for lambda fit")

    # N-scaling fit (chi ~ N^alpha) for comparison with Law IX convention
    ns_fit = []
    chis_fit = []
    for n in ns_sorted:
        mc = statistics.mean(chi_by_n[n])
        if mc > 0.011:
            ns_fit.append(n)
            chis_fit.append(mc)
    alpha, c_n, r2_n = fit_power_law(ns_fit, chis_fit)
    print(f"\n  N-scaling fit (chi ~ N^alpha, mean chi > floor):")
    print(f"    N values used: {ns_fit}")
    if alpha is not None:
        print(f"    alpha = {alpha:.2f}")
        print(f"    R^2 = {r2_n:.3f}")
        print(f"    Target: alpha ~ -2.81, R^2 >= 0.95")
        print(f"    Match: |alpha - (-2.81)| = {abs(alpha - (-2.81)):.2f}")
    else:
        print(f"    Insufficient data for N-scaling fit")

    print()
    return alpha, r2_n


# ============================================================
# Section 4: Radius-dependent exponent decomposition
# ============================================================

def section_4(cube):
    print("=" * 70)
    print("SECTION 4: Radius-Dependent Exponent Decomposition")
    print("=" * 70)

    # Group radii into bands
    all_radii = sorted(set(round(k[2], 3) for k in cube))

    # Per-radius mean chi at each N
    print(f"  Per-radius temporal exponent:")
    print(f"    {'radius':>7s}  {'p':>7s}  {'R^2':>6s}  {'chi@N=4':>8s}  {'chi@N=48':>8s}  {'class':>10s}")

    for r in all_radii:
        chi_by_n = defaultdict(list)
        for key, rec in cube.items():
            if abs(round(key[2], 3) - r) > 0.002:
                continue
            chi_by_n[key[0]].append(rec["chi_emp"])

        if len(chi_by_n) < 3:
            continue

        lambdas = []
        mean_chis = []
        chi_n4 = statistics.mean(chi_by_n.get(4, [0.01])) if 4 in chi_by_n else None
        chi_n48 = statistics.mean(chi_by_n.get(48, [0.01])) if 48 in chi_by_n else None
        for n in sorted(chi_by_n.keys()):
            lam = compute_lambda(n, r)
            mc = statistics.mean(chi_by_n[n])
            lambdas.append(lam)
            mean_chis.append(mc)

        p, c, r2 = fit_power_law(lambdas, mean_chis)
        if p is None:
            r_class = "floor-only"
            p_str = "   ---"
            r2_str = "  ---"
        else:
            r_class = "short-path" if abs(p) < 2.0 else "canonical" if abs(p) < 2.7 else "steep"
            p_str = f"{p:7.2f}"
            r2_str = f"{r2:6.3f}"

        chi4_str = f"{chi_n4:8.3f}" if chi_n4 is not None else "     ---"
        chi48_str = f"{chi_n48:8.3f}" if chi_n48 is not None else "     ---"
        print(f"    {r:7.3f}  {p_str}  {r2_str}  {chi4_str}  {chi48_str}  {r_class:>10s}")

    # Summary by radius class
    small_radii = [r for r in all_radii if r <= 0.15]
    mid_radii = [r for r in all_radii if 0.15 < r <= 0.60]
    large_radii = [r for r in all_radii if r > 0.60]

    print(f"\n  Radius bands:")
    print(f"    Small (r <= 0.15): {len(small_radii)} radii")
    print(f"    Mid   (0.15 < r <= 0.60): {len(mid_radii)} radii")
    print(f"    Large (r > 0.60): {len(large_radii)} radii")

    for label, rset in [("Small", small_radii), ("Mid", mid_radii), ("Large", large_radii)]:
        chi_by_n = defaultdict(list)
        for key, rec in cube.items():
            r_rounded = round(key[2], 3)
            if r_rounded in rset:
                chi_by_n[key[0]].append(rec["chi_emp"])
        lambdas = []
        mean_chis = []
        for n in sorted(chi_by_n.keys()):
            lam = compute_lambda(n, 0.50)
            mc = statistics.mean(chi_by_n[n])
            lambdas.append(lam)
            mean_chis.append(mc)
        p, c, r2 = fit_power_law(lambdas, mean_chis)
        if p is not None:
            print(f"    {label:5s} band: p = {p:.2f}, R^2 = {r2:.3f}")
        else:
            print(f"    {label:5s} band: insufficient above-floor data")

    print()


# ============================================================
# Section 5: Universal chi-floor
# ============================================================

def section_5(cube):
    print("=" * 70)
    print("SECTION 5: Universal Chi-Floor Verification")
    print("=" * 70)

    CHI_FLOOR = 0.02  # allow small tolerance

    for n in NS_ALL:
        cells = [(k, v) for k, v in cube.items() if k[0] == n]
        if not cells:
            continue
        at_floor = sum(1 for _, v in cells if v["chi_emp"] <= CHI_FLOOR)
        above = sum(1 for _, v in cells if v["chi_emp"] > CHI_FLOOR)
        chi_vals = [v["chi_emp"] for _, v in cells]
        chi_max = max(chi_vals)
        chi_mean = statistics.mean(chi_vals)
        pct = 100.0 * at_floor / len(cells) if cells else 0
        lam = compute_lambda(n, 0.50)
        print(f"    N={n:3d} (lambda={lam:.3f}): {pct:5.1f}% at floor, "
              f"max={chi_max:.3f}, mean={chi_mean:.4f}, above={above}")

    # Chi at floor vs lambda bins
    print(f"\n  Floor fraction by lambda bin:")
    bins = [(0, 0.05), (0.05, 0.1), (0.1, 0.3), (0.3, 1.0), (1.0, 3.0), (3.0, 100)]
    for lo, hi in bins:
        cells = [(k, v) for k, v in cube.items()
                 if lo <= compute_lambda(k[0], k[2]) < hi]
        if not cells:
            continue
        at_floor = sum(1 for _, v in cells if v["chi_emp"] <= CHI_FLOOR)
        pct = 100.0 * at_floor / len(cells)
        print(f"    lambda [{lo:.2f}, {hi:.2f}): {at_floor}/{len(cells)} ({pct:.1f}%) at floor")

    print()


# ============================================================
# Section 6: Floor-onset by radius
# ============================================================

def section_6(cube):
    print("=" * 70)
    print("SECTION 6: Chi-Floor Onset by Radius")
    print("=" * 70)

    CHI_FLOOR = 0.02
    all_radii = sorted(set(round(k[2], 3) for k in cube))

    print(f"  Floor onset N (first N where >90% of cells at floor):")
    print(f"    {'radius':>7s}  {'onset_N':>7s}  {'lambda@onset':>12s}")

    for r in all_radii:
        onset_n = None
        for n in NS_ALL:
            cells = [(k, v) for k, v in cube.items()
                     if k[0] == n and abs(round(k[2], 3) - r) < 0.002]
            if not cells:
                continue
            at_floor = sum(1 for _, v in cells if v["chi_emp"] <= CHI_FLOOR)
            pct = 100.0 * at_floor / len(cells)
            if pct >= 90.0:
                onset_n = n
                break
        if onset_n is not None:
            lam = compute_lambda(onset_n, r)
            print(f"    {r:7.3f}  N={onset_n:3d}     lambda={lam:.3f}")
        else:
            print(f"    {r:7.3f}  >128       ---")

    print()


# ============================================================
# Section 7: Temporal invisibility of boundary band
# ============================================================

def section_7(cube):
    print("=" * 70)
    print("SECTION 7: Temporal Invisibility of Boundary Band")
    print("=" * 70)

    # For N >= 48, compare chi across angle sectors
    print(f"  --- Mean chi by angle sector (N >= 48) ---")
    sectors = {
        "Inward (angle > 1)": [],
        "Outward (angle < -1)": [],
        "Boundary (-1 <= a <= 1)": [],
    }
    for key, rec in cube.items():
        if key[0] < 48:
            continue
        ang = key[1]
        chi = rec["chi_emp"]
        if ang > 1.0:
            sectors["Inward (angle > 1)"].append(chi)
        elif ang < -1.0:
            sectors["Outward (angle < -1)"].append(chi)
        else:
            sectors["Boundary (-1 <= a <= 1)"].append(chi)

    for label, vals in sectors.items():
        if vals:
            m = statistics.mean(vals)
            mx = max(vals)
            at_floor = sum(1 for v in vals if v <= 0.02)
            pct = 100.0 * at_floor / len(vals)
            print(f"    {label:30s}: mean={m:.4f}, max={mx:.3f}, "
                  f"floor={pct:.1f}%, n={len(vals)}")

    # Chi by mechanism type at N >= 48
    print(f"\n  --- Mean chi by mechanism (N >= 48) ---")
    mech_chi = defaultdict(list)
    for key, rec in cube.items():
        if key[0] < 48:
            continue
        mech_chi[rec["mechanism"]].append(rec["chi_emp"])
    for mech in sorted(mech_chi):
        vals = mech_chi[mech]
        m = statistics.mean(vals)
        mx = max(vals)
        print(f"    {mech:20s}: mean={m:.4f}, max={mx:.3f}, n={len(vals)}")

    # Chi by drift at N >= 48
    print(f"\n  --- Mean chi by drift (N >= 48) ---")
    drift_chi = defaultdict(list)
    for key, rec in cube.items():
        if key[0] < 48:
            continue
        drift_chi[key[3]].append(rec["chi_emp"])
    for d in sorted(drift_chi):
        vals = drift_chi[d]
        m = statistics.mean(vals)
        print(f"    drift={d:.3f}: mean={m:.4f}, n={len(vals)}")

    # Chi across radius bands at N >= 48
    print(f"\n  --- Mean chi by radius band (N >= 48) ---")
    for label, lo, hi in [("Small (r<=0.15)", 0, 0.155),
                           ("Mid (0.15<r<=0.60)", 0.155, 0.605),
                           ("Large (r>0.60)", 0.605, 1.0)]:
        vals = [v["chi_emp"] for k, v in cube.items()
                if k[0] >= 48 and lo <= k[2] < hi]
        if vals:
            m = statistics.mean(vals)
            mx = max(vals)
            print(f"    {label:25s}: mean={m:.4f}, max={mx:.3f}, n={len(vals)}")

    print()


# ============================================================
# Section 8: Strict temporal hierarchy
# ============================================================

def section_8(cube):
    print("=" * 70)
    print("SECTION 8: Strict Temporal Hierarchy Verification")
    print("=" * 70)

    # Mean chi per regime
    regime_chi = {"folded": [], "resonant": [], "planar": []}
    for key, rec in cube.items():
        n = key[0]
        lam = compute_lambda(n, key[2])
        chi = rec["chi_emp"]
        if lam > 3.0:
            regime_chi["folded"].append(chi)
        elif lam > 0.8:
            regime_chi["resonant"].append(chi)
        else:
            regime_chi["planar"].append(chi)

    print(f"  Regime statistics:")
    for regime in ["folded", "resonant", "planar"]:
        vals = regime_chi[regime]
        if vals:
            m = statistics.mean(vals)
            mx = max(vals)
            mn = min(vals)
            print(f"    {regime:10s}: mean={m:.4f}, min={mn:.3f}, max={mx:.3f}, n={len(vals)}")

    # Verify strict ordering
    folded_min = min(regime_chi["folded"]) if regime_chi["folded"] else 0
    resonant_max = max(regime_chi["resonant"]) if regime_chi["resonant"] else 0
    resonant_min = min(regime_chi["resonant"]) if regime_chi["resonant"] else 0
    planar_max = max(regime_chi["planar"]) if regime_chi["planar"] else 0

    # Monotone check: for every cell, chi(N_small) >= chi(N_large)?
    print(f"\n  Monotone hierarchy check (per-cell, same angle/radius/drift):")
    violations = 0
    total_pairs = 0
    violation_examples = []
    by_ard = defaultdict(dict)
    for key, rec in cube.items():
        by_ard[(key[1], key[2], key[3])][key[0]] = rec["chi_emp"]

    for (ang, r, d), n_chi in by_ard.items():
        ns = sorted(n_chi.keys())
        for i in range(len(ns)):
            for j in range(i + 1, len(ns)):
                n1, n2 = ns[i], ns[j]
                # n1 < n2, so lambda(n1) > lambda(n2), expect chi(n1) >= chi(n2)
                total_pairs += 1
                if n_chi[n1] < n_chi[n2] - 0.005:  # allow tiny tolerance
                    violations += 1
                    if len(violation_examples) < 5:
                        violation_examples.append(
                            f"    ang={ang}, r={r:.3f}, d={d:.3f}: "
                            f"chi(N={n1})={n_chi[n1]:.3f} < chi(N={n2})={n_chi[n2]:.3f}")

    pct_v = 100.0 * violations / total_pairs if total_pairs else 0
    print(f"    Total N-pairs tested: {total_pairs}")
    print(f"    Violations: {violations} ({pct_v:.2f}%)")
    if violation_examples:
        print(f"    Examples:")
        for ex in violation_examples:
            print(ex)

    print()
    return violations, total_pairs


# ============================================================
# Section 9: Temporal hierarchy invariance
# ============================================================

def section_9(cube):
    print("=" * 70)
    print("SECTION 9: Temporal Hierarchy Invariance")
    print("=" * 70)

    # Check: at every angle, chi(N=4) > chi(N=8) > chi(N=48)
    print(f"  --- Hierarchy holds at every angle? ---")
    angles = sorted(set(k[1] for k in cube))
    angle_fails = 0
    for ang in angles:
        chi_4 = [v["chi_emp"] for k, v in cube.items() if k[0] == 4 and k[1] == ang]
        chi_8 = [v["chi_emp"] for k, v in cube.items() if k[0] == 8 and k[1] == ang]
        chi_48 = [v["chi_emp"] for k, v in cube.items() if k[0] == 48 and k[1] == ang]
        if not chi_4 or not chi_8 or not chi_48:
            continue
        m4 = statistics.mean(chi_4)
        m8 = statistics.mean(chi_8)
        m48 = statistics.mean(chi_48)
        ok = m4 >= m8 >= m48
        if not ok:
            angle_fails += 1
            print(f"    FAIL angle={ang:+7.2f}: chi(4)={m4:.3f}, chi(8)={m8:.3f}, chi(48)={m48:.3f}")
    if angle_fails == 0:
        print(f"    ALL {len(angles)} angles: chi(N=4) >= chi(N=8) >= chi(N=48)  [CONFIRMED]")
    else:
        print(f"    {angle_fails}/{len(angles)} angles violate hierarchy")

    # Check: at every radius
    print(f"\n  --- Hierarchy holds at every radius? ---")
    radii = sorted(set(round(k[2], 3) for k in cube))
    radius_fails = 0
    for r in radii:
        chi_4 = [v["chi_emp"] for k, v in cube.items()
                 if k[0] == 4 and abs(round(k[2], 3) - r) < 0.002]
        chi_8 = [v["chi_emp"] for k, v in cube.items()
                 if k[0] == 8 and abs(round(k[2], 3) - r) < 0.002]
        chi_48 = [v["chi_emp"] for k, v in cube.items()
                  if k[0] == 48 and abs(round(k[2], 3) - r) < 0.002]
        if not chi_4 or not chi_8 or not chi_48:
            continue
        m4 = statistics.mean(chi_4)
        m8 = statistics.mean(chi_8)
        m48 = statistics.mean(chi_48)
        ok = m4 >= m8 >= m48
        if not ok:
            radius_fails += 1
            print(f"    FAIL r={r:.3f}: chi(4)={m4:.3f}, chi(8)={m8:.3f}, chi(48)={m48:.3f}")
    if radius_fails == 0:
        print(f"    ALL {len(radii)} radii: chi(N=4) >= chi(N=8) >= chi(N=48)  [CONFIRMED]")
    else:
        print(f"    {radius_fails}/{len(radii)} radii violate hierarchy")

    # Check: at every drift
    print(f"\n  --- Hierarchy holds at every drift? ---")
    drifts = sorted(set(k[3] for k in cube))
    drift_fails = 0
    for d in drifts:
        chi_4 = [v["chi_emp"] for k, v in cube.items() if k[0] == 4 and k[3] == d]
        chi_8 = [v["chi_emp"] for k, v in cube.items() if k[0] == 8 and k[3] == d]
        chi_48 = [v["chi_emp"] for k, v in cube.items() if k[0] == 48 and k[3] == d]
        if not chi_4 or not chi_8 or not chi_48:
            continue
        m4 = statistics.mean(chi_4)
        m8 = statistics.mean(chi_8)
        m48 = statistics.mean(chi_48)
        ok = m4 >= m8 >= m48
        if not ok:
            drift_fails += 1
            print(f"    FAIL d={d:.3f}: chi(4)={m4:.3f}, chi(8)={m8:.3f}, chi(48)={m48:.3f}")
    if drift_fails == 0:
        print(f"    ALL {len(drifts)} drifts: chi(N=4) >= chi(N=8) >= chi(N=48)  [CONFIRMED]")
    else:
        print(f"    {drift_fails}/{len(drifts)} drifts violate hierarchy")

    print()
    return angle_fails, radius_fails, drift_fails


# ============================================================
# Section 10: Cross-layer summary & verdict
# ============================================================

def section_10(cube, p_canon, r2_canon, p_ext, r2_ext, violations, total_pairs,
               angle_fails, radius_fails, drift_fails):
    print("=" * 70)
    print("SECTION 10: Cross-Layer Summary & Verdict")
    print("=" * 70)

    # Floor universality at N >= 48
    n48 = [(k, v) for k, v in cube.items() if k[0] >= 48]
    at_floor = sum(1 for _, v in n48 if v["chi_emp"] <= 0.02)
    pct_floor = 100.0 * at_floor / len(n48) if n48 else 0

    # Temporal invisibility: chi variance across angle sectors at N >= 48
    bl_chi = [v["chi_emp"] for k, v in cube.items()
              if k[0] >= 48 and -1 <= k[1] <= 1]
    bulk_chi = [v["chi_emp"] for k, v in cube.items()
                if k[0] >= 48 and abs(k[1]) > 1]
    bl_mean = statistics.mean(bl_chi) if bl_chi else 0
    bulk_mean = statistics.mean(bulk_chi) if bulk_chi else 0

    print(f"  Law XIV structural tests:")
    print(f"    Canonical exponent:    alpha = {p_canon:.2f}  (target: -2.33 +/- 0.50)")
    print(f"    Canonical R^2:         {r2_canon:.3f}  (target: >= 0.90)")
    print(f"    Extended exponent:     alpha = {p_ext:.2f}  (target: -2.81 +/- 0.50)")
    print(f"    Extended R^2:          {r2_ext:.3f}  (target: >= 0.90)")
    print(f"    Floor at N>=48:        {at_floor}/{len(n48)} ({pct_floor:.1f}%)")
    print(f"    BL mean chi (N>=48):   {bl_mean:.4f}")
    print(f"    Bulk mean chi (N>=48): {bulk_mean:.4f}")
    print(f"    Hierarchy violations:  {violations}/{total_pairs} ({100*violations/total_pairs:.2f}%)")
    print(f"    Angle invariance:      {angle_fails} failures")
    print(f"    Radius invariance:     {radius_fails} failures")
    print(f"    Drift invariance:      {drift_fails} failures")

    # Verdict
    print(f"\n  ===== LAW XIV CONSOLIDATION VERDICT =====")

    canon_ok = p_canon is not None and abs(p_canon - (-2.33)) <= 0.50 and r2_canon >= 0.90
    ext_ok = p_ext is not None and abs(p_ext - (-2.81)) <= 0.50 and r2_ext >= 0.90
    floor_ok = pct_floor >= 99.5
    hier_ok = violations / total_pairs < 0.02 if total_pairs else False
    invis_ok = abs(bl_mean - bulk_mean) < 0.005
    inv_ok = angle_fails == 0 and radius_fails == 0 and drift_fails == 0

    print(f"  Canonical exponent:      {'CONFIRMED' if canon_ok else 'PARTIAL (p=' + f'{p_canon:.2f}' + ')'}")
    print(f"  Extended exponent:       {'CONFIRMED' if ext_ok else 'PARTIAL (p=' + f'{p_ext:.2f}' + ')'}")
    print(f"  Universal chi-floor:     {'CONFIRMED' if floor_ok else 'PARTIAL'}")
    print(f"  Temporal hierarchy:      {'CONFIRMED' if hier_ok else 'PARTIAL'}")
    print(f"  Boundary-band invis.:    {'CONFIRMED' if invis_ok else 'PARTIAL'}")
    print(f"  Hierarchy invariance:    {'CONFIRMED' if inv_ok else 'PARTIAL'}")

    print()


# ============================================================
# Main
# ============================================================

def main():
    print("Law XIV Consolidation Sweep Analysis")
    print("=" * 70)
    print()

    cube = load_all_atlas()

    section_1(cube)
    p_canon, r2_canon = section_2(cube)
    p_ext, r2_ext = section_3(cube)
    section_4(cube)
    section_5(cube)
    section_6(cube)
    section_7(cube)
    violations, total_pairs = section_8(cube)
    angle_fails, radius_fails, drift_fails = section_9(cube)
    section_10(cube, p_canon or 0, r2_canon or 0, p_ext or 0, r2_ext or 0,
               violations, total_pairs, angle_fails, radius_fails, drift_fails)


if __name__ == "__main__":
    main()
