#!/usr/bin/env python
"""
check_global_phase_boundary.py
===============================
Global Phase Boundary analysis: constructs the IC/non-IC boundary as a
continuous geometric object across (N, angle, radius, drift), extracts
sub-boundaries, detects topological features (folds, cusps, multi-valued
regions, discontinuities), and computes boundary curvature.

Sections:
  1. Data loading and IC-onset extraction
  2. IC onset boundary curves per N (angle vs radius)
  3. Sub-boundary extraction (IC/OPBC, IC/OL, IC/CR, IC/DE)
  4. Boundary interpolation and smoothness
  5. Fold and multi-valued region detection
  6. Cusp detection (sharp corners in boundary)
  7. Discontinuity detection (jumps in onset angle)
  8. Drift perturbation of boundary
  9. Boundary topology changes with N
 10. Curvature and torsion computation
 11. Cross-N summary
"""

import json
import os
import glob
import math
import statistics
from collections import defaultdict, Counter

SWEEP_DIR = os.path.dirname(__file__)

NS     = [4, 8, 12, 20]
ANGLES = [-10, -5, -2, -1, 0, 1, 2, 5, 10]
RADII  = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 0.90, 0.95]
DRIFTS = [0.00, 0.01, 0.02, 0.05, 0.10]

MECHS = ["inward-collapse", "outward-PBC", "DECAY", "PBC-corner", "other-late"]
SHORT = {"inward-collapse": "IC", "outward-PBC": "OP", "DECAY": "DE",
         "PBC-corner": "CR", "other-late": "OL"}


# ============================================================
# Data loading
# ============================================================

def load_all_sweeps():
    """Load all 4D drift sweep files. Returns list of dicts."""
    records = []
    pattern = os.path.join(SWEEP_DIR, "n*_angle_*_r*_d*_sweep.json")
    for path in sorted(glob.glob(pattern)):
        with open(path) as f:
            data = json.load(f)
        records.append(data)
    return records


def build_cube(records):
    """Index records by (N, angle, radius, drift)."""
    cube = {}
    for rec in records:
        key = (rec["N"], rec["angle"], rec["radius"], rec["drift"])
        cube[key] = rec
    return cube


# ============================================================
# 1. IC onset extraction
# ============================================================

def extract_ic_onset(cube):
    """For each (N, radius, drift), find the onset angle where IC first
    appears when scanning from positive angles toward 0.
    Returns dict: (N, radius, drift) -> onset_angle (or None if no IC)."""
    onset = {}
    # Scan angles from most positive to most negative
    sorted_angles = sorted(ANGLES, reverse=True)

    for n in NS:
        for r in RADII:
            for d in DRIFTS:
                ic_angles = []
                non_ic_angles = []
                for ang in sorted_angles:
                    key = (n, ang, r, d)
                    if key not in cube:
                        continue
                    if cube[key]["mechanism"] == "inward-collapse":
                        ic_angles.append(ang)
                    else:
                        non_ic_angles.append(ang)

                if not ic_angles:
                    onset[(n, r, d)] = None  # no IC at this slice
                elif not non_ic_angles:
                    onset[(n, r, d)] = max(ANGLES) + 1  # all IC
                else:
                    # Onset = most positive angle that is IC, with a non-IC
                    # angle above it (or the most positive IC angle if it's
                    # the highest tested angle)
                    onset[(n, r, d)] = max(ic_angles)

    return onset


def check_ic_onset(cube):
    """Section 1: IC onset angles across the parameter space."""
    print("=" * 72)
    print("SECTION 1: IC ONSET ANGLE EXTRACTION")
    print("=" * 72)

    onset = extract_ic_onset(cube)

    for n in NS:
        print(f"\n  N={n}:")
        print(f"  {'radius':>8s}", end="")
        for d in DRIFTS:
            print(f"  d={d:.2f}", end="")
        print()

        for r in RADII:
            print(f"  {r:8.2f}", end="")
            for d in DRIFTS:
                val = onset.get((n, r, d))
                if val is None:
                    print(f"  {'none':>6s}", end="")
                elif val > max(ANGLES):
                    print(f"  {'all':>6s}", end="")
                else:
                    print(f"  {val:>+6d}", end="")
            print()

    return onset


# ============================================================
# 2. IC onset boundary curves per N
# ============================================================

def check_boundary_curves(cube, onset):
    """Section 2: Boundary curves -- onset angle vs radius for each N."""
    print("\n" + "=" * 72)
    print("SECTION 2: IC ONSET BOUNDARY CURVES")
    print("=" * 72)

    for n in NS:
        print(f"\n  N={n}, drift=0.00:")
        d = 0.0
        points = []
        for r in RADII:
            val = onset.get((n, r, d))
            if val is not None and val <= max(ANGLES):
                points.append((r, val))
                print(f"    r={r:.2f} -> onset={val:+d} deg")
            elif val is None:
                print(f"    r={r:.2f} -> no IC")
            else:
                print(f"    r={r:.2f} -> all IC")

        if len(points) >= 2:
            angles_list = [p[1] for p in points]
            r_range = [p[0] for p in points]
            span = max(angles_list) - min(angles_list)
            print(f"    Boundary span: {span} deg over r=[{min(r_range):.2f}, {max(r_range):.2f}]")
            if span == 0:
                print(f"    Boundary type: FLAT (constant onset angle)")
            else:
                print(f"    Boundary type: VARYING ({len(set(angles_list))} distinct onset angles)")
        else:
            print(f"    Boundary: insufficient points ({len(points)})")


# ============================================================
# 3. Sub-boundary extraction
# ============================================================

def check_sub_boundaries(cube):
    """Section 3: Sub-boundaries -- where IC meets each non-IC mechanism."""
    print("\n" + "=" * 72)
    print("SECTION 3: SUB-BOUNDARY EXTRACTION")
    print("=" * 72)

    # For each (N, radius, drift), identify what mechanism is adjacent
    # to IC on the positive-angle side
    sub_boundaries = defaultdict(lambda: Counter())

    for n in NS:
        for r in RADII:
            for d in DRIFTS:
                # Find boundary: last IC angle and first non-IC angle
                sorted_angles = sorted(ANGLES)
                mechs_at_angles = {}
                for ang in sorted_angles:
                    key = (n, ang, r, d)
                    if key in cube:
                        mechs_at_angles[ang] = cube[key]["mechanism"]

                # Find transitions from IC to non-IC
                prev_mech = None
                for ang in sorted_angles:
                    curr = mechs_at_angles.get(ang)
                    if prev_mech == "inward-collapse" and curr and curr != "inward-collapse":
                        sub_boundaries[n][curr] += 1
                    prev_mech = curr

    for n in NS:
        print(f"\n  N={n}: IC boundary neighbors (count of transitions IC -> X):")
        total = sum(sub_boundaries[n].values())
        if total == 0:
            print(f"    No IC boundaries found")
            continue
        for mech in MECHS:
            if mech == "inward-collapse":
                continue
            ct = sub_boundaries[n].get(mech, 0)
            pct = 100 * ct / total if total else 0
            print(f"    IC -> {SHORT.get(mech, mech):>2s}: {ct:4d} ({pct:5.1f}%)")


# ============================================================
# 4. Boundary interpolation and smoothness
# ============================================================

def check_smoothness(cube, onset):
    """Section 4: Smoothness -- first-differences of onset angle along radius."""
    print("\n" + "=" * 72)
    print("SECTION 4: BOUNDARY SMOOTHNESS (FIRST-DIFFERENCES)")
    print("=" * 72)

    for n in NS:
        print(f"\n  N={n}, drift=0.00:")
        d = 0.0
        onsets = []
        for r in RADII:
            val = onset.get((n, r, d))
            if val is not None and val <= max(ANGLES):
                onsets.append((r, val))

        if len(onsets) < 2:
            print(f"    Insufficient boundary points ({len(onsets)})")
            continue

        diffs = []
        for i in range(1, len(onsets)):
            r_prev, a_prev = onsets[i - 1]
            r_curr, a_curr = onsets[i]
            dr = r_curr - r_prev
            da = a_curr - a_prev
            grad = da / dr if dr > 0 else float('inf')
            diffs.append((r_prev, r_curr, da, grad))
            if da != 0:
                print(f"    r={r_prev:.2f}->{r_curr:.2f}: da={da:+d}, grad={grad:+.1f} deg/unit")

        total_jumps = sum(1 for _, _, da, _ in diffs if da != 0)
        max_jump = max(abs(da) for _, _, da, _ in diffs) if diffs else 0
        print(f"    Total jumps: {total_jumps}/{len(diffs)}, max jump: {max_jump} deg")

        if total_jumps == 0:
            print(f"    Smoothness: PERFECTLY FLAT")
        elif max_jump <= 2:
            print(f"    Smoothness: SMOOTH (max jump <= 2 deg)")
        elif max_jump <= 5:
            print(f"    Smoothness: MODERATE (max jump 3-5 deg)")
        else:
            print(f"    Smoothness: ROUGH (max jump > 5 deg)")


# ============================================================
# 5. Fold and multi-valued region detection
# ============================================================

def check_folds(cube, onset):
    """Section 5: Fold detection -- does the boundary curve double back?"""
    print("\n" + "=" * 72)
    print("SECTION 5: FOLD AND MULTI-VALUED REGION DETECTION")
    print("=" * 72)

    for n in NS:
        print(f"\n  N={n}:")
        for d in DRIFTS:
            onsets_at_d = []
            for r in RADII:
                val = onset.get((n, r, d))
                if val is not None and val <= max(ANGLES):
                    onsets_at_d.append((r, val))

            if len(onsets_at_d) < 3:
                continue

            # Check for non-monotonicity (fold): onset angle goes up then down
            angles_seq = [a for _, a in onsets_at_d]
            folds = 0
            for i in range(1, len(angles_seq) - 1):
                d1 = angles_seq[i] - angles_seq[i - 1]
                d2 = angles_seq[i + 1] - angles_seq[i]
                if d1 > 0 and d2 < 0:
                    folds += 1
                    r_fold = onsets_at_d[i][0]
                    print(f"    d={d:.2f}: FOLD at r={r_fold:.2f} "
                          f"(peak onset={angles_seq[i]:+d})")
                elif d1 < 0 and d2 > 0:
                    folds += 1
                    r_fold = onsets_at_d[i][0]
                    print(f"    d={d:.2f}: FOLD (valley) at r={r_fold:.2f} "
                          f"(min onset={angles_seq[i]:+d})")

            if folds == 0:
                # Check if monotonic or flat
                unique = set(angles_seq)
                if len(unique) == 1:
                    label = "FLAT"
                elif all(angles_seq[i] <= angles_seq[i+1]
                         for i in range(len(angles_seq)-1)):
                    label = "MONOTONE-INCREASING"
                elif all(angles_seq[i] >= angles_seq[i+1]
                         for i in range(len(angles_seq)-1)):
                    label = "MONOTONE-DECREASING"
                else:
                    label = "STEP-WISE (non-smooth but no fold)"
                print(f"    d={d:.2f}: No folds -- {label}")

    # Multi-valued detection: at a given (N, angle, drift), do both IC and
    # non-IC exist at different radii?
    print(f"\n  Multi-valued regions (IC and non-IC coexist at same angle):")
    for n in NS:
        mv_count = 0
        for d in DRIFTS:
            for ang in ANGLES:
                ic_radii = []
                non_ic_radii = []
                for r in RADII:
                    key = (n, ang, r, d)
                    if key not in cube:
                        continue
                    if cube[key]["mechanism"] == "inward-collapse":
                        ic_radii.append(r)
                    else:
                        non_ic_radii.append(r)
                # Multi-valued if IC radii are interspersed with non-IC
                if ic_radii and non_ic_radii:
                    if min(ic_radii) < max(non_ic_radii) and \
                       max(ic_radii) > min(non_ic_radii):
                        mv_count += 1
        print(f"    N={n}: {mv_count} multi-valued (angle, drift) slices")


# ============================================================
# 6. Cusp detection
# ============================================================

def check_cusps(cube, onset):
    """Section 6: Cusp detection -- sharp corners in the boundary curve."""
    print("\n" + "=" * 72)
    print("SECTION 6: CUSP DETECTION")
    print("=" * 72)

    for n in NS:
        print(f"\n  N={n}:")
        for d in [0.00, 0.05, 0.10]:
            onsets_at_d = []
            for r in RADII:
                val = onset.get((n, r, d))
                if val is not None and val <= max(ANGLES):
                    onsets_at_d.append((r, val))

            if len(onsets_at_d) < 3:
                continue

            # Cusp = large second-difference (abrupt change in gradient)
            cusps = []
            for i in range(1, len(onsets_at_d) - 1):
                r0, a0 = onsets_at_d[i - 1]
                r1, a1 = onsets_at_d[i]
                r2, a2 = onsets_at_d[i + 1]
                d2a = (a2 - a1) - (a1 - a0)
                if abs(d2a) >= 3:
                    cusps.append((r1, a1, d2a))

            if cusps:
                for r1, a1, d2a in cusps:
                    print(f"    d={d:.2f}: CUSP at r={r1:.2f} "
                          f"(onset={a1:+d}, d2a={d2a:+d})")
            else:
                print(f"    d={d:.2f}: No cusps (max |d2a| < 3)")


# ============================================================
# 7. Discontinuity detection
# ============================================================

def check_discontinuities(cube, onset):
    """Section 7: Discontinuities -- jumps in onset angle between
    adjacent radii."""
    print("\n" + "=" * 72)
    print("SECTION 7: DISCONTINUITY DETECTION")
    print("=" * 72)

    for n in NS:
        print(f"\n  N={n}:")
        for d in DRIFTS:
            onsets_at_d = []
            for r in RADII:
                val = onset.get((n, r, d))
                if val is not None and val <= max(ANGLES):
                    onsets_at_d.append((r, val))

            jumps = []
            for i in range(1, len(onsets_at_d)):
                r0, a0 = onsets_at_d[i - 1]
                r1, a1 = onsets_at_d[i]
                da = abs(a1 - a0)
                if da >= 5:
                    jumps.append((r0, r1, a0, a1, da))

            if jumps:
                for r0, r1, a0, a1, da in jumps:
                    print(f"    d={d:.2f}: JUMP r={r0:.2f}->{r1:.2f} "
                          f"onset={a0:+d}->{a1:+d} (|da|={da})")
            else:
                max_da = 0
                if len(onsets_at_d) >= 2:
                    max_da = max(abs(onsets_at_d[i][1] - onsets_at_d[i-1][1])
                                for i in range(1, len(onsets_at_d)))
                print(f"    d={d:.2f}: No discontinuities (max |da|={max_da})")


# ============================================================
# 8. Drift perturbation of boundary
# ============================================================

def check_drift_perturbation(cube, onset):
    """Section 8: How drift shifts the IC boundary."""
    print("\n" + "=" * 72)
    print("SECTION 8: DRIFT PERTURBATION OF BOUNDARY")
    print("=" * 72)

    for n in NS:
        print(f"\n  N={n}:")
        # Compare onset at drift=0 vs each other drift
        for r in RADII:
            base = onset.get((n, r, 0.0))
            if base is None or base > max(ANGLES):
                continue
            shifts = []
            for d in DRIFTS[1:]:
                val = onset.get((n, r, d))
                if val is not None and val <= max(ANGLES):
                    shift = val - base
                    shifts.append((d, shift))

            if shifts:
                shift_str = "  ".join(f"d={d:.2f}:{s:+d}" for d, s in shifts)
                print(f"    r={r:.2f} (base={base:+d}): {shift_str}")

    # Summary: mean onset shift per drift level
    print(f"\n  Mean onset shift (onset_d - onset_0) by N and drift:")
    print(f"  {'N':>4s}", end="")
    for d in DRIFTS[1:]:
        print(f"  d={d:.2f}", end="")
    print()

    for n in NS:
        print(f"  {n:4d}", end="")
        for d in DRIFTS[1:]:
            shifts = []
            for r in RADII:
                base = onset.get((n, r, 0.0))
                val = onset.get((n, r, d))
                if (base is not None and base <= max(ANGLES) and
                        val is not None and val <= max(ANGLES)):
                    shifts.append(val - base)
            if shifts:
                print(f"  {statistics.mean(shifts):>+6.2f}", end="")
            else:
                print(f"  {'n/a':>6s}", end="")
        print()


# ============================================================
# 9. Boundary topology changes with N
# ============================================================

def check_topology_vs_n(cube, onset):
    """Section 9: How boundary topology changes with N."""
    print("\n" + "=" * 72)
    print("SECTION 9: BOUNDARY TOPOLOGY CHANGES WITH N")
    print("=" * 72)

    for n in NS:
        d = 0.0
        # Count how many radii have IC, how many don't
        ic_count = 0
        non_ic_count = 0
        mixed_angles = 0
        dominant_non_ic = Counter()

        for r in RADII:
            has_ic = False
            has_non_ic = False
            for ang in ANGLES:
                key = (n, ang, r, d)
                if key not in cube:
                    continue
                if cube[key]["mechanism"] == "inward-collapse":
                    has_ic = True
                else:
                    has_non_ic = True
                    dominant_non_ic[cube[key]["mechanism"]] += 1
            if has_ic and has_non_ic:
                mixed_angles += 1
            elif has_ic:
                ic_count += 1
            elif has_non_ic:
                non_ic_count += 1

        # Boundary connectivity: is the IC region contiguous?
        # For each radius, find the angle range that is IC
        ic_ranges = []
        for r in RADII:
            ic_ang = []
            for ang in sorted(ANGLES):
                key = (n, ang, r, d)
                if key in cube and cube[key]["mechanism"] == "inward-collapse":
                    ic_ang.append(ang)
            if ic_ang:
                ic_ranges.append((r, min(ic_ang), max(ic_ang)))

        print(f"\n  N={n} (drift=0):")
        print(f"    Radii with only IC: {ic_count}")
        print(f"    Radii with only non-IC: {non_ic_count}")
        print(f"    Radii with mixed (boundary crosses): {mixed_angles}")
        print(f"    Non-IC neighbors: {dict(dominant_non_ic.most_common(5))}")

        if ic_ranges:
            min_ang_range = min(hi - lo for _, lo, hi in ic_ranges)
            max_ang_range = max(hi - lo for _, lo, hi in ic_ranges)
            print(f"    IC angle range: [{min_ang_range}, {max_ang_range}] deg")

        # Is boundary simply connected? Check if IC region has holes
        for r_idx in range(len(RADII)):
            r = RADII[r_idx]
            for ang_idx in range(len(ANGLES)):
                ang = ANGLES[ang_idx]
                key = (n, ang, r, d)
                if key not in cube:
                    continue
                if cube[key]["mechanism"] != "inward-collapse":
                    continue
                # Check if this IC cell is surrounded by non-IC on both sides
                # along radius
                if r_idx > 0 and r_idx < len(RADII) - 1:
                    key_lo = (n, ang, RADII[r_idx - 1], d)
                    key_hi = (n, ang, RADII[r_idx + 1], d)
                    if (key_lo in cube and key_hi in cube and
                            cube[key_lo]["mechanism"] != "inward-collapse" and
                            cube[key_hi]["mechanism"] != "inward-collapse"):
                        print(f"    HOLE: IC at r={r:.2f}, ang={ang:+d} "
                              f"surrounded by non-IC")


# ============================================================
# 10. Curvature computation
# ============================================================

def check_curvature(cube, onset):
    """Section 10: Boundary curvature in the (angle, radius) plane."""
    print("\n" + "=" * 72)
    print("SECTION 10: BOUNDARY CURVATURE AND TORSION")
    print("=" * 72)

    for n in NS:
        print(f"\n  N={n}:")
        for d in [0.00, 0.05, 0.10]:
            # Build boundary curve: (radius, onset_angle) pairs
            pts = []
            for r in RADII:
                val = onset.get((n, r, d))
                if val is not None and val <= max(ANGLES):
                    pts.append((r, float(val)))

            if len(pts) < 3:
                print(f"    d={d:.2f}: insufficient points for curvature")
                continue

            # Compute curvature at each interior point using 3-point formula
            # kappa = |x'y'' - y'x''| / (x'^2 + y'^2)^(3/2)
            curvatures = []
            for i in range(1, len(pts) - 1):
                x0, y0 = pts[i - 1]
                x1, y1 = pts[i]
                x2, y2 = pts[i + 1]
                dx1 = x1 - x0
                dy1 = y1 - y0
                dx2 = x2 - x1
                dy2 = y2 - y1
                # Second differences
                ddx = dx2 - dx1
                ddy = dy2 - dy1
                # Average first derivative
                dxa = (dx1 + dx2) / 2
                dya = (dy1 + dy2) / 2
                denom = (dxa**2 + dya**2)**1.5
                if denom > 1e-15:
                    kappa = abs(dxa * ddy - dya * ddx) / denom
                else:
                    kappa = 0.0
                curvatures.append((pts[i][0], kappa))

            if curvatures:
                max_k = max(k for _, k in curvatures)
                mean_k = statistics.mean(k for _, k in curvatures)
                max_r = [r for r, k in curvatures if k == max_k][0]
                print(f"    d={d:.2f}: mean kappa={mean_k:.3f}, "
                      f"max kappa={max_k:.3f} at r={max_r:.2f}")
                # Report high-curvature points
                for r, k in curvatures:
                    if k > 2 * mean_k and k > 0.1:
                        print(f"      HIGH CURVATURE: r={r:.2f}, kappa={k:.3f}")
            else:
                print(f"    d={d:.2f}: no interior points for curvature")

    # Torsion: how the boundary twists when drift is included as third dim
    print(f"\n  Torsion (boundary twist with drift):")
    for n in NS:
        torsion_vals = []
        for r in RADII:
            # Get onset angle at d=0, d=0.05, d=0.10
            a0 = onset.get((n, r, 0.0))
            a1 = onset.get((n, r, 0.05))
            a2 = onset.get((n, r, 0.10))
            if (a0 is not None and a0 <= max(ANGLES) and
                    a1 is not None and a1 <= max(ANGLES) and
                    a2 is not None and a2 <= max(ANGLES)):
                # Simple torsion measure: second difference in drift direction
                tau = (a2 - a1) - (a1 - a0)
                torsion_vals.append((r, tau))

        if torsion_vals:
            max_tau = max(abs(t) for _, t in torsion_vals)
            mean_tau = statistics.mean(abs(t) for _, t in torsion_vals)
            print(f"    N={n}: mean |tau|={mean_tau:.2f}, max |tau|={max_tau:.1f}")
            for r, tau in torsion_vals:
                if abs(tau) > 2:
                    print(f"      TWIST: r={r:.2f}, tau={tau:+d}")
        else:
            print(f"    N={n}: insufficient drift data for torsion")


# ============================================================
# 11. Cross-N summary
# ============================================================

def check_cross_n_summary(cube, onset):
    """Section 11: Cross-N summary of boundary properties."""
    print("\n" + "=" * 72)
    print("SECTION 11: CROSS-N BOUNDARY SUMMARY")
    print("=" * 72)

    print(f"\n  {'N':>4s}  {'IC%':>5s}  {'onset_range':>12s}  "
          f"{'n_folds':>7s}  {'n_cusps':>7s}  {'n_jumps':>7s}  "
          f"{'max_da':>6s}  {'topology':>12s}")

    for n in NS:
        d = 0.0
        # IC fraction
        total = 0
        ic_total = 0
        for r in RADII:
            for ang in ANGLES:
                key = (n, ang, r, d)
                if key in cube:
                    total += 1
                    if cube[key]["mechanism"] == "inward-collapse":
                        ic_total += 1
        ic_pct = 100 * ic_total / total if total else 0

        # Onset range
        onsets_at_d = []
        for r in RADII:
            val = onset.get((n, r, d))
            if val is not None and val <= max(ANGLES):
                onsets_at_d.append(val)

        if onsets_at_d:
            onset_range = f"[{min(onsets_at_d):+d},{max(onsets_at_d):+d}]"
        else:
            onset_range = "n/a"

        # Folds
        n_folds = 0
        angles_seq = []
        for r in RADII:
            val = onset.get((n, r, d))
            if val is not None and val <= max(ANGLES):
                angles_seq.append(val)
        for i in range(1, len(angles_seq) - 1):
            d1 = angles_seq[i] - angles_seq[i - 1]
            d2 = angles_seq[i + 1] - angles_seq[i]
            if (d1 > 0 and d2 < 0) or (d1 < 0 and d2 > 0):
                n_folds += 1

        # Cusps (|d2a| >= 3)
        n_cusps = 0
        for i in range(1, len(angles_seq) - 1):
            d2a = (angles_seq[i + 1] - angles_seq[i]) - \
                  (angles_seq[i] - angles_seq[i - 1])
            if abs(d2a) >= 3:
                n_cusps += 1

        # Jumps (|da| >= 5)
        n_jumps = 0
        max_da = 0
        for i in range(1, len(angles_seq)):
            da = abs(angles_seq[i] - angles_seq[i - 1])
            max_da = max(max_da, da)
            if da >= 5:
                n_jumps += 1

        # Topology
        if not onsets_at_d:
            topo = "no-boundary"
        elif len(set(onsets_at_d)) == 1:
            topo = "flat"
        elif n_folds == 0 and n_jumps == 0:
            topo = "smooth"
        elif n_folds > 0:
            topo = "folded"
        elif n_jumps > 0:
            topo = "piecewise"
        else:
            topo = "complex"

        print(f"  {n:4d}  {ic_pct:5.1f}  {onset_range:>12s}  "
              f"{n_folds:7d}  {n_cusps:7d}  {n_jumps:7d}  "
              f"{max_da:6d}  {topo:>12s}")

    # Mechanism at boundary
    print(f"\n  Mechanism at boundary (what IC transitions into):")
    for n in NS:
        boundary_mechs = Counter()
        for d in DRIFTS:
            for r in RADII:
                # Find the first non-IC angle scanning from negative to positive
                for ang_idx in range(len(ANGLES) - 1):
                    ang_lo = ANGLES[ang_idx]
                    ang_hi = ANGLES[ang_idx + 1]
                    key_lo = (n, ang_lo, r, d)
                    key_hi = (n, ang_hi, r, d)
                    if (key_lo in cube and key_hi in cube):
                        if (cube[key_lo]["mechanism"] == "inward-collapse" and
                                cube[key_hi]["mechanism"] != "inward-collapse"):
                            boundary_mechs[cube[key_hi]["mechanism"]] += 1
                        elif (cube[key_lo]["mechanism"] != "inward-collapse" and
                                cube[key_hi]["mechanism"] == "inward-collapse"):
                            boundary_mechs[cube[key_lo]["mechanism"]] += 1
        total = sum(boundary_mechs.values())
        if total:
            parts = []
            for mech, ct in boundary_mechs.most_common():
                parts.append(f"{SHORT.get(mech, mech)}={100*ct/total:.0f}%")
            print(f"    N={n}: {', '.join(parts)} ({total} transitions)")
        else:
            print(f"    N={n}: no transitions found")

    # Temporal gradient at boundary
    print(f"\n  Temporal gradient at boundary (chi just inside IC vs just outside):")
    for n in NS:
        chi_in = []
        chi_out = []
        for d in DRIFTS:
            for r in RADII:
                for ang_idx in range(len(ANGLES) - 1):
                    ang_lo = ANGLES[ang_idx]
                    ang_hi = ANGLES[ang_idx + 1]
                    key_lo = (n, ang_lo, r, d)
                    key_hi = (n, ang_hi, r, d)
                    if key_lo in cube and key_hi in cube:
                        if (cube[key_lo]["mechanism"] == "inward-collapse" and
                                cube[key_hi]["mechanism"] != "inward-collapse"):
                            chi_in.append(cube[key_lo]["chi_emp"])
                            chi_out.append(cube[key_hi]["chi_emp"])
        if chi_in:
            mean_in = statistics.mean(chi_in)
            mean_out = statistics.mean(chi_out)
            ratio = mean_out / mean_in if mean_in > 0 else float('inf')
            print(f"    N={n}: chi_IC={mean_in:.2f}, chi_nonIC={mean_out:.2f}, "
                  f"ratio={ratio:.1f}x")
        else:
            print(f"    N={n}: no boundary pairs for temporal gradient")


# ============================================================
# Main
# ============================================================

def main():
    print("Global Phase Boundary Analysis")
    print("=" * 72)

    records = load_all_sweeps()
    print(f"Loaded {len(records)} sweep records")

    cube = build_cube(records)
    print(f"Built 4D cube with {len(cube)} cells")
    print()

    onset = check_ic_onset(cube)
    check_boundary_curves(cube, onset)
    check_sub_boundaries(cube)
    check_smoothness(cube, onset)
    check_folds(cube, onset)
    check_cusps(cube, onset)
    check_discontinuities(cube, onset)
    check_drift_perturbation(cube, onset)
    check_topology_vs_n(cube, onset)
    check_curvature(cube, onset)
    check_cross_n_summary(cube, onset)


if __name__ == "__main__":
    main()
