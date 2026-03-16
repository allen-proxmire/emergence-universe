#!/usr/bin/env python
"""
check_4d_drift_robustness.py
==============================
4D robustness analysis: N x angle x radius x radial-drift.

Tests whether the key 3D manifold structures survive when a small
radial drift component is added to the initial velocity vector:
  - tangent discontinuity (0 -> +1 deg)
  - DECAY filament at (N=4, angle=0)
  - N=8 plateau
  - large-N flattening / near-binary partitioning

Data sources (in priority order):
  1. Per-point 4D files:  n{N}_angle_{label}_{rlabel}_d{DDD}_sweep.json
  2. Existing 3D files:   n{N}_angle_{label}_{rlabel}_sweep.json  (drift=0)
  3. Existing 2D grids:   n{N}_2d_angle_radius_grid.json          (drift=0)
"""

import json
import os
import glob
import re
from collections import Counter

SWEEP_DIR = os.path.dirname(__file__)

# Expected parameter ranges
NS = [4, 8, 12, 20]
ANGLES = [-10, -5, -2, -1, 0, 1, 2, 5, 10]
RADII_FRAC = [round(r / 100, 2) for r in range(5, 100, 5)]  # 0.05 .. 0.95
DRIFTS = [0.00, 0.01, 0.02, 0.05, 0.10]
BOX_PX = 400

MECHS = ["inward-collapse", "outward-PBC", "DECAY", "PBC-corner", "other-late"]
SHORT = {"inward-collapse": "IC", "outward-PBC": "OP", "DECAY": "DE",
         "PBC-corner": "CR", "other-late": "OL"}
CODE  = {"inward-collapse": "I", "outward-PBC": "O", "DECAY": "D",
         "PBC-corner": "C", "other-late": "."}


# ============================================================
# File loading helpers
# ============================================================

def angle_label(ang):
    if ang > 0:   return f"p{ang}deg"
    if ang < 0:   return f"m{abs(ang)}deg"
    return "0deg"

def radius_label(r_frac):
    return f"r{int(round(r_frac * 100)):03d}"

def drift_label(d):
    return f"d{int(round(d * 1000)):03d}"

def radius_from_d_px(d_px):
    return round(d_px / BOX_PX, 4)


def load_4d_file(n, ang, r_frac, drift):
    """Load a per-point 4D sweep file. Returns list of records or None."""
    fname = (f"n{n}_angle_{angle_label(ang)}_"
             f"{radius_label(r_frac)}_{drift_label(drift)}_sweep.json")
    path = os.path.join(SWEEP_DIR, fname)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict) and "results" in data:
        data = data["results"]
    if isinstance(data, dict):
        data = [data]
    return data


def load_3d_file(n, ang, r_frac):
    """Load a drift=0 3D sweep file (no drift suffix)."""
    fname = f"n{n}_angle_{angle_label(ang)}_{radius_label(r_frac)}_sweep.json"
    path = os.path.join(SWEEP_DIR, fname)
    if not os.path.exists(path):
        return None
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict) and "results" in data:
        data = data["results"]
    if isinstance(data, dict):
        data = [data]
    return data


def load_2d_grid(n):
    """Load an existing 2D angle x d_px grid. Returns dict or None."""
    path = os.path.join(SWEEP_DIR, f"n{n}_2d_angle_radius_grid.json")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        return json.load(f)


def mech_fracs(records):
    total = len(records)
    if total == 0:
        return {m: 0.0 for m in MECHS}
    counts = Counter(r["mechanism"] for r in records)
    return {m: 100.0 * counts.get(m, 0) / total for m in MECHS}


def dominant(fracs_dict):
    return max(fracs_dict, key=fracs_dict.get)


# ============================================================
# Build the 4D hypercube
# ============================================================

def build_hypercube():
    """
    Build hcube[n][angle][radius][drift] = {mech: frac, ...}
    Returns (hcube, loaded_ns, loaded_angles, loaded_radii, loaded_drifts).
    """
    hcube = {}
    all_radii = set()
    all_angles = set()
    all_ns = set()
    all_drifts = set()
    count_4d = 0
    count_3d = 0
    count_2d = 0

    # --- Pass 1: 4D sweep files (explicit drift) ---
    for n in NS:
        for ang in ANGLES:
            for r in RADII_FRAC:
                for d in DRIFTS:
                    records = load_4d_file(n, ang, r, d)
                    if records is not None:
                        fracs = mech_fracs(records)
                        (hcube.setdefault(n, {})
                               .setdefault(ang, {})
                               .setdefault(r, {})[d]) = fracs
                        all_ns.add(n); all_angles.add(ang)
                        all_radii.add(r); all_drifts.add(d)
                        count_4d += 1

    # --- Pass 2: scan for unexpected 4D files ---
    pattern = os.path.join(SWEEP_DIR, "n*_angle_*_r*_d*_sweep.json")
    for path in glob.glob(pattern):
        fname = os.path.basename(path)
        m = re.match(
            r"n(\d+)_angle_([pm]\d+deg|0deg)_r(\d+)_d(\d+)_sweep\.json",
            fname)
        if not m:
            continue
        n = int(m.group(1))
        ang_str = m.group(2)
        r_int = int(m.group(3))
        d_int = int(m.group(4))
        if ang_str == "0deg":
            ang = 0
        elif ang_str.startswith("p"):
            ang = int(ang_str[1:-3])
        else:
            ang = -int(ang_str[1:-3])
        r = round(r_int / 100, 2)
        d = round(d_int / 1000, 3)
        if (n in hcube and ang in hcube.get(n, {})
            and r in hcube.get(n, {}).get(ang, {})
            and d in hcube.get(n, {}).get(ang, {}).get(r, {})):
            continue
        with open(path) as f:
            data = json.load(f)
        if isinstance(data, dict) and "results" in data:
            data = data["results"]
        if isinstance(data, dict):
            data = [data]
        (hcube.setdefault(n, {})
               .setdefault(ang, {})
               .setdefault(r, {})[d]) = mech_fracs(data)
        all_ns.add(n); all_angles.add(ang)
        all_radii.add(r); all_drifts.add(d)

    # --- Pass 3: backfill drift=0 from 3D files ---
    for n in NS:
        for ang in ANGLES:
            for r in RADII_FRAC:
                if (n in hcube and ang in hcube.get(n, {})
                    and r in hcube.get(n, {}).get(ang, {})
                    and 0.0 in hcube.get(n, {}).get(ang, {}).get(r, {})):
                    continue
                records = load_3d_file(n, ang, r)
                if records is not None:
                    (hcube.setdefault(n, {})
                           .setdefault(ang, {})
                           .setdefault(r, {})[0.0]) = mech_fracs(records)
                    all_ns.add(n); all_angles.add(ang)
                    all_radii.add(r); all_drifts.add(0.0)
                    count_3d += 1

    # --- Pass 4: backfill drift=0 from 2D grids ---
    for n in NS:
        grid_data = load_2d_grid(n)
        if grid_data is None:
            continue
        for rec in grid_data["grid"]:
            ang = rec["angle"]
            d_px = rec["d_px"]
            r = radius_from_d_px(d_px)
            if (n in hcube and ang in hcube.get(n, {})
                and r in hcube.get(n, {}).get(ang, {})
                and 0.0 in hcube.get(n, {}).get(ang, {}).get(r, {})):
                continue
            fracs = {m: 0.0 for m in MECHS}
            fracs[rec["mechanism"]] = 100.0
            (hcube.setdefault(n, {})
                   .setdefault(ang, {})
                   .setdefault(r, {})[0.0]) = fracs
            all_ns.add(n); all_angles.add(ang)
            all_radii.add(r); all_drifts.add(0.0)
            count_2d += 1

    ns = sorted(all_ns)
    angles = sorted(all_angles)
    radii = sorted(all_radii)
    drifts = sorted(all_drifts)

    print(f"  Data loaded: {count_4d} from 4D files, "
          f"{count_3d} from 3D files, {count_2d} from 2D grids")
    print(f"  Coverage: {len(ns)} N, {len(angles)} angles, "
          f"{len(radii)} radii, {len(drifts)} drifts")
    total_cells = sum(
        sum(sum(len(hcube.get(n, {}).get(a, {}).get(r, {}))
                for r in radii)
            for a in angles)
        for n in ns
    )
    max_cells = len(ns) * len(angles) * len(radii) * len(drifts)
    pct = 100 * total_cells / max_cells if max_cells > 0 else 0
    print(f"  Filled cells: {total_cells} / {max_cells} ({pct:.1f}%)")
    print()
    return hcube, ns, angles, radii, drifts


# ============================================================
# Helper: get a fracs dict from the hypercube, or None
# ============================================================

def get_fracs(hcube, n, ang, r, d):
    return hcube.get(n, {}).get(ang, {}).get(r, {}).get(d)


# ============================================================
# Section 1: Per-N phase maps at each drift level
# ============================================================

def print_drift_phase_maps(hcube, ns, angles, radii, drifts):
    for n in ns:
        for d in drifts:
            # Check if any data exists at this (N, drift)
            has_data = False
            for ang in angles:
                for r in radii:
                    if get_fracs(hcube, n, ang, r, d) is not None:
                        has_data = True
                        break
                if has_data:
                    break
            if not has_data:
                continue

            print("=" * 90)
            print(f"  N = {n}   drift = {d:.3f}   "
                  f"PHASE MAP (angle x radius)")
            print(f"  Codes: I=IC  O=OPBC  D=DECAY  C=PBC-corner  "
                  f".=other-late  ?=missing")
            print("=" * 90)

            show_radii = radii if len(radii) <= 20 else radii[::2]
            hdr = "  angle "
            for r in show_radii:
                hdr += f" {r:.2f}"[1:]
            print(hdr)
            print("  " + "-" * (7 + 5 * len(show_radii)))

            for ang in angles:
                row = f"  {ang:+4d}  "
                for r in show_radii:
                    f = get_fracs(hcube, n, ang, r, d)
                    if f is None:
                        row += "   ? "
                    else:
                        dom = dominant(f)
                        row += f"   {CODE.get(dom, '?')} "
                    if r == show_radii[-1] and ang == 0:
                        row += "  <-- tangent"
                print(row)
            print()


# ============================================================
# Section 2: Mechanism fractions by drift (averaged over radii)
# ============================================================

def print_drift_profiles(hcube, ns, angles, radii, drifts):
    print("=" * 90)
    print("  MECHANISM FRACTIONS BY (N, drift) — averaged over all radii")
    print("=" * 90)
    print(f"  {'N':>3}  {'drift':>6}  {'ang':>4}  "
          f"{'IC':>6}  {'OP':>6}  {'DE':>6}  {'CR':>6}  {'OL':>6}  {'pts':>4}")
    print("  " + "-" * 56)

    for n in ns:
        for d in drifts:
            any_data = False
            for ang in angles:
                agg = {m: 0.0 for m in MECHS}
                count = 0
                for r in radii:
                    f = get_fracs(hcube, n, ang, r, d)
                    if f is not None:
                        for m in MECHS:
                            agg[m] += f[m]
                        count += 1
                if count == 0:
                    continue
                any_data = True
                for m in MECHS:
                    agg[m] /= count
                print(f"  {n:3d}  {d:6.3f}  {ang:+4d}  "
                      f"{agg['inward-collapse']:5.1f}%"
                      f"  {agg['outward-PBC']:5.1f}%"
                      f"  {agg['DECAY']:5.1f}%"
                      f"  {agg['PBC-corner']:5.1f}%"
                      f"  {agg['other-late']:5.1f}%"
                      f"  {count:4d}")
            if any_data:
                print()


# ============================================================
# Section 3: Tangent discontinuity survival under drift
# ============================================================

def check_tangent_discontinuity_vs_drift(hcube, ns, radii, drifts):
    print("=" * 90)
    print("  TANGENT DISCONTINUITY vs DRIFT: IC jump from 0 deg -> +1 deg")
    print("=" * 90)

    for n in ns:
        print(f"\n  N = {n}:")
        print(f"  {'drift':>6}  {'tested':>6}  {'jump>10':>8}  "
              f"{'jump~0':>7}  {'mean':>8}  {'verdict':>20}")
        print(f"  " + "-" * 60)

        for d in drifts:
            jumps = []
            for r in radii:
                f0 = get_fracs(hcube, n, 0, r, d)
                f1 = get_fracs(hcube, n, 1, r, d)
                if f0 is not None and f1 is not None:
                    jumps.append(f1["inward-collapse"] -
                                 f0["inward-collapse"])
            if not jumps:
                continue

            n_pos = sum(1 for j in jumps if j > 10)
            n_zero = sum(1 for j in jumps if abs(j) < 5)
            mn = sum(jumps) / len(jumps)

            if n_pos > 0.8 * len(jumps):
                verdict = "UNIVERSAL sheet"
            elif n_pos > 0.5 * len(jumps):
                verdict = "PARTIAL sheet"
            elif n_pos > 0.2 * len(jumps):
                verdict = "WEAK sheet"
            else:
                verdict = "DESTROYED"

            print(f"  {d:6.3f}  {len(jumps):6d}  "
                  f"{n_pos:5d}/{len(jumps):<3d}  "
                  f"{n_zero:4d}/{len(jumps):<3d}  "
                  f"{mn:+7.1f}pp  "
                  f"{verdict:>20}")


# ============================================================
# Section 4: DECAY survival under drift
# ============================================================

def check_decay_vs_drift(hcube, ns, radii, drifts):
    print("\n" + "=" * 90)
    print("  DECAY SURVIVAL vs DRIFT (N=4, angle=0)")
    print("=" * 90)

    for n in ns:
        for d in drifts:
            decay_cells = 0
            total = 0
            for r in radii:
                f = get_fracs(hcube, n, 0, r, d)
                if f is not None:
                    total += 1
                    if f["DECAY"] > 50:
                        decay_cells += 1
            if total == 0:
                continue
            if decay_cells > 0 or (n == 4 and d == 0.0):
                print(f"  N={n:2d}  drift={d:.3f}:  "
                      f"DECAY cells = {decay_cells}/{total} "
                      f"({100*decay_cells/total:.0f}%)")

    # Compact summary for N=4
    print()
    print("  N=4 DECAY survival trajectory:")
    for d in drifts:
        total = 0
        decay_cells = 0
        for r in radii:
            f = get_fracs(hcube, 4, 0, r, d)
            if f is not None:
                total += 1
                if f["DECAY"] > 50:
                    decay_cells += 1
        if total == 0:
            print(f"    drift={d:.3f}:  no data")
        else:
            bar = "#" * decay_cells + "." * (total - decay_cells)
            pct = 100 * decay_cells / total
            print(f"    drift={d:.3f}:  {decay_cells:2d}/{total} "
                  f"({pct:5.1f}%)  |{bar}|")

    # DECAY at non-zero angles under drift?
    print()
    print("  DECAY at non-zero angles (any N, any drift):")
    found = False
    for n in ns:
        for d in drifts:
            for ang in [-10, -5, -2, -1, 1, 2, 5, 10]:
                for r in radii:
                    f = get_fracs(hcube, n, ang, r, d)
                    if f is not None and f["DECAY"] > 50:
                        print(f"    N={n} ang={ang:+d} r={r:.2f} "
                              f"drift={d:.3f}: DECAY={f['DECAY']:.0f}%")
                        found = True
    if not found:
        print("    None found -- DECAY remains angle=0 exclusive")


# ============================================================
# Section 5: N=8 anomaly under drift
# ============================================================

def check_n8_anomaly_vs_drift(hcube, ns, angles, radii, drifts):
    print("\n" + "=" * 90)
    print("  N=8 ANOMALY vs DRIFT")
    print("=" * 90)

    # Compare anomaly footprint across drift levels
    print(f"\n  Anomaly footprint (OL>20% or CR>20% or IC<30% at ang>0):")
    print(f"  {'N':>3}  {'drift':>6}  {'anom/total':>12}  "
          f"{'pct':>5}  {'class':>12}")
    print(f"  " + "-" * 44)

    for n in ns:
        for d in drifts:
            anomaly = 0
            total = 0
            for ang in angles:
                for r in radii:
                    f = get_fracs(hcube, n, ang, r, d)
                    if f is None:
                        continue
                    total += 1
                    if (f["other-late"] > 20 or
                        f["PBC-corner"] > 20 or
                        (ang > 0 and f["inward-collapse"] < 30)):
                        anomaly += 1
            if total == 0:
                continue
            pct = 100 * anomaly / total
            if pct > 30:
                cls = "PLATEAU"
            elif pct > 10:
                cls = "RIDGE"
            elif pct > 3:
                cls = "LOCALIZED"
            else:
                cls = "ABSENT"
            print(f"  {n:3d}  {d:6.3f}  {anomaly:5d}/{total:<5d}"
                  f"  {pct:4.0f}%  {cls:>12}")
        print()

    # OL at angle=+1 across drift (key anomaly indicator)
    print(f"  Other-late at angle=+1 across drift:")
    print(f"  {'N':>3}  {'drift':>6}  {'mean OL':>8}  {'mean IC':>8}")
    print(f"  " + "-" * 30)
    for n in ns:
        for d in drifts:
            ols, ics = [], []
            for r in radii:
                f = get_fracs(hcube, n, 1, r, d)
                if f is not None:
                    ols.append(f["other-late"])
                    ics.append(f["inward-collapse"])
            if ols:
                print(f"  {n:3d}  {d:6.3f}  "
                      f"{sum(ols)/len(ols):7.1f}%  "
                      f"{sum(ics)/len(ics):7.1f}%")
        print()


# ============================================================
# Section 6: Large-N binary partitioning under drift
# ============================================================

def check_large_n_stability(hcube, ns, angles, radii, drifts):
    print("=" * 90)
    print("  LARGE-N BINARY PARTITIONING vs DRIFT")
    print("=" * 90)

    print(f"\n  Distinct dominant mechanisms per (N, drift) plane:")
    print(f"  {'N':>3}  {'drift':>6}  {'n_dom':>5}  "
          f"{'IC%':>5}  {'OP%':>5}  {'OL%':>5}  {'verdict':>16}")
    print(f"  " + "-" * 52)

    for n in ns:
        for d in drifts:
            dom_set = set()
            agg = {m: 0.0 for m in MECHS}
            count = 0
            for ang in angles:
                for r in radii:
                    f = get_fracs(hcube, n, ang, r, d)
                    if f is None:
                        continue
                    dom_set.add(dominant(f))
                    for m in MECHS:
                        agg[m] += f[m]
                    count += 1
            if count == 0:
                continue
            for m in MECHS:
                agg[m] /= count

            if len(dom_set) <= 2:
                verdict = "BINARY"
            elif len(dom_set) == 3:
                verdict = "NEAR-BINARY"
            else:
                verdict = "MULTI-MECHANISM"

            print(f"  {n:3d}  {d:6.3f}  {len(dom_set):5d}  "
                  f"{agg['inward-collapse']:4.0f}%"
                  f"  {agg['outward-PBC']:4.0f}%"
                  f"  {agg['other-late']:4.0f}%"
                  f"  {verdict:>16}")
        print()


# ============================================================
# Section 7: Fold / cusp creation under drift
# ============================================================

def check_drift_folds(hcube, ns, angles, radii, drifts):
    print("=" * 90)
    print("  FOLD / CUSP CREATION UNDER DRIFT")
    print("=" * 90)
    print(f"  {'N':>3}  {'drift':>6}  {'folds':>12}  {'cusps':>12}")
    print(f"  " + "-" * 38)

    for n in ns:
        for d in drifts:
            fold_count = 0
            cusp_count = 0
            total_r = 0

            for r in radii:
                dom_seq = []
                for ang in angles:
                    f = get_fracs(hcube, n, ang, r, d)
                    if f is None:
                        dom_seq.append(None)
                    else:
                        dom_seq.append(dominant(f))

                valid = [x for x in dom_seq if x is not None]
                if len(valid) < 3:
                    continue
                total_r += 1

                # Fold: A -> B -> A
                for i in range(len(valid) - 2):
                    if valid[i] != valid[i+1] and valid[i] == valid[i+2]:
                        fold_count += 1
                        break

                # Cusp: 3+ distinct mechanisms
                if len(set(valid)) >= 3:
                    cusp_count += 1

            if total_r == 0:
                continue
            print(f"  {n:3d}  {d:6.3f}  "
                  f"{fold_count:4d}/{total_r:<4d}  "
                  f"{cusp_count:4d}/{total_r:<4d}")
        print()


# ============================================================
# Section 8: Drift sensitivity — mechanism change rate
# ============================================================

def check_drift_sensitivity(hcube, ns, angles, radii, drifts):
    print("=" * 90)
    print("  DRIFT SENSITIVITY: mechanism changes from d=0.00 baseline")
    print("=" * 90)
    print(f"  {'N':>3}  {'drift':>6}  {'changed':>12}  {'pct':>5}  "
          f"{'IC shift':>9}  {'stability':>14}")
    print(f"  " + "-" * 54)

    for n in ns:
        for d in drifts:
            if d == 0.0:
                continue
            changes = 0
            total = 0
            ic_shifts = []
            for ang in angles:
                for r in radii:
                    f0 = get_fracs(hcube, n, ang, r, 0.0)
                    fd = get_fracs(hcube, n, ang, r, d)
                    if f0 is None or fd is None:
                        continue
                    total += 1
                    if dominant(f0) != dominant(fd):
                        changes += 1
                    ic_shifts.append(
                        fd["inward-collapse"] - f0["inward-collapse"])

            if total == 0:
                continue
            pct = 100 * changes / total
            mn_shift = sum(ic_shifts) / len(ic_shifts) if ic_shifts else 0

            if pct < 5:
                stab = "RIGID"
            elif pct < 15:
                stab = "STABLE"
            elif pct < 30:
                stab = "SENSITIVE"
            else:
                stab = "UNSTABLE"

            print(f"  {n:3d}  {d:6.3f}  "
                  f"{changes:5d}/{total:<5d}  "
                  f"{pct:4.0f}%  "
                  f"{mn_shift:+8.1f}pp  "
                  f"{stab:>14}")
        print()


# ============================================================
# Section 9: Cross-N drift summary table
# ============================================================

def print_summary(hcube, ns, angles, radii, drifts):
    print("=" * 90)
    print("  CROSS-N DRIFT ROBUSTNESS SUMMARY")
    print("=" * 90)
    print(f"  {'N':>3}  {'drift':>6}  {'boundary':>10}  "
          f"{'DECAY':>5}  {'discont':>8}  {'anom%':>6}  {'n_dom':>5}  "
          f"{'folds%':>6}")
    print(f"  " + "-" * 62)

    for n in ns:
        for d in drifts:
            # Check for data
            data_count = sum(
                1 for ang in angles for r in radii
                if get_fracs(hcube, n, ang, r, d) is not None)
            if data_count == 0:
                continue

            # Boundary: IC onset angles
            onsets = []
            for r in radii:
                for ang in angles:
                    f = get_fracs(hcube, n, ang, r, d)
                    if f is not None and f["inward-collapse"] > 50:
                        onsets.append(ang)
                        break

            if onsets:
                uniq = sorted(set(onsets))
                bnd = (f"{uniq[0]:+d}" if len(uniq) == 1
                       else f"[{min(uniq):+d}..{max(uniq):+d}]")
            else:
                bnd = "none"

            # DECAY count
            decay_n = sum(
                1 for ang in angles for r in radii
                if (get_fracs(hcube, n, ang, r, d) or {}).get("DECAY", 0) > 50)

            # Discontinuity at 0->+1
            jumps = []
            for r in radii:
                f0 = get_fracs(hcube, n, 0, r, d)
                f1 = get_fracs(hcube, n, 1, r, d)
                if f0 and f1:
                    jumps.append(f1["inward-collapse"] -
                                 f0["inward-collapse"])
            n_disc = sum(1 for j in jumps if j > 10)
            disc_s = f"{n_disc}/{len(jumps)}" if jumps else "?"

            # Anomaly footprint
            anom = 0
            anom_tot = 0
            for ang in angles:
                for r in radii:
                    f = get_fracs(hcube, n, ang, r, d)
                    if f is None:
                        continue
                    anom_tot += 1
                    if (f["other-late"] > 20 or
                        f["PBC-corner"] > 20 or
                        (ang > 0 and f["inward-collapse"] < 30)):
                        anom += 1
            anom_pct = f"{100*anom/anom_tot:.0f}%" if anom_tot else "?"

            # Distinct dominant mechanisms
            dom_set = set()
            for ang in angles:
                for r in radii:
                    f = get_fracs(hcube, n, ang, r, d)
                    if f:
                        dom_set.add(dominant(f))

            # Folds
            fold_count = 0
            total_r = 0
            for r in radii:
                valid = []
                for ang in angles:
                    f = get_fracs(hcube, n, ang, r, d)
                    if f:
                        valid.append(dominant(f))
                if len(valid) < 3:
                    continue
                total_r += 1
                for i in range(len(valid) - 2):
                    if valid[i] != valid[i+1] and valid[i] == valid[i+2]:
                        fold_count += 1
                        break
            fold_pct = (f"{100*fold_count/total_r:.0f}%"
                        if total_r > 0 else "?")

            print(f"  {n:3d}  {d:6.3f}  {bnd:>10}  "
                  f"{decay_n:>5}  {disc_s:>8}  {anom_pct:>6}  "
                  f"{len(dom_set):>5}  {fold_pct:>6}")
        print()


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 90)
    print("  4D DRIFT-ROBUSTNESS ANALYSIS: N x angle x radius x drift")
    print("=" * 90)
    print()

    hcube, ns, angles, radii, drifts = build_hypercube()

    if not hcube:
        print("  ERROR: No data loaded.")
        print(f"  Expected 4D pattern: "
              f"n{{N}}_angle_{{label}}_r{{RRR}}_d{{DDD}}_sweep.json")
        print(f"  Or 3D/2D fallbacks for drift=0 baseline.")
        return

    # Section 1: Phase maps per (N, drift)
    print_drift_phase_maps(hcube, ns, angles, radii, drifts)

    # Section 2: Mechanism fractions by drift
    print_drift_profiles(hcube, ns, angles, radii, drifts)

    # Section 3: Tangent discontinuity vs drift
    check_tangent_discontinuity_vs_drift(hcube, ns, radii, drifts)

    # Section 4: DECAY survival
    check_decay_vs_drift(hcube, ns, radii, drifts)

    # Section 5: N=8 anomaly vs drift
    check_n8_anomaly_vs_drift(hcube, ns, angles, radii, drifts)

    # Section 6: Large-N binary partitioning
    check_large_n_stability(hcube, ns, angles, radii, drifts)

    # Section 7: Fold/cusp creation
    check_drift_folds(hcube, ns, angles, radii, drifts)

    # Section 8: Drift sensitivity
    check_drift_sensitivity(hcube, ns, angles, radii, drifts)

    # Section 9: Cross-N summary
    print_summary(hcube, ns, angles, radii, drifts)

    print("=" * 90)


if __name__ == "__main__":
    main()
