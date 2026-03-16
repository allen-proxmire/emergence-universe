#!/usr/bin/env python
"""
check_3d_manifold.py
=====================
3D manifold analysis across N x angle x radius.

Maps the global phase-boundary structure separating:
  - collapse-dominated (IC) regions
  - PBC/OPBC-dominated regions
  - DECAY (if present)
  - other-late regions

Data sources (in priority order):
  1. Per-point sweep files:  n{N}_angle_{label}_{rlabel}_sweep.json
  2. Existing 2D grid files: n{N}_2d_angle_radius_grid.json
"""

import json
import os
import glob
import re
from collections import Counter

SWEEP_DIR = os.path.dirname(__file__)

# Expected parameter ranges
NS = [4, 6, 8, 10, 12, 16, 20]
ANGLES = [-10, -5, -2, -1, 0, 1, 2, 5, 10]
RADII_FRAC = [round(r / 100, 2) for r in range(5, 100, 5)]  # 0.05 .. 0.95
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
    """0.30 -> 'r030', 0.05 -> 'r005'"""
    return f"r{int(round(r_frac * 100)):03d}"

def radius_from_d_px(d_px):
    """Convert d_px to fractional radius (d_px / BOX_PX)."""
    return round(d_px / BOX_PX, 4)

def d_px_from_radius(r_frac):
    """Convert fractional radius back to d_px."""
    return int(round(r_frac * BOX_PX))


def load_3d_file(n, ang, r_frac):
    """Load a per-point sweep file. Returns list of records or None."""
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
    """Compute mechanism fractions (%) from a list of records."""
    total = len(records)
    if total == 0:
        return {m: 0.0 for m in MECHS}
    counts = Counter(r["mechanism"] for r in records)
    return {m: 100.0 * counts.get(m, 0) / total for m in MECHS}


def dominant(fracs_dict):
    """Return the mechanism with the highest fraction."""
    return max(fracs_dict, key=fracs_dict.get)


# ============================================================
# Build the 3D data cube
# ============================================================

def build_cube():
    """
    Build cube[n][angle][radius] = {mech: frac, ...}
    Returns (cube, loaded_ns, loaded_angles, loaded_radii).
    """
    cube = {}
    all_radii = set()
    all_angles = set()
    all_ns = set()

    # --- Pass 1: try per-point 3D sweep files ---
    count_3d = 0
    for n in NS:
        for ang in ANGLES:
            for r in RADII_FRAC:
                records = load_3d_file(n, ang, r)
                if records is not None:
                    cube.setdefault(n, {}).setdefault(ang, {})[r] = mech_fracs(records)
                    all_ns.add(n); all_angles.add(ang); all_radii.add(r)
                    count_3d += 1

    # --- Pass 2: backfill from existing 2D grids ---
    count_2d = 0
    for n in NS:
        grid_data = load_2d_grid(n)
        if grid_data is None:
            continue
        for rec in grid_data["grid"]:
            ang = rec["angle"]
            d_px = rec["d_px"]
            r = radius_from_d_px(d_px)
            # Only fill if not already loaded from 3D files
            if n in cube and ang in cube[n] and r in cube[n][ang]:
                continue
            # Single mechanism per grid cell → 100% for that mechanism
            fracs = {m: 0.0 for m in MECHS}
            fracs[rec["mechanism"]] = 100.0
            cube.setdefault(n, {}).setdefault(ang, {})[r] = fracs
            all_ns.add(n); all_angles.add(ang); all_radii.add(r)
            count_2d += 1

    # --- Also scan for any 3D files with N values we didn't expect ---
    pattern = os.path.join(SWEEP_DIR, "n*_angle_*_r*_sweep.json")
    for path in glob.glob(pattern):
        fname = os.path.basename(path)
        m = re.match(r"n(\d+)_angle_([pm]\d+deg|0deg)_r(\d+)_sweep\.json", fname)
        if not m:
            continue
        n = int(m.group(1))
        ang_str = m.group(2)
        r_int = int(m.group(3))
        # Decode angle
        if ang_str == "0deg":
            ang = 0
        elif ang_str.startswith("p"):
            ang = int(ang_str[1:-3])
        else:
            ang = -int(ang_str[1:-3])
        r = round(r_int / 100, 2)
        if n in cube and ang in cube.get(n, {}) and r in cube.get(n, {}).get(ang, {}):
            continue
        with open(path) as f:
            data = json.load(f)
        if isinstance(data, dict) and "results" in data:
            data = data["results"]
        if isinstance(data, dict):
            data = [data]
        cube.setdefault(n, {}).setdefault(ang, {})[r] = mech_fracs(data)
        all_ns.add(n); all_angles.add(ang); all_radii.add(r)

    ns = sorted(all_ns)
    angles = sorted(all_angles)
    radii = sorted(all_radii)

    print(f"  Data loaded: {count_3d} from 3D files, {count_2d} from 2D grids")
    print(f"  Coverage: {len(ns)} N-values, {len(angles)} angles, {len(radii)} radii")
    total_cells = sum(
        sum(len(cube[n].get(a, {})) for a in angles)
        for n in ns
    )
    max_cells = len(ns) * len(angles) * len(radii)
    print(f"  Filled cells: {total_cells} / {max_cells} "
          f"({100 * total_cells / max_cells:.1f}% if uniform grid)")
    print()
    return cube, ns, angles, radii


# ============================================================
# Section 1: Per-N ASCII phase maps (angle x radius)
# ============================================================

def print_phase_maps(cube, ns, angles, radii):
    for n in ns:
        if n not in cube:
            continue
        print("=" * 90)
        print(f"  N = {n}   PHASE MAP (angle x radius)")
        print(f"  Codes: I=IC  O=OPBC  D=DECAY  C=PBC-corner  .=other-late  ?=missing")
        print("=" * 90)

        # Header: radii (show subset to fit screen)
        # Show every other radius if there are too many
        show_radii = radii if len(radii) <= 20 else radii[::2]
        hdr = "  angle "
        for r in show_radii:
            hdr += f" {r:.2f}"[1:]  # trim leading 0
        print(hdr)
        print("  " + "-" * (7 + 5 * len(show_radii)))

        for ang in angles:
            row = f"  {ang:+4d}  "
            for r in show_radii:
                fracs = cube.get(n, {}).get(ang, {}).get(r)
                if fracs is None:
                    row += "   ? "
                else:
                    dom = dominant(fracs)
                    row += f"   {CODE.get(dom, '?')} "
            if ang == 0:
                row += "  <-- tangent"
            print(row)
        print()


# ============================================================
# Section 2: Mechanism fractions per N (aggregated over radii)
# ============================================================

def print_n_profiles(cube, ns, angles, radii):
    print("=" * 90)
    print("  MECHANISM FRACTIONS BY N (averaged over all radii)")
    print("=" * 90)
    print(f"  {'N':>3}  {'ang':>4}  {'IC':>6}  {'OP':>6}  {'DE':>6}  {'CR':>6}  {'OL':>6}  {'n_pts':>5}")
    print("  " + "-" * 50)

    for n in ns:
        if n not in cube:
            continue
        for ang in angles:
            agg = {m: 0.0 for m in MECHS}
            count = 0
            for r in radii:
                f = cube.get(n, {}).get(ang, {}).get(r)
                if f is not None:
                    for m in MECHS:
                        agg[m] += f[m]
                    count += 1
            if count == 0:
                continue
            for m in MECHS:
                agg[m] /= count
            print(f"  {n:3d}  {ang:+4d}  {agg['inward-collapse']:5.1f}%"
                  f"  {agg['outward-PBC']:5.1f}%  {agg['DECAY']:5.1f}%"
                  f"  {agg['PBC-corner']:5.1f}%  {agg['other-late']:5.1f}%  {count:5d}")
        print()


# ============================================================
# Section 3: Phase boundary extraction per N
# ============================================================

def extract_boundaries(cube, ns, angles, radii):
    """For each N and radius, find the angle where IC first appears (onset)."""
    print("=" * 90)
    print("  PHASE BOUNDARY: IC onset angle per (N, radius)")
    print("=" * 90)

    boundaries = {}  # boundaries[n] = [(radius, onset_angle), ...]

    for n in ns:
        if n not in cube:
            continue
        boundaries[n] = []
        for r in radii:
            onset = None
            for ang in angles:
                f = cube.get(n, {}).get(ang, {}).get(r)
                if f is not None and f["inward-collapse"] > 50:
                    onset = ang
                    break
            boundaries[n].append((r, onset))

        # Print compact boundary for this N
        print(f"\n  N = {n}:")
        print(f"  {'radius':>7}  {'onset':>7}  {'dom@0':>6}  {'dom@+1':>7}")
        print(f"  " + "-" * 32)
        for r, onset in boundaries[n]:
            f0 = cube.get(n, {}).get(0, {}).get(r)
            f1 = cube.get(n, {}).get(1, {}).get(r)
            dom0 = CODE.get(dominant(f0), "?") if f0 else "?"
            dom1 = CODE.get(dominant(f1), "?") if f1 else "?"
            onset_s = f"{onset:+3d}deg" if onset is not None else "  none"
            print(f"  {r:7.3f}  {onset_s:>7}  {dom0:>6}  {dom1:>7}")

    return boundaries


# ============================================================
# Section 4: Boundary geometry analysis
# ============================================================

def analyze_boundary_geometry(boundaries, ns):
    print("\n" + "=" * 90)
    print("  BOUNDARY GEOMETRY ANALYSIS")
    print("=" * 90)
    print(f"  {'N':>3}  {'n_pts':>5}  {'n_none':>6}  {'unique':>20}"
          f"  {'tilt':>10}  {'shape':>16}")
    print("  " + "-" * 72)

    for n in ns:
        if n not in boundaries:
            continue
        pairs = boundaries[n]
        valid = [(r, o) for r, o in pairs if o is not None]
        n_none = len(pairs) - len(valid)

        if not valid:
            print(f"  {n:3d}  {len(pairs):5d}  {n_none:6d}  {'(no IC anywhere)':>20}"
                  f"  {'---':>10}  {'NO BOUNDARY':>16}")
            continue

        onset_vals = [o for _, o in valid]
        unique = sorted(set(onset_vals))

        # Tilt: correlation of onset angle with radius
        if len(valid) >= 2:
            r_vals = [r for r, _ in valid]
            mean_r = sum(r_vals) / len(r_vals)
            mean_o = sum(onset_vals) / len(onset_vals)
            cov = sum((r - mean_r) * (o - mean_o)
                      for r, o in valid) / len(valid)
            var_r = sum((r - mean_r) ** 2 for r in r_vals) / len(r_vals)
            tilt = cov / var_r if var_r > 1e-12 else 0.0
        else:
            tilt = 0.0

        if len(unique) == 1:
            shape = "VERTICAL"
        elif abs(tilt) < 0.5:
            shape = "NEAR-VERTICAL"
        elif abs(tilt) < 3.0:
            shape = "SLIGHT TILT"
        elif len(unique) <= 3:
            shape = "PIECEWISE"
        else:
            shape = "CURVED"

        ustr = str(unique) if len(unique) <= 4 else f"[{min(unique)}..{max(unique)}]"
        print(f"  {n:3d}  {len(pairs):5d}  {n_none:6d}  {ustr:>20}"
              f"  {tilt:+9.3f}  {shape:>16}")

    # Cross-N boundary continuity
    print()
    print("  --- Cross-N boundary continuity ---")
    print(f"  {'N':>3}  {'mean onset':>11}  {'std onset':>10}  {'continuous?':>12}")
    print("  " + "-" * 40)
    prev_mean = None
    for n in ns:
        if n not in boundaries:
            continue
        valid = [o for _, o in boundaries[n] if o is not None]
        if not valid:
            print(f"  {n:3d}  {'---':>11}  {'---':>10}  {'---':>12}")
            prev_mean = None
            continue
        import statistics
        mn = statistics.mean(valid)
        sd = statistics.stdev(valid) if len(valid) > 1 else 0.0
        if prev_mean is not None:
            gap = abs(mn - prev_mean)
            cont = "YES" if gap < 3 else "JUMP"
        else:
            cont = "---"
        print(f"  {n:3d}  {mn:+10.2f}  {sd:9.2f}  {cont:>12}")
        prev_mean = mn


# ============================================================
# Section 5: Tangent discontinuity persistence
# ============================================================

def check_tangent_discontinuity(cube, ns, radii):
    print("\n" + "=" * 90)
    print("  TANGENT DISCONTINUITY: IC jump from 0 deg to +1 deg across all radii")
    print("=" * 90)

    for n in ns:
        if n not in cube:
            continue
        jumps = []
        for r in radii:
            f0 = cube.get(n, {}).get(0, {}).get(r)
            f1 = cube.get(n, {}).get(1, {}).get(r)
            if f0 is not None and f1 is not None:
                jump = f1["inward-collapse"] - f0["inward-collapse"]
                jumps.append((r, jump))

        if not jumps:
            print(f"\n  N={n}: no data at angle=0 and +1")
            continue

        jump_vals = [j for _, j in jumps]
        n_pos = sum(1 for j in jump_vals if j > 10)
        n_zero = sum(1 for j in jump_vals if abs(j) < 5)
        n_neg = sum(1 for j in jump_vals if j < -10)

        print(f"\n  N={n}: {len(jumps)} radii tested")
        print(f"    Jump > +10 pp: {n_pos}/{len(jumps)} "
              f"({100*n_pos/len(jumps):.0f}%)")
        print(f"    Jump ~ 0:      {n_zero}/{len(jumps)} "
              f"({100*n_zero/len(jumps):.0f}%)")
        print(f"    Jump < -10 pp: {n_neg}/{len(jumps)} "
              f"({100*n_neg/len(jumps):.0f}%)")

        if jump_vals:
            mn = sum(jump_vals) / len(jump_vals)
            mx = max(jump_vals)
            mi = min(jump_vals)
            print(f"    Mean jump: {mn:+.1f} pp  "
                  f"range: [{mi:+.1f}, {mx:+.1f}]")

            if n_pos > 0.8 * len(jumps):
                print(f"    -> UNIVERSAL discontinuity (vertical sheet in 3D)")
            elif n_pos > 0.5 * len(jumps):
                print(f"    -> PARTIAL discontinuity (sheet with holes)")
            else:
                print(f"    -> WEAK or ABSENT discontinuity")


# ============================================================
# Section 6: DECAY isolation in 3D
# ============================================================

def check_decay_isolation(cube, ns, angles, radii):
    print("\n" + "=" * 90)
    print("  DECAY ISOLATION IN 3D: where does DECAY appear?")
    print("=" * 90)

    decay_cells = []
    total_cells = 0
    for n in ns:
        for ang in angles:
            for r in radii:
                f = cube.get(n, {}).get(ang, {}).get(r)
                if f is None:
                    continue
                total_cells += 1
                if f["DECAY"] > 50:
                    decay_cells.append((n, ang, r, f["DECAY"]))

    print(f"\n  Total cells examined: {total_cells}")
    print(f"  DECAY-dominant cells (>50%): {len(decay_cells)}")

    if not decay_cells:
        print("  -> NO DECAY cells found in the 3D volume")
        return

    # Group by N
    by_n = {}
    for n, ang, r, dec in decay_cells:
        by_n.setdefault(n, []).append((ang, r, dec))

    for n in sorted(by_n):
        cells = by_n[n]
        angs_d = sorted(set(a for a, _, _ in cells))
        rads_d = sorted(set(r for _, r, _ in cells))
        print(f"\n  N={n}: {len(cells)} DECAY cells")
        print(f"    Angles:  {angs_d}")
        print(f"    Radii:   [{min(rads_d):.2f} .. {max(rads_d):.2f}] "
              f"({len(rads_d)} values)")
        if len(angs_d) == 1 and angs_d[0] == 0:
            print(f"    -> DECAY confined to angle=0 (1D filament in 3D)")
        elif len(angs_d) <= 3:
            print(f"    -> DECAY near-confined to angles {angs_d}")
        else:
            print(f"    -> DECAY extends across angles: "
                  f"possible 2D surface in 3D")

    # 3D isolation assessment
    n_vals = sorted(set(n for n, _, _, _ in decay_cells))
    if len(n_vals) == 1 and n_vals[0] == 4:
        a_vals = set(a for _, a, _, _ in decay_cells)
        if a_vals == {0}:
            print(f"\n  GLOBAL: DECAY exists ONLY at (N=4, angle=0, r=*)")
            print(f"  -> ISOLATED 1D FILAMENT in the 3D volume")
            print(f"     (a line, not a point -- extends along radius axis)")
        else:
            print(f"\n  GLOBAL: DECAY exists only at N=4")
    elif len(n_vals) == 1:
        print(f"\n  GLOBAL: DECAY exists only at N={n_vals[0]}")
    else:
        print(f"\n  GLOBAL: DECAY found at N = {n_vals}")


# ============================================================
# Section 7: N=8 anomaly structure
# ============================================================

def check_n8_anomaly(cube, ns, angles, radii):
    print("\n" + "=" * 90)
    print("  N=8 ANOMALY STRUCTURE: ridge, plateau, or bulge?")
    print("=" * 90)

    # Compare IC fraction at angle=+5 across N to see if N=8 forms a
    # dip (anomalous low IC) relative to neighbors
    test_angles = [1, 2, 5]
    for test_ang in test_angles:
        print(f"\n  IC fraction at angle={test_ang:+d} deg across N:")
        print(f"  {'N':>3}  {'mean IC':>8}  {'min IC':>7}  {'max IC':>7}  {'n_pts':>5}")
        print("  " + "-" * 36)
        profile = {}
        for n in ns:
            ics = []
            for r in radii:
                f = cube.get(n, {}).get(test_ang, {}).get(r)
                if f is not None:
                    ics.append(f["inward-collapse"])
            if ics:
                profile[n] = ics
                mn = sum(ics) / len(ics)
                print(f"  {n:3d}  {mn:7.1f}%  {min(ics):6.1f}%"
                      f"  {max(ics):6.1f}%  {len(ics):5d}")

    # Check other-late inflation at N=8 vs neighbors
    print(f"\n  Other-late fraction at angle=+1 deg across N:")
    print(f"  {'N':>3}  {'mean OL':>8}  {'mean CR':>8}")
    print("  " + "-" * 24)
    for n in ns:
        ols, crs = [], []
        for r in radii:
            f = cube.get(n, {}).get(1, {}).get(r)
            if f is not None:
                ols.append(f["other-late"])
                crs.append(f["PBC-corner"])
        if ols:
            print(f"  {n:3d}  {sum(ols)/len(ols):7.1f}%  {sum(crs)/len(crs):7.1f}%")

    # Anomaly footprint: for N=8, count cells where OL > 20% or IC < 30%
    if 8 in cube:
        anomaly_cells = 0
        total = 0
        for ang in angles:
            for r in radii:
                f = cube.get(8, {}).get(ang, {}).get(r)
                if f is None:
                    continue
                total += 1
                if (f["other-late"] > 20 or
                    f["PBC-corner"] > 20 or
                    (ang > 0 and f["inward-collapse"] < 30)):
                    anomaly_cells += 1
        if total > 0:
            print(f"\n  N=8 anomaly footprint: {anomaly_cells}/{total} cells "
                  f"({100*anomaly_cells/total:.0f}%)")
            if anomaly_cells / total > 0.3:
                print(f"  -> PLATEAU (anomaly covers >30% of angle-radius plane)")
            elif anomaly_cells / total > 0.1:
                print(f"  -> RIDGE (anomaly covers 10-30%)")
            else:
                print(f"  -> LOCALIZED (anomaly covers <10%)")


# ============================================================
# Section 8: Large-N flattening
# ============================================================

def check_large_n_flattening(cube, ns, angles, radii):
    print("\n" + "=" * 90)
    print("  LARGE-N FLATTENING: does the manifold simplify at large N?")
    print("=" * 90)

    print(f"\n  Mechanism diversity per (N, angle) slice — counted across radii:")
    print(f"  {'N':>3}  {'ang':>4}  {'n_mech':>6}  {'dominant':>9}  {'dom%':>5}  {'IC range':>9}")
    print("  " + "-" * 44)

    for n in ns:
        if n not in cube:
            continue
        for ang in angles:
            agg = {m: 0.0 for m in MECHS}
            count = 0
            ic_vals = []
            for r in radii:
                f = cube.get(n, {}).get(ang, {}).get(r)
                if f is None:
                    continue
                for m in MECHS:
                    agg[m] += f[m]
                count += 1
                ic_vals.append(f["inward-collapse"])
            if count == 0:
                continue
            for m in MECHS:
                agg[m] /= count
            n_active = sum(1 for m in MECHS if agg[m] > 5)
            dom = dominant(agg)
            ic_range = max(ic_vals) - min(ic_vals)
            print(f"  {n:3d}  {ang:+4d}  {n_active:6d}  {SHORT[dom]:>9}"
                  f"  {agg[dom]:4.0f}%  {ic_range:8.0f}")
        print()

    # Check if positive-angle IC converges across N
    print(f"\n  IC at angle=+5 deg, radius=0.50 across N:")
    for n in ns:
        r_test = min(radii, key=lambda r: abs(r - 0.50))
        f = cube.get(n, {}).get(5, {}).get(r_test)
        if f:
            print(f"    N={n:2d}: IC={f['inward-collapse']:5.1f}%  "
                  f"OP={f['outward-PBC']:5.1f}%  "
                  f"OL={f['other-late']:5.1f}%")


# ============================================================
# Section 9: Fold / cusp / multi-valued region detection
# ============================================================

def check_folds_and_cusps(cube, ns, angles, radii):
    print("\n" + "=" * 90)
    print("  FOLD / CUSP / MULTI-VALUED REGION DETECTION")
    print("=" * 90)

    for n in ns:
        if n not in cube:
            continue
        # For each radius, check if the dominant mechanism changes
        # non-monotonically as angle increases (indicates a fold)
        fold_count = 0
        cusp_count = 0
        total_r = 0

        for r in radii:
            dom_seq = []
            for ang in angles:
                f = cube.get(n, {}).get(ang, {}).get(r)
                if f is None:
                    dom_seq.append(None)
                else:
                    dom_seq.append(dominant(f))

            valid = [(a, d) for a, d in zip(angles, dom_seq) if d is not None]
            if len(valid) < 3:
                continue
            total_r += 1

            # Count mechanism transitions
            transitions = []
            for i in range(len(valid) - 1):
                if valid[i][1] != valid[i+1][1]:
                    transitions.append((valid[i][0], valid[i+1][0],
                                        valid[i][1], valid[i+1][1]))

            # A fold: mechanism A -> B -> A (returns to same mechanism)
            mechs_seen = [d for _, d in valid]
            for i in range(len(mechs_seen) - 2):
                if (mechs_seen[i] != mechs_seen[i+1] and
                    mechs_seen[i] == mechs_seen[i+2]):
                    fold_count += 1
                    break

            # A cusp: three or more distinct mechanisms in one radius slice
            unique_mechs = set(d for _, d in valid)
            if len(unique_mechs) >= 3:
                cusp_count += 1

        if total_r > 0:
            print(f"\n  N={n}: {total_r} radii examined")
            print(f"    Folds (A->B->A):  {fold_count}/{total_r} "
                  f"({100*fold_count/total_r:.0f}%)")
            print(f"    Cusps (3+ mechs): {cusp_count}/{total_r} "
                  f"({100*cusp_count/total_r:.0f}%)")


# ============================================================
# Section 10: Cross-N summary table
# ============================================================

def print_cross_n_summary(cube, boundaries, ns, angles, radii):
    print("\n" + "=" * 90)
    print("  CROSS-N SUMMARY TABLE")
    print("=" * 90)
    print(f"  {'N':>3}  {'boundary':>10}  {'shape':>14}  "
          f"{'DECAY':>6}  {'discont':>8}  {'n_mech':>7}  anomaly")
    print("  " + "-" * 72)

    for n in ns:
        if n not in cube:
            continue

        # Boundary shape
        if n in boundaries:
            valid = [o for _, o in boundaries[n] if o is not None]
            unique = sorted(set(valid)) if valid else []
            if not unique:
                bnd = "none"
                shape = "---"
            elif len(unique) == 1:
                bnd = f"{unique[0]:+d}deg"
                shape = "VERTICAL"
            else:
                bnd = f"[{min(unique):+d}..{max(unique):+d}]"
                shape = "PIECEWISE" if len(unique) <= 3 else "CURVED"
        else:
            bnd = "?"
            shape = "?"

        # DECAY count
        decay_n = 0
        for ang in angles:
            for r in radii:
                f = cube.get(n, {}).get(ang, {}).get(r)
                if f and f["DECAY"] > 50:
                    decay_n += 1

        # Discontinuity at 0->+1
        jumps = []
        for r in radii:
            f0 = cube.get(n, {}).get(0, {}).get(r)
            f1 = cube.get(n, {}).get(1, {}).get(r)
            if f0 and f1:
                jumps.append(f1["inward-collapse"] - f0["inward-collapse"])
        disc = f"{sum(j > 10 for j in jumps)}/{len(jumps)}" if jumps else "?"

        # Mechanism diversity: count distinct dominant mechanisms across all cells
        dom_set = set()
        for ang in angles:
            for r in radii:
                f = cube.get(n, {}).get(ang, {}).get(r)
                if f:
                    dom_set.add(dominant(f))
        mn_div = f"{len(dom_set):d}" if dom_set else "?"

        # Anomaly classification
        if n == 4 and decay_n > 0:
            anom = "tangent DECAY"
        elif n == 8:
            anom = "4k catastrophic"
        elif n in (12, 16) and any(
            cube.get(n, {}).get(a, {}).get(r, {}).get("other-late", 0) > 20
            for a in angles for r in radii
            if cube.get(n, {}).get(a, {}).get(r)):
            anom = "OL inflation"
        else:
            anom = "NORMAL"

        print(f"  {n:3d}  {bnd:>10}  {shape:>14}  "
              f"{decay_n:>6}  {disc:>8}  {mn_div:>7}  {anom}")


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 90)
    print("  3D MANIFOLD ANALYSIS: N x angle x radius")
    print("=" * 90)
    print()

    cube, ns, angles, radii = build_cube()

    if not cube:
        print("  ERROR: No data loaded. Check that sweep files exist.")
        print(f"  Expected pattern: n{{N}}_angle_{{label}}_r{{RRR}}_sweep.json")
        print(f"  Or 2D grids:      n{{N}}_2d_angle_radius_grid.json")
        return

    # Section 1: ASCII phase maps
    print_phase_maps(cube, ns, angles, radii)

    # Section 2: Mechanism fractions by N
    print_n_profiles(cube, ns, angles, radii)

    # Section 3: Phase boundary extraction
    boundaries = extract_boundaries(cube, ns, angles, radii)

    # Section 4: Boundary geometry
    analyze_boundary_geometry(boundaries, ns)

    # Section 5: Tangent discontinuity
    check_tangent_discontinuity(cube, ns, radii)

    # Section 6: DECAY isolation
    check_decay_isolation(cube, ns, angles, radii)

    # Section 7: N=8 anomaly
    check_n8_anomaly(cube, ns, angles, radii)

    # Section 8: Large-N flattening
    check_large_n_flattening(cube, ns, angles, radii)

    # Section 9: Folds and cusps
    check_folds_and_cusps(cube, ns, angles, radii)

    # Section 10: Cross-N summary
    print_cross_n_summary(cube, boundaries, ns, angles, radii)

    print()
    print("=" * 90)


if __name__ == "__main__":
    main()
