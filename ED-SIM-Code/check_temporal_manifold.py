#!/usr/bin/env python
"""
check_temporal_manifold.py
===========================
Temporal structure analysis of the ED phase manifold.

Loads existing 4D sweep files (N x angle x radius x drift) and
analyzes collapse time / orbit lifetime / DECAY timeout across the
parameter space.  Identifies:
  - temporal discontinuities at mechanism boundaries
  - slow manifolds (long-lived orbit regions)
  - metastable regions near the phase boundary
  - drift-dependent temporal smoothing or sharpening

Data sources:
  1. 4D drift sweep files: n{N}_angle_{label}_r{RRR}_d{DDD}_sweep.json
  2. Dedicated time files:  n{N}_angle_{label}_r{RRR}_d{DDD}_time.json
  3. 2D grid files:         n{N}_2d_angle_radius_grid.json
"""

import json
import os
import glob
import re
import statistics
from collections import defaultdict

SWEEP_DIR = os.path.dirname(__file__)

NS     = [4, 8, 12, 20]
ANGLES = [-10, -5, -2, -1, 0, 1, 2, 5, 10]
RADII  = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 0.90, 0.95]
DRIFTS = [0.00, 0.01, 0.02, 0.05, 0.10]

MECHS = ["inward-collapse", "outward-PBC", "DECAY", "PBC-corner", "other-late"]
SHORT = {"inward-collapse": "IC", "outward-PBC": "OP", "DECAY": "DE",
         "PBC-corner": "CR", "other-late": "OL"}

DECAY_TIMEOUT = 100.0   # chi_emp ceiling for DECAY


# ============================================================
# Data loading
# ============================================================

def angle_label(ang):
    if ang > 0:   return f"p{ang}deg"
    if ang < 0:   return f"m{abs(ang)}deg"
    return "0deg"

def radius_label(r):
    return f"r{int(round(r * 100)):03d}"

def drift_label(d):
    return f"d{int(round(d * 1000)):03d}"


def build_cube():
    """
    Build cube[n][angle][radius][drift] = {mechanism, chi_emp}
    Returns (cube, ns, angles, radii, drifts).
    """
    cube = {}
    all_ns, all_angles, all_radii, all_drifts = set(), set(), set(), set()
    count = 0

    # Load 4D sweep files and time files
    for suffix in ("_sweep.json", "_time.json"):
        pattern = os.path.join(SWEEP_DIR, f"n*_angle_*_r*_d*{suffix}")
        for path in glob.glob(pattern):
            fname = os.path.basename(path)
            m = re.match(
                r"n(\d+)_angle_([pm]\d+deg|0deg)_r(\d+)_d(\d+)"
                + re.escape(suffix), fname)
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
            drift = round(d_int / 1000, 3)

            # Skip if already loaded (sweep files take priority over time files)
            if (n in cube and ang in cube.get(n, {}) and
                r in cube.get(n, {}).get(ang, {}) and
                drift in cube.get(n, {}).get(ang, {}).get(r, {})):
                continue

            with open(path) as f:
                data = json.load(f)

            if isinstance(data, list):
                # Multiple records — take mechanism fractions and mean chi
                mechs_list = [rec["mechanism"] for rec in data]
                chis_list = [rec.get("chi_emp", 0) for rec in data]
                from collections import Counter
                mc = Counter(mechs_list)
                dom = mc.most_common(1)[0][0]
                chi = statistics.mean(chis_list) if chis_list else 0
            else:
                dom = data["mechanism"]
                chi = data.get("chi_emp", 0)

            (cube.setdefault(n, {})
                 .setdefault(ang, {})
                 .setdefault(r, {}))[drift] = {"mechanism": dom, "chi": chi}

            all_ns.add(n); all_angles.add(ang)
            all_radii.add(r); all_drifts.add(drift)
            count += 1

    ns = sorted(all_ns)
    angles = sorted(all_angles)
    radii = sorted(all_radii)
    drifts = sorted(all_drifts)

    print(f"  Data loaded: {count} cells")
    print(f"  N: {ns}")
    print(f"  Angles: {angles}")
    print(f"  Radii: {len(radii)} values [{min(radii):.2f}..{max(radii):.2f}]")
    print(f"  Drifts: {drifts}")
    print()
    return cube, ns, angles, radii, drifts


def get(cube, n, ang, r, drift):
    """Safely fetch a cell or None."""
    return cube.get(n, {}).get(ang, {}).get(r, {}).get(drift)


# ============================================================
# Section 1: Temporal heatmaps per N (angle x radius at drift=0)
# ============================================================

def print_temporal_heatmaps(cube, ns, angles, radii):
    print("=" * 90)
    print("  TEMPORAL HEATMAPS: chi_emp at drift=0 (angle x radius)")
    print("  Legend: chi values, colored by lifetime tier")
    print("=" * 90)

    for n in ns:
        print(f"\n  N = {n}   chi_emp (drift=0)")
        # Header
        hdr = "  angle "
        for r in radii:
            hdr += f"  {r:.2f}"[1:]
        print(hdr)
        print("  " + "-" * (7 + 6 * len(radii)))

        for ang in angles:
            row = f"  {ang:+4d}  "
            for r in radii:
                cell = get(cube, n, ang, r, 0.00)
                if cell is None:
                    row += "     ?"
                else:
                    chi = cell["chi"]
                    if chi >= DECAY_TIMEOUT:
                        row += "   DEC"
                    elif chi >= 50:
                        row += f"  {chi:4.0f}"
                    elif chi >= 10:
                        row += f"  {chi:4.1f}"
                    elif chi >= 1:
                        row += f"  {chi:4.1f}"
                    else:
                        row += f"  {chi:4.2f}"
            if ang == 0:
                row += "  <-- tangent"
            print(row)
        print()


# ============================================================
# Section 2: Temporal statistics by mechanism
# ============================================================

def print_mechanism_timing(cube, ns, angles, radii, drifts):
    print("=" * 90)
    print("  TEMPORAL STATISTICS BY MECHANISM (all drifts pooled)")
    print("=" * 90)

    for n in ns:
        by_mech = defaultdict(list)
        for ang in angles:
            for r in radii:
                for d in drifts:
                    cell = get(cube, n, ang, r, d)
                    if cell:
                        by_mech[cell["mechanism"]].append(cell["chi"])

        print(f"\n  N = {n}:")
        print(f"  {'mechanism':>18}  {'count':>5}  {'mean':>8}  "
              f"{'median':>8}  {'min':>8}  {'max':>8}  {'stdev':>8}")
        print(f"  " + "-" * 68)

        for m in MECHS:
            vals = by_mech.get(m, [])
            if not vals:
                continue
            mn = statistics.mean(vals)
            med = statistics.median(vals)
            sd = statistics.stdev(vals) if len(vals) > 1 else 0
            print(f"  {SHORT[m]:>18}  {len(vals):5d}  {mn:8.2f}  "
                  f"{med:8.2f}  {min(vals):8.2f}  {max(vals):8.2f}  {sd:8.2f}")


# ============================================================
# Section 3: Tangent discontinuity — temporal jump at 0 -> +1 deg
# ============================================================

def check_tangent_temporal_jump(cube, ns, radii, drifts):
    print("\n" + "=" * 90)
    print("  TANGENT TEMPORAL JUMP: chi at angle=0 vs angle=+1")
    print("=" * 90)

    for n in ns:
        print(f"\n  N = {n}:")
        print(f"  {'drift':>6}  {'n_pts':>5}  {'chi@0':>8}  {'chi@+1':>8}  "
              f"{'ratio':>8}  {'jump':>10}  verdict")
        print(f"  " + "-" * 62)

        for d in drifts:
            chis_0, chis_1 = [], []
            for r in radii:
                c0 = get(cube, n, ang=0, r=r, drift=d)
                c1 = get(cube, n, ang=1, r=r, drift=d)
                if c0 and c1:
                    chis_0.append(c0["chi"])
                    chis_1.append(c1["chi"])
            if not chis_0:
                continue

            mean0 = statistics.mean(chis_0)
            mean1 = statistics.mean(chis_1)
            ratio = mean0 / mean1 if mean1 > 0.01 else float('inf')
            jump = mean0 - mean1

            if abs(jump) > 20:
                verdict = "SHARP TEMPORAL JUMP"
            elif abs(jump) > 5:
                verdict = "MODERATE JUMP"
            elif ratio > 3 or ratio < 0.33:
                verdict = "SCALE CHANGE"
            else:
                verdict = "SMOOTH"

            print(f"  {d:6.3f}  {len(chis_0):5d}  {mean0:8.2f}  {mean1:8.2f}  "
                  f"{ratio:8.2f}  {jump:+9.2f}  {verdict}")


# ============================================================
# Section 4: DECAY temporal signature
# ============================================================

def check_decay_temporal(cube, ns, angles, radii, drifts):
    print("\n" + "=" * 90)
    print("  DECAY TEMPORAL SIGNATURE")
    print("=" * 90)

    # Collect chi values for DECAY cells vs their nearest non-DECAY neighbors
    for n in ns:
        decay_chis = []
        neighbor_chis = []   # chi of nearest non-DECAY cell (same N, same drift)

        for d in drifts:
            for r in radii:
                for ang in angles:
                    cell = get(cube, n, ang, r, d)
                    if not cell:
                        continue
                    if cell["mechanism"] == "DECAY":
                        decay_chis.append(cell["chi"])
                        # Find neighbor at angle+1
                        next_ang = ang + 1
                        if next_ang in angles or next_ang == 1:
                            cn = get(cube, n, next_ang, r, d)
                            if cn and cn["mechanism"] != "DECAY":
                                neighbor_chis.append(cn["chi"])

        if not decay_chis:
            print(f"\n  N={n}: no DECAY cells")
            continue

        print(f"\n  N={n}: {len(decay_chis)} DECAY cells")
        print(f"    DECAY chi:    mean={statistics.mean(decay_chis):.2f}  "
              f"all at timeout={sum(1 for c in decay_chis if c >= DECAY_TIMEOUT - 1)}"
              f"/{len(decay_chis)}")
        if neighbor_chis:
            print(f"    Neighbor chi:  mean={statistics.mean(neighbor_chis):.2f}  "
                  f"median={statistics.median(neighbor_chis):.2f}")
            ratio = statistics.mean(decay_chis) / statistics.mean(neighbor_chis)
            print(f"    DECAY/neighbor ratio: {ratio:.1f}x")
            if ratio > 5:
                print(f"    -> DECAY is a TEMPORAL CLIFF (sharp lifetime boundary)")
            elif ratio > 2:
                print(f"    -> DECAY is a TEMPORAL STEP")
            else:
                print(f"    -> DECAY has GRADUAL temporal transition")


# ============================================================
# Section 5: Slow manifolds — regions of anomalously long lifetime
# ============================================================

def find_slow_manifolds(cube, ns, angles, radii, drifts):
    print("\n" + "=" * 90)
    print("  SLOW MANIFOLDS: regions of anomalously long lifetime")
    print("=" * 90)

    for n in ns:
        # Collect all non-DECAY chi values for this N
        all_chi = []
        for ang in angles:
            for r in radii:
                for d in drifts:
                    cell = get(cube, n, ang, r, d)
                    if cell and cell["mechanism"] != "DECAY":
                        all_chi.append(cell["chi"])

        if len(all_chi) < 10:
            continue

        mean_chi = statistics.mean(all_chi)
        std_chi = statistics.stdev(all_chi)
        slow_thr = mean_chi + 2 * std_chi  # 2-sigma outliers

        print(f"\n  N={n}: mean chi = {mean_chi:.2f}, "
              f"std = {std_chi:.2f}, slow threshold = {slow_thr:.2f}")

        # Find slow cells
        slow_cells = []
        for ang in angles:
            for r in radii:
                for d in drifts:
                    cell = get(cube, n, ang, r, d)
                    if cell and cell["mechanism"] != "DECAY" and cell["chi"] > slow_thr:
                        slow_cells.append((ang, r, d, cell["chi"], cell["mechanism"]))

        print(f"  Slow cells (chi > {slow_thr:.1f}): {len(slow_cells)}/{len(all_chi)}")

        if slow_cells:
            # Group by mechanism
            by_mech = defaultdict(list)
            for ang, r, d, chi, mech in slow_cells:
                by_mech[mech].append((ang, r, d, chi))

            for mech in sorted(by_mech):
                cells = by_mech[mech]
                angs_s = sorted(set(a for a, _, _, _ in cells))
                drifts_s = sorted(set(d for _, _, d, _ in cells))
                chis_s = [c for _, _, _, c in cells]
                print(f"    {SHORT[mech]}: {len(cells)} cells, "
                      f"angles={angs_s}, drifts={drifts_s}, "
                      f"chi=[{min(chis_s):.1f}..{max(chis_s):.1f}]")

            # Check spatial clustering
            slow_angles = sorted(set(a for a, _, _, _, _ in slow_cells))
            slow_radii = sorted(set(r for _, r, _, _, _ in slow_cells))
            if len(slow_angles) <= 3:
                print(f"  -> Slow manifold LOCALIZED at angles {slow_angles}")
            elif len(slow_radii) <= 3:
                print(f"  -> Slow manifold LOCALIZED at radii {slow_radii}")
            else:
                print(f"  -> Slow manifold DISTRIBUTED across "
                      f"{len(slow_angles)} angles, {len(slow_radii)} radii")


# ============================================================
# Section 6: Temporal boundary sharpness
# ============================================================

def check_temporal_boundaries(cube, ns, angles, radii, drifts):
    print("\n" + "=" * 90)
    print("  TEMPORAL BOUNDARY SHARPNESS: max chi gradient along angle axis")
    print("=" * 90)

    for n in ns:
        print(f"\n  N={n}:")
        print(f"  {'drift':>6}  {'max_grad':>9}  {'at_angle':>9}  "
              f"{'at_r':>6}  {'chi_lo':>7}  {'chi_hi':>7}  verdict")
        print(f"  " + "-" * 56)

        for d in drifts:
            max_grad = 0
            max_loc = (0, 0)
            max_chis = (0, 0)

            for r in radii:
                prev_chi = None
                prev_ang = None
                for ang in ANGLES:
                    cell = get(cube, n, ang, r, d)
                    if cell is None:
                        prev_chi = None
                        continue
                    chi = cell["chi"]
                    if prev_chi is not None:
                        grad = abs(chi - prev_chi)
                        if grad > max_grad:
                            max_grad = grad
                            max_loc = (prev_ang, r)
                            max_chis = (prev_chi, chi)
                    prev_chi = chi
                    prev_ang = ang

            if max_grad > 30:
                verdict = "CLIFF"
            elif max_grad > 10:
                verdict = "STEEP"
            elif max_grad > 3:
                verdict = "MODERATE"
            else:
                verdict = "GENTLE"

            ang_at = f"{max_loc[0]:+d}->{max_loc[0]+1:+d}" if max_loc else "?"
            # Find the actual angle pair
            print(f"  {d:6.3f}  {max_grad:9.2f}  {max_loc[0]:+4d}deg  "
                  f"{max_loc[1]:5.2f}  {max_chis[0]:7.2f}  {max_chis[1]:7.2f}  {verdict}")


# ============================================================
# Section 7: N=8 plateau — temporal footprint
# ============================================================

def check_n8_temporal(cube, ns, angles, radii, drifts):
    print("\n" + "=" * 90)
    print("  N=8 PLATEAU TEMPORAL FOOTPRINT")
    print("=" * 90)

    # Compare mean chi at positive angles across N
    print(f"\n  Mean chi at positive angles (ang > 0) by N and drift:")
    print(f"  {'N':>3}  {'drift':>6}  {'mean_chi':>9}  {'n_pts':>5}  "
          f"{'IC_chi':>8}  {'OL_chi':>8}  {'IC_n':>4}  {'OL_n':>4}")
    print(f"  " + "-" * 60)

    for n in ns:
        for d in [0.00, 0.02, 0.10]:
            if d not in drifts:
                continue
            all_chi, ic_chi, ol_chi = [], [], []
            for ang in angles:
                if ang <= 0:
                    continue
                for r in radii:
                    cell = get(cube, n, ang, r, d)
                    if cell is None:
                        continue
                    all_chi.append(cell["chi"])
                    if cell["mechanism"] == "inward-collapse":
                        ic_chi.append(cell["chi"])
                    elif cell["mechanism"] == "other-late":
                        ol_chi.append(cell["chi"])

            if not all_chi:
                continue

            ic_m = statistics.mean(ic_chi) if ic_chi else 0
            ol_m = statistics.mean(ol_chi) if ol_chi else 0
            print(f"  {n:3d}  {d:6.3f}  {statistics.mean(all_chi):9.2f}  "
                  f"{len(all_chi):5d}  {ic_m:8.2f}  {ol_m:8.2f}  "
                  f"{len(ic_chi):4d}  {len(ol_chi):4d}")
        print()

    # N=8 vs neighbors at drift=0
    if 8 in ns and 4 in ns and 12 in ns:
        print(f"\n  N=8 vs neighbors at drift=0 (positive angles only):")
        for n in [4, 8, 12, 20]:
            chis = []
            for ang in angles:
                if ang <= 0:
                    continue
                for r in radii:
                    cell = get(cube, n, ang, r, 0.00)
                    if cell and cell["mechanism"] != "DECAY":
                        chis.append(cell["chi"])
            if chis:
                print(f"    N={n:2d}: mean={statistics.mean(chis):7.2f}  "
                      f"median={statistics.median(chis):7.2f}  "
                      f"max={max(chis):7.2f}  n={len(chis)}")


# ============================================================
# Section 8: Large-N temporal scaling
# ============================================================

def check_large_n_scaling(cube, ns, angles, radii, drifts):
    print("\n" + "=" * 90)
    print("  LARGE-N TEMPORAL SCALING")
    print("=" * 90)

    # Fixed angle (+5), fixed drift (0.00), vary N and radius
    print(f"\n  Chi at angle=+5 deg, drift=0.00 across N and radius:")
    print(f"  {'N':>3}  {'mean':>8}  {'med':>8}  {'min':>8}  {'max':>8}  {'n':>4}")
    print(f"  " + "-" * 42)

    for n in ns:
        chis = []
        for r in radii:
            cell = get(cube, n, 5, r, 0.00)
            if cell:
                chis.append(cell["chi"])
        if chis:
            print(f"  {n:3d}  {statistics.mean(chis):8.2f}  "
                  f"{statistics.median(chis):8.2f}  "
                  f"{min(chis):8.2f}  {max(chis):8.2f}  {len(chis):4d}")

    # How does chi scale with N for IC specifically?
    print(f"\n  IC-only chi at angle=+5 deg, drift=0.00:")
    print(f"  {'N':>3}  {'mean IC chi':>12}  {'n IC':>5}  "
          f"{'chi/N':>8}  {'chi/N^2':>8}")
    print(f"  " + "-" * 42)

    for n in ns:
        ic_chis = []
        for r in radii:
            cell = get(cube, n, 5, r, 0.00)
            if cell and cell["mechanism"] == "inward-collapse":
                ic_chis.append(cell["chi"])
        if ic_chis:
            mn = statistics.mean(ic_chis)
            print(f"  {n:3d}  {mn:12.2f}  {len(ic_chis):5d}  "
                  f"{mn/n:8.3f}  {mn/(n*n):8.4f}")

    # Drift effect on chi across N
    print(f"\n  Mean IC chi vs drift (angle=+5 deg):")
    print(f"  {'N':>3}", end="")
    for d in drifts:
        print(f"  d={d:.2f}", end="")
    print()
    print(f"  " + "-" * (4 + 8 * len(drifts)))

    for n in ns:
        row = f"  {n:3d}"
        for d in drifts:
            ic_chis = []
            for r in radii:
                cell = get(cube, n, 5, r, d)
                if cell and cell["mechanism"] == "inward-collapse":
                    ic_chis.append(cell["chi"])
            if ic_chis:
                row += f"  {statistics.mean(ic_chis):6.2f}"
            else:
                row += "     ---"
        print(row)


# ============================================================
# Section 9: Metastable regions near phase boundary
# ============================================================

def find_metastable_regions(cube, ns, angles, radii, drifts):
    print("\n" + "=" * 90)
    print("  METASTABLE REGIONS: high chi near mechanism transitions")
    print("=" * 90)

    for n in ns:
        meta_cells = []

        for d in drifts:
            for r in radii:
                for i, ang in enumerate(angles):
                    cell = get(cube, n, ang, r, d)
                    if cell is None:
                        continue

                    # Check if this cell borders a different mechanism
                    is_boundary = False
                    for neighbor_ang in [angles[i-1] if i > 0 else None,
                                         angles[i+1] if i < len(angles)-1 else None]:
                        if neighbor_ang is None:
                            continue
                        nc = get(cube, n, neighbor_ang, r, d)
                        if nc and nc["mechanism"] != cell["mechanism"]:
                            is_boundary = True
                            break

                    if is_boundary and cell["chi"] > 5 and cell["mechanism"] != "DECAY":
                        meta_cells.append(
                            (ang, r, d, cell["chi"], cell["mechanism"]))

        if not meta_cells:
            print(f"\n  N={n}: no metastable cells (chi>5 at boundaries)")
            continue

        chis_m = [c for _, _, _, c, _ in meta_cells]
        print(f"\n  N={n}: {len(meta_cells)} boundary cells with chi > 5")
        print(f"    chi range: [{min(chis_m):.2f}, {max(chis_m):.2f}]")
        print(f"    mean chi:  {statistics.mean(chis_m):.2f}")

        # Group by angle
        by_ang = defaultdict(list)
        for ang, r, d, chi, mech in meta_cells:
            by_ang[ang].append((r, d, chi, mech))

        for ang in sorted(by_ang):
            cells = by_ang[ang]
            chis_a = [c for _, _, c, _ in cells]
            mechs_a = set(m for _, _, _, m in cells)
            print(f"    angle={ang:+3d}: {len(cells)} cells, "
                  f"chi=[{min(chis_a):.1f}..{max(chis_a):.1f}], "
                  f"mechs={set(SHORT[m] for m in mechs_a)}")


# ============================================================
# Section 10: Drift smoothing of temporal transitions
# ============================================================

def check_drift_temporal_smoothing(cube, ns, radii, drifts):
    print("\n" + "=" * 90)
    print("  DRIFT TEMPORAL SMOOTHING: does drift smooth or sharpen chi transitions?")
    print("=" * 90)

    # For each N, compute the max angular chi gradient at each drift
    print(f"\n  Max angular chi-gradient vs drift:")
    print(f"  {'N':>3}  {'drift':>6}  {'max_grad':>9}  {'mean_grad':>10}  "
          f"{'chi_std':>8}  verdict")
    print(f"  " + "-" * 48)

    for n in ns:
        prev_max = None
        for d in drifts:
            grads = []
            chis_all = []

            for r in radii:
                prev_chi = None
                for ang in ANGLES:
                    cell = get(cube, n, ang, r, d)
                    if cell is None:
                        prev_chi = None
                        continue
                    chi = cell["chi"]
                    chis_all.append(chi)
                    if prev_chi is not None:
                        grads.append(abs(chi - prev_chi))
                    prev_chi = chi

            if not grads:
                continue

            max_g = max(grads)
            mean_g = statistics.mean(grads)
            chi_std = statistics.stdev(chis_all) if len(chis_all) > 1 else 0

            if prev_max is not None:
                if max_g < prev_max * 0.7:
                    verdict = "SMOOTHING"
                elif max_g > prev_max * 1.3:
                    verdict = "SHARPENING"
                else:
                    verdict = "STABLE"
            else:
                verdict = "---"

            print(f"  {n:3d}  {d:6.3f}  {max_g:9.2f}  {mean_g:10.2f}  "
                  f"{chi_std:8.2f}  {verdict}")
            prev_max = max_g
        print()


# ============================================================
# Section 11: Cross-N temporal summary
# ============================================================

def print_temporal_summary(cube, ns, angles, radii, drifts):
    print("=" * 90)
    print("  CROSS-N TEMPORAL SUMMARY")
    print("=" * 90)
    print(f"  {'N':>3}  {'drift':>6}  {'mean_chi':>9}  {'med_chi':>8}  "
          f"{'max_chi':>8}  {'DECAY_n':>7}  {'slow_n':>6}  "
          f"{'tangent_jump':>13}  {'max_grad':>9}")
    print(f"  " + "-" * 82)

    for n in ns:
        for d in [0.00, 0.02, 0.10]:
            if d not in drifts:
                continue

            all_chi = []
            decay_n = 0
            for ang in angles:
                for r in radii:
                    cell = get(cube, n, ang, r, d)
                    if cell:
                        all_chi.append(cell["chi"])
                        if cell["mechanism"] == "DECAY":
                            decay_n += 1

            if not all_chi:
                continue

            mean_c = statistics.mean(all_chi)
            med_c = statistics.median(all_chi)
            max_c = max(all_chi)

            # Slow cells (>2 sigma above non-DECAY mean)
            non_decay = [c for c in all_chi if c < DECAY_TIMEOUT - 1]
            if len(non_decay) > 2:
                thr = statistics.mean(non_decay) + 2 * statistics.stdev(non_decay)
                slow_n = sum(1 for c in non_decay if c > thr)
            else:
                slow_n = 0

            # Tangent jump
            chis_0 = [get(cube, n, 0, r, d)["chi"]
                       for r in radii if get(cube, n, 0, r, d)]
            chis_1 = [get(cube, n, 1, r, d)["chi"]
                       for r in radii if get(cube, n, 1, r, d)]
            if chis_0 and chis_1:
                t_jump = statistics.mean(chis_0) - statistics.mean(chis_1)
            else:
                t_jump = 0

            # Max gradient
            max_g = 0
            for r in radii:
                prev = None
                for ang in angles:
                    cell = get(cube, n, ang, r, d)
                    if cell:
                        if prev is not None:
                            max_g = max(max_g, abs(cell["chi"] - prev))
                        prev = cell["chi"]

            print(f"  {n:3d}  {d:6.3f}  {mean_c:9.2f}  {med_c:8.2f}  "
                  f"{max_c:8.2f}  {decay_n:7d}  {slow_n:6d}  "
                  f"{t_jump:+12.2f}  {max_g:9.2f}")
        print()


# ============================================================
# Main
# ============================================================

def main():
    print("=" * 90)
    print("  TEMPORAL MANIFOLD ANALYSIS: chi_emp across N x angle x radius x drift")
    print("=" * 90)
    print()

    cube, ns, angles, radii, drifts = build_cube()

    if not cube:
        print("  ERROR: No data loaded.")
        return

    # Section 1: Temporal heatmaps
    print_temporal_heatmaps(cube, ns, angles, radii)

    # Section 2: Mechanism timing statistics
    print_mechanism_timing(cube, ns, angles, radii, drifts)

    # Section 3: Tangent temporal jump
    check_tangent_temporal_jump(cube, ns, radii, drifts)

    # Section 4: DECAY temporal signature
    check_decay_temporal(cube, ns, angles, radii, drifts)

    # Section 5: Slow manifolds
    find_slow_manifolds(cube, ns, angles, radii, drifts)

    # Section 6: Temporal boundary sharpness
    check_temporal_boundaries(cube, ns, angles, radii, drifts)

    # Section 7: N=8 temporal footprint
    check_n8_temporal(cube, ns, angles, radii, drifts)

    # Section 8: Large-N temporal scaling
    check_large_n_scaling(cube, ns, angles, radii, drifts)

    # Section 9: Metastable regions
    find_metastable_regions(cube, ns, angles, radii, drifts)

    # Section 10: Drift temporal smoothing
    check_drift_temporal_smoothing(cube, ns, radii, drifts)

    # Section 11: Cross-N temporal summary
    print_temporal_summary(cube, ns, angles, radii, drifts)

    print("=" * 90)


if __name__ == "__main__":
    main()
