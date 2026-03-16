#!/usr/bin/env python
"""
check_2d_tangent_manifold.py
==============================
2D phase-diagram analysis of the tangent manifold in (angle, d_px) space
for N = 4, 8, 12, 20.
"""

import json
import os

SWEEP_DIR = os.path.dirname(__file__)
NS = [4, 8, 12, 20]

# Mechanism → single-char code for ASCII phase maps
CODE = {
    "inward-collapse": "I",
    "outward-PBC":     "O",
    "DECAY":           "D",
    "PBC-corner":      "C",
    "other-late":      ".",
}


def load_grid(n):
    path = os.path.join(SWEEP_DIR, f"n{n}_2d_angle_radius_grid.json")
    with open(path) as f:
        data = json.load(f)
    angles = data["angles"]
    d_px_list = data["d_px_list"]
    # Build lookup: (angle, d_px) -> mechanism
    lookup = {}
    for rec in data["grid"]:
        lookup[(rec["angle"], rec["d_px"])] = rec["mechanism"]
    return angles, d_px_list, lookup


def main():
    for n in NS:
        angles, d_px_list, lookup = load_grid(n)

        # ============================================================
        # Section 1: ASCII phase map
        # ============================================================
        print("=" * 80)
        print(f"  N = {n}   PHASE MAP (angle x d_px)")
        print(f"  Codes: I=IC  O=OPBC  D=DECAY  C=PBC-corner  .=other-late")
        print("=" * 80)

        # Header: angles
        hdr = "  d_px "
        for ang in angles:
            hdr += f" {ang:+3d}"
        print(hdr)
        print("  " + "-" * (6 + 4 * len(angles)))

        for d_px in d_px_list:
            row = f"  {d_px:4d}  "
            for ang in angles:
                mech = lookup.get((ang, d_px), "?")
                row += f"  {CODE.get(mech, '?')} "
            print(row)
        print()

        # ============================================================
        # Section 2: Mechanism fraction per angle column
        # ============================================================
        print(f"  --- Mechanism fractions by angle column ---")
        print(f"  {'ang':>4}  {'IC':>5}  {'OPBC':>5}  {'DEC':>5}  {'CRN':>5}  {'OL':>5}")
        print(f"  " + "-" * 36)
        for ang in angles:
            counts = {}
            total = 0
            for d_px in d_px_list:
                m = lookup.get((ang, d_px), "other-late")
                counts[m] = counts.get(m, 0) + 1
                total += 1
            row = f"  {ang:+4d}"
            for mech in ["inward-collapse", "outward-PBC", "DECAY",
                         "PBC-corner", "other-late"]:
                row += f"  {100 * counts.get(mech, 0) / total:5.1f}"
            print(row)
        print()

        # ============================================================
        # Section 3: Phase boundary detection
        # ============================================================
        # For each d_px row, find the angle where IC first exceeds 0
        # (the collapse onset) going from negative to positive
        print(f"  --- Phase boundary: IC onset angle per d_px ---")
        print(f"  {'d_px':>5}  {'IC onset':>9}  {'mech at 0':>10}  {'mech at +1':>11}")
        print(f"  " + "-" * 38)
        onsets = []
        for d_px in d_px_list:
            # Find first angle where IC appears
            onset = None
            for ang in angles:
                m = lookup.get((ang, d_px), "")
                if m == "inward-collapse":
                    onset = ang
                    break
            m_at_0 = CODE.get(lookup.get((0, d_px), ""), "?")
            m_at_p1 = CODE.get(lookup.get((1, d_px), ""), "?")
            onset_str = f"{onset:+3d}deg" if onset is not None else "  none"
            print(f"  {d_px:5d}  {onset_str:>9}  {m_at_0:>10}  {m_at_p1:>11}")
            onsets.append(onset)

        # Boundary statistics
        valid = [o for o in onsets if o is not None]
        if valid:
            from collections import Counter
            c = Counter(valid)
            print(f"\n  Boundary angle distribution:")
            for ang in sorted(c.keys()):
                pct = 100 * c[ang] / len(d_px_list)
                bar = "#" * int(pct / 2)
                print(f"    {ang:+3d}deg: {c[ang]:2d}/{len(d_px_list)} "
                      f"({pct:4.1f}%)  {bar}")
            none_count = len(d_px_list) - len(valid)
            if none_count:
                pct = 100 * none_count / len(d_px_list)
                print(f"    none:  {none_count:2d}/{len(d_px_list)} "
                      f"({pct:4.1f}%)  (no IC at any angle)")
        print()

        # ============================================================
        # Section 4: DECAY extent for N=4
        # ============================================================
        if n == 4:
            print(f"  --- N=4 DECAY extent in (angle, d_px) ---")
            decay_cells = [(ang, d_px) for (ang, d_px), m in lookup.items()
                           if m == "DECAY"]
            if decay_cells:
                d_angs = sorted(set(a for a, _ in decay_cells))
                d_rads = sorted(set(d for _, d in decay_cells))
                print(f"  DECAY present at angles: {d_angs}")
                print(f"  DECAY d_px range: [{min(d_rads)}, {max(d_rads)}]"
                      f"  ({len(d_rads)} of {len(d_px_list)} radii)")
                print(f"  Total DECAY cells: {len(decay_cells)} "
                      f"of {len(lookup)} ({100*len(decay_cells)/len(lookup):.1f}%)")
                # DECAY per d_px at angle=0
                print(f"\n  DECAY at angle=0 by d_px:")
                for d_px in d_px_list:
                    m = lookup.get((0, d_px), "")
                    if m == "DECAY":
                        print(f"    d_px={d_px:3d}: DECAY")
            else:
                print(f"  No DECAY cells found in the grid.")
            print()

        # ============================================================
        # Section 5: Boundary curvature assessment
        # ============================================================
        print(f"  --- Boundary geometry ---")
        # Check if boundary is vertical (same onset angle for all d_px)
        if valid:
            unique_onsets = set(valid)
            if len(unique_onsets) == 1:
                print(f"  VERTICAL boundary at angle = {valid[0]:+d}deg")
                print(f"  (identical onset across all {len(valid)} radii)")
            else:
                # Check for tilt: correlation between d_px and onset angle
                paired = [(d, o) for d, o in zip(d_px_list, onsets)
                          if o is not None]
                d_vals = [p[0] for p in paired]
                o_vals = [p[1] for p in paired]
                mean_d = sum(d_vals) / len(d_vals)
                mean_o = sum(o_vals) / len(o_vals)
                cov = sum((d - mean_d) * (o - mean_o)
                          for d, o in zip(d_vals, o_vals)) / len(d_vals)
                var_d = sum((d - mean_d) ** 2 for d in d_vals) / len(d_vals)
                slope = cov / var_d if var_d > 0 else 0
                print(f"  PIECEWISE boundary")
                print(f"  Onset angles: {sorted(unique_onsets)}")
                print(f"  Mean onset: {mean_o:+.1f}deg")
                print(f"  Tilt (d onset / d d_px): {slope:+.4f} deg/px")
                if abs(slope) < 0.001:
                    print(f"  -> NEAR-VERTICAL (negligible tilt)")
                elif abs(slope) < 0.01:
                    print(f"  -> SLIGHT TILT")
                else:
                    print(f"  -> SIGNIFICANT CURVATURE")
        print()
        print()

    # ================================================================
    # Section 6: Cross-N comparison
    # ================================================================
    print("=" * 80)
    print("  CROSS-N BOUNDARY COMPARISON")
    print("=" * 80)
    print(f"  {'N':>3}  {'boundary':>10}  {'tilt':>12}  {'DECAY cells':>12}  {'IC at +1':>9}")
    print("  " + "-" * 52)

    for n in NS:
        angles, d_px_list, lookup = load_grid(n)
        onsets = []
        for d_px in d_px_list:
            onset = None
            for ang in angles:
                if lookup.get((ang, d_px), "") == "inward-collapse":
                    onset = ang
                    break
            onsets.append(onset)
        valid = [o for o in onsets if o is not None]
        unique = sorted(set(valid)) if valid else []
        decay_n = sum(1 for m in lookup.values() if m == "DECAY")
        # IC fraction at +1 deg
        ic_at_1 = sum(1 for d in d_px_list
                      if lookup.get((1, d), "") == "inward-collapse")
        ic_pct = 100 * ic_at_1 / len(d_px_list)

        bnd = f"{unique}" if len(unique) <= 3 else f"[{min(unique)}..{max(unique)}]"
        print(f"  {n:3d}  {bnd:>10}  "
              f"{'vertical' if len(unique) <= 1 else 'piecewise':>12}  "
              f"{decay_n:>12}  {ic_pct:8.1f}%")

    print()
    print("=" * 80)


if __name__ == "__main__":
    main()
