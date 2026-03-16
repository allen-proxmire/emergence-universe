#!/usr/bin/env python
"""
check_tangent_manifold.py
==========================
Map the full angular structure of the tangent manifold from -10 to +10
degrees around the tangent gate (0 deg) for N = 4, 8, 12, 20.
"""

import json
import os

SWEEP_DIR = os.path.dirname(__file__)
NS = [4, 8, 12, 20]
ANGLES = list(range(-10, 11))
MECHS = ["inward-collapse", "outward-PBC", "DECAY", "PBC-corner", "other-late"]
SHORT = {"inward-collapse": "IC", "outward-PBC": "OPBC", "DECAY": "DEC",
         "PBC-corner": "CRN", "other-late": "OL"}


def load(n, ang):
    label = f"p{ang}deg" if ang > 0 else (f"m{abs(ang)}deg" if ang < 0 else "0deg")
    path = os.path.join(SWEEP_DIR, f"n{n}_angle_{label}_sweep.json")
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict) and "results" in data:
        data = data["results"]
    return data


def fracs(records):
    total = len(records)
    c = {}
    for r in records:
        m = r["mechanism"]
        c[m] = c.get(m, 0) + 1
    return {m: 100 * c.get(m, 0) / total for m in MECHS}


def main():
    # ================================================================
    # Section 1: Full mechanism table per N
    # ================================================================
    all_data = {}  # all_data[n][ang] = {mech: frac}

    for n in NS:
        print("=" * 78)
        print(f"  N = {n}   Tangent Manifold: angle = -10 .. +10 degrees")
        print("=" * 78)

        hdr = f"  {'ang':>4}"
        for m in MECHS:
            hdr += f"  {SHORT[m]:>5}"
        print(hdr)
        print("  " + "-" * 38)

        all_data[n] = {}
        for ang in ANGLES:
            f = fracs(load(n, ang))
            all_data[n][ang] = f
            row = f"  {ang:+4d}"
            for m in MECHS:
                row += f"  {f[m]:5.1f}"
            if ang == 0:
                row += "  <-- tangent"
            print(row)
        print()

    # ================================================================
    # Section 2: IC profile with ASCII spark line
    # ================================================================
    print("=" * 78)
    print("  IC(angle) PROFILES")
    print("=" * 78)
    for n in NS:
        print(f"\n  N={n}:")
        for ang in ANGLES:
            ic = all_data[n][ang]["inward-collapse"]
            bar = "#" * int(ic / 2)
            marker = " <--" if ang == 0 else ""
            print(f"    {ang:+3d}deg  {ic:5.1f}%  {bar}{marker}")

    # ================================================================
    # Section 3: Jump magnitude at 0 -> +1 deg
    # ================================================================
    print("\n" + "=" * 78)
    print("  CRITICAL DISCONTINUITY: IC jump from 0 deg to +1 deg")
    print("=" * 78)
    for n in NS:
        ic0 = all_data[n][0]["inward-collapse"]
        ic1 = all_data[n][1]["inward-collapse"]
        jump = ic1 - ic0
        print(f"  N={n:2d}:  IC(0°)={ic0:5.1f}%  IC(+1°)={ic1:5.1f}%  "
              f"jump={jump:+6.1f} pp")

    # ================================================================
    # Section 4: Derivative / slope around the tangent
    # ================================================================
    print("\n" + "=" * 78)
    print("  SLOPE ANALYSIS: dIC/dAngle (pp per degree)")
    print("=" * 78)
    print(f"  {'N':>3}  {'neg slope':>10}  {'pos slope':>10}  {'asymmetry':>10}  notes")
    print(f"       {'(-10..-1)':>10}  {'(+1..+10)':>10}  {'(pos/neg)':>10}")
    print("  " + "-" * 54)

    for n in NS:
        # Negative side: linear fit -10 to -1
        neg_ics = [all_data[n][a]["inward-collapse"] for a in range(-10, 0)]
        neg_angs = list(range(-10, 0))
        # Simple rise/run
        neg_slope = (neg_ics[-1] - neg_ics[0]) / (neg_angs[-1] - neg_angs[0])

        # Positive side: linear fit +1 to +10
        pos_ics = [all_data[n][a]["inward-collapse"] for a in range(1, 11)]
        pos_angs = list(range(1, 11))
        pos_slope = (pos_ics[-1] - pos_ics[0]) / (pos_angs[-1] - pos_angs[0])

        if abs(neg_slope) > 0.01:
            asym = abs(pos_slope / neg_slope)
        else:
            asym = float("inf")

        note = ""
        if all_data[n][0]["DECAY"] > 50:
            note = "DECAY at 0°"
        if asym > 10:
            note += " extreme asymmetry"

        print(f"  {n:3d}  {neg_slope:+9.2f}  {pos_slope:+9.2f}  {asym:9.1f}x  {note}")

    # ================================================================
    # Section 5: DECAY profile for N=4
    # ================================================================
    print("\n" + "=" * 78)
    print("  N=4 DECAY PROFILE (measure-zero test)")
    print("=" * 78)
    print(f"  {'ang':>4}  {'DECAY':>6}  {'IC':>5}  {'OPBC':>5}  {'CRN':>5}  bar")
    print("  " + "-" * 42)
    for ang in ANGLES:
        f = all_data[4][ang]
        dec = f["DECAY"]
        ic = f["inward-collapse"]
        opbc = f["outward-PBC"]
        crn = f["PBC-corner"]
        bar_d = "D" * int(dec / 2)
        bar_i = "I" * int(ic / 2)
        print(f"  {ang:+4d}  {dec:5.1f}%  {ic:5.1f}  {opbc:5.1f}  {crn:5.1f}  {bar_d}{bar_i}")

    # ================================================================
    # Section 6: Mechanism at 0 -> what replaces it on each side
    # ================================================================
    print("\n" + "=" * 78)
    print("  MECHANISM REPLACEMENT: what replaces the 0-deg dominant mechanism?")
    print("=" * 78)
    for n in NS:
        f0 = all_data[n][0]
        dom_mech = max(f0, key=f0.get)
        dom_frac = f0[dom_mech]
        fm1 = all_data[n][-1]
        fp1 = all_data[n][1]
        # What grew at -1 vs 0?
        neg_gains = {m: fm1[m] - f0[m] for m in MECHS if fm1[m] - f0[m] > 1}
        pos_gains = {m: fp1[m] - f0[m] for m in MECHS if fp1[m] - f0[m] > 1}
        print(f"\n  N={n}: dominant at 0° = {SHORT[dom_mech]} ({dom_frac:.1f}%)")
        if neg_gains:
            parts = ", ".join(f"{SHORT[m]}={v:+.1f}" for m, v in neg_gains.items())
            print(f"    at -1°: gains = {parts}")
        else:
            print(f"    at -1°: no significant change")
        if pos_gains:
            parts = ", ".join(f"{SHORT[m]}={v:+.1f}" for m, v in pos_gains.items())
            print(f"    at +1°: gains = {parts}")
        else:
            print(f"    at +1°: no significant change")

    # ================================================================
    # Section 7: Plateau detection on positive side
    # ================================================================
    print("\n" + "=" * 78)
    print("  POSITIVE-SIDE PLATEAU: does IC saturate?")
    print("=" * 78)
    for n in NS:
        ics_pos = [all_data[n][a]["inward-collapse"] for a in range(1, 11)]
        ic_at_1 = ics_pos[0]
        ic_at_10 = ics_pos[-1]
        ic_max = max(ics_pos)
        ic_min = min(ics_pos)
        span = ic_max - ic_min
        drift = ic_at_10 - ic_at_1
        if span < 10:
            shape = "PLATEAU"
        elif drift < -5:
            shape = "DECLINING"
        elif drift > 5:
            shape = "RISING"
        else:
            shape = "MIXED"
        print(f"  N={n:2d}: IC(+1)={ic_at_1:5.1f}  IC(+10)={ic_at_10:5.1f}  "
              f"range=[{ic_min:.1f},{ic_max:.1f}]  span={span:.1f}pp  -> {shape}")

    print("\n" + "=" * 78)


if __name__ == "__main__":
    main()
