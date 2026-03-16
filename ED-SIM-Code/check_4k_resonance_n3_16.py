#!/usr/bin/env python
"""
check_4k_resonance_n3_16.py
============================
Test the 4k resonance hypothesis across N = 3-16.  Prints a compact
table of IC fractions, complementarity sums, and anomaly indicators,
then a focused 4k family comparison (N = 4, 8, 12, 16).
"""

import json
import os

SWEEP_DIR = os.path.dirname(__file__)
NS = [3, 4, 5, 6, 7, 8, 9, 10, 12, 16]


def load_sweep(gate, n):
    path = os.path.join(SWEEP_DIR, f"n{n}_{gate}_sweep.json")
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict) and "results" in data:
        data = data["results"]
    return data


def frac(records, mech):
    return sum(1 for r in records if r["mechanism"] == mech) / len(records)


def main():
    # ---- Full table ----
    print("=" * 92)
    print("  4k Resonance Test — N = 3..16")
    print("=" * 92)
    print()
    hdr = (f"  {'N':>3}  {'IC(inw)':>8}  {'IC(out)':>8}  {'Sum':>8}  "
           f"{'Crn(inw)':>9}  {'OL(inw)':>8}  {'Crn(out)':>9}  {'OL(out)':>8}  Flag")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    for n in NS:
        inw = load_sweep("inward", n)
        out = load_sweep("outward", n)

        ic_i  = frac(inw, "inward-collapse") * 100
        ic_o  = frac(out, "inward-collapse") * 100
        s     = ic_i + ic_o
        crn_i = frac(inw, "PBC-corner") * 100
        ol_i  = frac(inw, "other-late") * 100
        crn_o = frac(out, "PBC-corner") * 100
        ol_o  = frac(out, "other-late") * 100

        flags = []
        if n == 3:           flags.append("LawV")
        if abs(s - 100) > 10: flags.append("LawII-break")
        if crn_i > 1:        flags.append("LawI(inw)")
        if crn_o > 1:        flags.append("LawI(out)")
        if ol_i > 10:        flags.append("OL(inw)")
        if ol_o > 10:        flags.append("OL(out)")
        if n % 4 == 0:       flags.append("4k")
        flag_str = "  ".join(flags)

        print(f"  {n:3d}  {ic_i:7.1f}%  {ic_o:7.1f}%  {s:7.1f}%  "
              f"{crn_i:8.1f}%  {ol_i:7.1f}%  {crn_o:8.1f}%  {ol_o:7.1f}%  {flag_str}")

    # ---- 4k family focus ----
    print()
    print("-" * 92)
    print("  4k Family Comparison (N = 4, 8, 12, 16)")
    print("-" * 92)
    print(f"  {'N':>3}  {'IC(inw)':>8}  {'IC(out)':>8}  {'Sum':>8}  "
          f"{'Crn(inw)':>9}  {'OL(inw)':>8}  {'DECAY(tan)':>11}  Anomaly type")
    print("  " + "-" * 82)

    for n in [4, 8, 12, 16]:
        inw = load_sweep("inward", n)
        out = load_sweep("outward", n)
        tan = load_sweep("tangent", n)

        ic_i    = frac(inw, "inward-collapse") * 100
        ic_o    = frac(out, "inward-collapse") * 100
        s       = ic_i + ic_o
        crn_i   = frac(inw, "PBC-corner") * 100
        ol_i    = frac(inw, "other-late") * 100
        decay_t = frac(tan, "DECAY") * 100

        if decay_t > 50:          atype = "tangent DECAY"
        elif crn_i > 10:          atype = "all-gate catastrophic"
        elif ol_i > 10:           atype = "OL inflation"
        elif abs(s - 100) > 5:    atype = "mild complementarity drift"
        else:                     atype = "NORMAL"

        print(f"  {n:3d}  {ic_i:7.1f}%  {ic_o:7.1f}%  {s:7.1f}%  "
              f"{crn_i:8.1f}%  {ol_i:7.1f}%  {decay_t:10.1f}%  {atype}")

    # ---- Trend within 4k ----
    print()
    print("-" * 92)
    print("  4k Anomaly Severity Trend")
    print("-" * 92)
    for n in [4, 8, 12, 16]:
        inw = load_sweep("inward", n)
        out = load_sweep("outward", n)
        s = (frac(inw, "inward-collapse") + frac(out, "inward-collapse")) * 100
        dev = abs(s - 100)
        ol_i = frac(inw, "other-late") * 100
        crn_i = frac(inw, "PBC-corner") * 100
        severity = dev + ol_i + crn_i  # composite anomaly score
        bar = "#" * int(severity / 2)
        print(f"  N={n:2d}: |sum-100|={dev:5.1f}  OL(inw)={ol_i:5.1f}  "
              f"Crn(inw)={crn_i:5.1f}  score={severity:5.1f}  {bar}")

    print()
    print("=" * 92)


if __name__ == "__main__":
    main()
