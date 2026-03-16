#!/usr/bin/env python
"""
check_4k_resonance_n3_12.py
============================
Test the 4k resonance hypothesis (N = 4, 8, 12, ...) by comparing IC
fractions, complementarity sums, PBC-corner intrusion, and other-late
inflation across N = 3-12.

Prints a compact table and flags any N that matches the N=8 anomaly
signature: suppressed IC, elevated PBC-corner/other-late under radial
gates, and broken complementarity.
"""

import json
import os

SWEEP_DIR = os.path.dirname(__file__)
NS = [3, 4, 5, 6, 7, 8, 9, 10, 12]


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
    print("=" * 86)
    print("  4k Resonance Test — N = 3..12")
    print("=" * 86)
    print()
    print(f"  {'N':>3}  {'IC(inw)':>8}  {'IC(out)':>8}  {'Sum':>8}  "
          f"{'Crn(inw)':>9}  {'OL(inw)':>8}  {'Crn(out)':>9}  {'OL(out)':>8}  Flag")
    print(f"  {'---':>3}  {'-------':>8}  {'-------':>8}  {'---':>8}  "
          f"{'---------':>9}  {'--------':>8}  {'---------':>9}  {'--------':>8}  ----")

    for n in NS:
        inw = load_sweep("inward", n)
        out = load_sweep("outward", n)

        ic_i = frac(inw, "inward-collapse") * 100
        ic_o = frac(out, "inward-collapse") * 100
        s    = ic_i + ic_o

        crn_i = frac(inw, "PBC-corner") * 100
        ol_i  = frac(inw, "other-late") * 100
        crn_o = frac(out, "PBC-corner") * 100
        ol_o  = frac(out, "other-late") * 100

        # Flag logic
        flags = []
        if n == 3:
            flags.append("Law V")
        if abs(s - 100) > 10:
            flags.append("LawII-break")
        if crn_i > 1:
            flags.append("LawI-break(inw)")
        if crn_o > 1:
            flags.append("LawI-break(out)")
        if ol_i > 10:
            flags.append("OL-high(inw)")
        if ol_o > 10:
            flags.append("OL-high(out)")
        if n % 4 == 0:
            flags.append("4k")

        flag_str = "  ".join(flags) if flags else ""

        print(f"  {n:3d}  {ic_i:7.1f}%  {ic_o:7.1f}%  {s:7.1f}%  "
              f"{crn_i:8.1f}%  {ol_i:7.1f}%  {crn_o:8.1f}%  {ol_o:7.1f}%  {flag_str}")

    # --- 4k summary ---
    print()
    print("-" * 86)
    print("  4k Family Comparison")
    print("-" * 86)
    for n in [4, 8, 12]:
        inw = load_sweep("inward", n)
        out = load_sweep("outward", n)
        tan = load_sweep("tangent", n)

        ic_i = frac(inw, "inward-collapse") * 100
        ic_o = frac(out, "inward-collapse") * 100
        s    = ic_i + ic_o
        crn_i = frac(inw, "PBC-corner") * 100
        ol_i  = frac(inw, "other-late") * 100
        decay_t = frac(tan, "DECAY") * 100

        print(f"  N={n:2d}:  IC(inw)={ic_i:5.1f}%  IC(out)={ic_o:5.1f}%  "
              f"sum={s:5.1f}%  corner(inw)={crn_i:4.1f}%  "
              f"OL(inw)={ol_i:4.1f}%  DECAY(tan)={decay_t:4.1f}%")

    print()
    print("  Verdict: does N=12 match the N=8 anomaly signature?")
    inw12 = load_sweep("inward", 12)
    s12 = (frac(inw12, "inward-collapse") +
           frac(load_sweep("outward", 12), "inward-collapse")) * 100
    crn12 = frac(inw12, "PBC-corner") * 100
    ol12  = frac(inw12, "other-late") * 100

    anomaly = abs(s12 - 100) > 10 or crn12 > 1 or ol12 > 10
    if anomaly:
        print("  --> YES: N=12 shows anomalous behavior (4k resonance supported)")
    else:
        print("  --> NO: N=12 behaves normally (4k resonance NOT confirmed at k=3)")
    print("=" * 86)


if __name__ == "__main__":
    main()
