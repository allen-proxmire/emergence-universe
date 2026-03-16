#!/usr/bin/env python
"""
check_large_n_drift_n3_32.py
==============================
Large-N drift test: do 4k anomalies weaken, stabilize, or re-amplify
as N increases through 20, 24, 28, 32?
"""

import json
import os

SWEEP_DIR = os.path.dirname(__file__)
NS = [3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24, 28, 32]


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
    print("=" * 100)
    print("  Large-N Drift Test — N = 3..32")
    print("=" * 100)
    print()
    print(f"  {'N':>3}  {'IC(inw)':>8}  {'IC(out)':>8}  {'Sum':>8}  "
          f"{'Crn(inw)':>9}  {'OL(inw)':>8}  {'Crn(out)':>9}  {'OL(out)':>8}  Flag")
    print("  " + "-" * 94)

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
        if n == 3:            flags.append("LawV")
        if abs(s - 100) > 10: flags.append("LawII-brk")
        if crn_i > 1:         flags.append("Crn(i)")
        if crn_o > 1:         flags.append("Crn(o)")
        if ol_i > 10:         flags.append("OL(i)")
        if ol_o > 10:         flags.append("OL(o)")
        if n % 4 == 0:        flags.append("4k")
        flag_str = "  ".join(flags)

        print(f"  {n:3d}  {ic_i:7.1f}%  {ic_o:7.1f}%  {s:7.1f}%  "
              f"{crn_i:8.1f}%  {ol_i:7.1f}%  {crn_o:8.1f}%  {ol_o:7.1f}%  {flag_str}")

    # ---- 4k family focus ----
    print()
    print("-" * 100)
    print("  4k Family Trajectory (k = 1..8)")
    print("-" * 100)
    print(f"  {'N':>3}  {'k':>2}  {'IC(inw)':>8}  {'IC(out)':>8}  {'Sum':>8}  "
          f"{'Crn(inw)':>9}  {'OL(inw)':>8}  {'DECAY(t)':>9}  {'Score':>6}  Anomaly type")
    print("  " + "-" * 92)

    for n in [4, 8, 12, 16, 20, 24, 28, 32]:
        k = n // 4
        inw = load_sweep("inward", n)
        out = load_sweep("outward", n)
        tan = load_sweep("tangent", n)

        ic_i    = frac(inw, "inward-collapse") * 100
        ic_o    = frac(out, "inward-collapse") * 100
        s       = ic_i + ic_o
        crn_i   = frac(inw, "PBC-corner") * 100
        ol_i    = frac(inw, "other-late") * 100
        decay_t = frac(tan, "DECAY") * 100
        dev     = abs(s - 100)
        score   = dev + ol_i + crn_i

        if decay_t > 50:       atype = "tangent DECAY"
        elif crn_i > 10:       atype = "catastrophic"
        elif ol_i > 10:
            if crn_i > 1:      atype = "OL + corner"
            else:              atype = "OL inflation"
        elif dev > 5:          atype = "complementarity drift"
        else:                  atype = "NORMAL"

        bar = "#" * int(score / 2)
        print(f"  {n:3d}  {k:2d}  {ic_i:7.1f}%  {ic_o:7.1f}%  {s:7.1f}%  "
              f"{crn_i:8.1f}%  {ol_i:7.1f}%  {decay_t:8.1f}%  {score:5.1f}  "
              f"{atype}  {bar}")

    # ---- Asymptotic convergence check ----
    print()
    print("-" * 100)
    print("  Asymptotic Convergence (N >= 16)")
    print("-" * 100)
    large_ns = [16, 20, 24, 28, 32]
    sums, ols_i, crns_i = [], [], []
    for n in large_ns:
        inw = load_sweep("inward", n)
        out = load_sweep("outward", n)
        s = (frac(inw, "inward-collapse") + frac(out, "inward-collapse")) * 100
        ol_i = frac(inw, "other-late") * 100
        crn_i = frac(inw, "PBC-corner") * 100
        sums.append(s)
        ols_i.append(ol_i)
        crns_i.append(crn_i)

    import statistics
    print(f"  Complementarity sum:  mean={statistics.mean(sums):.1f}%  "
          f"std={statistics.stdev(sums):.1f}%  range=[{min(sums):.1f}, {max(sums):.1f}]")
    print(f"  Other-late(inw):      mean={statistics.mean(ols_i):.1f}%  "
          f"std={statistics.stdev(ols_i):.1f}%  range=[{min(ols_i):.1f}, {max(ols_i):.1f}]")
    print(f"  PBC-corner(inw):      mean={statistics.mean(crns_i):.1f}%  "
          f"std={statistics.stdev(crns_i):.1f}%  range=[{min(crns_i):.1f}, {max(crns_i):.1f}]")

    # IC convergence
    print()
    print("  IC(inw) vs IC(out) convergence:")
    for n in large_ns:
        inw = load_sweep("inward", n)
        out = load_sweep("outward", n)
        ic_i = frac(inw, "inward-collapse") * 100
        ic_o = frac(out, "inward-collapse") * 100
        gap = ic_i - ic_o
        print(f"    N={n:2d}: IC(inw)={ic_i:5.1f}  IC(out)={ic_o:5.1f}  gap={gap:+5.1f}")

    print()
    print("=" * 100)


if __name__ == "__main__":
    main()
