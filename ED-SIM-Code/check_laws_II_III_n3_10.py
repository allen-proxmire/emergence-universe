#!/usr/bin/env python
"""
check_laws_II_III_n3_10.py
==========================
Maintenance script: test Law II (Radial Complementarity) and Law III
(Monotone Gate-Ordering) across the extended range N = 3-10.

Loads all inward/outward sweep files, computes IC fractions, prints a
compact table, and flags any violations.
"""

import json
import os

SWEEP_DIR = os.path.dirname(__file__)
NS = list(range(3, 11))


def load_sweep(gate, n):
    path = os.path.join(SWEEP_DIR, f"n{n}_{gate}_sweep.json")
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict) and "results" in data:
        data = data["results"]
    return data


def ic_frac(records):
    return sum(1 for r in records if r["mechanism"] == "inward-collapse") / len(records)


def main():
    print("=" * 66)
    print("  Law II / III check — N = 3..10")
    print("=" * 66)
    print()
    print(f"  {'N':>3}  {'IC(inw)':>8}  {'IC(out)':>8}  {'Sum':>8}  {'|Sum-100|':>9}")
    print(f"  {'---':>3}  {'-------':>8}  {'-------':>8}  {'---':>8}  {'---------':>9}")

    inw_seq = []
    out_seq = []

    for n in NS:
        fi = ic_frac(load_sweep("inward", n)) * 100
        fo = ic_frac(load_sweep("outward", n)) * 100
        s = fi + fo
        dev = abs(s - 100)
        inw_seq.append(fi)
        out_seq.append(fo)
        flag = "  <-- OUTLIER" if n == 3 else ("  *** DRIFT" if dev > 5 else "")
        print(f"  {n:3d}  {fi:7.1f}%  {fo:7.1f}%  {s:7.1f}%  {dev:8.1f}%{flag}")

    # --- Law III: monotonicity ---
    print()
    print("-" * 66)
    print("  Law III: Monotonicity Check")
    print("-" * 66)

    inw_diffs = [inw_seq[i+1] - inw_seq[i] for i in range(len(NS)-1)]
    out_diffs = [out_seq[i+1] - out_seq[i] for i in range(len(NS)-1)]

    print()
    print("  Inward IC step diffs (should all be < 0):")
    for i, d in enumerate(inw_diffs):
        tag = "OK" if d < 0 else "VIOLATION"
        print(f"    N={NS[i]}->{NS[i+1]}: {d:+6.1f} pp  [{tag}]")

    inw_ok = all(d < 0 for d in inw_diffs)
    print(f"  Inward strictly decreasing: {inw_ok}")

    print()
    print("  Outward IC step diffs (should all be > 0):")
    for i, d in enumerate(out_diffs):
        tag = "OK" if d > 0 else "VIOLATION"
        print(f"    N={NS[i]}->{NS[i+1]}: {d:+6.1f} pp  [{tag}]")

    out_ok = all(d > 0 for d in out_diffs)
    print(f"  Outward strictly increasing: {out_ok}")

    # --- Law II: complementarity for N >= 4 ---
    print()
    print("-" * 66)
    print("  Law II: Complementarity Summary (N >= 4)")
    print("-" * 66)
    sums_ge4 = [inw_seq[i] + out_seq[i] for i in range(len(NS)) if NS[i] >= 4]
    import statistics
    mean_s = statistics.mean(sums_ge4)
    max_dev = max(abs(s - 100) for s in sums_ge4)
    print(f"  Mean sum: {mean_s:.1f}%")
    print(f"  Max |deviation|: {max_dev:.1f}%")
    law2_ok = max_dev < 10
    print(f"  Law II holds (< 10% deviation): {law2_ok}")

    # --- Overall ---
    print()
    print("=" * 66)
    if law2_ok and inw_ok and out_ok:
        print("  ALL LAWS HOLD through N = 10.")
    else:
        issues = []
        if not law2_ok:
            issues.append("Law II violated (complementarity drift > 10%)")
        if not inw_ok:
            issues.append("Law III violated (inward IC not monotone decreasing)")
        if not out_ok:
            issues.append("Law III violated (outward IC not monotone increasing)")
        print("  VIOLATIONS FOUND:")
        for iss in issues:
            print(f"    - {iss}")
    print("=" * 66)


if __name__ == "__main__":
    main()
