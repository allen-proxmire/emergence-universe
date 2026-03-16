#!/usr/bin/env python
"""
Reproduce_This_Law — Law III: Monotone Gate-Ordering Law
========================================================
Verifies: Inward IC strictly decreases with N; outward IC strictly
increases; inward always exceeds outward.
"""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SWEEP_DIR = os.path.join(os.path.dirname(__file__), "..", "..")
EXPECTED_DIR = os.path.join(os.path.dirname(__file__), "expected_outputs")
PLOT_DIR = os.path.join(os.path.dirname(__file__), "plots")

NS = [3, 4, 5, 6, 7]


def load_sweep(gate, n):
    path = os.path.join(SWEEP_DIR, f"n{n}_{gate}_sweep.json")
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict) and "results" in data:
        data = data["results"]
    return data


def ic_frac(records):
    return sum(1 for r in records if r["mechanism"] == "inward-collapse") / len(records) * 100


def main():
    os.makedirs(PLOT_DIR, exist_ok=True)
    all_pass = True

    print("=" * 72)
    print("LAW III VERIFICATION: Monotone Gate-Ordering Law")
    print("=" * 72)

    inw_seq = []
    out_seq = []
    for n in NS:
        fi = ic_frac(load_sweep("inward", n))
        fo = ic_frac(load_sweep("outward", n))
        inw_seq.append(fi)
        out_seq.append(fo)
        spread = fi - fo
        print(f"  N={n}: IC(inward)={fi:5.1f}%  IC(outward)={fo:5.1f}%  spread={spread:5.1f} pp")

    # Check monotonicity
    inw_diffs = [inw_seq[i + 1] - inw_seq[i] for i in range(len(NS) - 1)]
    out_diffs = [out_seq[i + 1] - out_seq[i] for i in range(len(NS) - 1)]

    print("\n  Inward step diffs:", [f"{d:+.1f}" for d in inw_diffs])
    inw_mono = all(d < 0 for d in inw_diffs)
    print(f"  Inward strictly decreasing: {inw_mono}  "
          f"[{'PASS' if inw_mono else 'FAIL'}]")
    if not inw_mono:
        all_pass = False

    print(f"  Outward step diffs: {[f'{d:+.1f}' for d in out_diffs]}")
    out_mono = all(d > 0 for d in out_diffs)
    print(f"  Outward strictly increasing: {out_mono}  "
          f"[{'PASS' if out_mono else 'FAIL'}]")
    if not out_mono:
        all_pass = False

    # Check dominance
    dom = all(inw_seq[i] > out_seq[i] for i in range(len(NS)))
    print(f"\n  Inward dominates outward at all N: {dom}  "
          f"[{'PASS' if dom else 'FAIL'}]")
    if not dom:
        all_pass = False

    # Convergence
    spreads = [inw_seq[i] - out_seq[i] for i in range(len(NS))]
    print(f"\n  Spreads: {[f'{s:.1f}' for s in spreads]}")
    spread_diffs = [spreads[i + 1] - spreads[i] for i in range(len(NS) - 1)]
    narrowing = all(d < 0 for d in spread_diffs)
    print(f"  Spread narrowing: {narrowing}  (informational, not pass/fail)")

    # Compare with expected
    expected_path = os.path.join(EXPECTED_DIR, "metrics.json")
    if os.path.exists(expected_path):
        with open(expected_path) as f:
            expected = json.load(f)
        print("\n" + "-" * 72)
        print("EXPECTED-VALUE COMPARISON")
        print("-" * 72)
        for i, n in enumerate(NS):
            key = f"N{n}"
            if key in expected:
                exp_i = expected[key]["ic_inward"]
                exp_o = expected[key]["ic_outward"]
                act_i = round(inw_seq[i], 1)
                act_o = round(out_seq[i], 1)
                match = abs(exp_i - act_i) < 0.5 and abs(exp_o - act_o) < 0.5
                print(f"  {key}: inw exp={exp_i} act={act_i}  "
                      f"out exp={exp_o} act={act_o}  "
                      f"[{'OK' if match else 'DRIFT'}]")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.array(NS)
    ax.plot(x, inw_seq, "b^-", markersize=10, linewidth=2, label="IC(inward)")
    ax.plot(x, out_seq, "rs-", markersize=10, linewidth=2, label="IC(outward)")
    ax.fill_between(x, out_seq, inw_seq, alpha=0.15, color="gray", label="Spread")
    ax.set_xlabel("N (particle count)")
    ax.set_ylabel("IC fraction (%)")
    ax.set_title("Law III: Monotone Gate-Ordering")
    ax.legend()
    ax.set_xticks(NS)
    ax.set_ylim(0, 105)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plot_path = os.path.join(PLOT_DIR, "law_III_monotone_ordering.png")
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f"\n  Plot saved: {plot_path}")

    print("\n" + "=" * 72)
    if all_pass:
        print("LAW III VERIFICATION: PASS")
        print("  Monotone ordering confirmed: inward decreasing, outward increasing,")
        print("  inward dominates at all N.")
    else:
        print("LAW III VERIFICATION: FAIL")
    print("=" * 72)
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
