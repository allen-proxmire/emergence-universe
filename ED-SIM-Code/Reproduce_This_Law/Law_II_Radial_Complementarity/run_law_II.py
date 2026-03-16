#!/usr/bin/env python
"""
Reproduce_This_Law — Law II: Radial Complementarity Law
========================================================
Verifies: IC(inward) + IC(outward) ≈ 100% for N >= 4.
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
SUM_TOL = 5.0          # max deviation from 100% for N >= 4
MEAN_LOW = 98.0
MEAN_HIGH = 106.0
N3_OUTLIER_MIN = 115.0  # N=3 sum must exceed this (Law V interaction)


def load_sweep(gate, n):
    path = os.path.join(SWEEP_DIR, f"n{n}_{gate}_sweep.json")
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict) and "results" in data:
        data = data["results"]
    return data


def ic_fraction(records):
    return sum(1 for r in records if r["mechanism"] == "inward-collapse") / len(records)


def main():
    os.makedirs(PLOT_DIR, exist_ok=True)
    all_pass = True
    results = {}

    print("=" * 72)
    print("LAW II VERIFICATION: Radial Complementarity Law")
    print("=" * 72)

    sums = []
    ic_inw = []
    ic_out = []

    for n in NS:
        inw = load_sweep("inward", n)
        out = load_sweep("outward", n)
        fi = ic_fraction(inw) * 100
        fo = ic_fraction(out) * 100
        s = fi + fo
        ic_inw.append(fi)
        ic_out.append(fo)
        sums.append(s)

        results[f"N{n}"] = {
            "ic_inward": round(fi, 2),
            "ic_outward": round(fo, 2),
            "comp_sum": round(s, 2),
        }

        dev = abs(s - 100)
        if n >= 4:
            ok = dev < SUM_TOL
            tag = "PASS" if ok else "FAIL"
            if not ok:
                all_pass = False
        else:
            ok = s > N3_OUTLIER_MIN
            tag = "OUTLIER-OK" if ok else "OUTLIER-FAIL"
            if not ok:
                all_pass = False

        print(f"  N={n}: IC(inw)={fi:5.1f}%  IC(out)={fo:5.1f}%  "
              f"sum={s:6.1f}%  dev={dev:4.1f}%  [{tag}]")

    # Mean for N >= 4
    sums_ge4 = [sums[i] for i, n in enumerate(NS) if n >= 4]
    mean_s = np.mean(sums_ge4)
    mean_ok = MEAN_LOW <= mean_s <= MEAN_HIGH
    if not mean_ok:
        all_pass = False
    print(f"\n  Mean sum (N>=4): {mean_s:.2f}%  "
          f"[{'PASS' if mean_ok else 'FAIL'}, expected {MEAN_LOW}-{MEAN_HIGH}%]")

    # Compare with expected
    expected_path = os.path.join(EXPECTED_DIR, "metrics.json")
    if os.path.exists(expected_path):
        with open(expected_path) as f:
            expected = json.load(f)
        print("\n" + "-" * 72)
        print("EXPECTED-VALUE COMPARISON")
        print("-" * 72)
        for key in expected:
            exp_sum = expected[key]["comp_sum"]
            act_sum = results.get(key, {}).get("comp_sum", None)
            if act_sum is not None and abs(exp_sum - act_sum) > 0.5:
                print(f"  {key}: expected sum={exp_sum}  actual={act_sum}  DRIFT")

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(NS))
    w = 0.35
    bars_i = ax.bar(x - w / 2, ic_inw, w, label="IC(inward)", color="#2196F3")
    bars_o = ax.bar(x + w / 2, ic_out, w, label="IC(outward)", color="#FF9800")
    ax.axhline(100, color="red", linestyle="--", linewidth=1, label="100% reference")

    # Plot sum as line
    ax2 = ax.twinx()
    ax2.plot(x, sums, "ko-", linewidth=2, markersize=8, label="Sum")
    ax2.axhline(100, color="red", linestyle="--", linewidth=0.5)
    ax2.set_ylabel("Complementarity Sum (%)")
    ax2.set_ylim(90, 130)
    ax2.legend(loc="upper right")

    ax.set_xticks(x)
    ax.set_xticklabels([str(n) for n in NS])
    ax.set_xlabel("N")
    ax.set_ylabel("IC Fraction (%)")
    ax.set_title("Law II: Radial Complementarity — IC(inward) + IC(outward)")
    ax.legend(loc="upper left")
    plt.tight_layout()
    plot_path = os.path.join(PLOT_DIR, "law_II_complementarity.png")
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f"\n  Plot saved: {plot_path}")

    # Verdict
    print("\n" + "=" * 72)
    if all_pass:
        print("LAW II VERIFICATION: PASS")
        print(f"  Complementarity confirmed: mean sum = {mean_s:.1f}% for N >= 4.")
    else:
        print("LAW II VERIFICATION: FAIL")
    print("=" * 72)

    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
