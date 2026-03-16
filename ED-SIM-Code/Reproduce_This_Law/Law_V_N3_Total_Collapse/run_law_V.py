#!/usr/bin/env python
"""
Reproduce_This_Law — Law V: N=3 Inward Total Collapse Law
==========================================================
Verifies: Every diameter at N=3 inward produces inward-collapse (100% IC).
This is the unique (gate, N) pair achieving total collapse.
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

GATES = ["inward", "tangent", "outward"]
NS = [3, 4, 5, 6, 7]


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
    os.makedirs(PLOT_DIR, exist_ok=True)
    all_pass = True

    print("=" * 72)
    print("LAW V VERIFICATION: N=3 Inward Total Collapse Law")
    print("=" * 72)

    # Primary check: N=3 inward
    records = load_sweep("inward", 3)
    total = len(records)
    ic_count = sum(1 for r in records if r["mechanism"] == "inward-collapse")
    frac = ic_count / total

    print(f"\n  N=3 inward: {ic_count}/{total} IC = {frac * 100:.1f}%")

    if frac != 1.0:
        all_pass = False
        non_ic = [r for r in records if r["mechanism"] != "inward-collapse"]
        print(f"  FAIL: {len(non_ic)} non-IC events found!")
        for r in non_ic[:5]:
            chi_key = "chi" if "chi" in r else "chi_emp"
            print(f"    d_px={r.get('d_px', '?')} mechanism={r['mechanism']} "
                  f"chi={r.get(chi_key, '?')}")
    else:
        print("  PASS: 100% IC confirmed.")

    # Uniqueness check: scan all 15 cells
    print("\n  UNIQUENESS CHECK: no other cell at 100% IC")
    other_100 = []
    ic_grid = {}
    for gate in GATES:
        for n in NS:
            if gate == "inward" and n == 3:
                ic_grid[(gate, n)] = 100.0
                continue
            recs = load_sweep(gate, n)
            f = ic_frac(recs) * 100
            ic_grid[(gate, n)] = f
            if f == 100.0:
                other_100.append((gate, n))
                print(f"    FAIL: {gate} N={n} also has 100% IC!")
                all_pass = False

    if not other_100:
        print("    PASS: (inward, N=3) is the unique 100% IC cell.")

    # Compare with expected
    expected_path = os.path.join(EXPECTED_DIR, "metrics.json")
    if os.path.exists(expected_path):
        with open(expected_path) as f:
            expected = json.load(f)
        exp_count = expected.get("n3_inward_ic_count", 385)
        if ic_count != exp_count:
            print(f"\n  EXPECTED MISMATCH: expected {exp_count}, got {ic_count}")

    # Plot: 3x5 grid of IC fractions, highlighting the singular cell
    fig, ax = plt.subplots(figsize=(8, 4))
    data = np.zeros((len(GATES), len(NS)))
    for gi, gate in enumerate(GATES):
        for ni, n in enumerate(NS):
            data[gi, ni] = ic_grid[(gate, n)]

    im = ax.imshow(data, aspect="auto", cmap="RdYlGn", vmin=0, vmax=100)
    ax.set_xticks(range(len(NS)))
    ax.set_xticklabels([str(n) for n in NS])
    ax.set_yticks(range(len(GATES)))
    ax.set_yticklabels(GATES)
    ax.set_xlabel("N")
    ax.set_title("Law V: IC Fraction Across All 15 (gate, N) Cells")

    for i in range(len(GATES)):
        for j in range(len(NS)):
            v = data[i, j]
            color = "white" if v > 60 or v < 10 else "black"
            weight = "bold" if v == 100.0 else "normal"
            ax.text(j, i, f"{v:.1f}%", ha="center", va="center",
                    fontsize=9, color=color, fontweight=weight)

    # Highlight the singular cell
    from matplotlib.patches import Rectangle
    rect = Rectangle((-0.5, -0.5), 1, 1, linewidth=3, edgecolor="red",
                     facecolor="none")
    ax.add_patch(rect)

    plt.colorbar(im, ax=ax, label="IC fraction (%)")
    plt.tight_layout()
    plot_path = os.path.join(PLOT_DIR, "law_V_total_collapse.png")
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f"\n  Plot saved: {plot_path}")

    print("\n" + "=" * 72)
    if all_pass:
        print("LAW V VERIFICATION: PASS")
        print("  N=3 inward total collapse confirmed: 100% IC, unique cell.")
    else:
        print("LAW V VERIFICATION: FAIL")
    print("=" * 72)
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
