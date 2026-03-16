#!/usr/bin/env python
"""
Reproduce_This_Law — Law VII: Tangent DECAY Exclusivity Law
============================================================
Verifies: DECAY occurs if and only if gate=tangent AND N=4.
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
DECAY_MIN_AT_TAN4 = 0.50


def load_sweep(gate, n):
    path = os.path.join(SWEEP_DIR, f"n{n}_{gate}_sweep.json")
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict) and "results" in data:
        data = data["results"]
    return data


def main():
    os.makedirs(PLOT_DIR, exist_ok=True)
    all_pass = True

    print("=" * 72)
    print("LAW VII VERIFICATION: Tangent DECAY Exclusivity Law")
    print("=" * 72)

    decay_grid = np.zeros((len(GATES), len(NS)))
    frac_grid = np.zeros((len(GATES), len(NS)))
    nonzero_cells = []
    total_zero = 0
    total_decay = 0

    for gi, gate in enumerate(GATES):
        for ni, n in enumerate(NS):
            recs = load_sweep(gate, n)
            total = len(recs)
            dc = sum(1 for r in recs if r["mechanism"] == "DECAY")
            df = dc / total
            decay_grid[gi, ni] = dc
            frac_grid[gi, ni] = df * 100

            if dc > 0:
                nonzero_cells.append((gate, n, dc, df))
                total_decay += dc
            else:
                total_zero += total

            status = "DECAY" if dc > 0 else "clean"
            print(f"  {gate:>7} N={n}: DECAY={dc:3d}/{total}  "
                  f"({df * 100:5.1f}%)  [{status}]")

    # --- Criterion 1: exactly one nonzero cell ---
    print(f"\n  Nonzero DECAY cells: {len(nonzero_cells)}")
    if len(nonzero_cells) != 1:
        all_pass = False
        print(f"  FAIL: expected 1, found {len(nonzero_cells)}")
        for g, n, c, f in nonzero_cells:
            print(f"    ({g}, N={n}): {c} events ({f * 100:.1f}%)")
    else:
        print(f"  PASS: exactly 1 nonzero cell")

    # --- Criterion 2: that cell is (tangent, N=4) ---
    if nonzero_cells:
        g0, n0, c0, f0 = nonzero_cells[0]
        is_tan4 = (g0 == "tangent" and n0 == 4)
        print(f"\n  Nonzero cell: ({g0}, N={n0})  "
              f"[{'PASS' if is_tan4 else 'FAIL'}]")
        if not is_tan4:
            all_pass = False

        # --- Criterion 3: DECAY fraction >= 50% ---
        if is_tan4:
            if f0 < DECAY_MIN_AT_TAN4:
                all_pass = False
                print(f"  FAIL: DECAY frac {f0:.3f} < {DECAY_MIN_AT_TAN4}")
            else:
                print(f"  PASS: DECAY frac = {f0:.3f} >= {DECAY_MIN_AT_TAN4}")

    print(f"\n  Total zero-DECAY points: {total_zero}")
    print(f"  Total DECAY points: {total_decay}")

    # Compare with expected
    expected_path = os.path.join(EXPECTED_DIR, "metrics.json")
    if os.path.exists(expected_path):
        with open(expected_path) as f:
            expected = json.load(f)
        print("\n" + "-" * 72)
        print("EXPECTED-VALUE COMPARISON")
        print("-" * 72)
        for key, exp_val in expected.get("decay_counts", {}).items():
            parts = key.split("_")
            gate = parts[0]
            n = int(parts[1][1:])
            gi = GATES.index(gate)
            ni = NS.index(n)
            actual = int(decay_grid[gi, ni])
            if actual != exp_val:
                print(f"  {key}: expected={exp_val} actual={actual}  MISMATCH")
                all_pass = False

    # Plot: 3x5 heatmap of DECAY fraction
    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.imshow(frac_grid, aspect="auto", cmap="Reds", vmin=0,
                   vmax=max(80, frac_grid.max()))
    ax.set_xticks(range(len(NS)))
    ax.set_xticklabels([str(n) for n in NS])
    ax.set_yticks(range(len(GATES)))
    ax.set_yticklabels(GATES)
    ax.set_xlabel("N")
    ax.set_title("Law VII: DECAY Fraction Across All 15 (gate, N) Cells")

    for i in range(len(GATES)):
        for j in range(len(NS)):
            v = frac_grid[i, j]
            color = "white" if v > 40 else "black"
            txt = f"{v:.1f}%" if v > 0 else "0"
            weight = "bold" if v > 0 else "normal"
            ax.text(j, i, txt, ha="center", va="center",
                    fontsize=9, color=color, fontweight=weight)

    # Highlight the singular cell
    from matplotlib.patches import Rectangle
    tan_idx = GATES.index("tangent")
    n4_idx = NS.index(4)
    rect = Rectangle((n4_idx - 0.5, tan_idx - 0.5), 1, 1,
                      linewidth=3, edgecolor="gold", facecolor="none")
    ax.add_patch(rect)

    plt.colorbar(im, ax=ax, label="DECAY fraction (%)")
    plt.tight_layout()
    plot_path = os.path.join(PLOT_DIR, "law_VII_decay_heatmap.png")
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f"\n  Plot saved: {plot_path}")

    print("\n" + "=" * 72)
    if all_pass:
        print("LAW VII VERIFICATION: PASS")
        print("  DECAY exclusivity confirmed: only (tangent, N=4).")
    else:
        print("LAW VII VERIFICATION: FAIL")
    print("=" * 72)
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
