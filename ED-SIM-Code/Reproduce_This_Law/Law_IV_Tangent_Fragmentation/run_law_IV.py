#!/usr/bin/env python
"""
Reproduce_This_Law — Law IV: Tangent Fragmentation Law
======================================================
Verifies: Tangent IC shatters into 7-10 bands for N >= 5; radial gates
stay at 1-3 bands; PBC-corner is tangent-exclusive.
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
TANGENT_MIN_BANDS_N5PLUS = 7
RADIAL_MAX_BANDS = 3


def load_sweep(gate, n):
    path = os.path.join(SWEEP_DIR, f"n{n}_{gate}_sweep.json")
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict) and "results" in data:
        data = data["results"]
    return data


def count_ic_bands(records):
    """Count contiguous runs of inward-collapse mechanism."""
    bands = 0
    in_ic = False
    max_width = 0
    current_width = 0
    for r in records:
        if r["mechanism"] == "inward-collapse":
            if not in_ic:
                bands += 1
                in_ic = True
                current_width = 1
            else:
                current_width += 1
        else:
            if in_ic:
                max_width = max(max_width, current_width)
                in_ic = False
                current_width = 0
    if in_ic:
        max_width = max(max_width, current_width)
    return bands, max_width


def main():
    os.makedirs(PLOT_DIR, exist_ok=True)
    all_pass = True

    print("=" * 72)
    print("LAW IV VERIFICATION: Tangent Fragmentation Law")
    print("=" * 72)

    band_counts = {}
    max_widths = {}
    corner_counts = {}

    for gate in GATES:
        band_counts[gate] = {}
        max_widths[gate] = {}
        corner_counts[gate] = {}
        print(f"\n  --- {gate.upper()} GATE ---")
        for n in NS:
            records = load_sweep(gate, n)
            bands, mw = count_ic_bands(records)
            corners = sum(1 for r in records if r["mechanism"] == "PBC-corner")
            band_counts[gate][n] = bands
            max_widths[gate][n] = mw
            corner_counts[gate][n] = corners
            print(f"  N={n}: IC bands={bands:2d}  max_width={mw:3d}  "
                  f"PBC-corner={corners}")

    # --- Criterion 1: Tangent bands >= 7 for N >= 5 ---
    print("\n  CRITERION 1: Tangent IC bands >= 7 for N >= 5")
    for n in [5, 6, 7]:
        b = band_counts["tangent"][n]
        ok = b >= TANGENT_MIN_BANDS_N5PLUS
        print(f"    N={n}: bands={b}  [{'PASS' if ok else 'FAIL'}]")
        if not ok:
            all_pass = False

    # --- Criterion 2: Radial bands <= 3 for all N ---
    print("\n  CRITERION 2: Radial IC bands <= 3 for all N")
    for gate in ["inward", "outward"]:
        for n in NS:
            b = band_counts[gate][n]
            ok = b <= RADIAL_MAX_BANDS
            print(f"    {gate} N={n}: bands={b}  [{'PASS' if ok else 'FAIL'}]")
            if not ok:
                all_pass = False

    # --- Criterion 3: PBC-corner == 0 for radial ---
    print("\n  CRITERION 3: PBC-corner == 0 for radial gates")
    for gate in ["inward", "outward"]:
        for n in NS:
            c = corner_counts[gate][n]
            ok = c == 0
            print(f"    {gate} N={n}: PBC-corner={c}  [{'PASS' if ok else 'FAIL'}]")
            if not ok:
                all_pass = False

    # --- Criterion 4: PBC-corner > 0 for tangent N >= 4 ---
    print("\n  CRITERION 4: PBC-corner > 0 for tangent N >= 4")
    for n in [4, 5, 6, 7]:
        c = corner_counts["tangent"][n]
        ok = c > 0
        print(f"    tangent N={n}: PBC-corner={c}  [{'PASS' if ok else 'FAIL'}]")
        if not ok:
            all_pass = False

    # Compare with expected
    expected_path = os.path.join(EXPECTED_DIR, "metrics.json")
    if os.path.exists(expected_path):
        with open(expected_path) as f:
            expected = json.load(f)
        print("\n" + "-" * 72)
        print("EXPECTED-VALUE COMPARISON")
        print("-" * 72)
        for gate in GATES:
            for n in NS:
                key = f"{gate}_N{n}"
                if key in expected:
                    exp_b = expected[key]["ic_bands"]
                    act_b = band_counts[gate][n]
                    if exp_b != act_b:
                        print(f"  {key}: bands expected={exp_b} actual={act_b}  MISMATCH")

    # Plot
    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(NS))
    w = 0.25
    colors = {"inward": "#2196F3", "tangent": "#4CAF50", "outward": "#FF9800"}
    for gi, gate in enumerate(GATES):
        vals = [band_counts[gate][n] for n in NS]
        ax.bar(x + gi * w - w, vals, w, label=gate, color=colors[gate])
    ax.axhline(TANGENT_MIN_BANDS_N5PLUS, color="red", linestyle="--",
               linewidth=1, label=f"Tangent threshold ({TANGENT_MIN_BANDS_N5PLUS})")
    ax.axhline(RADIAL_MAX_BANDS, color="blue", linestyle=":",
               linewidth=1, label=f"Radial max ({RADIAL_MAX_BANDS})")
    ax.set_xticks(x)
    ax.set_xticklabels([str(n) for n in NS])
    ax.set_xlabel("N")
    ax.set_ylabel("# IC Bands")
    ax.set_title("Law IV: Tangent Fragmentation — IC Band Counts")
    ax.legend(fontsize=8)
    plt.tight_layout()
    plot_path = os.path.join(PLOT_DIR, "law_IV_band_counts.png")
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f"\n  Plot saved: {plot_path}")

    print("\n" + "=" * 72)
    if all_pass:
        print("LAW IV VERIFICATION: PASS")
        print("  Tangent fragmentation confirmed; radial gates contiguous;")
        print("  PBC-corner is tangent-exclusive.")
    else:
        print("LAW IV VERIFICATION: FAIL")
    print("=" * 72)
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
