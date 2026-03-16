#!/usr/bin/env python
"""
Reproduce_This_Law — Law VI: N=4 Rotational Degeneracy Law
===========================================================
Verifies: N=4 tangent has 0% IC, ~77% DECAY. Radial gates at N=4 are normal.
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

DECAY_MIN = 0.70
RADIAL_IC_MIN = 0.20


def load_sweep(gate, n):
    path = os.path.join(SWEEP_DIR, f"n{n}_{gate}_sweep.json")
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict) and "results" in data:
        data = data["results"]
    return data


def mech_dist(records):
    counts = {}
    for r in records:
        m = r["mechanism"]
        counts[m] = counts.get(m, 0) + 1
    return counts


def main():
    os.makedirs(PLOT_DIR, exist_ok=True)
    all_pass = True

    print("=" * 72)
    print("LAW VI VERIFICATION: N=4 Rotational Degeneracy Law")
    print("=" * 72)

    # --- N=4 tangent: the anomalous cell ---
    tan = load_sweep("tangent", 4)
    total_t = len(tan)
    mc_t = mech_dist(tan)

    ic_t = mc_t.get("inward-collapse", 0)
    decay_t = mc_t.get("DECAY", 0)
    decay_f = decay_t / total_t

    print(f"\n  N=4 TANGENT (anomalous cell):")
    print(f"    Total records: {total_t}")
    for m, c in sorted(mc_t.items(), key=lambda x: -x[1]):
        print(f"    {m:>20}: {c:3d} ({100 * c / total_t:5.1f}%)")

    # Criterion 1: zero IC
    if ic_t != 0:
        all_pass = False
        print(f"    FAIL: IC count = {ic_t} (expected 0)")
    else:
        print(f"    PASS: IC count = 0")

    # Criterion 2: DECAY >= 70%
    if decay_f < DECAY_MIN:
        all_pass = False
        print(f"    FAIL: DECAY frac = {decay_f:.3f} (expected >= {DECAY_MIN})")
    else:
        print(f"    PASS: DECAY frac = {decay_f:.3f}")

    # --- N=4 radial gates: normal behavior ---
    for gate in ["inward", "outward"]:
        recs = load_sweep(gate, 4)
        total = len(recs)
        mc = mech_dist(recs)
        ic = mc.get("inward-collapse", 0)
        decay = mc.get("DECAY", 0)
        ic_f = ic / total

        print(f"\n  N=4 {gate.upper()} (control):")
        print(f"    IC = {ic}/{total} ({ic_f * 100:.1f}%)  DECAY = {decay}")

        if ic_f < RADIAL_IC_MIN:
            all_pass = False
            print(f"    FAIL: IC frac {ic_f:.3f} < {RADIAL_IC_MIN}")
        else:
            print(f"    PASS: IC frac >= {RADIAL_IC_MIN}")

        if decay != 0:
            all_pass = False
            print(f"    FAIL: DECAY = {decay} (expected 0)")
        else:
            print(f"    PASS: DECAY = 0")

    # Compare with expected
    expected_path = os.path.join(EXPECTED_DIR, "metrics.json")
    if os.path.exists(expected_path):
        with open(expected_path) as f:
            expected = json.load(f)
        print("\n" + "-" * 72)
        print("EXPECTED-VALUE COMPARISON")
        print("-" * 72)
        exp_ic = expected.get("n4_tangent_ic_count", 0)
        exp_decay_f = expected.get("n4_tangent_decay_frac", 0.774)
        if ic_t != exp_ic:
            print(f"  IC count: expected={exp_ic} actual={ic_t}")
        if abs(decay_f - exp_decay_f) > 0.02:
            print(f"  DECAY frac: expected={exp_decay_f:.3f} actual={decay_f:.3f}")

    # Plot: three pie charts side by side
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    gate_labels = ["inward", "tangent", "outward"]
    for i, gate in enumerate(gate_labels):
        recs = load_sweep(gate, 4)
        mc = mech_dist(recs)
        labels = list(mc.keys())
        sizes = list(mc.values())
        colors_map = {
            "inward-collapse": "#2196F3",
            "outward-PBC": "#FF9800",
            "DECAY": "#f44336",
            "PBC-corner": "#9C27B0",
            "other-late": "#9E9E9E",
        }
        colors = [colors_map.get(l, "#CCCCCC") for l in labels]
        axes[i].pie(sizes, labels=labels, colors=colors, autopct="%1.1f%%",
                    textprops={"fontsize": 7})
        title = f"N=4 {gate}"
        if gate == "tangent":
            title += " (ANOMALY)"
        axes[i].set_title(title, fontweight="bold" if gate == "tangent" else "normal")

    fig.suptitle("Law VI: N=4 Rotational Degeneracy", fontsize=12)
    plt.tight_layout()
    plot_path = os.path.join(PLOT_DIR, "law_VI_n4_anomaly.png")
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f"\n  Plot saved: {plot_path}")

    print("\n" + "=" * 72)
    if all_pass:
        print("LAW VI VERIFICATION: PASS")
        print("  N=4 rotational degeneracy confirmed: 0% IC, 77% DECAY")
        print("  under tangent gate; normal IC under radial gates.")
    else:
        print("LAW VI VERIFICATION: FAIL")
    print("=" * 72)
    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
