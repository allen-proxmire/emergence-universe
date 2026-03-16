#!/usr/bin/env python
"""
Reproduce_This_Law — Law I: Binary Partition Law
=================================================
Verifies: Under radial gates (inward/outward), mechanism space contains
only {inward-collapse, outward-PBC, other-late}. Zero DECAY, zero PBC-corner.
"""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SWEEP_DIR = os.path.join(os.path.dirname(__file__), "..", "..")
EXPECTED_DIR = os.path.join(os.path.dirname(__file__), "expected_outputs")
PLOT_DIR = os.path.join(os.path.dirname(__file__), "plots")

GATES = ["inward", "outward"]
NS = [3, 4, 5, 6, 7]
ALLOWED_MECHS = {"inward-collapse", "outward-PBC", "other-late"}
OL_THRESHOLD = 0.06  # max other-late fraction

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_sweep(gate, n):
    path = os.path.join(SWEEP_DIR, f"n{n}_{gate}_sweep.json")
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict) and "results" in data:
        data = data["results"]
    return data


def count_mechanisms(records):
    counts = {}
    for r in records:
        m = r["mechanism"]
        counts[m] = counts.get(m, 0) + 1
    return counts


# ---------------------------------------------------------------------------
# Main verification
# ---------------------------------------------------------------------------

def main():
    os.makedirs(PLOT_DIR, exist_ok=True)

    all_pass = True
    results = {}
    heatmap_data = {mech: np.zeros((len(GATES), len(NS))) for mech in
                    ["inward-collapse", "outward-PBC", "other-late", "DECAY", "PBC-corner"]}

    print("=" * 72)
    print("LAW I VERIFICATION: Binary Partition Law")
    print("=" * 72)

    for gi, gate in enumerate(GATES):
        for ni, n in enumerate(NS):
            records = load_sweep(gate, n)
            total = len(records)
            mc = count_mechanisms(records)

            ic = mc.get("inward-collapse", 0)
            opbc = mc.get("outward-PBC", 0)
            ol = mc.get("other-late", 0)
            decay = mc.get("DECAY", 0)
            corner = mc.get("PBC-corner", 0)

            ic_f = ic / total
            opbc_f = opbc / total
            ol_f = ol / total

            # Fill heatmap
            for mech in heatmap_data:
                heatmap_data[mech][gi, ni] = mc.get(mech, 0) / total * 100

            # Check criteria
            cell_pass = True
            issues = []

            if decay > 0:
                cell_pass = False
                issues.append(f"DECAY={decay}")
            if corner > 0:
                cell_pass = False
                issues.append(f"PBC-corner={corner}")
            if ol_f > OL_THRESHOLD:
                cell_pass = False
                issues.append(f"other-late={ol_f:.3f} > {OL_THRESHOLD}")

            unexpected = set(mc.keys()) - ALLOWED_MECHS
            if unexpected:
                cell_pass = False
                issues.append(f"unexpected mechanisms: {unexpected}")

            status = "PASS" if cell_pass else "FAIL"
            if not cell_pass:
                all_pass = False

            print(f"\n  {gate:>7} N={n}: IC={ic_f:.3f}  OPBC={opbc_f:.3f}  "
                  f"OL={ol_f:.3f}  DECAY={decay}  Corner={corner}  [{status}]")
            if issues:
                for iss in issues:
                    print(f"    !! {iss}")

            results[f"{gate}_N{n}"] = {
                "ic_frac": round(ic_f, 4),
                "opbc_frac": round(opbc_f, 4),
                "other_late_frac": round(ol_f, 4),
                "decay_count": decay,
                "pbc_corner_count": corner,
                "pass": cell_pass,
            }

    # -----------------------------------------------------------------------
    # Compare with expected metrics
    # -----------------------------------------------------------------------
    expected_path = os.path.join(EXPECTED_DIR, "metrics.json")
    if os.path.exists(expected_path):
        with open(expected_path) as f:
            expected = json.load(f)
        print("\n" + "-" * 72)
        print("EXPECTED-VALUE COMPARISON")
        print("-" * 72)
        for key in expected:
            exp = expected[key]
            act = results.get(key, {})
            if exp.get("decay_count") != act.get("decay_count"):
                print(f"  {key}: DECAY expected={exp['decay_count']} "
                      f"actual={act.get('decay_count')} MISMATCH")
                all_pass = False
            if exp.get("pbc_corner_count") != act.get("pbc_corner_count"):
                print(f"  {key}: PBC-corner expected={exp['pbc_corner_count']} "
                      f"actual={act.get('pbc_corner_count')} MISMATCH")
                all_pass = False

    # -----------------------------------------------------------------------
    # Generate plot
    # -----------------------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(14, 3.5))
    for idx, mech in enumerate(["inward-collapse", "outward-PBC", "other-late"]):
        ax = axes[idx]
        data = heatmap_data[mech]
        im = ax.imshow(data, aspect="auto", cmap="YlOrRd", vmin=0, vmax=100)
        ax.set_xticks(range(len(NS)))
        ax.set_xticklabels([str(n) for n in NS])
        ax.set_yticks(range(len(GATES)))
        ax.set_yticklabels(GATES)
        ax.set_xlabel("N")
        ax.set_title(mech)
        for i in range(len(GATES)):
            for j in range(len(NS)):
                ax.text(j, i, f"{data[i, j]:.1f}%", ha="center", va="center",
                        fontsize=8, color="black" if data[i, j] < 60 else "white")
    plt.colorbar(im, ax=axes[-1], label="fraction (%)")
    fig.suptitle("Law I: Binary Partition — Radial Gate Mechanism Fractions", fontsize=12)
    plt.tight_layout()
    plot_path = os.path.join(PLOT_DIR, "law_I_mechanism_heatmap.png")
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f"\n  Plot saved: {plot_path}")

    # -----------------------------------------------------------------------
    # Verdict
    # -----------------------------------------------------------------------
    print("\n" + "=" * 72)
    if all_pass:
        print("LAW I VERIFICATION: PASS")
        print("  Binary partition confirmed: zero DECAY, zero PBC-corner")
        print("  across all 10 radial (gate, N) pairs.")
    else:
        print("LAW I VERIFICATION: FAIL")
        print("  One or more criteria violated. See details above.")
    print("=" * 72)

    sys.exit(0 if all_pass else 1)


if __name__ == "__main__":
    main()
