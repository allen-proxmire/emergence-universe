#!/usr/bin/env python
"""
check_gate_stability.py
========================
Gate-stability analysis: how sensitive are mechanism distributions to
small angular perturbations of the initial velocity direction?

All files (baselines + perturbations) use the init_from_arrays path
with gamma_gate='tangent' for uniform classification.
Baselines: n{N}_{gate}_0deg_sweep.json  (angle = -90, 0, +90)
Perturbed: n{N}_{gate}_{p|m}{1|2|5}deg_sweep.json
"""

import json
import os
import sys

SWEEP_DIR = os.path.dirname(__file__)
NS = [4, 8, 12, 20]
GATES = ["inward", "tangent", "outward"]
GATE_ANGLES = {"inward": -90, "tangent": 0, "outward": 90}
OFFSETS = [-5, -2, -1, 0, 1, 2, 5]
MECHS = ["inward-collapse", "outward-PBC", "DECAY", "PBC-corner", "other-late"]
MECH_SHORT = {"inward-collapse": "IC", "outward-PBC": "OPBC", "DECAY": "DEC",
              "PBC-corner": "CRN", "other-late": "OL"}


def load(gate, n, offset):
    if offset == 0:
        fname = f"n{n}_{gate}_0deg_sweep.json"
    else:
        label = f"p{offset}deg" if offset > 0 else f"m{abs(offset)}deg"
        fname = f"n{n}_{gate}_{label}_sweep.json"
    path = os.path.join(SWEEP_DIR, fname)
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict) and "results" in data:
        data = data["results"]
    return data


def mech_fracs(records):
    total = len(records)
    counts = {}
    for r in records:
        m = r["mechanism"]
        counts[m] = counts.get(m, 0) + 1
    return {m: 100 * counts.get(m, 0) / total for m in MECHS}


def main():
    # =========================================================
    # Section 1: Full sensitivity table per gate
    # =========================================================
    for gate in GATES:
        print("=" * 90)
        print(f"  GATE: {gate.upper()} (base angle = {GATE_ANGLES[gate]}°)")
        print("=" * 90)

        header = f"  {'N':>3}  {'off':>4}  {'angle':>6}"
        for m in MECHS:
            header += f"  {MECH_SHORT[m]:>5}"
        header += "  delta_IC"
        print(header)
        print("  " + "-" * 82)

        for n in NS:
            base_fracs = mech_fracs(load(gate, n, 0))
            for off in OFFSETS:
                ang = GATE_ANGLES[gate] + off
                fracs = mech_fracs(load(gate, n, off))
                delta_ic = fracs["inward-collapse"] - base_fracs["inward-collapse"]

                row = f"  {n:3d}  {off:+4d}  {ang:+5d}°"
                for m in MECHS:
                    row += f"  {fracs[m]:5.1f}"
                row += f"  {delta_ic:+6.1f}"
                if off == 0:
                    row += "  <-- BASE"
                print(row)
            print()

    # =========================================================
    # Section 2: IC sensitivity summary (delta per degree)
    # =========================================================
    print("=" * 90)
    print("  IC SENSITIVITY SUMMARY: |delta_IC| per degree of perturbation")
    print("=" * 90)
    print(f"  {'N':>3}  {'gate':>8}  {'base IC':>8}  "
          f"{'@+/-1°':>7}  {'@+/-2°':>7}  {'@+/-5°':>7}  "
          f"{'slope':>7}  stability")
    print("  " + "-" * 76)

    for n in NS:
        for gate in GATES:
            base_ic = mech_fracs(load(gate, n, 0))["inward-collapse"]
            deltas = {}
            for off in [1, 2, 5]:
                ic_p = mech_fracs(load(gate, n, off))["inward-collapse"]
                ic_m = mech_fracs(load(gate, n, -off))["inward-collapse"]
                deltas[off] = (abs(ic_p - base_ic) + abs(ic_m - base_ic)) / 2

            # Slope: average IC change per degree
            slope = deltas[5] / 5

            if slope < 0.5:
                stab = "STABLE"
            elif slope < 3:
                stab = "MODERATE"
            elif slope < 10:
                stab = "SENSITIVE"
            else:
                stab = "FRAGILE"

            print(f"  {n:3d}  {gate:>8}  {base_ic:7.1f}%  "
                  f"{deltas[1]:6.1f}  {deltas[2]:6.1f}  {deltas[5]:6.1f}  "
                  f"{slope:6.2f}  {stab}")
        print()

    # =========================================================
    # Section 3: Key law-specific checks
    # =========================================================
    print("=" * 90)
    print("  LAW-SPECIFIC STABILITY CHECKS")
    print("=" * 90)

    # Law VI: N=4 tangent DECAY
    print("\n  Law VI — N=4 tangent DECAY under perturbation:")
    for off in OFFSETS:
        fracs = mech_fracs(load("tangent", 4, off))
        ang = off
        print(f"    offset={off:+2d}° (angle={ang:+3d}°): "
              f"IC={fracs['inward-collapse']:5.1f}%  "
              f"DECAY={fracs['DECAY']:5.1f}%  "
              f"OPBC={fracs['outward-PBC']:5.1f}%")

    # N=8 anomaly
    print("\n  N=8 anomaly — inward gate under perturbation:")
    for off in OFFSETS:
        fracs = mech_fracs(load("inward", 8, off))
        ang = -90 + off
        print(f"    offset={off:+2d}° (angle={ang:+4d}°): "
              f"IC={fracs['inward-collapse']:5.1f}%  "
              f"CRN={fracs['PBC-corner']:5.1f}%  "
              f"OL={fracs['other-late']:5.1f}%")

    # Outward gate: all N show 100% IC — check if it breaks
    print("\n  Outward gate universality (IC=100% for all N, all offsets?):")
    for n in NS:
        all_100 = True
        for off in OFFSETS:
            fracs = mech_fracs(load("outward", n, off))
            if fracs["inward-collapse"] < 100.0:
                all_100 = False
                ang = 90 + off
                print(f"    N={n} offset={off:+2d}° (angle={ang:+4d}°): IC={fracs['inward-collapse']:.1f}% !!!")
        if all_100:
            print(f"    N={n}: all offsets IC=100.0%  [STABLE]")

    # Tangent asymmetry: +offset vs -offset
    print("\n  Tangent gate asymmetry (IC at +offset vs -offset):")
    for n in NS:
        print(f"    N={n}:", end="")
        for off in [1, 2, 5]:
            ic_p = mech_fracs(load("tangent", n, off))["inward-collapse"]
            ic_m = mech_fracs(load("tangent", n, -off))["inward-collapse"]
            print(f"  +{off}°={ic_p:5.1f}  -{off}°={ic_m:5.1f}  gap={ic_p-ic_m:+6.1f}", end="")
        print()

    print("\n" + "=" * 90)


if __name__ == "__main__":
    main()
