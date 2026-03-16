#!/usr/bin/env python
"""
generate_4d_drift_sweeps.py
============================
Generate sweep files for the 4D drift-robustness experiment:
    N x angle x radius x radial-drift

For each (N, angle, radius, drift) combination, constructs a velocity
vector v = cos(angle)*that + (sin(angle) + drift)*(-rhat), normalizes,
runs the simulation via init_from_arrays, classifies the micro-event,
and saves the result as a JSON file.

Filenames:  n{N}_angle_{label}_r{RRR}_d{DDD}_sweep.json
"""

import json
import math
import os
import sys
import time as time_mod

import numpy as np

# --- Engine imports ---
sys.path.insert(0, os.path.dirname(__file__))
from event_lattice import RingParams, RingState
from event_update import apply_pbc
from micro_event_operator import detect_micro_event


# ============================================================
# Constants (pixel-to-engine conversion, matching prior sweeps)
# ============================================================

BOX_PX      = 400
DT          = 1.0        # pixel-space dt
MERGE_THR_PX = 23.5      # pixel-space merge threshold
T_DECAY     = 100
SCALE       = 1.0 / BOX_PX

# Engine-space parameters
DT_ENGINE       = DT * SCALE           # 0.0025
MERGE_THR_ENGINE = MERGE_THR_PX * SCALE  # 0.05875


# ============================================================
# Parameter grid
# ============================================================

NS     = [4, 8, 12, 20]
ANGLES = [-10, -5, -2, -1, 0, 1, 2, 5, 10]
RADII  = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 0.90, 0.95]
DRIFTS = [0.00, 0.01, 0.02, 0.05, 0.10]

SWEEP_DIR = os.path.dirname(__file__)


# ============================================================
# Filename helpers
# ============================================================

def angle_label(ang):
    if ang > 0:   return f"p{ang}deg"
    if ang < 0:   return f"m{abs(ang)}deg"
    return "0deg"

def radius_label(r_frac):
    return f"r{int(round(r_frac * 100)):03d}"

def drift_label(d):
    return f"d{int(round(d * 1000)):03d}"

def output_path(n, ang, r_frac, drift):
    fname = (f"n{n}_angle_{angle_label(ang)}_"
             f"{radius_label(r_frac)}_{drift_label(drift)}_sweep.json")
    return os.path.join(SWEEP_DIR, fname)


# ============================================================
# Single simulation
# ============================================================

def run_one(ncl, gamma_deg, r_frac, drift):
    """
    Run a single simulation at (N, angle, radius, drift).

    Returns a dict with the sweep record, or None on error.
    """
    d_px = int(round(r_frac * BOX_PX))
    if d_px < 1:
        d_px = 1

    sin_pi_n = math.sin(math.pi / ncl)
    diameter_engine = (d_px / sin_pi_n) * SCALE

    # Build RingParams — use gamma_gate='tangent' as label
    # (actual velocity direction is set via init_from_arrays)
    params = RingParams(
        N=ncl,
        diameter=diameter_engine,
        gamma_gate="tangent",
        box_size=1.0,
        dt=DT_ENGINE,
        merge_thr=MERGE_THR_ENGINE,
        speed=1.0,
        max_steps=10_000,
    )

    # Positions: equispaced on ring centered at (0.5, 0.5)
    radius_engine = diameter_engine / 2.0
    center = np.array([0.5, 0.5])
    theta = np.linspace(0, 2 * np.pi, ncl, endpoint=False)
    positions = np.column_stack([
        center[0] + radius_engine * np.cos(theta),
        center[1] + radius_engine * np.sin(theta),
    ])
    positions = apply_pbc(positions, 1.0)

    # Radial and tangent unit vectors
    radial = positions - center
    rn = np.maximum(np.linalg.norm(radial, axis=1, keepdims=True), 1e-15)
    rhat = radial / rn
    that = np.column_stack([-rhat[:, 1], rhat[:, 0]])

    # Velocity: angle component + drift component, then normalize
    gamma_rad = math.radians(gamma_deg)
    # Convention: positive angle -> inward radial component
    # v_base = cos(gamma) * tangent + sin(gamma) * (-radial)
    # drift adds extra inward radial: drift * (-radial)
    v_raw = (math.cos(gamma_rad) * that
             + (math.sin(gamma_rad) + drift) * (-rhat))

    # Normalize each velocity to unit speed
    v_norms = np.linalg.norm(v_raw, axis=1, keepdims=True)
    v_norms = np.maximum(v_norms, 1e-15)
    vel = v_raw / v_norms

    # Initialize and run
    ring = RingState(params)
    ring.init_from_arrays(positions, vel)
    me = detect_micro_event(ring)

    # Compute chi_emp (scaled collapse time)
    chi_val = round(ring.time * BOX_PX / (DT * T_DECAY), 3)
    if chi_val < 0.01:
        chi_val = 0.01

    return {
        "N":         ncl,
        "angle":     gamma_deg,
        "d_px":      d_px,
        "radius":    r_frac,
        "drift":     drift,
        "mechanism": me.mechanism,
        "chi_emp":   chi_val,
        "d_min":     round(me.d_min_at_collapse, 6),
        "dR":        round(me.delta_R_at_collapse, 6),
    }


# ============================================================
# Main generator
# ============================================================

def main():
    total = len(NS) * len(ANGLES) * len(RADII) * len(DRIFTS)
    print(f"4D Drift Sweep Generator")
    print(f"========================")
    print(f"  N values:  {NS}")
    print(f"  Angles:    {ANGLES}")
    print(f"  Radii:     {RADII}")
    print(f"  Drifts:    {DRIFTS}")
    print(f"  Total combinations: {total}")
    print()

    generated = 0
    skipped = 0
    errors = 0
    extended = []   # sweeps that hit max_steps

    t0 = time_mod.time()

    for ni, n in enumerate(NS):
        for ai, ang in enumerate(ANGLES):
            for ri, r in enumerate(RADII):
                for di, drift in enumerate(DRIFTS):
                    path = output_path(n, ang, r, drift)

                    # Skip if already exists
                    if os.path.exists(path):
                        skipped += 1
                        continue

                    try:
                        result = run_one(n, ang, r, drift)
                    except Exception as e:
                        print(f"  ERROR: N={n} ang={ang:+d} "
                              f"r={r:.2f} d={drift:.3f}: {e}")
                        errors += 1
                        continue

                    # Check for extended runtime (DECAY = hit max_steps)
                    if result["mechanism"] == "DECAY":
                        extended.append(
                            f"N={n} ang={ang:+d} r={r:.2f} d={drift:.3f}")

                    # Save
                    with open(path, "w") as f:
                        json.dump(result, f, indent=2)

                    generated += 1

        # Progress per N
        elapsed = time_mod.time() - t0
        done = (ni + 1) * len(ANGLES) * len(RADII) * len(DRIFTS)
        pct = 100 * done / total
        print(f"  N={n:2d} complete  "
              f"({done}/{total}, {pct:.0f}%, "
              f"{elapsed:.1f}s elapsed)")

    elapsed = time_mod.time() - t0

    # Summary
    print()
    print(f"=" * 60)
    print(f"  GENERATION SUMMARY")
    print(f"=" * 60)
    print(f"  Total combinations:  {total}")
    print(f"  Generated (new):     {generated}")
    print(f"  Skipped (existing):  {skipped}")
    print(f"  Errors:              {errors}")
    print(f"  Extended runtime:    {len(extended)}")
    print(f"  Wall time:           {elapsed:.1f}s")

    if errors > 0:
        print(f"\n  WARNING: {errors} combinations failed!")

    missing = total - generated - skipped
    if missing > 0 and missing != errors:
        print(f"\n  Missing combinations: {missing}")

    if extended:
        print(f"\n  Extended runtime (DECAY) cases ({len(extended)}):")
        for desc in extended[:20]:
            print(f"    {desc}")
        if len(extended) > 20:
            print(f"    ... and {len(extended) - 20} more")

    print()
    print(f"  Files written to: {SWEEP_DIR}")
    print(f"  Run check_4d_drift_robustness.py to analyze results.")


if __name__ == "__main__":
    main()
