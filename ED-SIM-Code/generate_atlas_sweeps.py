#!/usr/bin/env python
"""
generate_atlas_sweeps.py
=========================
Generate mechanism-atlas sweep files for the ED sub-mechanism taxonomy.

For each (N, angle, radius, drift) combination, runs the simulation
step-by-step and captures per-trajectory diagnostics:
  - PBC boundary crossing count and angles
  - corner approach count and sequence
  - winding number of center of mass
  - circumradius extrema and trend
  - min approach distance trajectory
  - multi-stage mechanism transitions
  - periodicity score

Filenames: n{N}_angle_{label}_r{RRR}_d{DDD}_atlas.json
"""

import json
import math
import os
import sys
import time as time_mod
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from event_lattice import RingParams, RingState
from event_update import (
    apply_pbc, ballistic_step, circumradius, d_min, delta_R, pbc_proximity,
)
from micro_event_operator import detect_micro_event, classify_mechanism

# ============================================================
# Constants
# ============================================================

BOX_PX       = 400
DT           = 1.0
MERGE_THR_PX = 23.5
T_DECAY      = 100
SCALE        = 1.0 / BOX_PX

DT_ENGINE        = DT * SCALE
MERGE_THR_ENGINE = MERGE_THR_PX * SCALE

CORNER_APPROACH_THR = 0.08   # fraction of box_size for corner approach
EDGE_APPROACH_THR   = 0.05   # fraction of box_size for edge approach
PERIODICITY_WINDOW  = 200    # steps to check for periodicity

CORNERS = np.array([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=float)

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

def radius_label(r):
    return f"r{int(round(r * 100)):03d}"

def drift_label(d):
    return f"d{int(round(d * 1000)):03d}"

def output_path(n, ang, r, drift):
    fname = (f"n{n}_angle_{angle_label(ang)}_"
             f"{radius_label(r)}_{drift_label(drift)}_atlas.json")
    return os.path.join(SWEEP_DIR, fname)


# ============================================================
# Instrumented simulation
# ============================================================

def run_instrumented(ncl, gamma_deg, r_frac, drift):
    """
    Run a single simulation with step-by-step trajectory diagnostics.
    Returns a dict with mechanism + atlas diagnostics.
    """
    d_px = max(1, int(round(r_frac * BOX_PX)))
    sin_pi_n = math.sin(math.pi / ncl)
    diameter_engine = (d_px / sin_pi_n) * SCALE
    box = 1.0

    params = RingParams(
        N=ncl, diameter=diameter_engine, gamma_gate="tangent",
        box_size=box, dt=DT_ENGINE, merge_thr=MERGE_THR_ENGINE,
        speed=1.0, max_steps=10_000,
    )

    # --- init positions + velocities ---
    radius_engine = diameter_engine / 2.0
    center = np.array([0.5, 0.5])
    theta = np.linspace(0, 2 * np.pi, ncl, endpoint=False)
    positions = np.column_stack([
        center[0] + radius_engine * np.cos(theta),
        center[1] + radius_engine * np.sin(theta),
    ])
    positions = apply_pbc(positions, box)

    radial = positions - center
    rn = np.maximum(np.linalg.norm(radial, axis=1, keepdims=True), 1e-15)
    rhat = radial / rn
    that = np.column_stack([-rhat[:, 1], rhat[:, 0]])

    gamma_rad = math.radians(gamma_deg)
    v_raw = (math.cos(gamma_rad) * that
             + (math.sin(gamma_rad) + drift) * (-rhat))
    v_norms = np.maximum(np.linalg.norm(v_raw, axis=1, keepdims=True), 1e-15)
    vel = v_raw / v_norms

    ring = RingState(params)
    ring.init_from_arrays(positions, vel)

    # --- Trajectory diagnostics ---
    n_pbc_crossings = 0
    crossing_angles = []         # velocity angle at each crossing
    corner_approaches = [0, 0, 0, 0]   # count per corner index
    corner_sequence = []         # temporal sequence of corner indices
    n_edge_approaches = 0
    com_unwrapped = center.copy()  # unwrapped center of mass tracking

    R_history = []
    dmin_history = []
    mechanism_snapshots = []     # (step_fraction, mechanism_label)

    R_init = circumradius(ring.positions)
    R_max = R_init
    R_min = R_init

    # Milestone fractions for multi-stage check
    milestones = {0.25: False, 0.50: False, 0.75: False}

    max_steps = params.max_steps
    step = 0

    while step < max_steps and not ring.collapsed:
        old_pos = ring.positions.copy()
        old_com = np.mean(old_pos, axis=0)

        # Step the simulation
        ring.step()
        step += 1

        new_pos = ring.positions
        new_com = np.mean(new_pos, axis=0)

        # --- PBC crossing detection ---
        # unwrapped = old + v*dt
        unwrapped = old_pos + vel * DT_ENGINE
        diff = new_pos - unwrapped
        # Crossings: |diff| > box/2 in either dimension
        crossed_mask = np.any(np.abs(diff) > box / 2, axis=1)
        n_crossed = int(np.sum(crossed_mask))
        if n_crossed > 0:
            n_pbc_crossings += n_crossed
            # Record crossing angles from velocity
            for i in np.where(crossed_mask)[0]:
                angle_deg = math.degrees(math.atan2(vel[i, 1], vel[i, 0]))
                crossing_angles.append(round(angle_deg, 1))

        # --- Unwrapped COM tracking for winding number ---
        com_delta = new_com - old_com
        # If COM jumped by > box/2, it wrapped
        for dim in range(2):
            if com_delta[dim] > box / 2:
                com_delta[dim] -= box
            elif com_delta[dim] < -box / 2:
                com_delta[dim] += box
        com_unwrapped += com_delta

        # --- Corner and edge proximity ---
        corner_thr = CORNER_APPROACH_THR * box
        edge_thr = EDGE_APPROACH_THR * box
        edge_dist, corner_dist = pbc_proximity(new_pos, box)

        if corner_dist < corner_thr:
            # Find which corner
            for ci, corner in enumerate(CORNERS):
                dists_to_corner = np.linalg.norm(
                    new_pos - corner[np.newaxis, :], axis=1)
                # Also check PBC-wrapped distances
                for wrap_x in [-box, 0, box]:
                    for wrap_y in [-box, 0, box]:
                        wrapped_corner = corner + np.array([wrap_x, wrap_y])
                        d2c = np.linalg.norm(
                            new_pos - wrapped_corner[np.newaxis, :], axis=1)
                        if np.min(d2c) < corner_thr:
                            corner_approaches[ci] += 1
                            if (not corner_sequence or
                                    corner_sequence[-1] != ci):
                                corner_sequence.append(ci)
                            break
                    else:
                        continue
                    break

        if edge_dist < edge_thr:
            n_edge_approaches += 1

        # --- R tracking (every 50 steps for efficiency) ---
        if step % 50 == 0 or step <= 10:
            R_now = circumradius(new_pos)
            R_history.append(R_now)
            R_max = max(R_max, R_now)
            R_min = min(R_min, R_now)

            dmin_val, _ = d_min(new_pos, box)
            dmin_history.append(dmin_val)

        # --- Multi-stage mechanism snapshot at milestones ---
        frac = step / max_steps
        for mf in milestones:
            if not milestones[mf] and frac >= mf:
                milestones[mf] = True
                R_now = circumradius(new_pos)
                cumul_dR = R_now - R_init
                e_d, c_d = pbc_proximity(new_pos, box)
                mech = classify_mechanism(
                    collapsed=ring.collapsed,
                    delta_R_val=cumul_dR,
                    pbc_edge_dist=e_d,
                    pbc_corner_dist=c_d,
                    chi_emp=ring.time,
                    box_size=box,
                )
                mechanism_snapshots.append((round(mf, 2), mech.value))

    # --- Final classification ---
    me = detect_micro_event(ring)
    mechanism_snapshots.append((1.0, me.mechanism))

    # --- Winding number ---
    total_displacement = com_unwrapped - center
    winding_x = round(total_displacement[0] / box, 3)
    winding_y = round(total_displacement[1] / box, 3)

    # --- Periodicity score ---
    # Autocorrelation of R_history at lag = PERIODICITY_WINDOW
    periodicity_score = 0.0
    if len(R_history) > PERIODICITY_WINDOW + 10:
        R_arr = np.array(R_history)
        R_mean = np.mean(R_arr)
        R_var = np.var(R_arr)
        if R_var > 1e-15:
            lag = min(PERIODICITY_WINDOW, len(R_arr) // 3)
            shifted = R_arr[lag:]
            original = R_arr[:len(shifted)]
            cov = np.mean((original - R_mean) * (shifted - R_mean))
            periodicity_score = round(float(cov / R_var), 4)

    # --- Bounce-angle histogram (binned to 45-degree sectors) ---
    angle_bins = [0] * 8  # 8 bins of 45 degrees each
    for a in crossing_angles:
        bin_idx = int(((a + 180) % 360) / 45)
        if bin_idx >= 8:
            bin_idx = 7
        angle_bins[bin_idx] += 1

    # --- Multi-stage detection ---
    stage_mechs = [m for _, m in mechanism_snapshots]
    unique_stages = []
    for m in stage_mechs:
        if not unique_stages or unique_stages[-1] != m:
            unique_stages.append(m)
    is_multi_stage = len(unique_stages) > 1

    # --- Tangent crossing ---
    # Did R ever cross back through R_init?
    tangent_crossings = 0
    if len(R_history) > 2:
        above = R_history[0] > R_init
        for R_val in R_history[1:]:
            now_above = R_val > R_init
            if now_above != above:
                tangent_crossings += 1
                above = now_above

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
        "atlas": {
            "n_pbc_crossings":    n_pbc_crossings,
            "n_edge_approaches":  n_edge_approaches,
            "corner_approaches":  corner_approaches,
            "corner_sequence":    corner_sequence[:50],  # cap for JSON size
            "winding_number":     [winding_x, winding_y],
            "bounce_angle_bins":  angle_bins,
            "R_init":             round(R_init, 6),
            "R_max":              round(R_max, 6),
            "R_min":              round(R_min, 6),
            "R_range":            round(R_max - R_min, 6),
            "tangent_crossings":  tangent_crossings,
            "is_multi_stage":     is_multi_stage,
            "stage_sequence":     unique_stages,
            "periodicity_score":  periodicity_score,
            "n_steps":            step,
            "min_dmin_seen":      round(min(dmin_history), 6) if dmin_history else 0,
        },
    }


# ============================================================
# Main
# ============================================================

def main():
    total = len(NS) * len(ANGLES) * len(RADII) * len(DRIFTS)
    print(f"Mechanism Atlas Generator")
    print(f"=========================")
    print(f"  N: {NS}   Angles: {ANGLES}")
    print(f"  Radii: {RADII}")
    print(f"  Drifts: {DRIFTS}")
    print(f"  Total: {total}")
    print()

    generated = 0
    skipped = 0
    errors = 0
    t0 = time_mod.time()

    for ni, n in enumerate(NS):
        for ang in ANGLES:
            for r in RADII:
                for drift in DRIFTS:
                    path = output_path(n, ang, r, drift)
                    if os.path.exists(path):
                        skipped += 1
                        continue
                    try:
                        result = run_instrumented(n, ang, r, drift)
                    except Exception as e:
                        print(f"  ERROR: N={n} a={ang:+d} r={r:.2f} "
                              f"d={drift:.3f}: {e}")
                        errors += 1
                        continue

                    with open(path, "w") as f:
                        json.dump(result, f, indent=2)
                    generated += 1

        elapsed = time_mod.time() - t0
        done = (ni + 1) * len(ANGLES) * len(RADII) * len(DRIFTS)
        print(f"  N={n:2d} complete  ({done}/{total}, "
              f"{100*done/total:.0f}%, {elapsed:.1f}s)")

    elapsed = time_mod.time() - t0
    print(f"\n{'='*60}")
    print(f"  Generated: {generated}  Skipped: {skipped}  "
          f"Errors: {errors}  Time: {elapsed:.1f}s")


if __name__ == "__main__":
    main()
