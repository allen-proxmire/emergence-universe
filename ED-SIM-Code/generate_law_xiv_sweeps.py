#!/usr/bin/env python
"""
generate_law_xiv_sweeps.py
============================
Law XIV Consolidation Sweep: generate sweep + atlas files to test the
collapse-time scaling law chi ~ lambda^p, the temporal hierarchy, and
the temporal invisibility of the boundary band.

Extensions beyond prior coverage:
    A. Non-canonical radii (fine-grained for radius-modulation analysis):
       0.075, 0.175, 0.225, 0.325, 0.375, 0.425, 0.475,
       0.575, 0.625, 0.725, 0.775
    B. Non-canonical angles (temporal coverage):
       -12, -6, +6  (others already exist from prior sweeps)
    C. Non-canonical drifts: 0.003, 0.007, 0.013, 0.017, 0.030
       (all exist from prior sweeps; included for cross-products)

All tested at all N = 4,8,12,20,24,28,32,40,48,56,64,80,96,128.
Skips any (N, angle, radius, drift) combination where both sweep and atlas
files already exist.
"""

import json
import math
import os
import sys
import time as time_mod

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from event_lattice import RingParams, RingState
from event_update import (
    apply_pbc, ballistic_step, circumradius, d_min, delta_R, pbc_proximity,
)
from micro_event_operator import detect_micro_event, classify_mechanism


# ============================================================
# Constants (identical to all prior generators)
# ============================================================

BOX_PX       = 400
DT           = 1.0
MERGE_THR_PX = 23.5
T_DECAY      = 100
SCALE        = 1.0 / BOX_PX

DT_ENGINE        = DT * SCALE
MERGE_THR_ENGINE = MERGE_THR_PX * SCALE

CORNER_APPROACH_THR = 0.08
EDGE_APPROACH_THR   = 0.05
PERIODICITY_WINDOW  = 200
CORNERS = np.array([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=float)


# ============================================================
# Parameter grid
# ============================================================

NS_ALL = [4, 8, 12, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 128]

# Canonical values from all prior sweeps
ANGLES_CANONICAL = [-10, -5, -2, -1, 0, 1, 2, 5, 10]
RADII_CANONICAL  = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 0.90, 0.95]
DRIFTS_CANONICAL = [0.00, 0.01, 0.02]

# Prior non-canonical (from Law IX + Law X sweeps)
RADII_PRIOR_NEW  = [0.025, 0.125, 0.275, 0.525, 0.675, 0.825]
ANGLES_PRIOR_NEW = [-7, -3, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 3, 7, 12, 15]
DRIFTS_PRIOR_NEW = [0.003, 0.007, 0.013, 0.017, 0.030]

# Law XIV extensions: fine-grained radii for temporal radius-modulation analysis
RADII_XIV_NEW = [0.075, 0.175, 0.225, 0.325, 0.375, 0.425, 0.475,
                 0.575, 0.625, 0.725, 0.775]

# Law XIV new angles: -12, -6, +6 (not in any prior sweep)
ANGLES_XIV_NEW = [-12, -6, 6]

# Drifts: no truly new drifts; reuse prior non-canonical for cross-products
DRIFTS_XIV_NEW = DRIFTS_PRIOR_NEW

ALL_RADII  = sorted(set(RADII_CANONICAL + RADII_PRIOR_NEW + RADII_XIV_NEW))
ALL_ANGLES = sorted(set(ANGLES_CANONICAL + ANGLES_PRIOR_NEW + ANGLES_XIV_NEW))
ALL_DRIFTS = sorted(set(DRIFTS_CANONICAL + DRIFTS_PRIOR_NEW))

SWEEP_DIR = os.path.dirname(__file__)


# ============================================================
# Build extended grid
# ============================================================

def build_extended_grid():
    """
    Return list of (N, angle, radius, drift) tuples that need generation.
    Includes all cross-combinations involving at least one Law XIV new parameter.
    """
    cells = set()

    # A: New radii x all N x canonical angles x canonical drifts
    for n in NS_ALL:
        for a in ANGLES_CANONICAL:
            for r in RADII_XIV_NEW:
                for d in DRIFTS_CANONICAL:
                    cells.add((n, a, r, d))

    # B: New radii x all N x ALL angles x canonical drifts
    for n in NS_ALL:
        for a in ALL_ANGLES:
            for r in RADII_XIV_NEW:
                for d in DRIFTS_CANONICAL:
                    cells.add((n, a, r, d))

    # C: New angles x all N x ALL radii x canonical drifts
    for n in NS_ALL:
        for a in ANGLES_XIV_NEW:
            for r in ALL_RADII:
                for d in DRIFTS_CANONICAL:
                    cells.add((n, a, r, d))

    # D: New radii x all N x canonical angles x non-canonical drifts
    for n in NS_ALL:
        for a in ANGLES_CANONICAL:
            for r in RADII_XIV_NEW:
                for d in DRIFTS_XIV_NEW:
                    cells.add((n, a, r, d))

    # E: New angles x all N x canonical radii x non-canonical drifts
    for n in NS_ALL:
        for a in ANGLES_XIV_NEW:
            for r in RADII_CANONICAL:
                for d in DRIFTS_XIV_NEW:
                    cells.add((n, a, r, d))

    return sorted(cells)


# ============================================================
# Filename helpers
# ============================================================

def angle_label(ang):
    if ang == 0:
        return "0deg"
    sign = "p" if ang > 0 else "m"
    val = abs(ang)
    if val == int(val):
        return f"{sign}{int(val)}deg"
    val_str = f"{val:.2f}".rstrip('0').rstrip('.')
    val_str = val_str.replace('.', 'p')
    return f"{sign}{val_str}deg"


def radius_label(r):
    return f"r{int(round(r * 100)):03d}"


def drift_label(d):
    return f"d{int(round(d * 1000)):03d}"


def output_path(n, ang, r, drift, suffix):
    fname = (f"n{n}_angle_{angle_label(ang)}_"
             f"{radius_label(r)}_{drift_label(drift)}_{suffix}.json")
    return os.path.join(SWEEP_DIR, fname)


# ============================================================
# Instrumented simulation (produces sweep + atlas in one pass)
# ============================================================

def run_instrumented(ncl, gamma_deg, r_frac, drift):
    d_px = max(1, int(round(r_frac * BOX_PX)))
    sin_pi_n = math.sin(math.pi / ncl)
    diameter_engine = (d_px / sin_pi_n) * SCALE
    box = 1.0

    params = RingParams(
        N=ncl, diameter=diameter_engine, gamma_gate="tangent",
        box_size=box, dt=DT_ENGINE, merge_thr=MERGE_THR_ENGINE,
        speed=1.0, max_steps=10_000,
    )

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

    n_pbc_crossings = 0
    crossing_angles = []
    corner_approaches = [0, 0, 0, 0]
    corner_sequence = []
    n_edge_approaches = 0
    com_unwrapped = center.copy()

    R_history = []
    dmin_history = []
    mechanism_snapshots = []

    R_init = circumradius(ring.positions)
    R_max = R_init
    R_min = R_init

    milestones = {0.25: False, 0.50: False, 0.75: False}
    max_steps = params.max_steps
    step = 0

    while step < max_steps and not ring.collapsed:
        old_pos = ring.positions.copy()
        old_com = np.mean(old_pos, axis=0)

        ring.step()
        step += 1

        new_pos = ring.positions
        new_com = np.mean(new_pos, axis=0)

        unwrapped = old_pos + vel * DT_ENGINE
        diff = new_pos - unwrapped
        crossed_mask = np.any(np.abs(diff) > box / 2, axis=1)
        n_crossed = int(np.sum(crossed_mask))
        if n_crossed > 0:
            n_pbc_crossings += n_crossed
            for i in np.where(crossed_mask)[0]:
                a_deg = math.degrees(math.atan2(vel[i, 1], vel[i, 0]))
                crossing_angles.append(round(a_deg, 1))

        com_delta = new_com - old_com
        for dim in range(2):
            if com_delta[dim] > box / 2:
                com_delta[dim] -= box
            elif com_delta[dim] < -box / 2:
                com_delta[dim] += box
        com_unwrapped += com_delta

        corner_thr = CORNER_APPROACH_THR * box
        edge_thr = EDGE_APPROACH_THR * box
        edge_dist, corner_dist = pbc_proximity(new_pos, box)

        if corner_dist < corner_thr:
            for ci, corner in enumerate(CORNERS):
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

        if step % 50 == 0 or step <= 10:
            R_now = circumradius(new_pos)
            R_history.append(R_now)
            R_max = max(R_max, R_now)
            R_min = min(R_min, R_now)
            dmin_val, _ = d_min(new_pos, box)
            dmin_history.append(dmin_val)

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

    me = detect_micro_event(ring)
    mechanism_snapshots.append((1.0, me.mechanism))

    total_displacement = com_unwrapped - center
    winding_x = round(total_displacement[0] / box, 3)
    winding_y = round(total_displacement[1] / box, 3)

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

    angle_bins = [0] * 8
    for a in crossing_angles:
        bin_idx = int(((a + 180) % 360) / 45)
        if bin_idx >= 8:
            bin_idx = 7
        angle_bins[bin_idx] += 1

    stage_mechs = [m for _, m in mechanism_snapshots]
    unique_stages = []
    for m in stage_mechs:
        if not unique_stages or unique_stages[-1] != m:
            unique_stages.append(m)
    is_multi_stage = len(unique_stages) > 1

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
        "N": ncl, "angle": gamma_deg, "d_px": d_px,
        "radius": r_frac, "drift": drift,
        "mechanism": me.mechanism,
        "chi_emp": chi_val,
        "d_min": round(me.d_min_at_collapse, 6),
        "dR": round(me.delta_R_at_collapse, 6),
        "atlas": {
            "n_pbc_crossings":   n_pbc_crossings,
            "n_edge_approaches": n_edge_approaches,
            "corner_approaches": corner_approaches,
            "corner_sequence":   corner_sequence[:50],
            "winding_number":    [winding_x, winding_y],
            "bounce_angle_bins": angle_bins,
            "R_init":            round(R_init, 6),
            "R_max":             round(R_max, 6),
            "R_min":             round(R_min, 6),
            "R_range":           round(R_max - R_min, 6),
            "tangent_crossings": tangent_crossings,
            "is_multi_stage":    is_multi_stage,
            "stage_sequence":    unique_stages,
            "periodicity_score": periodicity_score,
            "n_steps":           step,
            "min_dmin_seen":     round(min(dmin_history), 6) if dmin_history else 0,
        },
    }


# ============================================================
# Main
# ============================================================

def main():
    grid = build_extended_grid()
    total = len(grid)

    print("Law XIV Consolidation Sweep Generator")
    print("=" * 60)
    print(f"  N values:          {NS_ALL}")
    print(f"  New radii (XIV):   {RADII_XIV_NEW}")
    print(f"  New angles (XIV):  {ANGLES_XIV_NEW}")
    print(f"  Total grid cells:  {total}")
    print()

    generated_sweep = 0
    generated_atlas = 0
    skipped = 0
    errors = 0
    t0 = time_mod.time()

    current_n = None
    n_start = None

    for idx, (n, ang, r, drift) in enumerate(grid):
        if n != current_n:
            if current_n is not None:
                n_elapsed = time_mod.time() - n_start
                total_elapsed = time_mod.time() - t0
                rate = idx / total_elapsed if total_elapsed > 0 else 1
                remaining = total - idx
                eta = remaining / rate if rate > 0 else 0
                print(f"  N={current_n:3d} complete  ({idx}/{total}, "
                      f"{100*idx/total:.0f}%, {n_elapsed:.0f}s this N, "
                      f"~{eta:.0f}s remaining)")
            current_n = n
            n_start = time_mod.time()

        sweep_path = output_path(n, ang, r, drift, "sweep")
        atlas_path = output_path(n, ang, r, drift, "atlas")

        if os.path.exists(sweep_path) and os.path.exists(atlas_path):
            skipped += 1
            continue

        try:
            result = run_instrumented(n, ang, r, drift)
        except Exception as e:
            print(f"  ERROR: N={n} ang={ang} r={r:.3f} d={drift:.3f}: {e}")
            errors += 1
            continue

        if not os.path.exists(sweep_path):
            sweep_rec = {k: result[k] for k in
                         ["N", "angle", "d_px", "radius", "drift",
                          "mechanism", "chi_emp", "d_min", "dR"]}
            with open(sweep_path, "w") as f:
                json.dump(sweep_rec, f, indent=2)
            generated_sweep += 1

        if not os.path.exists(atlas_path):
            with open(atlas_path, "w") as f:
                json.dump(result, f, indent=2)
            generated_atlas += 1

    if current_n is not None:
        n_elapsed = time_mod.time() - n_start
        print(f"  N={current_n:3d} complete  ({total}/{total}, "
              f"100%, {n_elapsed:.0f}s this N)")

    elapsed = time_mod.time() - t0

    print()
    print("=" * 60)
    print("  LAW XIV CONSOLIDATION SWEEP SUMMARY")
    print("=" * 60)
    print(f"  Total grid cells:      {total}")
    print(f"  Sweep files generated: {generated_sweep}")
    print(f"  Atlas files generated: {generated_atlas}")
    print(f"  Skipped (existing):    {skipped}")
    print(f"  Errors:                {errors}")
    print(f"  Wall time:             {elapsed:.1f}s")

    if errors > 0:
        print(f"\n  WARNING: {errors} combinations failed!")

    print(f"\n  Files written to: {SWEEP_DIR}")
    print(f"  Run check_law_xiv_consolidation.py to analyze results.")


if __name__ == "__main__":
    main()
