#!/usr/bin/env python
"""
generate_consolidation_sweeps.py
=================================
Post-Law VIII Consolidation Sweep: generate sweep + atlas files for
the extended parameter space.

Three extensions beyond the canonical 4D grid:
  A. Higher N:   N = 24, 28, 32
  B. Finer angles near tangent:  -1, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1
     (0 is already in the canonical grid; we add 0.25-step resolution)
  C. Finer drifts:  drift = 0.005, 0.015, 0.025

Generates BOTH sweep.json and atlas.json for each new (N, angle, radius, drift).
Reuses the exact same simulation engine and conventions as
generate_4d_drift_sweeps.py and generate_atlas_sweeps.py.

Filenames follow the same n{N}_angle_{label}_r{RRR}_d{DDD}_{type}.json scheme,
with fractional angles encoded as e.g. "p0p25deg" for +0.25 and "m0p75deg" for -0.75.
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
# Constants (identical to canonical generators)
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
# Extended parameter grid
# ============================================================

# Extension A: Higher N
NS_HIGH = [24, 28, 32]

# Extension B: Finer angles near tangent (canonical angles that are new)
ANGLES_FINE = [-0.75, -0.5, -0.25, 0.25, 0.5, 0.75]

# Extension C: Finer drifts
DRIFTS_FINE = [0.005, 0.015, 0.025]

# Canonical grids (for combinations)
NS_CANONICAL     = [4, 8, 12, 20]
ANGLES_CANONICAL = [-10, -5, -2, -1, 0, 1, 2, 5, 10]
RADII            = [0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 0.90, 0.95]
DRIFTS_CANONICAL = [0.00, 0.01, 0.02, 0.05, 0.10]

SWEEP_DIR = os.path.dirname(__file__)


# ============================================================
# Build the list of all extended (N, angle, radius, drift) tuples
# ============================================================

def build_extended_grid():
    """
    Return set of (N, angle, radius, drift) tuples that are NEW—
    i.e., not already covered by the canonical 4x9x11x5 grid.
    """
    canonical = set()
    for n in NS_CANONICAL:
        for a in ANGLES_CANONICAL:
            for r in RADII:
                for d in DRIFTS_CANONICAL:
                    canonical.add((n, a, r, d))

    extended = set()

    # A: Higher N with all canonical angles + drifts
    for n in NS_HIGH:
        for a in ANGLES_CANONICAL:
            for r in RADII:
                for d in DRIFTS_CANONICAL:
                    extended.add((n, a, r, d))

    # B: Finer angles at ALL N (canonical + high) with canonical drifts
    all_ns = NS_CANONICAL + NS_HIGH
    for n in all_ns:
        for a in ANGLES_FINE:
            for r in RADII:
                for d in DRIFTS_CANONICAL:
                    extended.add((n, a, r, d))

    # C: Finer drifts at ALL N with ALL angles (canonical + fine)
    all_angles = ANGLES_CANONICAL + ANGLES_FINE
    for n in all_ns:
        for a in all_angles:
            for r in RADII:
                for d in DRIFTS_FINE:
                    extended.add((n, a, r, d))

    # Remove anything that's already in the canonical grid
    extended -= canonical

    return sorted(extended)


# ============================================================
# Filename helpers (extended for fractional angles)
# ============================================================

def angle_label(ang):
    """Encode angle into filename-safe label.
    Integer angles: p5deg, m2deg, 0deg
    Fractional angles: p0p25deg, m0p75deg
    """
    if ang == 0:
        return "0deg"
    sign = "p" if ang > 0 else "m"
    val = abs(ang)
    if val == int(val):
        return f"{sign}{int(val)}deg"
    # Fractional: replace '.' with 'p'
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
# Simulation: sweep only
# ============================================================

def run_sweep(ncl, gamma_deg, r_frac, drift):
    """Run a single simulation, return sweep record."""
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
    me = detect_micro_event(ring)

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
    }


# ============================================================
# Simulation: instrumented (atlas)
# ============================================================

def run_instrumented(ncl, gamma_deg, r_frac, drift):
    """Run simulation with step-by-step diagnostics. Return sweep+atlas record."""
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

    # --- Trajectory diagnostics ---
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

        # PBC crossing detection
        unwrapped = old_pos + vel * DT_ENGINE
        diff = new_pos - unwrapped
        crossed_mask = np.any(np.abs(diff) > box / 2, axis=1)
        n_crossed = int(np.sum(crossed_mask))
        if n_crossed > 0:
            n_pbc_crossings += n_crossed
            for i in np.where(crossed_mask)[0]:
                a_deg = math.degrees(math.atan2(vel[i, 1], vel[i, 0]))
                crossing_angles.append(round(a_deg, 1))

        # Unwrapped COM
        com_delta = new_com - old_com
        for dim in range(2):
            if com_delta[dim] > box / 2:
                com_delta[dim] -= box
            elif com_delta[dim] < -box / 2:
                com_delta[dim] += box
        com_unwrapped += com_delta

        # Corner and edge proximity
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

        # R tracking
        if step % 50 == 0 or step <= 10:
            R_now = circumradius(new_pos)
            R_history.append(R_now)
            R_max = max(R_max, R_now)
            R_min = min(R_min, R_now)
            dmin_val, _ = d_min(new_pos, box)
            dmin_history.append(dmin_val)

        # Multi-stage snapshots
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

    # Final classification
    me = detect_micro_event(ring)
    mechanism_snapshots.append((1.0, me.mechanism))

    # Winding number
    total_displacement = com_unwrapped - center
    winding_x = round(total_displacement[0] / box, 3)
    winding_y = round(total_displacement[1] / box, 3)

    # Periodicity score
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

    # Bounce-angle histogram
    angle_bins = [0] * 8
    for a in crossing_angles:
        bin_idx = int(((a + 180) % 360) / 45)
        if bin_idx >= 8:
            bin_idx = 7
        angle_bins[bin_idx] += 1

    # Multi-stage detection
    stage_mechs = [m for _, m in mechanism_snapshots]
    unique_stages = []
    for m in stage_mechs:
        if not unique_stages or unique_stages[-1] != m:
            unique_stages.append(m)
    is_multi_stage = len(unique_stages) > 1

    # Tangent crossings
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

    print("Post-Law VIII Consolidation Sweep Generator")
    print("=" * 60)
    print(f"  Higher N:       {NS_HIGH}")
    print(f"  Fine angles:    {ANGLES_FINE}")
    print(f"  Fine drifts:    {DRIFTS_FINE}")
    print(f"  Radii:          {RADII}")
    print(f"  Total NEW combinations: {total}")
    print()

    generated_sweep = 0
    generated_atlas = 0
    skipped = 0
    errors = 0
    t0 = time_mod.time()

    for idx, (n, ang, r, drift) in enumerate(grid):
        sweep_path = output_path(n, ang, r, drift, "sweep")
        atlas_path = output_path(n, ang, r, drift, "atlas")

        sweep_exists = os.path.exists(sweep_path)
        atlas_exists = os.path.exists(atlas_path)

        if sweep_exists and atlas_exists:
            skipped += 1
            continue

        try:
            # Generate atlas (which includes sweep data)
            result = run_instrumented(n, ang, r, drift)
        except Exception as e:
            print(f"  ERROR: N={n} ang={ang:+.2f} r={r:.2f} d={drift:.3f}: {e}")
            errors += 1
            continue

        # Write sweep file
        if not sweep_exists:
            sweep_rec = {k: result[k] for k in
                         ["N", "angle", "d_px", "radius", "drift",
                          "mechanism", "chi_emp", "d_min", "dR"]}
            with open(sweep_path, "w") as f:
                json.dump(sweep_rec, f, indent=2)
            generated_sweep += 1

        # Write atlas file
        if not atlas_exists:
            with open(atlas_path, "w") as f:
                json.dump(result, f, indent=2)
            generated_atlas += 1

        # Progress every 500
        if (idx + 1) % 500 == 0:
            elapsed = time_mod.time() - t0
            rate = (idx + 1) / elapsed if elapsed > 0 else 0
            eta = (total - idx - 1) / rate if rate > 0 else 0
            print(f"  {idx+1}/{total} ({100*(idx+1)/total:.0f}%) "
                  f"  {elapsed:.0f}s elapsed, ~{eta:.0f}s remaining")

    elapsed = time_mod.time() - t0

    print()
    print("=" * 60)
    print("  CONSOLIDATION SWEEP GENERATION SUMMARY")
    print("=" * 60)
    print(f"  Total new combinations:  {total}")
    print(f"  Sweep files generated:   {generated_sweep}")
    print(f"  Atlas files generated:   {generated_atlas}")
    print(f"  Skipped (existing):      {skipped}")
    print(f"  Errors:                  {errors}")
    print(f"  Wall time:               {elapsed:.1f}s")

    if errors > 0:
        print(f"\n  WARNING: {errors} combinations failed!")

    # Breakdown by extension type
    n_high = sum(1 for (n, a, r, d) in grid if n in NS_HIGH)
    n_fine_a = sum(1 for (n, a, r, d) in grid if a in ANGLES_FINE)
    n_fine_d = sum(1 for (n, a, r, d) in grid if d in DRIFTS_FINE)
    print(f"\n  Extension breakdown:")
    print(f"    Higher N (24/28/32):   {n_high} cells")
    print(f"    Fine angles:           {n_fine_a} cells")
    print(f"    Fine drifts:           {n_fine_d} cells")
    print(f"    (categories overlap)")

    print(f"\n  Files written to: {SWEEP_DIR}")
    print(f"  Run check_consolidation_sweep.py to analyze results.")


if __name__ == "__main__":
    main()
