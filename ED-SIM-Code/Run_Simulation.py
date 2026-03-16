"""
Run_Simulation.py
=================
Top-level runner script for 2D Event Density lattice simulations.

Ties together ED_Lattice, ED_Update_Rule, and ED_Visualization into three
ready-to-run demo scenarios that map directly onto the four cosmological
phases described in ED-05.5:

    SCENARIO A – ED Inflation → Structure Formation  (standard regime)
    SCENARIO B – Two-Body Interaction                (competition / merger)
    SCENARIO C – Black-Hole Saturation               (mobility-weighted)
    SCENARIO D – Noisy Universe                      (Langevin / stochastic ED)

HOW TO RUN
----------
Run all three scenarios and save output to the current directory:

    python Run_Simulation.py

Run a single named scenario:

    python Run_Simulation.py --scenario A
    python Run_Simulation.py --scenario B
    python Run_Simulation.py --scenario C

Override grid size and step count:

    python Run_Simulation.py --scenario A --size 128 --steps 1000

Save animation GIFs (requires Pillow):

    python Run_Simulation.py --scenario A --gif

Suppress all plots (useful in headless / CI environments):

    python Run_Simulation.py --no-show

OUTPUT FILES
------------
For each scenario the following are written to --outdir (default: "output/"):

    <scenario>_overview_t{step}.png    Three-panel snapshot at selected steps
    <scenario>_history.png             Four-panel coarse-grained time series
    <scenario>_phase_diagram.png       p̂(t) vs G(t) trajectory
    <scenario>_snapshots.png           Horizontal strip of field snapshots
    <scenario>_animation.gif           (only with --gif flag)

PARAMETER DEFAULTS
------------------
Scenario A: alpha=0.05, beta=0.25, gamma=0.5, dt=0.1,  size=64, steps=600
Scenario B: alpha=0.05, beta=0.25, gamma=0.5, dt=0.1,  size=64, steps=400
Scenario C: alpha=0.02, beta=0.25, gamma=0.4, dt=0.05, size=64, steps=300
            (mobility-weighted, absorbing boundary; --mobility-exp controls M(rho)^n)
Scenario D: alpha=0.03, beta=0.20, gamma=0.5, dt=0.05, size=64, steps=500
            noise_scale=0.02 (Langevin fluctuations seed structure), periodic boundary
            mobility mode  --  M(rho)^n weights diffusion; --mobility-exp controls n
            --noise-amp overrides noise_scale (default 0.02)
"""

from __future__ import annotations

import argparse
import os
import sys
import warnings
from pathlib import Path
from typing import Optional

import matplotlib
import numpy as np

# ---------------------------------------------------------------------------
# Local imports
# ---------------------------------------------------------------------------
from ED_Lattice import (
    EDLattice,
    EDParams,
    make_black_hole_lattice,
    make_inflation_lattice,
    make_structure_lattice,
)
from ED_Visualization import (
    animate_field,
    plot_history,
    plot_overview,
    plot_phase_diagram,
    plot_snapshot_strip,
    save_figure,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _ensure_dir(path: str) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def _print_header(title: str) -> None:
    line = "=" * 60
    print(f"\n{line}\n  {title}\n{line}")


def _print_stats(lat: EDLattice, label: str = "") -> None:
    s = lat.stats
    tag = f"[{label}] " if label else ""
    print(
        f"  {tag}step={lat.step_count:5d}  t={lat.time:7.2f}  "
        f"p_hat={s['p_hat']:.4f}  G={s['G']:.5f}  "
        f"L={s['L']:.2f}  phase={s['phase']}"
    )


# ---------------------------------------------------------------------------
# SCENARIO A – ED Inflation → Structure Formation
# ---------------------------------------------------------------------------

def run_scenario_A(
    size:          int  = 64,
    steps:         int  = 600,
    outdir:        str  = "output",
    save_gif:      bool = False,
    show:          bool = True,
) -> EDLattice:
    """
    Scenario A: high-ED nearly-uniform initial condition.

    Traces the full sequence:
        Phase 1 (inflation)  →  Phase 2 (structure seeds)
                             →  Phase 3 (structure formation)
                             →  Phase 4 (late thinning)

    Uses:
        make_inflation_lattice()  →  standard ed_step  →  periodic BC
    """
    _print_header("SCENARIO A -- ED Inflation -> Structure Formation")

    lat = make_inflation_lattice(
        size=size,
        noise_scale=0.02,
        seed=42,
    )
    print(f"  Initial state: {lat}")
    print(f"  Running {steps} steps ...")

    # Snapshot at 0, 25%, 50%, 75%, 100% of run
    snap_interval = max(1, steps // 4)
    lat.run(
        steps=steps,
        record_every=5,
        snapshot_every=snap_interval,
        verbose=True,
        verbose_every=max(1, steps // 6),
    )
    print(f"  Final state:   {lat}")
    print(f"  Structures detected: {lat.structure_count()}")

    # ---- Save plots ----
    d = _ensure_dir(outdir)

    # Overview at final state
    fig = plot_overview(lat, suptitle="Scenario A: Inflation → Structure Formation")
    save_figure(fig, str(d / "A_overview_final.png"))

    # History
    fig = plot_history(lat.history, title="Scenario A: Coarse-grained History")
    save_figure(fig, str(d / "A_history.png"))

    # Phase diagram
    fig = plot_phase_diagram(lat.history, title="Scenario A: Phase-Space Trajectory")
    save_figure(fig, str(d / "A_phase_diagram.png"))

    # Snapshot strip
    fig = plot_snapshot_strip(lat, title="Scenario A: Field Snapshots")
    save_figure(fig, str(d / "A_snapshots.png"))

    # Animation
    if save_gif:
        anim = animate_field(lat, title="Scenario A: ED Inflation → Structure", interval=120)
        gif_path = str(d / "A_animation.gif")
        anim.save(gif_path, writer="pillow", fps=8)
        print(f"Saved: {gif_path}")

    if show:
        import matplotlib.pyplot as plt
        plt.show()

    return lat


# ---------------------------------------------------------------------------
# SCENARIO B – Two-Body Interaction
# ---------------------------------------------------------------------------

def run_scenario_B(
    size:     int  = 64,
    steps:    int  = 400,
    outdir:   str  = "output",
    save_gif: bool = False,
    show:     bool = True,
) -> EDLattice:
    """
    Scenario B: two competing high-ED Gaussian blobs.

    Studies:
      - Competition / merging of two overdense pockets.
      - Effect of the relational (−α·p^γ) drain on adjacent structures.
      - Emergence of a single dominant structure vs. mutual suppression.

    Uses:
        init_two_body()  →  standard ed_step  →  periodic BC
    """
    _print_header("SCENARIO B -- Two-Body Interaction")

    params = EDParams(
        alpha=0.05,
        beta=0.25,
        gamma=0.5,
        dt=0.1,
        boundary="periodic",
        mode="standard",
    )
    lat = EDLattice(rows=size, cols=size, params=params, seed=99)
    lat.init_two_body(separation=0.4, amplitude=0.9, width=0.08, background=0.1)

    print(f"  Initial state: {lat}")
    print(f"  Running {steps} steps ...")

    snap_interval = max(1, steps // 5)
    lat.run(
        steps=steps,
        record_every=4,
        snapshot_every=snap_interval,
        verbose=True,
        verbose_every=max(1, steps // 6),
    )
    print(f"  Final state:   {lat}")
    print(f"  Structures detected: {lat.structure_count()}")

    d = _ensure_dir(outdir)

    fig = plot_overview(lat, suptitle="Scenario B: Two-Body Interaction")
    save_figure(fig, str(d / "B_overview_final.png"))

    fig = plot_history(lat.history, title="Scenario B: Coarse-grained History")
    save_figure(fig, str(d / "B_history.png"))

    fig = plot_phase_diagram(lat.history, title="Scenario B: Phase-Space Trajectory")
    save_figure(fig, str(d / "B_phase_diagram.png"))

    fig = plot_snapshot_strip(lat, title="Scenario B: Two-Body Field Snapshots")
    save_figure(fig, str(d / "B_snapshots.png"))

    if save_gif:
        anim = animate_field(lat, title="Scenario B: Two-Body Interaction", interval=130)
        gif_path = str(d / "B_animation.gif")
        anim.save(gif_path, writer="pillow", fps=7)
        print(f"Saved: {gif_path}")

    if show:
        import matplotlib.pyplot as plt
        plt.show()

    return lat


# ---------------------------------------------------------------------------
# SCENARIO C – Black-Hole Saturation (mobility-weighted)
# ---------------------------------------------------------------------------

def run_scenario_C(
    size:         int   = 64,
    steps:        int   = 300,
    outdir:       str   = "output",
    save_gif:     bool  = False,
    show:         bool  = True,
    mobility_exp: float = 1.0,
) -> EDLattice:
    """
    Scenario C: saturated central patch (black-hole analogue).

    Demonstrates:
      - Horizon-like freezing of the saturated core (M(ρ) → 0 as ρ → ρ_max).
      - Absorbing boundary: edge cells fixed at p_min (ED horizon decoupling).
      - Mobility-weighted diffusion from ED-12.5 §1.

    Uses:
        make_black_hole_lattice()  →  ed_step_mobility  →  absorbing BC
    """
    _print_header("SCENARIO C -- Black-Hole Saturation (mobility-weighted)")

    print(f"  mobility_exp={mobility_exp}")
    lat = make_black_hole_lattice(
        size=size,
        radius=0.14,
        seed=11,
        alpha=0.02,
        beta=0.25,
        gamma=0.4,
        dt=0.05,
        mobility_exp=mobility_exp,
    )

    print(f"  Initial state: {lat}")
    print(f"  Running {steps} steps ...")

    snap_interval = max(1, steps // 5)
    lat.run(
        steps=steps,
        record_every=3,
        snapshot_every=snap_interval,
        verbose=True,
        verbose_every=max(1, steps // 6),
    )
    print(f"  Final state:   {lat}")

    d = _ensure_dir(outdir)

    fig = plot_overview(lat, suptitle="Scenario C: Black-Hole Saturation")
    save_figure(fig, str(d / "C_overview_final.png"))

    fig = plot_history(lat.history, title="Scenario C: Coarse-grained History")
    save_figure(fig, str(d / "C_history.png"))

    fig = plot_phase_diagram(lat.history, title="Scenario C: Phase-Space Trajectory")
    save_figure(fig, str(d / "C_phase_diagram.png"))

    fig = plot_snapshot_strip(lat, title="Scenario C: Black-Hole Field Snapshots")
    save_figure(fig, str(d / "C_snapshots.png"))

    if save_gif:
        anim = animate_field(lat, title="Scenario C: Black-Hole Saturation", interval=110)
        gif_path = str(d / "C_animation.gif")
        anim.save(gif_path, writer="pillow", fps=9)
        print(f"Saved: {gif_path}")

    if show:
        import matplotlib.pyplot as plt
        plt.show()

    return lat


# ---------------------------------------------------------------------------
# SCENARIO D – Noisy Universe  (Langevin / stochastic ED)
# ---------------------------------------------------------------------------

def run_scenario_D(
    size:         int   = 64,
    steps:        int   = 500,
    outdir:       str   = "output",
    save_gif:     bool  = False,
    show:         bool  = True,
    mobility_exp: float = 1.0,
    noise_amp:    float = 0.02,
) -> EDLattice:
    """
    Scenario D: stochastic Langevin ED — the noisy universe.

    Adds per-site Gaussian noise eta ~ N(0, noise_scale) to delta_p at every
    step, implementing the Langevin extension of the ED update rule:

        delta_p = beta * L  -  alpha * p^gamma  +  eta(x, t)

    Physical interpretation:
      - eta models thermal / quantum fluctuations in the ED field.
      - Noise continuously seeds new overdensities that compete via the
        relational term — producing a richer, more spatially diverse set
        of structures than the deterministic scenarios.
      - Directly analogous to the role of inflation-era quantum fluctuations
        in seeding large-scale cosmic structure.

    Starting from a disordered mid-ED state (init_random_noise), the noise
    sustains structure formation well into the late-thinning phase, preventing
    the clean heat-death seen in Scenarios A and B.

    Uses:
        init_random_noise(lo=0.3, hi=0.7)  ->  ed_step_mobility with noise_scale=noise_amp
        mobility mode  ->  periodic BC  ->  seed=77
        noise_amp (default 0.02) is overridden by --noise-amp CLI flag
    """
    _print_header("SCENARIO D -- Noisy Universe (Langevin / stochastic ED)")

    params = EDParams(
        alpha=0.03,
        beta=0.20,
        gamma=0.5,
        dt=0.05,
        boundary="periodic",
        mode="mobility",
        noise_scale=noise_amp,
        mobility_exp=mobility_exp,
    )
    lat = EDLattice(rows=size, cols=size, params=params, seed=77)
    lat.init_random_noise(lo=0.3, hi=0.7)

    print(f"  Initial state: {lat}")
    print(f"  noise_scale={params.noise_scale}  mobility_exp={params.mobility_exp}")
    print(f"  Running {steps} steps ...")

    snap_interval = max(1, steps // 5)
    lat.run(
        steps=steps,
        record_every=5,
        snapshot_every=snap_interval,
        verbose=True,
        verbose_every=max(1, steps // 6),
    )
    print(f"  Final state:   {lat}")
    print(f"  Structures detected: {lat.structure_count()}")

    d = _ensure_dir(outdir)

    fig = plot_overview(lat, suptitle="Scenario D: Noisy Universe")
    save_figure(fig, str(d / "D_overview_final.png"))

    fig = plot_history(lat.history, title="Scenario D: Coarse-grained History")
    save_figure(fig, str(d / "D_history.png"))

    fig = plot_phase_diagram(lat.history, title="Scenario D: Phase-Space Trajectory")
    save_figure(fig, str(d / "D_phase_diagram.png"))

    fig = plot_snapshot_strip(lat, title="Scenario D: Noisy Universe Snapshots")
    save_figure(fig, str(d / "D_snapshots.png"))

    if save_gif:
        anim = animate_field(lat, title="Scenario D: Noisy Universe", interval=110)
        gif_path = str(d / "D_animation.gif")
        anim.save(gif_path, writer="pillow", fps=9)
        print(f"Saved: {gif_path}")

    if show:
        import matplotlib.pyplot as plt
        plt.show()

    return lat


# ---------------------------------------------------------------------------
# Comparative summary
# ---------------------------------------------------------------------------

def run_comparative_summary(
    size:         int   = 64,
    steps:        int   = 400,
    outdir:       str   = "output",
    show:         bool  = True,
    mobility_exp: float = 1.0,
) -> None:
    """
    Run all four scenarios and print a side-by-side summary table.

    Useful for quickly comparing final-state observables across scenarios.
    """
    import matplotlib.pyplot as plt

    _print_header("COMPARATIVE SUMMARY -- All Four Scenarios")

    results: dict[str, dict] = {}

    for label, fn, kw in [
        ("A", run_scenario_A, dict(size=size, steps=steps, outdir=outdir, show=False)),
        ("B", run_scenario_B, dict(size=size, steps=steps, outdir=outdir, show=False)),
        ("C", run_scenario_C, dict(size=size, steps=min(steps, 300), outdir=outdir, show=False, mobility_exp=mobility_exp)),
        ("D", run_scenario_D, dict(size=size, steps=min(steps, 500), outdir=outdir, show=False, mobility_exp=mobility_exp)),
    ]:
        lat = fn(**kw)
        s   = lat.stats
        results[label] = {
            "steps":      lat.step_count,
            "p_hat":      s["p_hat"],
            "G":          s["G"],
            "L":          s["L"] if np.isfinite(s["L"]) else float("inf"),
            "phase":      s["phase"],
            "structures": lat.structure_count(),
        }

    # ---- Summary table ----
    print("\n" + "=" * 70)
    print(f"  {'':4s}  {'steps':>6s}  {'p_hat':>7s}  {'G':>8s}  {'L':>7s}  {'structures':>11s}  phase")
    print("-" * 70)
    for lbl, r in results.items():
        L_str = f"{r['L']:.2f}" if np.isfinite(r["L"]) else "  inf"
        print(
            f"  {lbl:4s}  {r['steps']:>6d}  {r['p_hat']:>7.4f}  "
            f"{r['G']:>8.5f}  {L_str:>7s}  {r['structures']:>11d}  {r['phase']}"
        )
    print("=" * 70 + "\n")

    if show:
        plt.show()


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="Run_Simulation.py",
        description="2D Event Density lattice simulation runner.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument(
        "--scenario", "-s",
        choices=["A", "B", "C", "D", "all", "compare"],
        default="all",
        help=(
            "Which scenario to run: "
            "A (inflation->structure), B (two-body), C (black-hole), "
            "D (noisy universe), "
            "all (A+B+C+D sequentially), compare (all + summary table). "
            "Default: all."
        ),
    )
    p.add_argument(
        "--size", type=int, default=64,
        help="Lattice side length in cells (square grid).  Default: 64.",
    )
    p.add_argument(
        "--steps", type=int, default=600,
        help="Number of time steps for Scenario A/B (C uses min(steps,300)).  Default: 600.",
    )
    p.add_argument(
        "--outdir", type=str, default="output",
        help="Directory for saved PNG / GIF files.  Default: output/.",
    )
    p.add_argument(
        "--mobility-exp", type=float, default=1.0,
        dest="mobility_exp",
        help=(
            "Exponent n for mobility M(rho) = ((rho_max-rho)/rho_max)^n "
            "used in Scenario C (mode=mobility) and Scenario D.  "
            "Default: 1.0 (linear, standard ED-12.5).  "
            "Larger values give a sharper horizon freeze near saturation."
        ),
    )
    p.add_argument(
        "--noise-amp", type=float, default=0.02,
        dest="noise_amp",
        help=(
            "Std-dev of per-site Langevin noise eta ~ N(0, noise_amp) added to "
            "delta_p at each step in Scenario D.  Default: 0.02.  "
            "Larger values seed structure more aggressively."
        ),
    )
    p.add_argument(
        "--gif", action="store_true",
        help="Also save an animated GIF for each scenario (requires Pillow).",
    )
    p.add_argument(
        "--no-show", dest="show", action="store_false",
        help="Do not call plt.show() (useful in headless / CI environments).",
    )
    p.set_defaults(show=True)
    return p.parse_args(argv)


def main(argv=None) -> None:
    args = _parse_args(argv)

    # Use a non-interactive backend when --no-show is requested
    if not args.show:
        matplotlib.use("Agg")

    # Shared keyword args for each scenario runner
    common = dict(
        size=args.size,
        outdir=args.outdir,
        save_gif=args.gif,
        show=args.show,
    )

    scenario = args.scenario.lower()

    if scenario in ("a", "all", "compare"):
        run_scenario_A(steps=args.steps, **common)

    if scenario in ("b", "all", "compare"):
        run_scenario_B(steps=args.steps, **common)

    if scenario in ("c", "all", "compare"):
        run_scenario_C(steps=min(args.steps, 300), mobility_exp=args.mobility_exp, **common)

    if scenario in ("d", "all", "compare"):
        run_scenario_D(steps=min(args.steps, 500), mobility_exp=args.mobility_exp,
                       noise_amp=args.noise_amp, **common)

    if scenario == "compare":
        run_comparative_summary(
            size=args.size,
            steps=args.steps,
            outdir=args.outdir,
            show=args.show,
            mobility_exp=args.mobility_exp,
        )


# ---------------------------------------------------------------------------
# Script entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()
