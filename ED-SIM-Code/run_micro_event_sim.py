"""
run_micro_event_sim.py
======================
Top-level CLI runner for the ED MicroEvent engine.  Analogous to
Run_Simulation.py (which runs 2D lattice scenarios), this script runs
ring-particle scenarios that validate the micro-event operator.

SCENARIOS  (from the ED Micro-Event Operator paper)
----------------------------------------------------
    H       : Golden inward collapse (canonical single ring, N=3)
    triads  : Three independent inward rings (compositionality)
    matrix  : 3x3 gamma-pair compositionality matrix
    scaling : Sweep N=3..12 at fixed d and gamma (scale-invariance)
    sweep   : Parameter sweep over (d, gamma_gate) → regime map

HOW TO RUN
----------
    python run_micro_event_sim.py                       # all scenarios
    python run_micro_event_sim.py --scenario H          # single scenario
    python run_micro_event_sim.py --scenario sweep      # regime map
    python run_micro_event_sim.py --no-show             # headless

OUTPUT
------
For each scenario, writes to --outdir (default: "me_output/"):
    <scenario>_events.json     MicroEvent records
    <scenario>_*.png           Diagnostic plots
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List

import matplotlib
import numpy as np

from event_lattice import (
    RingParams,
    RingState,
    make_inward_ring,
    make_outward_ring,
    make_tangent_ring,
    make_triad,
)
from micro_event_operator import (
    MicroEvent,
    detect_micro_event,
    save_micro_events,
)
from event_visualization import (
    plot_collapse_scatter,
    plot_decision_snapshot,
    plot_mechanism_pie,
    plot_observables,
    plot_regime_map,
    plot_ring,
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


def _print_event(event: MicroEvent, label: str = "") -> None:
    tag = f"[{label}] " if label else ""
    print(
        f"  {tag}mechanism={event.mechanism:<18s}  "
        f"chi_emp={event.chi_emp:.4f}  d_min={event.d_min_at_collapse:.6f}  "
        f"dR={event.delta_R_at_collapse:.6f}  N={event.N}  d={event.diameter}"
    )


# ---------------------------------------------------------------------------
# SCENARIO H — Golden Inward Collapse
# ---------------------------------------------------------------------------

def run_scenario_H(
    outdir: str = "me_output",
    show: bool = True,
) -> MicroEvent:
    """
    Canonical inward collapse: single ring, N=3, inward gamma gate.

    Expected: mechanism = inward-collapse, chi_emp > 0, delta_R < 0.
    """
    _print_header("SCENARIO H -- Golden Inward Collapse")

    ring = make_inward_ring(N=3, diameter=0.4, dt=0.001, max_steps=10_000)
    print(f"  Initial: {ring}")

    event = detect_micro_event(ring)
    _print_event(event)

    d = _ensure_dir(outdir)
    save_micro_events([event], str(d / "H_events.json"))

    fig = plot_ring(ring, title="Scenario H: Golden Inward Collapse (final state)")
    save_figure(fig, str(d / "H_ring_final.png"))

    fig = plot_observables(ring.history, title="Scenario H: Observables")
    save_figure(fig, str(d / "H_observables.png"))

    fig = plot_decision_snapshot(event, title="Scenario H: Decision Snapshot")
    save_figure(fig, str(d / "H_decision.png"))

    if show:
        import matplotlib.pyplot as plt
        plt.show()

    return event


# ---------------------------------------------------------------------------
# SCENARIO TRIADS — Compositionality
# ---------------------------------------------------------------------------

def run_scenario_triads(
    outdir: str = "me_output",
    show: bool = True,
) -> List[MicroEvent]:
    """
    Three independent inward rings (same params).

    Expected: 3 micro-events, all inward-collapse, identical chi_emp.
    """
    _print_header("SCENARIO TRIADS -- Compositionality")

    configs = [
        {"N": 3, "diameter": 0.3, "gamma_gate": "inward"},
        {"N": 3, "diameter": 0.3, "gamma_gate": "inward"},
        {"N": 3, "diameter": 0.3, "gamma_gate": "inward"},
    ]
    rings = make_triad(configs, dt=0.001, max_steps=10_000)

    events = []
    for i, ring in enumerate(rings):
        ev = detect_micro_event(ring)
        events.append(ev)
        _print_event(ev, label=f"Ring {i}")

    # Check compositionality: chi_emp should be identical
    chis = [e.chi_emp for e in events]
    spread = max(chis) - min(chis)
    print(f"  chi_emp spread: {spread:.6f}  {'PASS' if spread < 1e-10 else 'WARN'}")

    d = _ensure_dir(outdir)
    save_micro_events(events, str(d / "triads_events.json"))

    fig = plot_mechanism_pie(events, title="Triads: Mechanism Distribution")
    save_figure(fig, str(d / "triads_pie.png"))

    if show:
        import matplotlib.pyplot as plt
        plt.show()

    return events


# ---------------------------------------------------------------------------
# SCENARIO MATRIX — 3x3 Gamma-Pair Compositionality
# ---------------------------------------------------------------------------

def run_scenario_matrix(
    outdir: str = "me_output",
    show: bool = True,
) -> List[MicroEvent]:
    """
    3x3 matrix of gamma-gate combinations (inward, tangent, outward).

    Tests all single-ring gamma gates.  Expected:
        inward  → inward-collapse
        tangent → DECAY (or other-late)
        outward → outward-PBC (or PBC-corner)
    """
    _print_header("SCENARIO MATRIX -- 3x3 Gamma-Gate Combinations")

    gates = ["inward", "tangent", "outward"]
    events = []

    for gate in gates:
        params = RingParams(
            N=3, diameter=0.4, gamma_gate=gate,
            dt=0.001, max_steps=15_000,
        )
        ring = RingState(params).init_ring()
        ev = detect_micro_event(ring)
        events.append(ev)
        _print_event(ev, label=gate)

    d = _ensure_dir(outdir)
    save_micro_events(events, str(d / "matrix_events.json"))

    fig = plot_mechanism_pie(events, title="Matrix: Mechanism Distribution")
    save_figure(fig, str(d / "matrix_pie.png"))

    fig = plot_regime_map(events, x_key="diameter", y_key="gamma_gate",
                          title="Matrix: Regime Map")
    save_figure(fig, str(d / "matrix_regime.png"))

    if show:
        import matplotlib.pyplot as plt
        plt.show()

    return events


# ---------------------------------------------------------------------------
# SCENARIO SCALING — N=3..12
# ---------------------------------------------------------------------------

def run_scenario_scaling(
    outdir: str = "me_output",
    show: bool = True,
) -> List[MicroEvent]:
    """
    Scale-invariance test: sweep N=3..12 with fixed diameter and inward gate.

    Expected: mechanism = inward-collapse for all N, consistent chi_emp scaling.
    """
    _print_header("SCENARIO SCALING -- N=3..12 Sweep")

    events = []
    for N in range(3, 13):
        ring = make_inward_ring(N=N, diameter=0.4, dt=0.001, max_steps=10_000)
        ev = detect_micro_event(ring)
        events.append(ev)
        _print_event(ev, label=f"N={N}")

    d = _ensure_dir(outdir)
    save_micro_events(events, str(d / "scaling_events.json"))

    fig = plot_collapse_scatter(events, title="Scaling: chi_emp vs N (inward)")
    save_figure(fig, str(d / "scaling_scatter.png"))

    if show:
        import matplotlib.pyplot as plt
        plt.show()

    return events


# ---------------------------------------------------------------------------
# SCENARIO SWEEP — Parameter Regime Map
# ---------------------------------------------------------------------------

def run_scenario_sweep(
    outdir: str = "me_output",
    show: bool = True,
) -> List[MicroEvent]:
    """
    Parameter sweep over diameter and gamma_gate to reconstruct the ED manifold
    regime map.

    Sweeps d ∈ [0.1, 0.9] × gamma ∈ {inward, tangent, outward}.
    """
    _print_header("SCENARIO SWEEP -- Regime Map")

    diameters = np.linspace(0.1, 0.9, 9)
    gates = ["inward", "tangent", "outward"]
    events = []

    for gate in gates:
        for diam in diameters:
            ring = RingState(RingParams(
                N=3, diameter=float(diam), gamma_gate=gate,
                dt=0.001, max_steps=15_000,
            )).init_ring()
            ev = detect_micro_event(ring)
            events.append(ev)
            _print_event(ev, label=f"{gate} d={diam:.2f}")

    d = _ensure_dir(outdir)
    save_micro_events(events, str(d / "sweep_events.json"))

    fig = plot_regime_map(events, x_key="diameter", y_key="gamma_gate",
                          title="Sweep: Regime Map (d × gamma)")
    save_figure(fig, str(d / "sweep_regime.png"))

    fig = plot_mechanism_pie(events, title="Sweep: Mechanism Distribution")
    save_figure(fig, str(d / "sweep_pie.png"))

    fig = plot_collapse_scatter(events, title="Sweep: chi_emp vs Diameter")
    save_figure(fig, str(d / "sweep_scatter.png"))

    if show:
        import matplotlib.pyplot as plt
        plt.show()

    return events


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="run_micro_event_sim.py",
        description="ED MicroEvent engine — scenario runner.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--scenario", "-s",
        choices=["H", "triads", "matrix", "scaling", "sweep", "all"],
        default="all",
        help="Which scenario to run.  Default: all.",
    )
    p.add_argument(
        "--outdir", type=str, default="me_output",
        help="Output directory for JSON and PNG files.  Default: me_output/.",
    )
    p.add_argument(
        "--no-show", dest="show", action="store_false",
        help="Suppress plt.show() (headless mode).",
    )
    p.set_defaults(show=True)
    return p.parse_args(argv)


def main(argv=None) -> None:
    args = _parse_args(argv)

    if not args.show:
        matplotlib.use("Agg")

    common = dict(outdir=args.outdir, show=args.show)
    scenario = args.scenario.lower()

    if scenario in ("h", "all"):
        run_scenario_H(**common)

    if scenario in ("triads", "all"):
        run_scenario_triads(**common)

    if scenario in ("matrix", "all"):
        run_scenario_matrix(**common)

    if scenario in ("scaling", "all"):
        run_scenario_scaling(**common)

    if scenario in ("sweep", "all"):
        run_scenario_sweep(**common)


if __name__ == "__main__":
    main()
