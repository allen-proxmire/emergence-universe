"""
micro_event_operator.py
=======================
The diagnostic core of the ED MicroEvent engine: collapse detection,
micro-event snapshot capture, mechanism classification, and JSON
serialization.

This module implements the micro-event operator as formalized in
ED-Arch-00 (The ED Micro-Event Operator).  It is read-only with respect
to the simulation — it observes, never modifies.

MECHANISM LABELS
----------------
The operator assigns one of five labels via a deterministic decision tree:

    INWARD_COLLAPSE : Monotone contraction, dR < 0, no PBC interaction.
    OUTWARD_PBC     : Ring expanding (dR > 0), collapse via PBC re-entry.
    PBC_CORNER      : Particle near domain corner, asymmetric PBC collapse.
    OTHER_LATE      : Late collapse (chi > 3), no clear radial/boundary driver.
    DECAY           : No collapse within simulation horizon.

DECISION TREE  (Section 5.6)
-----------------------------
    1. No collapse within horizon  →  DECAY
    2. dR < 0                      →  INWARD_COLLAPSE
    3. dR > 0:
       a. near corner              →  PBC_CORNER
       b. near edge                →  OUTWARD_PBC
    4. otherwise                   →  OTHER_LATE

USAGE
-----
    from event_lattice import make_inward_ring
    from micro_event_operator import detect_micro_event

    ring = make_inward_ring(N=3, diameter=0.4)
    event = detect_micro_event(ring)
    print(event.mechanism, event.chi_emp)
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, asdict
from enum import Enum
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
from event_update import circumradius, d_min, delta_R, pbc_proximity
from event_lattice import RingState


# ---------------------------------------------------------------------------
# Mechanism enum
# ---------------------------------------------------------------------------

class Mechanism(str, Enum):
    """The five mechanism labels of the ED micro-event operator."""
    INWARD_COLLAPSE = "inward-collapse"
    OUTWARD_PBC     = "outward-PBC"
    PBC_CORNER      = "PBC-corner"
    OTHER_LATE      = "other-late"
    DECAY           = "DECAY"


# ---------------------------------------------------------------------------
# MicroEvent dataclass
# ---------------------------------------------------------------------------

@dataclass
class MicroEvent:
    """
    The atomic unit of ED dynamics — a complete geometric and dynamical
    snapshot at the moment of collapse (or at horizon for DECAY).

    Attributes
    ----------
    collapse_step          : Step at which collapse occurred (-1 for DECAY).
    chi_emp                : Empirical collapse time (simulation time at collapse).
    R_at_collapse          : Circumradius at collapse moment.
    K                      : Ring constant N/d.
    N                      : Number of particles.
    diameter               : Initial ring diameter.
    gamma_gate             : Active gamma gate ("inward", "tangent", "outward").
    mechanism              : Classified mechanism label.
    d_min_at_collapse      : d_min value at collapse.
    delta_R_at_collapse    : Circumradius trend at collapse.
    pbc_edge_dist          : Min edge distance at collapse.
    pbc_corner_dist        : Min corner distance at collapse.
    positions              : (N, 2) particle positions at collapse.
    velocities             : (N, 2) particle velocities at collapse.
    """
    collapse_step:       int
    chi_emp:             float
    R_at_collapse:       float
    K:                   float
    N:                   int
    diameter:            float
    gamma_gate:          str
    mechanism:           str
    d_min_at_collapse:   float
    delta_R_at_collapse: float
    pbc_edge_dist:       float
    pbc_corner_dist:     float
    positions:           list   # serializable as nested list
    velocities:          list   # serializable as nested list
    chi_pred:            float = 0.0
    mechanism_pred:      str   = ""
    chi_delta_pred_emp:  float = 0.0


# ---------------------------------------------------------------------------
# Classification decision tree
# ---------------------------------------------------------------------------

# Thresholds for classification
_CORNER_THR = 0.05   # near-corner if min corner dist < this fraction of box_size
_EDGE_THR   = 0.05   # near-edge if min edge dist < this fraction of box_size
_LATE_CHI   = 3.0    # "late" collapse threshold on chi_emp
_DR_ZERO    = 1e-6   # delta_R magnitude below this is treated as ~0


def classify_mechanism(
    collapsed: bool,
    delta_R_val: float,
    pbc_edge_dist: float,
    pbc_corner_dist: float,
    chi_emp: float,
    box_size: float,
) -> Mechanism:
    """
    Deterministic decision tree for mechanism classification.

    Implements Section 5.6 of the ED Micro-Event Operator paper:

        1. No collapse within horizon  →  DECAY
        2. dR < 0                      →  INWARD_COLLAPSE
        3. dR > 0:
           a. corner proximity         →  PBC_CORNER
           b. edge proximity           →  OUTWARD_PBC
        4. Otherwise                   →  OTHER_LATE

    Parameters
    ----------
    collapsed       : Whether the ring collapsed within the horizon.
    delta_R_val     : Circumradius trend at collapse (R_current - R_previous).
    pbc_edge_dist   : Min distance from any particle to nearest PBC edge.
    pbc_corner_dist : Min distance from any particle to nearest PBC corner.
    chi_emp         : Empirical collapse time.
    box_size        : PBC domain side length (for normalizing proximity).

    Returns
    -------
    Mechanism enum value.
    """
    # Step 1: no collapse → DECAY
    if not collapsed:
        return Mechanism.DECAY

    # Step 2: dR < 0 → inward collapse
    if delta_R_val < -_DR_ZERO:
        return Mechanism.INWARD_COLLAPSE

    # Step 3: dR > 0 → outward / PBC mechanisms
    if delta_R_val > _DR_ZERO:
        corner_thr = _CORNER_THR * box_size
        edge_thr = _EDGE_THR * box_size

        # 3a: near corner
        if pbc_corner_dist < corner_thr:
            return Mechanism.PBC_CORNER

        # 3b: near edge
        if pbc_edge_dist < edge_thr:
            return Mechanism.OUTWARD_PBC

    # Step 4: fallback → other-late
    return Mechanism.OTHER_LATE


# ---------------------------------------------------------------------------
# Main operator function
# ---------------------------------------------------------------------------

def detect_micro_event(ring: RingState) -> MicroEvent:
    """
    Run the micro-event operator on a ring.

    If the ring has not yet been evolved, this function evolves it to
    collapse or the simulation horizon.  It then captures the geometric
    snapshot and classifies the mechanism.

    Parameters
    ----------
    ring : An initialized RingState (via init_ring() or init_from_arrays()).

    Returns
    -------
    MicroEvent with all fields populated.
    """
    # Evolve if not already done
    if ring.step_count == 0 and not ring.collapsed:
        ring.run_until_collapse()

    # Capture snapshot
    obs = ring.current_observables()

    # Use cumulative delta_R (R_now - R_initial) for robust classification.
    # The single-step delta_R can be near-zero at the exact collapse moment
    # even for clearly inward or outward trajectories.
    cumulative_dR = ring.cumulative_delta_R

    mechanism = classify_mechanism(
        collapsed=ring.collapsed,
        delta_R_val=cumulative_dR,
        pbc_edge_dist=obs["pbc_edge_dist"],
        pbc_corner_dist=obs["pbc_corner_dist"],
        chi_emp=ring.time,
        box_size=ring.params.box_size,
    )


    # --- Analytic predictions ---
    from ed_analytic_law import analytic_chi, analytic_mechanism

    chi_pred = analytic_chi(
        N=ring.params.N,
        d=ring.params.diameter,
        gamma=ring.params.gamma_gate,
    )

    mechanism_pred = analytic_mechanism(
        N=ring.params.N,
        d=ring.params.diameter,
        gamma=ring.params.gamma_gate,
    )

    chi_delta_pred_emp = ring.time - chi_pred

    return MicroEvent(
        collapse_step=ring.collapse_step,
        chi_emp=ring.time,
        R_at_collapse=obs["R"],
        K=ring.params.K,
        N=ring.params.N,
        diameter=ring.params.diameter,
        gamma_gate=ring.params.gamma_gate,
        mechanism=mechanism.value,
        d_min_at_collapse=obs["d_min"],
        delta_R_at_collapse=cumulative_dR,
        pbc_edge_dist=obs["pbc_edge_dist"],
        pbc_corner_dist=obs["pbc_corner_dist"],
        positions=ring.positions.tolist(),
        velocities=ring.velocities.tolist(),
        chi_pred=chi_pred,
        mechanism_pred=mechanism_pred,
        chi_delta_pred_emp=chi_delta_pred_emp,
    )


# ---------------------------------------------------------------------------
# JSON serialization
# ---------------------------------------------------------------------------

def micro_event_to_dict(event: MicroEvent) -> dict:
    """Convert a MicroEvent to a JSON-serializable dict."""
    d = asdict(event)
    return d


def micro_event_from_dict(d: dict) -> MicroEvent:
    """Reconstruct a MicroEvent from a dict."""
    return MicroEvent(**d)


def save_micro_events(events: List[MicroEvent], path: str) -> None:
    """
    Save a list of MicroEvents to a JSON file.

    Parameters
    ----------
    events : List of MicroEvent objects.
    path   : Output file path.
    """
    data = [micro_event_to_dict(e) for e in events]
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w") as f:
        json.dump(data, f, indent=2)
    print(f"Saved {len(events)} micro-event(s) to {path}")


def load_micro_events(path: str) -> List[MicroEvent]:
    """
    Load a list of MicroEvents from a JSON file.

    Parameters
    ----------
    path : Input file path.

    Returns
    -------
    List of MicroEvent objects.
    """
    with open(path, "r") as f:
        data = json.load(f)
    return [micro_event_from_dict(d) for d in data]
