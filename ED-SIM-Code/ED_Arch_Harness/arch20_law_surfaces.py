"""
arch20_law_surfaces.py
======================
Analytic law surface predictions and manifold construction for the ED-Arch
series.  This module is the ground-truth oracle — it contains no simulation
code.  Every function computes a closed-form prediction from initial geometry
(d, ncl, γ) using the results of ED-Arch-20.

LAW SURFACES
------------
The ED manifold comprises three non-interpolating γ-indexed sheets:

    γ = +1  (inward-collapse)
        χ_N(d) = (d − MERGE_THR) / (T_DECAY · K)
        K = 2·sin(π / ncl)
        Collapse is strictly linear for all d > MERGE_THR.

    γ = −1  (outward-PBC)
        χ_N(d) = (d_far − d) / (T_DECAY · K)
        Bilateral symmetry: |∂χ/∂d|_{γ=−1} = ∂χ/∂d|_{γ=+1}

    γ = 0   (tangential saddle)
        Contains DECAY windows bounded by (d_lower, d_upper).
        Drift law: d_sw = √(d_init² + (t · DT)²)

ANGULAR-RATE LAW
----------------
    ω_N(d) = K / d        (angular rate)
    ω_N(d) · d = K        (invariant product)
    |dd/dt|  = K           (constant closing rate)

UNIT CONVENTIONS
----------------
The Arch papers use pixel units:
    d         in pixels
    MERGE_THR = 23.5 px
    T_DECAY   = 100 steps
    box       = 400 × 400 px

The MicroEvent engine uses continuous coordinates:
    box_size  = 1.0  (default)
    merge_thr = 0.01 (default)

The to_engine_params / from_engine_chi functions handle translation.

USAGE
-----
    from arch20_law_surfaces import chi_inward, closing_rate, deviation

    K   = closing_rate(ncl=6)
    chi = chi_inward(d=50.0, ncl=6)
    dev = deviation(chi_empirical=0.265, chi_predicted=chi)
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Allow imports from the parent ED Simulations directory
# ---------------------------------------------------------------------------
_PARENT = str(Path(__file__).resolve().parent.parent)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

from event_lattice import RingParams


# ---------------------------------------------------------------------------
# Constants (Arch pixel-scale conventions)
# ---------------------------------------------------------------------------

MERGE_THR_PX: float = 23.5     # universal merge threshold (pixels)
T_DECAY:      int   = 100      # normalization horizon (steps)
DT:           float = 1.0      # step duration (pixels / step)
BOX_PX:       int   = 400      # default domain side length (pixels)


# ═══════════════════════════════════════════════════════════════════════════
# Core rate law
# ═══════════════════════════════════════════════════════════════════════════

def closing_rate(ncl: int) -> float:
    """
    Constant closing rate for a regular ncl-gon.

    K = 2 · sin(π / ncl)

    This is the central geometric quantity.  It appears in the denominator
    of every collapse-time prediction and defines the angular-rate invariant.

    Parameters
    ----------
    ncl : Number of cores participating in the polygonal ring.
          RING topology: ncl = N.
          CSAT topology: ncl = N − 1.

    Returns
    -------
    K : Closing rate (pixels / step).
    """
    return 2.0 * math.sin(math.pi / ncl)


def effective_radius(d: float, ncl: int) -> float:
    """
    Effective radius of the collapse polygon.

    R_eff = d / (2 · sin(π / ncl)) = d / K

    For a regular ncl-gon with adjacent pairwise separation d, this is the
    circumradius of the polygon.

    Parameters
    ----------
    d   : Initial minimum pairwise separation (pixels).
    ncl : Number of cores in the polygon.

    Returns
    -------
    R_eff : Effective radius (pixels).
    """
    return d / closing_rate(ncl)


def angular_rate(d: float, ncl: int) -> float:
    """
    Angular rate of collapse.

    ω_N(d) = K / d

    Invariant: ω_N(d) · d = K for all d and ncl.

    Parameters
    ----------
    d   : Pairwise separation (pixels).
    ncl : Number of cores.

    Returns
    -------
    omega : Angular rate (radians / step).
    """
    return closing_rate(ncl) / d


# ═══════════════════════════════════════════════════════════════════════════
# γ = +1 sheet (inward collapse)
# ═══════════════════════════════════════════════════════════════════════════

def chi_inward(d: float, ncl: int) -> float:
    """
    Normalized inward collapse time.

    χ_N(d) = (d − MERGE_THR) / (T_DECAY · K)

    Strictly linear in d for all d > MERGE_THR.  The denominator
    T_DECAY · K = 100 · 2·sin(π/ncl) is the effective length scale
    that unifies all ncl values onto a single rate law.

    Parameters
    ----------
    d   : Initial pairwise separation (pixels), must be > MERGE_THR.
    ncl : Number of cores.

    Returns
    -------
    chi_N : Normalized collapse time (dimensionless).
    """
    K = closing_rate(ncl)
    return (d - MERGE_THR_PX) / (T_DECAY * K)


def chi_inward_steps(d: float, ncl: int) -> float:
    """
    Raw step count for inward collapse.

    steps = χ_N · T_DECAY = (d − MERGE_THR) / K

    Parameters
    ----------
    d   : Initial pairwise separation (pixels).
    ncl : Number of cores.

    Returns
    -------
    steps : Collapse step count.
    """
    return chi_inward(d, ncl) * T_DECAY


# ═══════════════════════════════════════════════════════════════════════════
# γ = −1 sheet (outward-PBC)
# ═══════════════════════════════════════════════════════════════════════════

def chi_outward(d: float, ncl: int, d_far: float) -> float:
    """
    Normalized outward-PBC collapse time.

    χ_N(d) = (d_far − d) / (T_DECAY · K)

    Particles fly outward and PBC-wrap; d_far is the configuration-dependent
    distance at which the wrap produces a merge event.

    Bilateral symmetry with inward sheet:
        |∂χ/∂d|_{γ=−1} = ∂χ/∂d|_{γ=+1} = 1 / (T_DECAY · K)

    The slopes are equal in magnitude, opposite in sign — the outward sheet
    is a mirror of the inward sheet reflected about d_far.

    Parameters
    ----------
    d     : Initial pairwise separation (pixels).
    ncl   : Number of cores.
    d_far : Configuration-dependent far-boundary distance (pixels).
            For RING N=6: d_far ≈ BOX_PX / 2 = 200.

    Returns
    -------
    chi_N : Normalized collapse time.
    """
    K = closing_rate(ncl)
    return (d_far - d) / (T_DECAY * K)


def chi_zero_d(ncl: int, d_far: float) -> float:
    """
    Diameter at which outward collapse time reaches zero.

    Setting χ_N = 0 in the outward formula:
        0 = (d_far − d) / (T_DECAY · K)  →  d = d_far

    Parameters
    ----------
    ncl   : Number of cores.
    d_far : Far-boundary distance (pixels).

    Returns
    -------
    d_zero : Critical diameter (pixels) where χ → 0.
    """
    return d_far


# ═══════════════════════════════════════════════════════════════════════════
# γ = 0 sheet (tangential saddle)
# ═══════════════════════════════════════════════════════════════════════════

# Known DECAY window boundaries from Arch-19 Section 7 and Arch-20 Section 4.
# Keys are (ncl, config_type) → (d_lower, d_upper) in pixels.
# np.inf encodes an open-ended upper boundary (axis-immunity).
_DECAY_WINDOWS = {
    (6, "RING"): (39.0, 61.0),
    (5, "CSAT"): (MERGE_THR_PX, 31.0),
    (3, "CSAT"): (MERGE_THR_PX, 61.0),
    (4, "RING"): (25.5, np.inf),
}


def decay_window(ncl: int, config_type: str = "RING") -> Tuple[float, float]:
    """
    DECAY window boundaries for the tangential sheet.

    Within (d_lower, d_upper), γ = 0 trajectories do not collapse within
    the simulation horizon — the mechanism is DECAY.

    Outside this window:
        d < d_lower : collapse via inward-spiral (lower lobe)
        d > d_upper : collapse via PBC-corner or satellite-PBC (upper lobe)

    Known boundaries are taken from Arch-19 Section 7 (pixel-precise empirical
    measurements).  For configurations not in the lookup table, a heuristic
    estimate is provided based on the closing rate.

    Parameters
    ----------
    ncl         : Number of cores.
    config_type : "RING" or "CSAT".  Topology affects window bounds.

    Returns
    -------
    (d_lower, d_upper) : Window boundaries in pixels.
    """
    key = (ncl, config_type.upper())
    if key in _DECAY_WINDOWS:
        return _DECAY_WINDOWS[key]

    # Heuristic for unknown configurations:
    # Lower bound: near MERGE_THR for CSAT, somewhat above for RING
    # Upper bound: scales with box geometry and K
    K = closing_rate(ncl)
    if config_type.upper() == "CSAT":
        d_lower = MERGE_THR_PX
    else:
        d_lower = MERGE_THR_PX + 2.0 / K
    d_upper = MERGE_THR_PX + T_DECAY * K * 0.4
    return (d_lower, d_upper)


def pythagorean_drift(d_init: float, t: float, dt: float = DT) -> float:
    """
    Effective diameter after tangential ballistic drift under γ = 0.

    d_sw = √(d_init² + (t · dt)²)

    Under tangential motion, particles drift perpendicular to the radial
    direction.  The pairwise separation follows a Pythagorean law because
    the radial component is constant (d_init) while the tangential
    displacement grows linearly (t · dt).

    Parameters
    ----------
    d_init : Initial pairwise separation (pixels).
    t      : Elapsed time (steps).
    dt     : Step duration (pixels / step).

    Returns
    -------
    d_sw : Effective diameter at time t (pixels).
    """
    return math.sqrt(d_init ** 2 + (t * dt) ** 2)


# ═══════════════════════════════════════════════════════════════════════════
# Manifold sampling
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class LawSurfacePoint:
    """One point on the predicted law surface."""
    d:                  float
    ncl:                int
    gamma:              int       # +1, 0, or -1
    chi_predicted:      float
    mechanism_predicted: str


def sample_law_surface(
    gamma: int,
    ncl_range: List[int],
    d_range: np.ndarray,
    d_far: Optional[float] = None,
) -> List[LawSurfacePoint]:
    """
    Sample the predicted law surface for a given γ-sheet.

    Sweeps (d, ncl) and returns predicted χ and mechanism at each point.

    For γ = +1:
        d > MERGE_THR → inward-collapse with χ = chi_inward(d, ncl)
        d ≤ MERGE_THR → skipped (pre-merged, no meaningful collapse)

    For γ = −1:
        All d → outward-PBC with χ = chi_outward(d, ncl, d_far)
        Requires d_far.

    For γ = 0:
        d in DECAY window → DECAY (no collapse, χ = inf)
        d < d_lower        → inward-collapse (lower lobe)
        d > d_upper        → PBC-corner (upper lobe)

    Parameters
    ----------
    gamma     : Sheet selector: +1, 0, or −1.
    ncl_range : List of ncl values to sweep.
    d_range   : Array of d values to sweep (pixels).
    d_far     : Far-boundary distance for γ = −1 (required for outward sheet).

    Returns
    -------
    points : List of LawSurfacePoint.
    """
    if gamma == -1 and d_far is None:
        raise ValueError("d_far is required for the outward (γ = −1) sheet.")

    points: List[LawSurfacePoint] = []

    for ncl in ncl_range:
        for d in d_range:
            d_val = float(d)

            if gamma == 1:
                if d_val <= MERGE_THR_PX:
                    continue
                chi = chi_inward(d_val, ncl)
                mech = "inward-collapse"

            elif gamma == -1:
                chi = chi_outward(d_val, ncl, d_far)
                mech = "outward-PBC"

            elif gamma == 0:
                d_lo, d_hi = decay_window(ncl)
                if d_lo <= d_val <= d_hi:
                    chi = np.inf
                    mech = "DECAY"
                elif d_val < d_lo:
                    # Lower lobe: estimate chi from adjusted inward formula
                    # using the effective remaining distance to MERGE_THR
                    chi = chi_inward(d_val, ncl)
                    mech = "inward-collapse"
                else:
                    # Upper lobe: PBC-corner collapse at late times
                    # Use Pythagorean drift to estimate when d_sw crosses
                    # the upper boundary and PBC interaction begins
                    K = closing_rate(ncl)
                    chi = (d_val - d_hi) / (T_DECAY * K) + 3.0
                    mech = "PBC-corner"
            else:
                raise ValueError(f"gamma must be +1, 0, or −1, got {gamma}")

            points.append(LawSurfacePoint(
                d=d_val,
                ncl=ncl,
                gamma=gamma,
                chi_predicted=chi,
                mechanism_predicted=mech,
            ))

    return points


# ═══════════════════════════════════════════════════════════════════════════
# Deviation metric
# ═══════════════════════════════════════════════════════════════════════════

def deviation(chi_empirical: float, chi_predicted: float) -> float:
    """
    Central comparison metric: dev = χ_emp − χ_pred.

    Pass criterion for most invariants: dev = 0.

    Parameters
    ----------
    chi_empirical : Measured collapse time from MicroEvent engine.
    chi_predicted : Analytic prediction from law surface.

    Returns
    -------
    dev : Signed deviation.
    """
    return chi_empirical - chi_predicted


# ═══════════════════════════════════════════════════════════════════════════
# Unit conversion: Arch pixel-scale ↔ MicroEvent engine
# ═══════════════════════════════════════════════════════════════════════════

_GAMMA_MAP = {+1: "inward", 0: "tangent", -1: "outward"}


def to_engine_params(
    d_px: float,
    ncl: int,
    gamma: int,
    box_px: int = BOX_PX,
    config_type: str = "RING",
) -> RingParams:
    """
    Convert Arch pixel-scale parameters to a MicroEvent engine RingParams.

    Scaling convention:
        scale = 1.0 / box_px          (engine units per pixel)
        diameter_engine = d_px / (sin(π/ncl) · box_px)
        merge_thr_engine = MERGE_THR_PX / box_px
        dt_engine = DT / box_px        (so displacement/step matches)

    The engine's ``diameter`` parameter is the ring diameter (circumscribed
    circle diameter), not the pairwise separation.  For a regular ncl-gon:
        d_pairwise = diameter · sin(π / ncl)
    so:
        diameter = d_pairwise / sin(π / ncl)

    Parameters
    ----------
    d_px        : Initial pairwise separation in pixels.
    ncl         : Number of cores.
    gamma       : Sheet selector (+1, 0, −1).
    box_px      : Domain side length in pixels.
    config_type : "RING" (N = ncl) or "CSAT" (N = ncl + 1).

    Returns
    -------
    RingParams ready for RingState construction.
    """
    scale = 1.0 / box_px

    # Engine diameter from pairwise separation
    sin_pi_n = math.sin(math.pi / ncl)
    diameter_engine = (d_px / sin_pi_n) * scale

    # Particle count depends on topology
    N = ncl if config_type.upper() == "RING" else ncl + 1

    # Map gamma integer to engine gate string
    gamma_gate = _GAMMA_MAP[gamma]

    return RingParams(
        N=N,
        diameter=diameter_engine,
        gamma_gate=gamma_gate,
        box_size=1.0,
        dt=DT * scale,
        merge_thr=MERGE_THR_PX * scale,
        speed=1.0,
        max_steps=10_000,
    )


def from_engine_chi(
    chi_engine: float,
    box_px: int = BOX_PX,
) -> float:
    """
    Convert MicroEvent engine collapse time to Arch normalized χ_N.

    The engine reports collapse time as ring.time = step_count · dt_engine.
    With dt_engine = DT / box_px, step_count = chi_engine · box_px / DT.
    Normalizing: χ_N = step_count / T_DECAY = chi_engine · box_px / (DT · T_DECAY).

    Parameters
    ----------
    chi_engine : Collapse time from engine (ring.time at collapse).
    box_px     : Domain side length used for scaling.

    Returns
    -------
    chi_N : Normalized collapse time in Arch units.
    """
    return chi_engine * box_px / (DT * T_DECAY)
