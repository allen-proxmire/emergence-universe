"""
ed_analytic_law.py
==================
Analytic predictions for the ED MicroEvent operator.

This module encodes:
    - analytic_chi(N, d, gamma)
    - analytic_mechanism(N, d, gamma)

It is the single source of truth for analytic laws used by the
micro_event_operator.
"""

from __future__ import annotations
from typing import Literal


GammaGate = Literal["inward", "tangent", "outward"]


# ---------------------------------------------------------------------------
# Analytic chi law
# ---------------------------------------------------------------------------

def analytic_chi(N: int, d: float, gamma: GammaGate) -> float:
    """
    Analytic prediction for collapse time chi as a function of:
        - N      : number of particles
        - d      : ring diameter
        - gamma  : gamma gate ("inward", "tangent", "outward")

    Returns
    -------
    chi_pred : float
        Predicted collapse time.

    NOTE: This is currently a placeholder. Replace the body with the
    actual ED-Arch analytic law once specified.
    """
    # def analytic_chi(N: int, d: float, gamma: GammaGate) -> float:
    """
    Analytic prediction for collapse time chi.

    Test 1 hypothesis:
        For inward collapse: chi_pred = alpha * d
    """
    if gamma == "inward":
        # Inward-collapse analytic law with N-dependence:
        #   alpha(N) = 0.4701 + 0.0719 / N
        #   chi_inward(N, d) = alpha(N) * d
        alpha_N = 0.4701 + 0.0719 / N
        return alpha_N * d

    elif gamma == "tangent":
        # Placeholder until tangent law is derived
        return 10.0

    elif gamma == "outward":
        # Outward-PBC analytic law (provisional, N=3-based):
        #   chi_outward(d) = 0.5720 - 0.5 * d
        return 0.5720 - 0.5 * d

    else:
        raise ValueError(f"Unknown gamma gate: {gamma}")



# ---------------------------------------------------------------------------
# Analytic mechanism prediction
# ---------------------------------------------------------------------------

def analytic_mechanism(N: int, d: float, gamma: GammaGate) -> str:
    """
    Analytic mechanism prediction based on ED-Arch regime boundaries.
    """

    if gamma == "inward":
        return "inward-collapse"

    elif gamma == "tangent":
        return "DECAY"   # placeholder until tangent law is derived

    elif gamma == "outward":
        return "outward-PBC"   # placeholder until outward law is derived

    else:
        raise ValueError(f"Unknown gamma gate: {gamma}")


