"""
ED_Update_Rule.py
=================
Minimal 2D lattice update rule derived from the ED Compositional Rule.

SOURCE DERIVATION
-----------------
The update rule is extracted from the cosmological specialization of the
ED compositional rule (ED-05.5, §3.1 and §3.4):

    p(A∪B) = p(A) + p(B)
             − α ∫_{A∩B}   p^γ        dμ       [relational penalty]
             − β ∫_{A∪B}   |∇p|²      dμ       [gradient penalty]
             − ζ ∫_{∂(A∪B)} h(|∇p|)  dS        [boundary term]

with functional choices for the cosmological/intermediate regime:
    f(p) = p^γ,  0 < γ < 1   (concave; drives competition + stabilization)
    g(u) = u²                 (quadratic; drives exponential smoothing)
    h(u) = u / (1 + u)        (saturating; finite-cost horizons)

Taking the gradient flow of the cost functional implied by the
compositional rule (i.e., ∂p/∂t = −δH/δp) yields the continuous PDE:

    ∂p/∂t = β · ∇²p  −  α · p^γ

This is the leading-order equation for the intermediate regime:
  * ED not near saturation (M(ρ) ≈ 1, so nonlinear mobility is inactive)
  * Gradients small but nonzero (boundary term sub-dominant in bulk)

The ED Analogues paper (ED-12.5) gives the general diffusion law as:
    ∂_t ρ = ∇·(M(ρ) ∇ρ)
with M(ρ) → 0 as ρ → ρ_max for the saturation (black-hole) regime.
The mobility-weighted variant below recovers this.

2D LATTICE DISCRETIZATION
--------------------------
State:   p[i, j] ∈ [p_min, p_max]  for a grid of shape (N, N)

The 5-point Laplacian (von Neumann neighbourhood):
    L[i,j] = p[i+1,j] + p[i-1,j] + p[i,j+1] + p[i,j-1] − 4·p[i,j]

THE MINIMAL UPDATE RULE  (Standard regime):
    Δp[i,j] = β · L[i,j]  −  α · p[i,j]^γ
    p_new[i,j] = clamp(p[i,j] + dt · Δp[i,j], p_min, p_max)

THE MOBILITY-WEIGHTED UPDATE RULE  (Saturation-aware regime):
    Flux from neighbour k → (i,j):
        J_k = M( (p[i,j] + p[k]) / 2 ) · (p[k] − p[i,j])
    where  M(ρ) = ((ρ_max − ρ) / ρ_max)^n   (n = mobility_exp; stiffens as ρ → ρ_max)

    Δp[i,j] = Σ_{k ∈ neighbours} J_k  −  α · p[i,j]^γ
    p_new[i,j] = clamp(p[i,j] + dt · Δp[i,j], p_min, p_max)

STOCHASTIC (LANGEVIN) EXTENSION:
    Both update rules accept an optional noise term drawn from N(0, noise_scale):
        Δp[i,j] += η[i,j],   η ~ N(0, noise_scale)
    This models thermal / quantum fluctuations in the ED field and is the basis
    for Scenario D (noisy universe).  Pass noise_scale=0.0 (default) to recover
    the fully deterministic dynamics.

BOUNDARY CONDITIONS
-------------------
Three options are supported:
  "periodic"   : wrap-around (eliminates boundary term; useful for bulk dynamics)
  "absorbing"  : edge cells fixed at p_min (models horizon-like decoupling surface)
  "reflecting" : zero-flux Neumann condition (closed system, conserves total ED)

PHYSICAL INTERPRETATION OF EACH TERM
--------------------------------------
  β · ∇²p          Gradient penalty → diffusion → inflation-like smoothing.
                    High-gradient regions decay; homogeneity scale L(t) ~ 1/G(t) grows.

  −α · p^γ          Relational competition → concave drain.
                    At moderate p: drain is strong → suppresses overlap.
                    At high p:    drain saturates (p^γ grows sublinearly) →
                                  overdense pockets persist as metastable structures
                                  (ED analogue of gravitational stabilization).

  Mobility M(ρ)     Stiffness near ρ_max → horizon-like freezing.
                    Reproduces the black-hole saturation case from ED-12.5.

PARAMETER GUIDE
---------------
  alpha  (α)   Relational penalty strength.    Default: 0.05
               Larger α → faster competition; stronger structure suppression.

  beta   (β)   Gradient smoothing strength.    Default: 0.25
               Larger β → faster diffusion; quicker inflation-like phase.

  gamma  (γ)   Concavity exponent, 0 < γ < 1. Default: 0.5
               Closer to 0 → stronger saturation at high ED (more stable pockets).
               Closer to 1 → more linear drain (weaker structure).

  dt           Time step.                      Default: 0.1
               Must satisfy dt · 4β ≤ 1 for numerical stability (CFL-like).

  p_min        Floor on ED density.            Default: 0.01
  p_max        Ceiling on ED density.          Default: 1.0

EXPECTED DYNAMICAL PHASES  (maps to ED-05.5 cosmological history)
------------------------------------------------------------------
  Phase 1 – ED Inflation:
    High p everywhere, small gradients. β term dominates.
    G(t) = ⟨|∇p|⟩ decays exponentially; homogeneity scale grows.

  Phase 2 – Residual Gradients / Structure Seeds:
    Smoothing weakens as G → 0. Relational term (−α·p^γ) takes over.
    Small δp variations amplify; positive feedback loop begins.

  Phase 3 – Structure Formation:
    High-ED pockets stabilize via concavity saturation.
    Low-ED background thins: dp̂/dt ≈ −λ · G(t)².

  Phase 4 – Late-Time Thinning (Heat Death):
    Gradients vanish; structures dissolve; boundary term dominates.
    p̂ → p_min;  G → 0.

USAGE EXAMPLE
-------------
    import numpy as np
    from ED_Update_Rule import ed_step, ed_step_mobility

    p = np.random.uniform(0.8, 1.0, size=(64, 64))   # high-ED initial state
    for t in range(1000):
        p = ed_step(p, alpha=0.05, beta=0.25, gamma=0.5, dt=0.1)
"""

import numpy as np


# ---------------------------------------------------------------------------
# Core helper: Laplacian with boundary handling
# ---------------------------------------------------------------------------

def _laplacian(p: np.ndarray, boundary: str = "periodic") -> np.ndarray:
    """
    5-point discrete Laplacian on a 2D array.

    Parameters
    ----------
    p        : 2D array of ED density values, shape (N, M).
    boundary : "periodic" | "absorbing" | "reflecting"

    Returns
    -------
    L : 2D array of the same shape as p.
    """
    if boundary == "periodic":
        L = (
            np.roll(p, -1, axis=0) + np.roll(p, +1, axis=0)
            + np.roll(p, -1, axis=1) + np.roll(p, +1, axis=1)
            - 4.0 * p
        )
    elif boundary == "absorbing":
        p_pad = np.pad(p, 1, mode="constant", constant_values=0.0)
        L = (
            p_pad[2:, 1:-1] + p_pad[:-2, 1:-1]
            + p_pad[1:-1, 2:] + p_pad[1:-1, :-2]
            - 4.0 * p
        )
    elif boundary == "reflecting":
        p_pad = np.pad(p, 1, mode="edge")
        L = (
            p_pad[2:, 1:-1] + p_pad[:-2, 1:-1]
            + p_pad[1:-1, 2:] + p_pad[1:-1, :-2]
            - 4.0 * p
        )
    else:
        raise ValueError(f"Unknown boundary condition: '{boundary}'")
    return L


# ---------------------------------------------------------------------------
# Relational penalty term
# ---------------------------------------------------------------------------

def _relational(p: np.ndarray, gamma: float) -> np.ndarray:
    """
    Local relational competition term: p^γ.

    The concavity (0 < γ < 1) ensures:
      - Small overdensities are amplified (structure seeds survive).
      - Large overdensities saturate and stabilize (metastable pockets).
    """
    return p ** gamma


# ---------------------------------------------------------------------------
# Mobility function M(rho) for the saturation-aware variant
# ---------------------------------------------------------------------------

def _mobility(p: np.ndarray, p_max: float, mobility_exp: float = 1.0) -> np.ndarray:
    """
    Mobility M(rho) = ((rho_max - rho) / rho_max) ** mobility_exp.

    mobility_exp = 1.0 (default): linear falloff — standard ED-12.5 case.
    mobility_exp > 1.0: faster freeze near saturation (sharper horizon).
    mobility_exp < 1.0: slower freeze — mobility persists longer at high rho.

    Approaches 1 at low ED (free diffusion) and 0 at rho_max (frozen horizon).
    """
    return np.clip((p_max - p) / p_max, 0.0, 1.0) ** mobility_exp


# ---------------------------------------------------------------------------
# STANDARD UPDATE STEP  (minimal / default)
# ---------------------------------------------------------------------------

def ed_step(
    p: np.ndarray,
    alpha:       float = 0.05,
    beta:        float = 0.25,
    gamma:       float = 0.5,
    dt:          float = 0.1,
    p_min:       float = 0.01,
    p_max:       float = 1.0,
    boundary:    str   = "periodic",
    noise_scale: float = 0.0,
    rng: "np.random.Generator | None" = None,
) -> np.ndarray:
    """
    One time step of the minimal 2D ED update rule.

    Implements:
        delta_p[i,j] = beta * L[i,j]  -  alpha * p[i,j]^gamma
                       [+ N(0, noise_scale)  if noise_scale > 0]
        p_new        = clamp(p + dt * delta_p, p_min, p_max)

    Parameters
    ----------
    p           : Current ED density field, 2D array.
    alpha       : Relational penalty coefficient (competition strength).
    beta        : Gradient penalty coefficient (diffusion strength).
    gamma       : Concavity exponent for relational term, 0 < gamma < 1.
    dt          : Time step size.
    p_min       : Minimum allowed density (floor).
    p_max       : Maximum allowed density (ceiling / saturation).
    boundary    : Boundary condition: "periodic" | "absorbing" | "reflecting".
    noise_scale : Std dev of per-site Gaussian noise added to delta_p each step.
                  0.0 (default) gives fully deterministic dynamics.
                  > 0 gives Langevin / stochastic ED (Scenario D: noisy universe).
    rng         : numpy Generator for reproducible noise.  If None and
                  noise_scale > 0, a fresh Generator is used each call.

    Returns
    -------
    p_new : Updated ED density field, same shape as p.

    Notes
    -----
    Stability condition (CFL-like): dt * 4 * beta <= 1.
    For defaults (beta=0.25, dt=0.1): 4 * 0.25 * 0.1 = 0.1  -- stable.
    """
    L  = _laplacian(p, boundary=boundary)
    R  = _relational(p, gamma)

    delta_p = beta * L - alpha * R

    if noise_scale > 0.0:
        _rng = rng if rng is not None else np.random.default_rng()
        delta_p = delta_p + _rng.normal(0.0, noise_scale, p.shape)

    return np.clip(p + dt * delta_p, p_min, p_max)


# ---------------------------------------------------------------------------
# MOBILITY-WEIGHTED UPDATE STEP  (saturation-aware)
# ---------------------------------------------------------------------------

def ed_step_mobility(
    p: np.ndarray,
    alpha:        float = 0.05,
    gamma:        float = 0.5,
    dt:           float = 0.1,
    p_min:        float = 0.01,
    p_max:        float = 1.0,
    boundary:     str   = "periodic",
    mobility_exp: float = 1.0,
    noise_scale:  float = 0.0,
    rng: "np.random.Generator | None" = None,
) -> np.ndarray:
    """
    One time step using the mobility-weighted diffusion form.

    Implements the ED diffusion law from ED-12.5:
        d_t rho = grad·(M(rho) grad_rho)  -  alpha * rho^gamma

    Discretized with edge-centred mobility:
        J_k = M( (p[i,j] + p_k) / 2 ) * (p_k - p[i,j])
        delta_p = sum_k J_k  -  alpha * p^gamma

    Mobility: M(rho) = ((rho_max - rho) / rho_max) ** mobility_exp
      mobility_exp=1.0 (default): linear — standard ED-12.5 case.
      mobility_exp>1.0: sharper horizon freeze near saturation.

    This naturally recovers:
      - Free diffusion at low ED  (M ~ 1)
      - Frozen / horizon-like behaviour at ED saturation  (M -> 0)

    Parameters
    ----------
    p            : Current ED density field, 2D array.
    alpha        : Relational penalty coefficient.
    gamma        : Concavity exponent, 0 < gamma < 1.
    dt           : Time step size.
    p_min        : Density floor.
    p_max        : Density ceiling (saturation point rho_max).
    boundary     : "periodic" | "absorbing" | "reflecting".
    mobility_exp : Exponent n in M(rho) = ((rho_max-rho)/rho_max)^n.
                   Default 1.0 (linear).  Increase for sharper horizon freezing.
    noise_scale  : Std dev of per-site Gaussian noise added to delta_p.
                   0.0 (default) = deterministic.
    rng          : numpy Generator for reproducible noise.

    Returns
    -------
    p_new : Updated ED density field.
    """
    if boundary == "periodic":
        p_up    = np.roll(p, -1, axis=0)
        p_down  = np.roll(p, +1, axis=0)
        p_right = np.roll(p, -1, axis=1)
        p_left  = np.roll(p, +1, axis=1)
    elif boundary in ("absorbing", "reflecting"):
        mode = "constant" if boundary == "absorbing" else "edge"
        kw   = {"constant_values": 0.0} if boundary == "absorbing" else {}
        p_pad   = np.pad(p, 1, mode=mode, **kw)
        p_up    = p_pad[2:,  1:-1]
        p_down  = p_pad[:-2, 1:-1]
        p_right = p_pad[1:-1, 2:]
        p_left  = p_pad[1:-1, :-2]
    else:
        raise ValueError(f"Unknown boundary condition: '{boundary}'")

    # Edge-centred mobility: M evaluated at the average of the two endpoints
    flux_up    = _mobility(0.5 * (p + p_up),    p_max, mobility_exp) * (p_up    - p)
    flux_down  = _mobility(0.5 * (p + p_down),  p_max, mobility_exp) * (p_down  - p)
    flux_right = _mobility(0.5 * (p + p_right), p_max, mobility_exp) * (p_right - p)
    flux_left  = _mobility(0.5 * (p + p_left),  p_max, mobility_exp) * (p_left  - p)

    diffusion  = flux_up + flux_down + flux_right + flux_left
    relational = _relational(p, gamma)

    delta_p = diffusion - alpha * relational

    if noise_scale > 0.0:
        _rng = rng if rng is not None else np.random.default_rng()
        delta_p = delta_p + _rng.normal(0.0, noise_scale, p.shape)

    return np.clip(p + dt * delta_p, p_min, p_max)


# ---------------------------------------------------------------------------
# Diagnostic helpers
# ---------------------------------------------------------------------------

def gradient_magnitude(p: np.ndarray, boundary: str = "periodic") -> np.ndarray:
    """
    Compute |grad p| at every lattice site using central differences.

    G(t) = mean(|grad p|) is the coarse-grained gradient magnitude used in
    ED-05.5 to track the inflation / structure-formation / thinning phases.
    """
    if boundary == "periodic":
        dp_dx = (np.roll(p, -1, axis=1) - np.roll(p, +1, axis=1)) / 2.0
        dp_dy = (np.roll(p, -1, axis=0) - np.roll(p, +1, axis=0)) / 2.0
    else:
        dp_dx = np.gradient(p, axis=1)
        dp_dy = np.gradient(p, axis=0)
    return np.sqrt(dp_dx**2 + dp_dy**2)


def coarse_grained_stats(p: np.ndarray, boundary: str = "periodic") -> dict:
    """
    Return the coarse-grained observables tracked in ED-05.5:

        p_hat : average ED density  p_hat(t)  = mean(p)
        G     : mean gradient magnitude  G(t) = mean(|grad p|)
        L     : homogeneity scale proxy  L(t) ~ 1 / G(t)
        a     : ED scale factor proxy    a(t) proportional to L(t)

    Returns
    -------
    dict with keys: "p_hat", "G", "L", "a"
    """
    p_hat = float(np.mean(p))
    G     = float(np.mean(gradient_magnitude(p, boundary=boundary)))
    L     = 1.0 / G if G > 1e-12 else float("inf")
    return {"p_hat": p_hat, "G": G, "L": L, "a": L}
