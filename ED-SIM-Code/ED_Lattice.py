"""
ED_Lattice.py
=============
2D Event Density lattice: grid creation, initialization, and time evolution.

Wraps the update rules in ED_Update_Rule.py inside an EDLattice class that
manages the full lifecycle of a simulation:

    grid creation  →  initial condition  →  time evolution  →  history

QUICK START
-----------
    from ED_Lattice import EDLattice

    # Create a 64×64 lattice with default cosmological parameters
    lattice = EDLattice(rows=64, cols=64)

    # Initialize in the "inflation" regime: high-ED, nearly uniform
    lattice.init_uniform_high(noise_scale=0.02)

    # Evolve for 500 steps, recording stats every step
    lattice.run(steps=500, record_every=1)

    # Inspect the coarse-grained history
    import pandas as pd
    df = pd.DataFrame(lattice.history)
    print(df[["step", "time", "p_hat", "G", "L", "phase"]].tail(10))

    # Access the current field
    print(lattice.p.shape, lattice.p.mean())

INITIAL CONDITIONS
------------------
    init_uniform_high(noise_scale)    High-ED, nearly uniform  (ED inflation start)
    init_gaussian_seeds(n, amplitude) High background + planted Gaussian overdensities
    init_random_noise(lo, hi)         Uniform random field over [lo, hi]
    init_black_hole_patch(radius)     Saturated central disk in low-ED background
    init_two_body(separation)         Two high-ED Gaussian blobs (interaction test)
    init_from_array(arr)              Load an arbitrary pre-built field

EVOLUTION MODES
---------------
    step()            Single step (standard or mobility-weighted)
    run(steps)        Multi-step loop with optional history recording
    run_until(cond)   Evolve until a callable condition returns True

HISTORY & DIAGNOSTICS
----------------------
Each recorded snapshot stores:
    step      : integer step counter
    time      : accumulated simulation time (step * dt)
    p_hat     : mean ED density  p̂(t) = mean(p)
    G         : mean gradient magnitude  G(t) = mean(|∇p|)
    L         : homogeneity scale proxy  L(t) ≈ 1/G(t)
    a         : ED scale factor proxy  a(t) ∝ L(t)
    p_min_val : min(p) on the grid
    p_max_val : max(p) on the grid
    p_std     : std(p) — spread / structure amplitude
    phase     : detected dynamical phase label (string)

Snapshots of the full field p can be stored at chosen intervals via the
`snapshot_every` argument to run().
"""

from __future__ import annotations

import copy
import warnings
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple

import numpy as np

from ED_Update_Rule import (
    coarse_grained_stats,
    ed_step,
    ed_step_mobility,
    gradient_magnitude,
)


# ---------------------------------------------------------------------------
# Parameter container
# ---------------------------------------------------------------------------

@dataclass
class EDParams:
    """
    All tunable parameters for an ED lattice simulation.

    Attributes
    ----------
    alpha        : Relational penalty strength (competition / structure drain).
    beta         : Gradient smoothing strength (diffusion coefficient).
    gamma        : Concavity exponent for relational term, strictly 0 < gamma < 1.
    dt           : Time step.  Stability requires  dt * 4 * beta <= 1.
    p_min        : Hard floor on ED density (vacuum level).
    p_max        : Hard ceiling on ED density (saturation / rho_max).
    boundary     : "periodic" | "absorbing" | "reflecting"
    mode         : "standard" uses ed_step; "mobility" uses ed_step_mobility.
    noise_scale  : Std dev of per-site Gaussian noise added each step.
                   0.0 (default) = fully deterministic.
                   > 0 = Langevin dynamics (Scenario D: noisy universe).
    mobility_exp : Exponent n in M(rho) = ((rho_max-rho)/rho_max)^n.
                   Only used in mode="mobility".  Default 1.0 (linear, ED-12.5).
                   Larger values give a sharper horizon freeze.
    """
    alpha:        float = 0.05
    beta:         float = 0.25
    gamma:        float = 0.5
    dt:           float = 0.1
    p_min:        float = 0.01
    p_max:        float = 1.0
    boundary:     str   = "periodic"
    mode:         str   = "standard"   # "standard" | "mobility"
    noise_scale:  float = 0.0
    mobility_exp: float = 1.0

    def __post_init__(self) -> None:
        if not (0.0 < self.gamma < 1.0):
            raise ValueError(f"gamma must satisfy 0 < gamma < 1, got {self.gamma}")
        if self.dt * 4.0 * self.beta > 1.0:
            warnings.warn(
                f"CFL condition violated: dt*4*beta = {self.dt * 4 * self.beta:.3f} > 1. "
                "Consider reducing dt or beta to avoid numerical instability.",
                RuntimeWarning,
                stacklevel=2,
            )
        if self.mode not in ("standard", "mobility"):
            raise ValueError(f"mode must be 'standard' or 'mobility', got '{self.mode}'")
        if self.boundary not in ("periodic", "absorbing", "reflecting"):
            raise ValueError(f"boundary must be 'periodic', 'absorbing', or 'reflecting'")
        if self.noise_scale < 0.0:
            raise ValueError(f"noise_scale must be >= 0, got {self.noise_scale}")
        if self.mobility_exp <= 0.0:
            raise ValueError(f"mobility_exp must be > 0, got {self.mobility_exp}")


# ---------------------------------------------------------------------------
# Phase detection
# ---------------------------------------------------------------------------

# Thresholds for coarse-grained phase labelling.
# Phase boundaries are approximate; see ED-05.5 §4-8 for derivations.
_PHASE_G_HIGH  = 0.05   # above this G → still in inflation / smoothing
_PHASE_G_LOW   = 0.005  # below this G → late-time thinning or heat death
_PHASE_P_LOW   = 0.1    # p_hat below this → late-time / heat-death regime


def _detect_phase(p_hat: float, G: float) -> str:
    """
    Heuristic label for the current dynamical phase.

    Phase 1 – ED Inflation      : G large, p_hat high
    Phase 2 – Structure Seeds   : G moderate, p_hat high-moderate
    Phase 3 – Structure Formation: G low but nonzero, p_hat moderate
    Phase 4 – Late Thinning     : G very small, p_hat low
    """
    if G >= _PHASE_G_HIGH:
        return "1-inflation"
    if G >= _PHASE_G_LOW and p_hat >= _PHASE_P_LOW:
        return "2-structure_seeds"
    if G >= _PHASE_G_LOW and p_hat < _PHASE_P_LOW:
        return "4-late_thinning"
    if G < _PHASE_G_LOW and p_hat >= _PHASE_P_LOW:
        return "3-structure_formation"
    return "4-late_thinning"


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class EDLattice:
    """
    2D Event Density lattice simulation.

    Manages the ED density field p[i,j], drives time evolution via the
    update rules in ED_Update_Rule.py, and accumulates a diagnostic history.

    Parameters
    ----------
    rows   : Number of lattice rows.
    cols   : Number of lattice columns (defaults to rows for a square grid).
    params : EDParams instance.  Keyword arguments are forwarded to EDParams
             when params is None.
    seed   : Random seed for reproducible initial conditions.

    Examples
    --------
    >>> lattice = EDLattice(64, params=EDParams(alpha=0.05, beta=0.25))
    >>> lattice.init_uniform_high(noise_scale=0.02)
    >>> lattice.run(steps=200, record_every=10, snapshot_every=50)
    >>> lattice.history[-1]["phase"]
    '3-structure_formation'
    """

    def __init__(
        self,
        rows: int = 64,
        cols: Optional[int] = None,
        params: Optional[EDParams] = None,
        seed: Optional[int] = None,
        **param_kwargs,
    ) -> None:
        self.rows   = rows
        self.cols   = cols if cols is not None else rows
        self.params = params if params is not None else EDParams(**param_kwargs)
        self.rng    = np.random.default_rng(seed)

        # State
        self.p: np.ndarray = np.full(
            (self.rows, self.cols), self.params.p_max, dtype=np.float64
        )

        # Time tracking
        self.step_count: int   = 0
        self.time:       float = 0.0

        # History: list of dicts, one per recorded step
        self.history: List[dict] = []

        # Optional field snapshots: list of (step, array) tuples
        self.snapshots: List[Tuple[int, np.ndarray]] = []

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def shape(self) -> Tuple[int, int]:
        """Grid shape (rows, cols)."""
        return (self.rows, self.cols)

    @property
    def total_ed(self) -> float:
        """Sum of p over the entire grid — total Event Density."""
        return float(self.p.sum())

    @property
    def stats(self) -> dict:
        """
        Current coarse-grained statistics (not recorded to history).
        Returns p_hat, G, L, a, p_min_val, p_max_val, p_std, phase.
        """
        cg = coarse_grained_stats(self.p, boundary=self.params.boundary)
        return {
            **cg,
            "p_min_val": float(self.p.min()),
            "p_max_val": float(self.p.max()),
            "p_std":     float(self.p.std()),
            "phase":     _detect_phase(cg["p_hat"], cg["G"]),
        }

    # ------------------------------------------------------------------
    # Initial condition methods
    # ------------------------------------------------------------------

    def init_uniform_high(self, noise_scale: float = 0.02) -> "EDLattice":
        """
        High-ED, nearly uniform initial condition.

        Models the early-universe / ED-inflation starting state described in
        ED-05.5 §2.3: p ≈ p_max everywhere with small nonzero gradients.

        Parameters
        ----------
        noise_scale : Amplitude of random perturbations (relative to p_max).
                      Provides the 'small but nonzero initial gradients'
                      G(0) ≪ p_max required for later structure formation.
        """
        noise = self.rng.uniform(-noise_scale, noise_scale, size=self.shape)
        self.p = np.clip(
            self.params.p_max + noise * self.params.p_max,
            self.params.p_min,
            self.params.p_max,
        )
        self._reset_clock()
        return self

    def init_gaussian_seeds(
        self,
        n_seeds: int = 4,
        amplitude: float = 0.15,
        width: float = 0.06,
        background: float = 0.85,
    ) -> "EDLattice":
        """
        High-ED background with planted Gaussian overdensities.

        Seeds high-ED pockets at random locations to study structure
        formation from well-defined initial overdensities (δp ≪ p̂).

        Parameters
        ----------
        n_seeds    : Number of overdense seed regions.
        amplitude  : Peak δp above background (fraction of p_max).
        width      : Gaussian sigma as a fraction of grid size.
        background : Uniform background level (fraction of p_max).
        """
        p0 = np.full(self.shape, background * self.params.p_max)
        sigma_r = width * self.rows
        sigma_c = width * self.cols

        rows_c = self.rng.uniform(0.1, 0.9, n_seeds) * self.rows
        cols_c = self.rng.uniform(0.1, 0.9, n_seeds) * self.cols

        r_idx, c_idx = np.indices(self.shape, dtype=float)

        for rc, cc in zip(rows_c, cols_c):
            # Periodic-aware distance
            dr = np.abs(r_idx - rc)
            dc = np.abs(c_idx - cc)
            dr = np.minimum(dr, self.rows - dr)
            dc = np.minimum(dc, self.cols - dc)
            gauss = amplitude * self.params.p_max * np.exp(
                -0.5 * ((dr / sigma_r) ** 2 + (dc / sigma_c) ** 2)
            )
            p0 += gauss

        self.p = np.clip(p0, self.params.p_min, self.params.p_max)
        self._reset_clock()
        return self

    def init_random_noise(
        self, lo: float = 0.3, hi: float = 0.7
    ) -> "EDLattice":
        """
        Fully random uniform field over [lo, hi] (fractions of p_max).

        Useful for stress-testing or observing structure formation from a
        maximally disordered state.

        Parameters
        ----------
        lo : Lower bound as fraction of p_max.
        hi : Upper bound as fraction of p_max.
        """
        self.p = self.rng.uniform(
            lo * self.params.p_max,
            hi * self.params.p_max,
            size=self.shape,
        )
        self._reset_clock()
        return self

    def init_black_hole_patch(
        self,
        radius: float = 0.15,
        inner_value: Optional[float] = None,
        outer_value: float = 0.2,
    ) -> "EDLattice":
        """
        Saturated (or near-saturated) central disk in a low-ED background.

        Initializes the saturation regime described in ED-12.5 §1 (black hole
        analogue): a compact region with ρ ≈ ρ_max surrounded by low ED.
        The mobility-weighted update rule (mode='mobility') is recommended.

        Parameters
        ----------
        radius      : Disk radius as fraction of min(rows, cols).
        inner_value : ED inside the disk (defaults to p_max).
        outer_value : ED outside the disk (fraction of p_max).
        """
        inner = inner_value if inner_value is not None else self.params.p_max
        self.p = np.full(self.shape, outer_value * self.params.p_max)

        cr = self.rows / 2.0
        cc = self.cols / 2.0
        r_px = radius * min(self.rows, self.cols)

        r_idx, c_idx = np.indices(self.shape, dtype=float)
        dist = np.sqrt((r_idx - cr) ** 2 + (c_idx - cc) ** 2)
        self.p[dist <= r_px] = inner

        self.p = np.clip(self.p, self.params.p_min, self.params.p_max)
        self._reset_clock()
        return self

    def init_two_body(
        self,
        separation: float = 0.4,
        amplitude: float = 0.9,
        width: float = 0.08,
        background: float = 0.1,
    ) -> "EDLattice":
        """
        Two symmetric Gaussian high-ED blobs.

        Useful for studying interaction / competition between two high-ED
        pockets — the ED analogue of two gravitating bodies.

        Parameters
        ----------
        separation : Centre-to-centre distance as fraction of grid width.
        amplitude  : Peak value of each blob as fraction of p_max.
        width      : Gaussian sigma as fraction of grid size.
        background : Background level as fraction of p_max.
        """
        p0 = np.full(self.shape, background * self.params.p_max)
        sigma = width * min(self.rows, self.cols)
        r_idx, c_idx = np.indices(self.shape, dtype=float)

        cr = self.rows / 2.0
        half_sep = (separation / 2.0) * self.cols

        for sign in (-1.0, +1.0):
            cc = self.cols / 2.0 + sign * half_sep
            dr = r_idx - cr
            dc = c_idx - cc
            gauss = amplitude * self.params.p_max * np.exp(
                -0.5 * (dr ** 2 + dc ** 2) / sigma ** 2
            )
            p0 += gauss

        self.p = np.clip(p0, self.params.p_min, self.params.p_max)
        self._reset_clock()
        return self

    def init_from_array(self, arr: np.ndarray) -> "EDLattice":
        """
        Load an arbitrary pre-built density field.

        Parameters
        ----------
        arr : 2D numpy array. Must match the lattice shape (rows, cols).
              Values are clipped to [p_min, p_max].
        """
        if arr.shape != self.shape:
            raise ValueError(
                f"Array shape {arr.shape} does not match lattice shape {self.shape}."
            )
        self.p = np.clip(arr.astype(np.float64), self.params.p_min, self.params.p_max)
        self._reset_clock()
        return self

    # ------------------------------------------------------------------
    # Evolution methods
    # ------------------------------------------------------------------

    def step(self) -> "EDLattice":
        """
        Advance the lattice by exactly one time step.

        Uses ed_step (standard) or ed_step_mobility (saturation-aware)
        depending on params.mode.

        Returns self for chaining.
        """
        pr = self.params

        if pr.mode == "standard":
            self.p = ed_step(
                self.p,
                alpha=pr.alpha,
                beta=pr.beta,
                gamma=pr.gamma,
                dt=pr.dt,
                p_min=pr.p_min,
                p_max=pr.p_max,
                boundary=pr.boundary,
                noise_scale=pr.noise_scale,
                rng=self.rng,
            )
        else:  # "mobility"
            self.p = ed_step_mobility(
                self.p,
                alpha=pr.alpha,
                gamma=pr.gamma,
                dt=pr.dt,
                p_min=pr.p_min,
                p_max=pr.p_max,
                boundary=pr.boundary,
                mobility_exp=pr.mobility_exp,
                noise_scale=pr.noise_scale,
                rng=self.rng,
            )

        self.step_count += 1
        self.time       += pr.dt
        return self

    def run(
        self,
        steps: int,
        record_every: int = 1,
        snapshot_every: Optional[int] = None,
        verbose: bool = False,
        verbose_every: int = 50,
    ) -> "EDLattice":
        """
        Evolve the lattice for a fixed number of steps.

        Parameters
        ----------
        steps          : Total number of update steps to perform.
        record_every   : Record a history entry every this many steps.
                         Set to 0 to disable history recording.
        snapshot_every : Store a copy of the full field every this many steps.
                         None disables snapshots.
        verbose        : Print progress to stdout.
        verbose_every  : Print interval (in steps) when verbose=True.

        Returns
        -------
        self  (for chaining)
        """
        for i in range(steps):
            self.step()

            if record_every > 0 and (self.step_count % record_every == 0):
                self._record()

            if snapshot_every is not None and (self.step_count % snapshot_every == 0):
                self.snapshots.append((self.step_count, self.p.copy()))

            if verbose and (i % verbose_every == 0 or i == steps - 1):
                s = self.stats
                print(
                    f"  step={self.step_count:6d}  t={self.time:8.2f}"
                    f"  p_hat={s['p_hat']:.4f}  G={s['G']:.5f}"
                    f"  L={s['L']:.2f}  phase={s['phase']}"
                )

        return self

    def run_until(
        self,
        condition: Callable[["EDLattice"], bool],
        max_steps: int = 100_000,
        record_every: int = 1,
        snapshot_every: Optional[int] = None,
        verbose: bool = False,
        verbose_every: int = 100,
    ) -> "EDLattice":
        """
        Evolve until a user-supplied condition becomes True.

        Parameters
        ----------
        condition   : Callable(lattice) → bool.  Evolution stops when True.
                      Example: lambda lat: lat.stats["G"] < 1e-4
        max_steps   : Hard ceiling on iterations (safety valve).
        record_every, snapshot_every, verbose, verbose_every : same as run().

        Returns
        -------
        self
        """
        for _ in range(max_steps):
            if condition(self):
                break
            self.step()

            if record_every > 0 and (self.step_count % record_every == 0):
                self._record()

            if snapshot_every is not None and (self.step_count % snapshot_every == 0):
                self.snapshots.append((self.step_count, self.p.copy()))

            if verbose and (self.step_count % verbose_every == 0):
                s = self.stats
                print(
                    f"  step={self.step_count:6d}  t={self.time:8.2f}"
                    f"  p_hat={s['p_hat']:.4f}  G={s['G']:.5f}"
                    f"  phase={s['phase']}"
                )
        else:
            warnings.warn(
                f"run_until reached max_steps={max_steps} without condition being met.",
                RuntimeWarning,
                stacklevel=2,
            )
        return self

    # ------------------------------------------------------------------
    # Analysis helpers
    # ------------------------------------------------------------------

    def structure_mask(self, threshold: Optional[float] = None) -> np.ndarray:
        """
        Boolean mask of 'high-ED' sites (ED pockets / proto-structures).

        Parameters
        ----------
        threshold : Sites above this value are labelled as structures.
                    Defaults to mean(p) + 0.5 * std(p).

        Returns
        -------
        mask : Boolean array of shape (rows, cols).
        """
        if threshold is None:
            threshold = self.p.mean() + 0.5 * self.p.std()
        return self.p >= threshold

    def structure_count(self, threshold: Optional[float] = None) -> int:
        """
        Count disconnected high-ED pockets using connected-component labelling.

        Requires scipy.  Falls back to a simple threshold count if unavailable.

        Parameters
        ----------
        threshold : Passed to structure_mask().

        Returns
        -------
        n_structures : Number of distinct connected overdense regions.
        """
        mask = self.structure_mask(threshold)
        try:
            from scipy.ndimage import label as nd_label
            _, n = nd_label(mask)
            return int(n)
        except ImportError:
            # Fallback: just count True cells (not connected components)
            warnings.warn(
                "scipy not available; structure_count returns pixel count, "
                "not connected component count.",
                ImportWarning,
                stacklevel=2,
            )
            return int(mask.sum())

    def clone(self) -> "EDLattice":
        """Return a deep copy of this lattice (state, history, snapshots)."""
        return copy.deepcopy(self)

    def reset(self) -> "EDLattice":
        """
        Reset state to p_max everywhere, clear history and snapshots.
        Does NOT change params or the random seed state.
        """
        self.p = np.full(self.shape, self.params.p_max, dtype=np.float64)
        self._reset_clock()
        return self

    # ------------------------------------------------------------------
    # Snapshot utilities
    # ------------------------------------------------------------------

    def get_snapshot(self, step: int) -> Optional[np.ndarray]:
        """
        Retrieve a stored field snapshot by step number.

        Returns None if no snapshot was recorded at that step.
        """
        for s, arr in self.snapshots:
            if s == step:
                return arr
        return None

    def snapshot_steps(self) -> List[int]:
        """Return the list of step numbers at which snapshots were stored."""
        return [s for s, _ in self.snapshots]

    # ------------------------------------------------------------------
    # History utilities
    # ------------------------------------------------------------------

    def history_array(self, key: str) -> np.ndarray:
        """
        Extract a single scalar series from history as a 1D numpy array.

        Parameters
        ----------
        key : Any key in a history entry dict, e.g. "p_hat", "G", "L", "a".

        Returns
        -------
        arr : 1D float array of length len(history).
        """
        return np.array([entry[key] for entry in self.history], dtype=float)

    def current_phase(self) -> str:
        """Return the phase label for the current state."""
        s = self.stats
        return _detect_phase(s["p_hat"], s["G"])

    # ------------------------------------------------------------------
    # String representation
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        s = self.stats
        return (
            f"EDLattice(shape={self.shape}, "
            f"step={self.step_count}, t={self.time:.2f}, "
            f"p_hat={s['p_hat']:.4f}, G={s['G']:.5f}, "
            f"phase='{s['phase']}', "
            f"mode='{self.params.mode}', boundary='{self.params.boundary}')"
        )

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _record(self) -> None:
        """Append a history entry for the current state."""
        cg = coarse_grained_stats(self.p, boundary=self.params.boundary)
        entry = {
            "step":      self.step_count,
            "time":      self.time,
            "p_hat":     cg["p_hat"],
            "G":         cg["G"],
            "L":         cg["L"],
            "a":         cg["a"],
            "p_min_val": float(self.p.min()),
            "p_max_val": float(self.p.max()),
            "p_std":     float(self.p.std()),
            "phase":     _detect_phase(cg["p_hat"], cg["G"]),
        }
        self.history.append(entry)

    def _reset_clock(self) -> None:
        """Reset step counter, time, history, and snapshots."""
        self.step_count = 0
        self.time       = 0.0
        self.history    = []
        self.snapshots  = []


# ---------------------------------------------------------------------------
# Convenience factory functions
# ---------------------------------------------------------------------------

def make_inflation_lattice(
    size: int = 64,
    noise_scale: float = 0.02,
    seed: Optional[int] = 42,
    **param_kwargs,
) -> EDLattice:
    """
    Create a lattice pre-initialized for the ED inflation phase.

    High-ED (≈ p_max), small random perturbations, periodic boundary.
    Default parameters tuned for the inflation → structure formation sequence.

    Parameters
    ----------
    size        : Grid side length (square grid).
    noise_scale : Perturbation amplitude relative to p_max.
    seed        : Random seed.
    **param_kwargs : Forwarded to EDParams (override any defaults).

    Returns
    -------
    EDLattice ready to run.
    """
    params = EDParams(boundary="periodic", mode="standard", **param_kwargs)
    lat = EDLattice(rows=size, cols=size, params=params, seed=seed)
    lat.init_uniform_high(noise_scale=noise_scale)
    return lat


def make_structure_lattice(
    size: int = 64,
    n_seeds: int = 5,
    seed: Optional[int] = 42,
    **param_kwargs,
) -> EDLattice:
    """
    Create a lattice pre-initialized with planted structure seeds.

    Moderate background ED with n Gaussian overdensities — designed to
    fast-forward past the inflation phase directly to structure formation.

    Parameters
    ----------
    size     : Grid side length.
    n_seeds  : Number of Gaussian overdensity seeds to plant.
    seed     : Random seed.
    **param_kwargs : Forwarded to EDParams.

    Returns
    -------
    EDLattice ready to run.
    """
    params = EDParams(boundary="periodic", mode="standard", **param_kwargs)
    lat = EDLattice(rows=size, cols=size, params=params, seed=seed)
    lat.init_gaussian_seeds(n_seeds=n_seeds)
    return lat


def make_black_hole_lattice(
    size: int = 64,
    radius: float = 0.12,
    seed: Optional[int] = 42,
    **param_kwargs,
) -> EDLattice:
    """
    Create a lattice with a central saturated (black-hole-like) patch.

    Uses mobility-weighted evolution by default so that M(ρ) → 0 near ρ_max
    correctly freezes the saturated core (ED-12.5 §1).

    Parameters
    ----------
    size   : Grid side length.
    radius : Saturated patch radius as fraction of grid size.
    seed   : Random seed.
    **param_kwargs : Forwarded to EDParams (mode defaults to "mobility").

    Returns
    -------
    EDLattice ready to run.
    """
    kw = {"boundary": "absorbing", "mode": "mobility"}
    kw.update(param_kwargs)
    params = EDParams(**kw)
    lat = EDLattice(rows=size, cols=size, params=params, seed=seed)
    lat.init_black_hole_patch(radius=radius)
    return lat
