"""
event_lattice.py
================
Ring state management, initialization, and time evolution for the ED
MicroEvent engine.  Analogous to ED_Lattice.py (which manages a 2D PDE
lattice), this module manages N-particle rings evolving ballistically
inside a PBC domain.

QUICK START
-----------
    from event_lattice import make_inward_ring

    ring = make_inward_ring(N=3, diameter=0.4)
    ring.run_until_collapse()

    if ring.collapsed:
        print(f"Collapsed at step {ring.collapse_step}, t={ring.time:.4f}")
    else:
        print("DECAY: no collapse within horizon")

CLASSES
-------
    RingParams   : Dataclass of all tunable simulation parameters.
    RingState    : Manages particle state, evolution loop, and history.

FACTORIES
---------
    make_inward_ring(N, diameter, ...)
    make_tangent_ring(N, diameter, ...)
    make_outward_ring(N, diameter, ...)
    make_triad(configs)
"""

from __future__ import annotations

import copy
import warnings
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

from event_update import (
    angular_spread,
    apply_pbc,
    ballistic_step,
    circumradius,
    d_min,
    delta_R,
    init_ring,
    pbc_proximity,
)


# ---------------------------------------------------------------------------
# Parameter container
# ---------------------------------------------------------------------------

@dataclass
class RingParams:
    """
    All tunable parameters for a ring-particle simulation.

    Attributes
    ----------
    N          : Number of particles on the ring (>= 2).
    diameter   : Diameter of the initial ring.
    gamma_gate : Velocity mode: "inward", "tangent", or "outward".
    box_size   : Side length of the square PBC domain.
    dt         : Time step for ballistic evolution.
    speed      : Magnitude of each particle's initial velocity.
    merge_thr  : MERGE_THR — collapse threshold on d_min.
    max_steps  : Simulation horizon (max steps before declaring DECAY).
    seed       : Random seed for reproducible perturbation experiments.
    """
    N:          int   = 3
    diameter:   float = 0.4
    gamma_gate: str   = "inward"
    box_size:   float = 1.0
    dt:         float = 0.001
    speed:      float = 1.0
    merge_thr:  float = 0.01
    max_steps:  int   = 10_000
    seed:       int | None = None

    def __post_init__(self) -> None:
        if self.N < 2:
            raise ValueError(f"N must be >= 2, got {self.N}")
        if self.gamma_gate not in ("inward", "tangent", "outward"):
            raise ValueError(
                f"gamma_gate must be 'inward', 'tangent', or 'outward', "
                f"got '{self.gamma_gate}'"
            )
        if self.diameter <= 0:
            raise ValueError(f"diameter must be > 0, got {self.diameter}")
        if self.box_size <= 0:
            raise ValueError(f"box_size must be > 0, got {self.box_size}")
        if self.merge_thr <= 0:
            raise ValueError(f"merge_thr must be > 0, got {self.merge_thr}")

    @property
    def K(self) -> float:
        """Ring constant K = N / diameter."""
        return self.N / self.diameter


# ---------------------------------------------------------------------------
# Ring state
# ---------------------------------------------------------------------------

class RingState:
    """
    Manages the full lifecycle of a single ring simulation:
        initialization  →  ballistic evolution  →  collapse detection  →  history

    Parameters
    ----------
    params : RingParams instance.

    Examples
    --------
    >>> ring = RingState(RingParams(N=3, diameter=0.4, gamma_gate="inward"))
    >>> ring.init_ring()
    >>> ring.run_until_collapse()
    >>> ring.collapsed
    True
    """

    def __init__(self, params: RingParams) -> None:
        self.params = params

        # Particle state
        self.positions:  np.ndarray = np.empty((0, 2))
        self.velocities: np.ndarray = np.empty((0, 2))

        # Time tracking
        self.step_count: int   = 0
        self.time:       float = 0.0

        # Collapse state
        self._collapsed:    bool = False
        self._collapse_step: int = -1

        # History: list of dicts, one per recorded step
        self.history: List[dict] = []

        # Circumradius tracking (for delta_R computation)
        self._R_prev:    float = 0.0
        self._R_initial: float = 0.0

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    @property
    def collapsed(self) -> bool:
        """Whether the ring has collapsed (d_min <= merge_thr)."""
        return self._collapsed

    @property
    def collapse_step(self) -> int:
        """Step at which collapse occurred (-1 if no collapse)."""
        return self._collapse_step

    @property
    def N(self) -> int:
        """Number of particles."""
        return self.params.N

    # ------------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------------

    def init_ring(self) -> "RingState":
        """
        Initialize particles equispaced on a circle with gamma-gate velocities.
        Resets time, history, and collapse state.

        Returns self for chaining.
        """
        self.positions, self.velocities = init_ring(
            N=self.params.N,
            diameter=self.params.diameter,
            gamma_gate=self.params.gamma_gate,
            box_size=self.params.box_size,
            speed=self.params.speed,
            seed=self.params.seed,
        )
        self._reset()
        R0 = circumradius(self.positions)
        self._R_prev = R0
        self._R_initial = R0
        return self

    def init_from_arrays(
        self,
        positions: np.ndarray,
        velocities: np.ndarray,
    ) -> "RingState":
        """
        Load arbitrary particle state (for pause-resume tests).

        Parameters
        ----------
        positions  : (N, 2) array.
        velocities : (N, 2) array.

        Returns self for chaining.
        """
        if positions.shape != (self.params.N, 2):
            raise ValueError(
                f"positions shape {positions.shape} does not match (N={self.params.N}, 2)"
            )
        if velocities.shape != (self.params.N, 2):
            raise ValueError(
                f"velocities shape {velocities.shape} does not match (N={self.params.N}, 2)"
            )
        self.positions = positions.copy().astype(np.float64)
        self.velocities = velocities.copy().astype(np.float64)
        self._reset()
        R0 = circumradius(self.positions)
        self._R_prev = R0
        self._R_initial = R0
        return self

    # ------------------------------------------------------------------
    # Evolution
    # ------------------------------------------------------------------

    def step(self) -> "RingState":
        """
        Advance by one ballistic timestep.  Checks for collapse after
        the step.  Returns self for chaining.
        """
        if self._collapsed:
            return self

        # Track R before step for delta_R
        self._R_prev = circumradius(self.positions)

        self.positions = ballistic_step(
            self.positions, self.velocities, self.params.dt
        )
        self.positions = apply_pbc(self.positions, self.params.box_size)

        self.step_count += 1
        self.time += self.params.dt

        # Check collapse
        dmin_val, _ = d_min(self.positions, self.params.box_size)
        if dmin_val <= self.params.merge_thr:
            self._collapsed = True
            self._collapse_step = self.step_count

        return self

    def run(
        self,
        steps: int,
        record_every: int = 1,
    ) -> "RingState":
        """
        Evolve for a fixed number of steps (or until collapse).

        Parameters
        ----------
        steps        : Number of steps to attempt.
        record_every : Record a history entry every this many steps.

        Returns self for chaining.
        """
        for _ in range(steps):
            if self._collapsed:
                break
            self.step()
            if record_every > 0 and (self.step_count % record_every == 0):
                self._record()
        return self

    def run_until_collapse(
        self,
        record_every: int = 1,
    ) -> "RingState":
        """
        Evolve until collapse or the simulation horizon (max_steps).

        Parameters
        ----------
        record_every : Record a history entry every this many steps.

        Returns self for chaining.
        """
        return self.run(steps=self.params.max_steps, record_every=record_every)

    # ------------------------------------------------------------------
    # History and diagnostics
    # ------------------------------------------------------------------

    @property
    def cumulative_delta_R(self) -> float:
        """R(t_current) - R(t=0): overall circumradius trend since init."""
        return circumradius(self.positions) - self._R_initial

    def current_observables(self) -> dict:
        """
        Compute the full set of geometric observables at the current state.

        Returns
        -------
        dict with keys: step, time, d_min, d_min_pair, R, delta_R,
                        angular_spread, pbc_edge_dist, pbc_corner_dist
        """
        dmin_val, dmin_pair = d_min(self.positions, self.params.box_size)
        R_now = circumradius(self.positions)
        dR = delta_R(R_now, self._R_prev)
        ang_spread = angular_spread(self.positions)
        edge_dist, corner_dist = pbc_proximity(self.positions, self.params.box_size)

        return {
            "step":            self.step_count,
            "time":            self.time,
            "d_min":           dmin_val,
            "d_min_pair":      dmin_pair,
            "R":               R_now,
            "delta_R":         dR,
            "angular_spread":  ang_spread,
            "pbc_edge_dist":   edge_dist,
            "pbc_corner_dist": corner_dist,
        }

    def clone(self) -> "RingState":
        """Return a deep copy of this ring state."""
        return copy.deepcopy(self)

    # ------------------------------------------------------------------
    # String representation
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        status = "COLLAPSED" if self._collapsed else "running"
        return (
            f"RingState(N={self.params.N}, d={self.params.diameter}, "
            f"gamma={self.params.gamma_gate}, step={self.step_count}, "
            f"t={self.time:.4f}, {status})"
        )

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _record(self) -> None:
        """Append a history entry for the current state."""
        obs = self.current_observables()
        self.history.append(obs)
        # Update R_prev for next delta_R
        self._R_prev = obs["R"]

    def _reset(self) -> None:
        """Reset time, history, and collapse state."""
        self.step_count = 0
        self.time = 0.0
        self._collapsed = False
        self._collapse_step = -1
        self.history = []


# ---------------------------------------------------------------------------
# Factory functions
# ---------------------------------------------------------------------------

def make_inward_ring(
    N: int = 3,
    diameter: float = 0.4,
    box_size: float = 1.0,
    dt: float = 0.001,
    speed: float = 1.0,
    merge_thr: float = 0.01,
    max_steps: int = 10_000,
    seed: int | None = None,
) -> RingState:
    """Create and initialize a ring with inward gamma gate."""
    params = RingParams(
        N=N, diameter=diameter, gamma_gate="inward",
        box_size=box_size, dt=dt, speed=speed,
        merge_thr=merge_thr, max_steps=max_steps, seed=seed,
    )
    return RingState(params).init_ring()


def make_tangent_ring(
    N: int = 3,
    diameter: float = 0.4,
    box_size: float = 1.0,
    dt: float = 0.001,
    speed: float = 1.0,
    merge_thr: float = 0.01,
    max_steps: int = 10_000,
    seed: int | None = None,
) -> RingState:
    """Create and initialize a ring with tangent gamma gate."""
    params = RingParams(
        N=N, diameter=diameter, gamma_gate="tangent",
        box_size=box_size, dt=dt, speed=speed,
        merge_thr=merge_thr, max_steps=max_steps, seed=seed,
    )
    return RingState(params).init_ring()


def make_outward_ring(
    N: int = 3,
    diameter: float = 0.4,
    box_size: float = 1.0,
    dt: float = 0.001,
    speed: float = 1.0,
    merge_thr: float = 0.01,
    max_steps: int = 10_000,
    seed: int | None = None,
) -> RingState:
    """Create and initialize a ring with outward gamma gate."""
    params = RingParams(
        N=N, diameter=diameter, gamma_gate="outward",
        box_size=box_size, dt=dt, speed=speed,
        merge_thr=merge_thr, max_steps=max_steps, seed=seed,
    )
    return RingState(params).init_ring()


def make_triad(
    configs: List[dict],
    **shared_kwargs,
) -> List[RingState]:
    """
    Create multiple independent rings (for compositionality tests).

    Parameters
    ----------
    configs       : List of dicts, each with keys for RingParams
                    (e.g. {"N": 3, "diameter": 0.3, "gamma_gate": "inward"}).
    **shared_kwargs: Shared parameters applied to all rings (e.g. box_size, dt).

    Returns
    -------
    rings : List of initialized RingState objects.
    """
    rings = []
    for cfg in configs:
        merged = {**shared_kwargs, **cfg}
        params = RingParams(**merged)
        ring = RingState(params).init_ring()
        rings.append(ring)
    return rings
