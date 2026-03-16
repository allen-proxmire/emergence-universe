"""
event_update.py
===============
Pure functions for single-step ballistic particle evolution and geometric
observables on a ring of N particles in a 2D periodic boundary condition (PBC)
domain.

This module is the geometric foundation of the ED MicroEvent engine.  All
functions are stateless (numpy in, numpy out) and contain no classification
or physics logic.

SOURCE
------
Derived from the ED Micro-Event Operator paper (ED-Arch-00):
    - d_min(t) = min_{i!=j} ||x_i(t) - x_j(t)||_PBC
    - R(t) = (1/N) * sum_i ||x_i(t) - c(t)||
    - Angular structure, PBC proximity, and boundary geometry

GAMMA GATES
-----------
Three velocity modes determine how particles move from their initial ring:
    "inward"  : velocities point radially inward  (toward centroid)
    "tangent" : velocities point tangentially      (along the ring)
    "outward" : velocities point radially outward  (away from centroid)

USAGE
-----
    from event_update import init_ring, ballistic_step, apply_pbc, d_min

    pos, vel = init_ring(N=3, diameter=0.4, gamma_gate="inward", box_size=1.0)
    for t in range(1000):
        pos = ballistic_step(pos, vel, dt=0.001)
        pos = apply_pbc(pos, box_size=1.0)
        dmin_val, pair = d_min(pos, box_size=1.0)
"""

from __future__ import annotations

from typing import Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Ballistic evolution
# ---------------------------------------------------------------------------

def ballistic_step(
    positions: np.ndarray,
    velocities: np.ndarray,
    dt: float,
) -> np.ndarray:
    """
    Advance particle positions by one ballistic timestep: x_new = x + v * dt.

    Parameters
    ----------
    positions  : (N, 2) array of current particle positions.
    velocities : (N, 2) array of particle velocities.
    dt         : Time step size.

    Returns
    -------
    new_positions : (N, 2) array.
    """
    return positions + velocities * dt


def apply_pbc(positions: np.ndarray, box_size: float) -> np.ndarray:
    """
    Wrap particle positions into the periodic domain [0, box_size)^2.

    Parameters
    ----------
    positions : (N, 2) array.
    box_size  : Side length of the square PBC domain.

    Returns
    -------
    wrapped : (N, 2) array with all coordinates in [0, box_size).
    """
    return positions % box_size


# ---------------------------------------------------------------------------
# PBC-corrected distances
# ---------------------------------------------------------------------------

def pairwise_distances_pbc(
    positions: np.ndarray,
    box_size: float,
) -> np.ndarray:
    """
    Full PBC-corrected pairwise distance matrix.

    Parameters
    ----------
    positions : (N, 2) array.
    box_size  : PBC domain side length.

    Returns
    -------
    dist : (N, N) symmetric distance matrix with zeros on the diagonal.
    """
    # (N, 1, 2) - (1, N, 2) → (N, N, 2)
    delta = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
    # Minimum image convention
    delta = delta - box_size * np.round(delta / box_size)
    return np.sqrt((delta ** 2).sum(axis=-1))


def d_min(
    positions: np.ndarray,
    box_size: float,
) -> Tuple[float, Tuple[int, int]]:
    """
    Minimum PBC-corrected pairwise distance and the index pair that achieves it.

    Parameters
    ----------
    positions : (N, 2) array.
    box_size  : PBC domain side length.

    Returns
    -------
    dmin_val : Scalar minimum distance.
    pair     : (i, j) indices of the closest particle pair.
    """
    dist = pairwise_distances_pbc(positions, box_size)
    # Mask diagonal
    np.fill_diagonal(dist, np.inf)
    flat_idx = int(np.argmin(dist))
    N = positions.shape[0]
    i, j = divmod(flat_idx, N)
    return float(dist[i, j]), (i, j)


# ---------------------------------------------------------------------------
# Circumradius and radial trends
# ---------------------------------------------------------------------------

def circumradius(positions: np.ndarray) -> float:
    """
    Circumradius R(t) = (1/N) * sum_i ||x_i - c||,  where c is the centroid.

    Parameters
    ----------
    positions : (N, 2) array.

    Returns
    -------
    R : Scalar circumradius.
    """
    centroid = positions.mean(axis=0)
    return float(np.mean(np.linalg.norm(positions - centroid, axis=1)))


def delta_R(R_current: float, R_previous: float) -> float:
    """
    Circumradius trend.  Sign encodes the dynamical regime:
        < 0 : inward contraction
        ~ 0 : tangent drift
        > 0 : outward expansion

    Parameters
    ----------
    R_current  : R(t).
    R_previous : R(t - dt).

    Returns
    -------
    dR : R_current - R_previous.
    """
    return R_current - R_previous


# ---------------------------------------------------------------------------
# Angular structure
# ---------------------------------------------------------------------------

def angular_positions(positions: np.ndarray) -> np.ndarray:
    """
    Angular positions of particles relative to the centroid (radians, [0, 2pi)).

    Parameters
    ----------
    positions : (N, 2) array.

    Returns
    -------
    angles : (N,) array of angles in [0, 2*pi).
    """
    centroid = positions.mean(axis=0)
    dx = positions[:, 0] - centroid[0]
    dy = positions[:, 1] - centroid[1]
    return np.arctan2(dy, dx) % (2.0 * np.pi)


def angular_spread(positions: np.ndarray) -> float:
    """
    Standard deviation of angular positions (circular std).

    A stable angular configuration (DECAY) has low angular spread over time.
    Tangent drift produces slowly increasing angular spread.

    Parameters
    ----------
    positions : (N, 2) array.

    Returns
    -------
    spread : Circular standard deviation in radians.
    """
    angles = angular_positions(positions)
    # Circular mean direction
    mean_sin = np.mean(np.sin(angles))
    mean_cos = np.mean(np.cos(angles))
    R_bar = np.sqrt(mean_sin ** 2 + mean_cos ** 2)
    # Circular std = sqrt(-2 * ln(R_bar)), clamped for numerical safety
    R_bar = np.clip(R_bar, 1e-15, 1.0)
    return float(np.sqrt(-2.0 * np.log(R_bar)))


# ---------------------------------------------------------------------------
# PBC boundary proximity
# ---------------------------------------------------------------------------

def pbc_proximity(
    positions: np.ndarray,
    box_size: float,
) -> Tuple[float, float]:
    """
    Distance from closest particle to nearest PBC edge and nearest PBC corner.

    Parameters
    ----------
    positions : (N, 2) array (assumed already wrapped into [0, box_size)).
    box_size  : PBC domain side length.

    Returns
    -------
    min_edge_dist   : Minimum distance from any particle to the nearest domain edge.
    min_corner_dist : Minimum distance from any particle to the nearest domain corner.
    """
    half = box_size / 2.0

    # Distance to nearest edge in each axis: min(x, box_size - x)
    dx_edge = np.minimum(positions[:, 0], box_size - positions[:, 0])
    dy_edge = np.minimum(positions[:, 1], box_size - positions[:, 1])
    min_edge_dist = float(np.min(np.minimum(dx_edge, dy_edge)))

    # Distance to nearest corner: corners at (0,0), (0,L), (L,0), (L,L)
    corners = np.array([
        [0.0, 0.0],
        [0.0, box_size],
        [box_size, 0.0],
        [box_size, box_size],
    ])
    # (N, 1, 2) - (1, 4, 2) → (N, 4)
    diffs = positions[:, np.newaxis, :] - corners[np.newaxis, :, :]
    # PBC-aware: use minimum image for corner distances too
    diffs = diffs - box_size * np.round(diffs / box_size)
    corner_dists = np.sqrt((diffs ** 2).sum(axis=-1))
    min_corner_dist = float(corner_dists.min())

    return min_edge_dist, min_corner_dist


# ---------------------------------------------------------------------------
# Ring initialization
# ---------------------------------------------------------------------------

def init_ring(
    N: int,
    diameter: float,
    gamma_gate: str,
    box_size: float = 1.0,
    speed: float = 1.0,
    seed: int | None = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Initialize N particles equispaced on a circle of given diameter, centered
    in the PBC domain, with velocities set by the gamma gate.

    Parameters
    ----------
    N          : Number of particles (>= 2).
    diameter   : Diameter of the initial ring.
    gamma_gate : One of "inward", "tangent", "outward".
    box_size   : Side length of the square PBC domain.
    speed      : Magnitude of each particle's velocity.
    seed       : Random seed (unused for deterministic init, reserved for
                 perturbation experiments).

    Returns
    -------
    positions  : (N, 2) array of initial particle positions.
    velocities : (N, 2) array of initial particle velocities.
    """
    if N < 2:
        raise ValueError(f"N must be >= 2, got {N}")
    if gamma_gate not in ("inward", "tangent", "outward"):
        raise ValueError(
            f"gamma_gate must be 'inward', 'tangent', or 'outward', got '{gamma_gate}'"
        )

    radius = diameter / 2.0
    center = np.array([box_size / 2.0, box_size / 2.0])

    # Equispaced angles
    angles = np.linspace(0, 2.0 * np.pi, N, endpoint=False)

    # Positions on the ring
    positions = np.column_stack([
        center[0] + radius * np.cos(angles),
        center[1] + radius * np.sin(angles),
    ])

    # Unit radial vectors (outward from center)
    radial = positions - center
    radial_norm = np.linalg.norm(radial, axis=1, keepdims=True)
    radial_norm = np.maximum(radial_norm, 1e-15)  # avoid division by zero
    radial_hat = radial / radial_norm

    # Unit tangent vectors (90 deg CCW from radial)
    tangent_hat = np.column_stack([-radial_hat[:, 1], radial_hat[:, 0]])

    if gamma_gate == "inward":
        velocities = -speed * radial_hat
    elif gamma_gate == "tangent":
        velocities = speed * tangent_hat
    elif gamma_gate == "outward":
        velocities = speed * radial_hat
    else:
        raise ValueError(f"Unknown gamma_gate: '{gamma_gate}'")

    # Wrap into PBC domain
    positions = apply_pbc(positions, box_size)

    return positions, velocities
