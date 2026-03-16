"""
ED_Visualization.py
===================
Plotting and animation utilities for 2D Event Density lattice simulations.

Works standalone (accepts raw numpy arrays) or integrated with EDLattice
objects from ED_Lattice.py.  Default boundary condition throughout is
``"periodic"``, matching the EDLattice / EDParams default.

AVAILABLE FUNCTIONS
-------------------
Field plots  (accept a 2D array *p* or an EDLattice):
    plot_field(p, ...)           Heat-map of ED density ρ(x,y)
    plot_gradient(p, ...)        Heat-map of |∇ρ|(x,y)
    plot_overview(p, ...)        Three-panel: density, gradient, stats text box

History plots  (accept lattice.history or a list of dicts):
    plot_history(history, ...)   4-panel time series: p̂, G, L, a
    plot_phase_diagram(history, ...) Phase-space trajectory p̂(t) vs G(t)

Snapshot strip  (accepts lattice.snapshots or an EDLattice):
    plot_snapshot_strip(snapshots, ...) Row of field images at recorded steps

Animation:
    animate_field(frames, ...)   Animated density map from a list of arrays
                                 or an EDLattice object with stored snapshots

Utility:
    save_figure(fig, path, ...)  Save any Figure to a file

QUICK START
-----------
    import matplotlib.pyplot as plt
    from ED_Lattice import make_inflation_lattice
    from ED_Visualization import plot_overview, plot_history, animate_field

    # Run a simulation
    lat = make_inflation_lattice(size=64)
    lat.run(steps=500, record_every=5, snapshot_every=50)

    # Inspect current state (density + gradient + stats)
    plot_overview(lat)
    plt.show()

    # Plot the coarse-grained history
    plot_history(lat.history)
    plt.show()

    # Animate the stored snapshots
    anim = animate_field(lat)
    plt.show()

    # Save to GIF
    anim.save("ed_evolution.gif", writer="pillow", fps=10)

DESIGN NOTES
------------
* All functions that call into ED_Update_Rule (e.g. gradient_magnitude,
  coarse_grained_stats) default to ``boundary="periodic"`` when not
  otherwise specified — consistent with the EDParams default.

* EDLattice objects are accepted via duck-typing (no hard import) to avoid
  circular dependencies.  Functions check for .p, .history, .snapshots,
  and .params attributes as needed.

* All plot functions return matplotlib Figure objects so callers can
  further customise axes, add text, or call fig.savefig() themselves.

* Pass ``ax=...`` to overlay on an existing axes where supported.

* Animations require a matplotlib GUI backend or must be saved to a file
  (e.g. anim.save("out.gif", writer="pillow")).
"""

from __future__ import annotations

import warnings
from typing import Any, List, Optional, Sequence, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.animation import FuncAnimation
from matplotlib.figure import Figure

from ED_Update_Rule import coarse_grained_stats, gradient_magnitude


# ---------------------------------------------------------------------------
# Module-level constants: phase-boundary thresholds
# (must mirror _PHASE_G_HIGH / _PHASE_G_LOW / _PHASE_P_LOW in ED_Lattice.py)
# ---------------------------------------------------------------------------

_G_HIGH = 0.05    # G above this → inflation / smoothing phase
_G_LOW  = 0.005   # G below this → late-time thinning or structure formation
_P_LOW  = 0.10    # p_hat below this → late-time / heat-death regime

_PHASE_BG_COLORS = {
    "1-inflation":           "#4c72b0",
    "2-structure_seeds":     "#55a868",
    "3-structure_formation": "#c44e52",
    "4-late_thinning":       "#8172b2",
}

_SERIES_COLORS = {
    "p_hat": "#4c72b0",
    "G":     "#c44e52",
    "L":     "#55a868",
    "a":     "#8172b2",
    "p_std": "#ccb974",
}

_SERIES_YLABELS = {
    "p_hat":     r"$\hat{\rho}$(t)  (mean density)",
    "G":         r"G(t)  (mean |∇ρ|)",
    "L":         r"L(t) ≈ 1/G  (homogeneity scale)",
    "a":         r"a(t) ∝ L(t)  (scale factor proxy)",
    "p_std":     r"σ(ρ)(t)  (density std)",
    "p_min_val": r"ρ_min(t)",
    "p_max_val": r"ρ_max(t)",
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _extract_p(obj: Any) -> np.ndarray:
    """Return the 2D density array from either a raw array or an EDLattice."""
    if isinstance(obj, np.ndarray):
        return obj
    if hasattr(obj, "p") and isinstance(obj.p, np.ndarray):
        return obj.p
    raise TypeError(
        f"Expected a 2D numpy array or EDLattice-like object with a .p attribute, "
        f"got {type(obj).__name__}."
    )


def _extract_boundary(obj: Any, default: str = "periodic") -> str:
    """Read boundary condition from an EDLattice-like object, or return default."""
    try:
        return obj.params.boundary
    except AttributeError:
        return default


def _build_phase_segments(
    history: List[dict],
    t: np.ndarray,
) -> List[Tuple[float, float, str]]:
    """
    Build a list of (t_start, t_end, phase_label) segments from history.
    Returns an empty list if 'phase' is not in the history entries.
    """
    if not history or "phase" not in history[0]:
        return []

    segments: List[Tuple[float, float, str]] = []
    current_phase = history[0]["phase"]
    seg_start     = t[0]

    for i, h in enumerate(history[1:], 1):
        ph = h.get("phase", current_phase)
        if ph != current_phase:
            segments.append((seg_start, t[i], current_phase))
            current_phase = ph
            seg_start     = t[i]

    segments.append((seg_start, t[-1], current_phase))
    return segments


# ---------------------------------------------------------------------------
# FIELD PLOTS
# ---------------------------------------------------------------------------

def plot_field(
    p_or_lattice: Any,
    ax: Optional[plt.Axes] = None,
    cmap: str = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    title: Optional[str] = None,
    colorbar: bool = True,
    boundary: str = "periodic",
) -> Figure:
    """
    Heat-map of the 2D ED density field ρ(x, y).

    Parameters
    ----------
    p_or_lattice : 2D numpy array OR EDLattice object.
    ax           : Optional existing Axes to draw on.  A new figure is created
                   if None.
    cmap         : Matplotlib colormap name.  Default ``"viridis"``.
    vmin, vmax   : Colour-scale limits.  Defaults to [array.min(), array.max()].
    title        : Axes title.  Defaults to ``"ED Density Field  ρ(x, y)"``.
    colorbar     : Whether to attach a colorbar.
    boundary     : Boundary condition label shown in the plot footer.
                   Overridden automatically if an EDLattice is passed.
                   Default ``"periodic"``.

    Returns
    -------
    fig : matplotlib Figure.
    """
    p        = _extract_p(p_or_lattice)
    boundary = _extract_boundary(p_or_lattice, default=boundary)

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    else:
        fig = ax.get_figure()

    v0 = vmin if vmin is not None else float(p.min())
    v1 = vmax if vmax is not None else float(p.max())

    im = ax.imshow(
        p,
        origin="lower",
        cmap=cmap,
        vmin=v0,
        vmax=v1,
        interpolation="nearest",
        aspect="equal",
    )
    ax.set_xlabel("x (lattice units)")
    ax.set_ylabel("y (lattice units)")
    ax.set_title(title or "ED Density Field  ρ(x, y)")
    ax.text(
        0.01, 0.01,
        f"BC: {boundary}  |  shape: {p.shape[0]}×{p.shape[1]}",
        transform=ax.transAxes,
        fontsize=7,
        color="white",
        va="bottom",
    )

    if colorbar:
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label("ρ  (ED density)")

    return fig


def plot_gradient(
    p_or_lattice: Any,
    ax: Optional[plt.Axes] = None,
    cmap: str = "plasma",
    vmin: float = 0.0,
    vmax: Optional[float] = None,
    title: Optional[str] = None,
    colorbar: bool = True,
    boundary: str = "periodic",
) -> Figure:
    """
    Heat-map of the gradient magnitude |∇ρ|(x, y).

    Uses central differences via ``gradient_magnitude()`` from ED_Update_Rule.
    Defaults to periodic boundary conditions, matching the EDParams default.

    Parameters
    ----------
    p_or_lattice : 2D numpy array OR EDLattice object.
    ax           : Optional existing Axes.
    cmap         : Matplotlib colormap name.  Default ``"plasma"``.
    vmin         : Lower colour-scale limit.  Default 0.
    vmax         : Upper colour-scale limit.  Defaults to 95th percentile of |∇ρ|.
    title        : Axes title.
    colorbar     : Whether to attach a colorbar.
    boundary     : Boundary condition for gradient computation.
                   Overridden automatically if an EDLattice is passed.
                   Default ``"periodic"``.

    Returns
    -------
    fig : matplotlib Figure.
    """
    p        = _extract_p(p_or_lattice)
    boundary = _extract_boundary(p_or_lattice, default=boundary)

    G_field = gradient_magnitude(p, boundary=boundary)

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    else:
        fig = ax.get_figure()

    v1 = vmax if vmax is not None else float(np.percentile(G_field, 95))

    im = ax.imshow(
        G_field,
        origin="lower",
        cmap=cmap,
        vmin=vmin,
        vmax=max(v1, 1e-10),
        interpolation="nearest",
        aspect="equal",
    )
    ax.set_xlabel("x (lattice units)")
    ax.set_ylabel("y (lattice units)")
    ax.set_title(title or "Gradient Magnitude  |∇ρ|(x, y)")
    ax.text(
        0.01, 0.01,
        f"BC: {boundary}  |  mean |∇ρ| = {G_field.mean():.5f}",
        transform=ax.transAxes,
        fontsize=7,
        color="white",
        va="bottom",
    )

    if colorbar:
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label("|∇ρ|")

    return fig


def plot_overview(
    p_or_lattice: Any,
    boundary: str = "periodic",
    suptitle: Optional[str] = None,
) -> Figure:
    """
    Three-panel overview: density field | gradient field | coarse-grained stats.

    Parameters
    ----------
    p_or_lattice : 2D numpy array OR EDLattice object.
                   If an EDLattice, step count, time, and phase are shown
                   in the stats panel.
    boundary     : Default ``"periodic"``.  Overridden by
                   lattice.params.boundary if an EDLattice is passed.
    suptitle     : Figure super-title.  Defaults to ``"ED Lattice Overview"``.

    Returns
    -------
    fig : matplotlib Figure.
    """
    p        = _extract_p(p_or_lattice)
    boundary = _extract_boundary(p_or_lattice, default=boundary)

    stats = coarse_grained_stats(p, boundary=boundary)

    step_info  = ""
    phase_info = ""
    if hasattr(p_or_lattice, "step_count"):
        step_info  = f"step = {p_or_lattice.step_count}   t = {p_or_lattice.time:.2f}"
    if hasattr(p_or_lattice, "current_phase"):
        phase_info = f"phase: {p_or_lattice.current_phase()}"

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.subplots_adjust(wspace=0.38)

    plot_field(p, ax=axes[0], boundary=boundary, colorbar=True)
    plot_gradient(p, ax=axes[1], boundary=boundary, colorbar=True)

    # ---- Stats text panel ----
    ax_txt = axes[2]
    ax_txt.axis("off")

    L_str = f"{stats['L']:.3f}" if np.isfinite(stats["L"]) else "∞"
    a_str = f"{stats['a']:.3f}" if np.isfinite(stats["a"]) else "∞"

    rows: List[Tuple[str, str]] = []
    if step_info:
        rows.append(("", step_info))
    if phase_info:
        rows.append(("", phase_info))
    rows += [
        ("p̂  (mean density)",     f"{stats['p_hat']:.5f}"),
        ("G  (mean |∇ρ|)",        f"{stats['G']:.5f}"),
        ("L  (homogeneity scale)", L_str),
        ("a  (scale factor)",      a_str),
        ("ρ_min",                  f"{float(p.min()):.5f}"),
        ("ρ_max",                  f"{float(p.max()):.5f}"),
        ("σ(ρ)",                   f"{float(p.std()):.5f}"),
        ("shape",                  f"{p.shape[0]} × {p.shape[1]}"),
        ("boundary",               boundary),
    ]

    y = 0.97
    dy = 0.088
    for label, value in rows:
        if label:
            ax_txt.text(0.04, y, label + ":", fontsize=10, va="top", color="#555555")
            ax_txt.text(0.60, y, value, fontsize=10, va="top", fontweight="bold")
        else:
            ax_txt.text(0.04, y, value, fontsize=11, va="top",
                        fontweight="bold", color="#222222")
        y -= dy

    ax_txt.set_title("Coarse-grained Statistics", fontsize=11)
    ax_txt.set_xlim(0, 1)
    ax_txt.set_ylim(0, 1)

    fig.suptitle(suptitle or "ED Lattice Overview", fontsize=13, y=1.01)
    return fig


# ---------------------------------------------------------------------------
# HISTORY PLOTS
# ---------------------------------------------------------------------------

def plot_history(
    history: List[dict],
    keys: Optional[Sequence[str]] = None,
    figsize: Tuple[float, float] = (12, 9),
    shade_phases: bool = True,
    title: Optional[str] = None,
) -> Figure:
    """
    Multi-panel time series of coarse-grained ED observables.

    Default panels: p̂(t), G(t), L(t), a(t).

    Parameters
    ----------
    history      : List of dicts from ``EDLattice.history``.  Each entry must
                   contain ``"time"`` and all keys in `keys`.
    keys         : Which history keys to plot, one panel each.
                   Default ``["p_hat", "G", "L", "a"]``.
    figsize      : Figure (width, height) in inches.
    shade_phases : If True, shade each panel background by the dominant
                   dynamical phase (uses the ``"phase"`` key if present).
    title        : Figure super-title.  Default ``"ED Lattice History"``.

    Returns
    -------
    fig : matplotlib Figure.
    """
    if not history:
        warnings.warn(
            "plot_history: history list is empty.", UserWarning, stacklevel=2
        )
        fig, _ = plt.subplots()
        return fig

    keys     = list(keys) if keys is not None else ["p_hat", "G", "L", "a"]
    n_panels = len(keys)
    t        = np.array([h["time"] for h in history])

    fig, axes = plt.subplots(n_panels, 1, figsize=figsize, sharex=True)
    if n_panels == 1:
        axes = [axes]
    fig.subplots_adjust(hspace=0.08)

    # Phase segments for background shading
    segments = _build_phase_segments(history, t) if shade_phases else []

    for ax, key in zip(axes, keys):
        vals  = np.array([h.get(key, np.nan) for h in history], dtype=float)
        color = _SERIES_COLORS.get(key, "#333333")

        ax.plot(t, vals, color=color, lw=1.6, label=key)
        ax.set_ylabel(_SERIES_YLABELS.get(key, key), fontsize=9)
        ax.grid(True, alpha=0.3, linestyle="--")

        for ts, te, ph in segments:
            bg = _PHASE_BG_COLORS.get(ph, "grey")
            ax.axvspan(ts, te, alpha=0.10, color=bg, lw=0)

    axes[-1].set_xlabel("Simulation time  t", fontsize=10)

    # Legend for phase shading
    if segments:
        unique_phases = list(dict.fromkeys(ph for _, _, ph in segments))
        patches = [
            mpatches.Patch(
                color=_PHASE_BG_COLORS.get(ph, "grey"),
                alpha=0.45,
                label=ph.replace("-", "  –  "),
            )
            for ph in unique_phases
        ]
        axes[0].legend(
            handles=patches,
            loc="upper right",
            fontsize=8,
            framealpha=0.85,
            title="Dynamical phase",
        )

    fig.suptitle(title or "ED Lattice History", fontsize=12, y=1.01)
    return fig


def plot_phase_diagram(
    history: List[dict],
    ax: Optional[plt.Axes] = None,
    title: Optional[str] = None,
    colorbar: bool = True,
) -> Figure:
    """
    Phase-space trajectory: p̂(t) vs G(t), coloured by simulation time.

    Overlays the three phase-boundary lines (G_HIGH, G_LOW, P_LOW) that
    separate the four ED dynamical phases (ED-05.5 §4-8).

    Parameters
    ----------
    history  : List of dicts from ``EDLattice.history``.
    ax       : Optional existing Axes to draw on.
    title    : Plot title.  Default ``"ED Phase-Space Trajectory"``.
    colorbar : Whether to attach a colour-bar indicating simulation time.

    Returns
    -------
    fig : matplotlib Figure.
    """
    if not history:
        warnings.warn(
            "plot_phase_diagram: history list is empty.", UserWarning, stacklevel=2
        )
        fig, _ = plt.subplots()
        return fig

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 6))
    else:
        fig = ax.get_figure()

    p_hat_arr = np.array([h["p_hat"] for h in history])
    G_arr     = np.array([h["G"]     for h in history])
    t_arr     = np.array([h["time"]  for h in history])

    # Trajectory line (faint) then scatter coloured by time
    ax.plot(p_hat_arr, G_arr, color="grey", lw=0.6, alpha=0.35, zorder=2)
    sc = ax.scatter(
        p_hat_arr, G_arr, c=t_arr,
        cmap="cool", s=14, zorder=3, alpha=0.8,
    )

    # Phase-boundary lines
    x_lo = max(float(p_hat_arr.min()) - 0.05, 0.0)
    x_hi = min(float(p_hat_arr.max()) + 0.05, 1.0)
    ax.set_xlim(x_lo, x_hi)

    ax.axhline(
        _G_HIGH, color=_PHASE_BG_COLORS["1-inflation"], lw=1.1, ls="--", alpha=0.7,
        label=f"G = {_G_HIGH}  (inflation boundary)",
    )
    ax.axhline(
        _G_LOW, color=_PHASE_BG_COLORS["4-late_thinning"], lw=1.1, ls="--", alpha=0.7,
        label=f"G = {_G_LOW}  (late-time boundary)",
    )
    ax.axvline(
        _P_LOW, color=_PHASE_BG_COLORS["3-structure_formation"], lw=1.1, ls=":", alpha=0.7,
        label=f"p̂ = {_P_LOW}  (thinning boundary)",
    )

    # Phase region labels
    ax.text(x_hi - 0.01, _G_HIGH * 1.05, "inflation",
            ha="right", va="bottom", fontsize=8, color=_PHASE_BG_COLORS["1-inflation"], alpha=0.8)
    ax.text(x_hi - 0.01, (_G_HIGH + _G_LOW) / 2, "structure seeds",
            ha="right", va="center", fontsize=8, color=_PHASE_BG_COLORS["2-structure_seeds"], alpha=0.8)
    ax.text(x_hi - 0.01, _G_LOW * 0.6, "structure / late thinning",
            ha="right", va="top", fontsize=8, color=_PHASE_BG_COLORS["3-structure_formation"], alpha=0.8)

    ax.set_xlabel(r"$\hat{\rho}$(t)  (mean density)", fontsize=10)
    ax.set_ylabel("G(t)  (mean |∇ρ|)",               fontsize=10)
    ax.set_title(title or "ED Phase-Space Trajectory")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.85)
    ax.grid(True, alpha=0.3, linestyle="--")

    if colorbar:
        cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label("Simulation time  t")

    return fig


# ---------------------------------------------------------------------------
# SNAPSHOT STRIP
# ---------------------------------------------------------------------------

def plot_snapshot_strip(
    snapshots: Union[List[Tuple[int, np.ndarray]], Any],
    max_panels: int = 8,
    cmap: str = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    figsize_per_panel: Tuple[float, float] = (3.2, 3.0),
    title: Optional[str] = None,
    boundary: str = "periodic",
) -> Figure:
    """
    Horizontal strip of density-field snapshots at recorded time steps.

    Parameters
    ----------
    snapshots         : List of ``(step, array)`` tuples
                        (as stored in ``EDLattice.snapshots``), OR an
                        EDLattice object with a ``snapshots`` attribute.
    max_panels        : Maximum panels to display.  If there are more
                        snapshots than this, evenly-spaced ones are selected.
    cmap              : Colormap.  Default ``"viridis"``.
    vmin, vmax        : Global colour-scale limits shared across all panels.
                        Default to [global min, global max] of all frames.
    figsize_per_panel : ``(width, height)`` in inches for each sub-panel.
    title             : Figure super-title.  Default ``"ED Field Snapshots"``.
    boundary          : Boundary condition label.  Overridden automatically
                        if an EDLattice is passed.  Default ``"periodic"``.

    Returns
    -------
    fig : matplotlib Figure.
    """
    # Accept an EDLattice object
    if hasattr(snapshots, "snapshots"):
        boundary  = _extract_boundary(snapshots, default=boundary)
        snapshots = snapshots.snapshots

    if not snapshots:
        warnings.warn(
            "plot_snapshot_strip: no snapshots available.  "
            "Call EDLattice.run(..., snapshot_every=N) first.",
            UserWarning,
            stacklevel=2,
        )
        fig, ax = plt.subplots()
        ax.text(
            0.5, 0.5, "No snapshots recorded.",
            ha="center", va="center", transform=ax.transAxes, fontsize=12,
        )
        ax.axis("off")
        return fig

    # Subsample if needed
    if len(snapshots) > max_panels:
        idx       = np.round(np.linspace(0, len(snapshots) - 1, max_panels)).astype(int)
        snapshots = [snapshots[i] for i in idx]

    n    = len(snapshots)
    w, h = figsize_per_panel
    fig, axes = plt.subplots(1, n, figsize=(w * n, h + 0.7))
    if n == 1:
        axes = [axes]
    fig.subplots_adjust(wspace=0.04)

    all_vals = np.concatenate([arr.ravel() for _, arr in snapshots])
    v0 = vmin if vmin is not None else float(all_vals.min())
    v1 = vmax if vmax is not None else float(all_vals.max())

    im_last = None
    for ax, (step, arr) in zip(axes, snapshots):
        im_last = ax.imshow(
            arr,
            origin="lower",
            cmap=cmap,
            vmin=v0,
            vmax=v1,
            interpolation="nearest",
            aspect="equal",
        )
        ax.set_title(f"step {step}", fontsize=9)
        ax.axis("off")

    if im_last is not None:
        fig.colorbar(im_last, ax=axes[-1], fraction=0.07, pad=0.04, label="ρ")

    fig.suptitle(title or "ED Field Snapshots", fontsize=11, y=1.02)
    return fig


# ---------------------------------------------------------------------------
# ANIMATION
# ---------------------------------------------------------------------------

def animate_field(
    frames: Union[List[np.ndarray], Any],
    interval: int = 100,
    cmap: str = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    title: Optional[str] = None,
    show_step: bool = True,
    boundary: str = "periodic",
    figsize: Tuple[float, float] = (6, 5),
    repeat: bool = True,
) -> FuncAnimation:
    """
    Create a matplotlib FuncAnimation of the ED density field over time.

    Parameters
    ----------
    frames    : One of:
                  • A list of 2D numpy arrays — one per animation frame.
                  • An EDLattice object — uses ``lattice.snapshots`` as frames.
                    A warning is issued and a single-frame animation is returned
                    if no snapshots are stored.
    interval  : Delay between frames in milliseconds.  Default 100.
    cmap      : Colormap.  Default ``"viridis"``.
    vmin, vmax: Colour-scale limits.  Default to [min, max] of the first frame.
    title     : Figure title (static).  Default ``"ED Density Field"``.
    show_step : If True, overlay a step / time label on each frame.
    boundary  : Boundary condition label for the footer annotation.
                Overridden automatically if an EDLattice is passed.
                Default ``"periodic"``.
    figsize   : Figure size in inches.
    repeat    : Whether the animation loops.  Default True.

    Returns
    -------
    anim : ``matplotlib.animation.FuncAnimation``.
           Call ``plt.show()`` for interactive display, or
           ``anim.save("out.gif", writer="pillow", fps=10)`` to export.

    Notes
    -----
    To save as a GIF::

        anim = animate_field(lat)
        anim.save("ed_evolution.gif", writer="pillow", fps=10)

    To save as an MP4 (requires ffmpeg)::

        anim.save("ed_evolution.mp4", writer="ffmpeg", fps=15)
    """
    step_labels: List[str] = []

    # ---- Unpack an EDLattice ----
    if hasattr(frames, "snapshots"):
        lattice  = frames
        boundary = _extract_boundary(lattice, default=boundary)

        if lattice.snapshots:
            step_labels  = [
                f"step {s}  t={s * lattice.params.dt:.1f}"
                for s, _ in lattice.snapshots
            ]
            frame_arrays = [arr for _, arr in lattice.snapshots]
        else:
            warnings.warn(
                "animate_field: no snapshots stored on the lattice.  "
                "Run with snapshot_every=N to cache frames.  "
                "Returning single-frame animation of the current state.",
                UserWarning,
                stacklevel=2,
            )
            frame_arrays = [lattice.p.copy()]
            step_labels  = [f"step {lattice.step_count}  t={lattice.time:.1f}"]
    else:
        frame_arrays = list(frames)
        step_labels  = [f"frame {i}" for i in range(len(frame_arrays))]

    if not frame_arrays:
        raise ValueError("animate_field: frame list is empty.")

    # ---- Colour scale ----
    first = frame_arrays[0]
    v0    = vmin if vmin is not None else float(first.min())
    v1    = vmax if vmax is not None else float(first.max())
    v1    = max(v1, v0 + 1e-10)   # guard against flat fields

    # ---- Figure ----
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        first,
        origin="lower",
        cmap=cmap,
        vmin=v0,
        vmax=v1,
        interpolation="nearest",
        aspect="equal",
        animated=True,
    )
    ax.set_xlabel("x (lattice units)")
    ax.set_ylabel("y (lattice units)")
    ax.set_title(title or "ED Density Field  ρ(x, y)")

    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("ρ  (ED density)")

    # Step label overlay
    step_text = None
    if show_step and step_labels:
        step_text = ax.text(
            0.02, 0.97, step_labels[0],
            transform=ax.transAxes,
            fontsize=9,
            va="top",
            color="white",
            bbox=dict(boxstyle="round,pad=0.25", fc="black", alpha=0.40),
            animated=True,
        )

    bc_text = ax.text(
        0.01, 0.01, f"BC: {boundary}",
        transform=ax.transAxes,
        fontsize=7,
        color="white",
        va="bottom",
        animated=True,
    )

    def _update(i: int):
        im.set_data(frame_arrays[i])
        artists = [im, bc_text]
        if step_text is not None:
            lbl = step_labels[i] if i < len(step_labels) else f"frame {i}"
            step_text.set_text(lbl)
            artists.append(step_text)
        return tuple(artists)

    anim = FuncAnimation(
        fig,
        _update,
        frames=len(frame_arrays),
        interval=interval,
        blit=True,
        repeat=repeat,
    )
    return anim


# ---------------------------------------------------------------------------
# Utility: save any figure
# ---------------------------------------------------------------------------

def save_figure(fig: Figure, path: str, dpi: int = 150, **kwargs) -> None:
    """
    Save a matplotlib Figure to a file.

    Parameters
    ----------
    fig    : The Figure to save.
    path   : Output file path, e.g. ``"density_field.png"`` or ``"history.pdf"``.
    dpi    : Resolution in dots per inch.  Default 150.
    **kwargs : Forwarded to ``fig.savefig()``.

    Examples
    --------
    >>> fig = plot_overview(lat)
    >>> save_figure(fig, "overview.png", dpi=200)
    """
    fig.savefig(path, dpi=dpi, bbox_inches="tight", **kwargs)
    print(f"Saved: {path}")
