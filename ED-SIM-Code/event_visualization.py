"""
event_visualization.py
======================
Plotting and animation utilities for the ED MicroEvent engine.
Analogous to ED_Visualization.py (which handles 2D lattice fields),
this module visualizes ring-particle dynamics and micro-event diagnostics.

AVAILABLE FUNCTIONS
-------------------
Ring geometry:
    plot_ring(ring_state, ...)           Scatter of particles on PBC domain
    plot_ring_evolution(ring_state, ...) Snapshot strip of ring at recorded steps

Time-series:
    plot_observables(history, ...)       Multi-panel: d_min(t), R(t), delta_R(t)

Mechanism diagnostics:
    plot_mechanism_pie(events, ...)      Pie chart of mechanism label distribution
    plot_collapse_scatter(events, ...)   chi_emp vs diameter, colored by mechanism
    plot_regime_map(events, ...)         2D heatmap of mechanisms over (d, gamma)

Decision-tree snapshot:
    plot_decision_snapshot(event, ...)   Annotated geometric quantities + branch

Utility:
    save_figure(fig, path, dpi=150)

USAGE
-----
    from event_lattice import make_inward_ring
    from micro_event_operator import detect_micro_event
    from event_visualization import plot_ring, plot_observables

    ring = make_inward_ring(N=3, diameter=0.4)
    ring.run_until_collapse()
    plot_ring(ring)

    import matplotlib.pyplot as plt
    plt.show()
"""

from __future__ import annotations

import warnings
from collections import Counter
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.figure import Figure


# ---------------------------------------------------------------------------
# Mechanism color palette
# ---------------------------------------------------------------------------

MECHANISM_COLORS: Dict[str, str] = {
    "inward-collapse": "#2166ac",   # blue
    "outward-PBC":     "#b2182b",   # red
    "PBC-corner":      "#d6604d",   # salmon
    "other-late":      "#8073ac",   # purple
    "DECAY":           "#4dac26",   # green
}

_SERIES_COLORS = {
    "d_min":          "#2166ac",
    "R":              "#b2182b",
    "delta_R":        "#4dac26",
    "angular_spread": "#8073ac",
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _extract_ring(obj: Any):
    """Duck-type: return (positions, box_size) from a RingState or arrays."""
    if hasattr(obj, "positions") and hasattr(obj, "params"):
        return obj.positions, obj.params.box_size
    raise TypeError(f"Expected a RingState-like object, got {type(obj).__name__}")


def _ev_get(event: Any, attr: str, default=None):
    """Get attribute from a MicroEvent object or dict."""
    if hasattr(event, attr):
        return getattr(event, attr)
    if isinstance(event, dict):
        return event.get(attr, default)
    return default


# ---------------------------------------------------------------------------
# RING GEOMETRY PLOTS
# ---------------------------------------------------------------------------

def plot_ring(
    ring_state: Any,
    ax: Optional[plt.Axes] = None,
    title: Optional[str] = None,
    show_velocities: bool = True,
    arrow_scale: float = 0.05,
) -> Figure:
    """
    Scatter plot of particle positions on the PBC domain.

    Parameters
    ----------
    ring_state      : RingState object.
    ax              : Optional existing Axes.
    title           : Plot title.
    show_velocities : Draw velocity arrows.
    arrow_scale     : Length scaling for velocity arrows.

    Returns
    -------
    fig : matplotlib Figure.
    """
    pos, box = _extract_ring(ring_state)

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.get_figure()

    # Domain box
    ax.set_xlim(0, box)
    ax.set_ylim(0, box)
    ax.set_aspect("equal")
    ax.add_patch(plt.Rectangle((0, 0), box, box, fill=False, ec="grey", lw=1.5, ls="--"))

    # Particles
    ax.scatter(pos[:, 0], pos[:, 1], s=80, c="#2166ac", edgecolors="black", zorder=5)

    # Velocity arrows
    if show_velocities and hasattr(ring_state, "velocities"):
        vel = ring_state.velocities
        for i in range(pos.shape[0]):
            ax.annotate(
                "",
                xy=(pos[i, 0] + vel[i, 0] * arrow_scale,
                    pos[i, 1] + vel[i, 1] * arrow_scale),
                xytext=(pos[i, 0], pos[i, 1]),
                arrowprops=dict(arrowstyle="->", color="#b2182b", lw=1.5),
                zorder=6,
            )

    # Labels
    for i in range(pos.shape[0]):
        ax.text(
            pos[i, 0], pos[i, 1] + 0.015 * box,
            str(i), ha="center", va="bottom", fontsize=8, color="#333",
        )

    status = ""
    if hasattr(ring_state, "collapsed"):
        status = f"  |  {'COLLAPSED' if ring_state.collapsed else 'running'}"
    if hasattr(ring_state, "params"):
        p = ring_state.params
        ax.set_title(
            title or f"Ring  N={p.N}  d={p.diameter}  γ={p.gamma_gate}{status}"
        )
    else:
        ax.set_title(title or "Ring Geometry")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True, alpha=0.2)
    return fig


def plot_ring_evolution(
    ring_state: Any,
    max_panels: int = 6,
    figsize_per_panel: Tuple[float, float] = (3.5, 3.5),
    title: Optional[str] = None,
) -> Figure:
    """
    Snapshot strip of the ring at recorded history steps.

    Parameters
    ----------
    ring_state        : RingState with a populated history.
    max_panels        : Max snapshots to show.
    figsize_per_panel : Size per panel.
    title             : Figure super-title.

    Returns
    -------
    fig : matplotlib Figure.
    """
    history = getattr(ring_state, "history", [])
    if not history:
        warnings.warn("plot_ring_evolution: no history recorded.", UserWarning, stacklevel=2)
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No history.", ha="center", va="center", transform=ax.transAxes)
        ax.axis("off")
        return fig

    # Subsample
    indices = np.round(np.linspace(0, len(history) - 1, min(max_panels, len(history)))).astype(int)
    selected = [history[i] for i in indices]

    n = len(selected)
    w, h = figsize_per_panel
    fig, axes = plt.subplots(1, n, figsize=(w * n, h))
    if n == 1:
        axes = [axes]

    box = ring_state.params.box_size if hasattr(ring_state, "params") else 1.0

    for ax, entry in zip(axes, selected):
        ax.set_xlim(0, box)
        ax.set_ylim(0, box)
        ax.set_aspect("equal")
        ax.set_title(f"step {entry['step']}\nt={entry['time']:.3f}", fontsize=8)
        ax.tick_params(labelsize=6)

        # We don't store full positions in history — show d_min text instead
        ax.text(
            0.5, 0.5,
            f"d_min={entry['d_min']:.4f}\nR={entry['R']:.4f}\ndR={entry['delta_R']:.4f}",
            ha="center", va="center", transform=ax.transAxes, fontsize=8,
        )
        ax.grid(True, alpha=0.2)

    fig.suptitle(title or "Ring Evolution", fontsize=11, y=1.02)
    fig.subplots_adjust(wspace=0.3)
    return fig


# ---------------------------------------------------------------------------
# TIME-SERIES PLOTS
# ---------------------------------------------------------------------------

def plot_observables(
    history: List[dict],
    keys: Optional[Sequence[str]] = None,
    figsize: Tuple[float, float] = (10, 8),
    title: Optional[str] = None,
) -> Figure:
    """
    Multi-panel time series of ring geometric observables.

    Default panels: d_min(t), R(t), delta_R(t).

    Parameters
    ----------
    history : List of dicts from RingState.history.
    keys    : Which keys to plot.  Default ["d_min", "R", "delta_R"].
    figsize : Figure size.
    title   : Super-title.

    Returns
    -------
    fig : matplotlib Figure.
    """
    if not history:
        warnings.warn("plot_observables: empty history.", UserWarning, stacklevel=2)
        fig, _ = plt.subplots()
        return fig

    keys = list(keys) if keys is not None else ["d_min", "R", "delta_R"]
    n = len(keys)
    t = np.array([h["time"] for h in history])

    fig, axes = plt.subplots(n, 1, figsize=figsize, sharex=True)
    if n == 1:
        axes = [axes]
    fig.subplots_adjust(hspace=0.08)

    ylabels = {
        "d_min":          r"$d_{\min}(t)$",
        "R":              r"$R(t)$  (circumradius)",
        "delta_R":        r"$\Delta R(t)$",
        "angular_spread": r"angular spread (rad)",
        "pbc_edge_dist":  r"PBC edge dist",
        "pbc_corner_dist": r"PBC corner dist",
    }

    for ax, key in zip(axes, keys):
        vals = np.array([h.get(key, np.nan) for h in history], dtype=float)
        color = _SERIES_COLORS.get(key, "#333")
        ax.plot(t, vals, color=color, lw=1.5)
        ax.set_ylabel(ylabels.get(key, key), fontsize=9)
        ax.grid(True, alpha=0.3, ls="--")

    axes[-1].set_xlabel("Simulation time  t", fontsize=10)
    fig.suptitle(title or "Ring Observables", fontsize=12, y=1.01)
    return fig


# ---------------------------------------------------------------------------
# MECHANISM DIAGNOSTICS
# ---------------------------------------------------------------------------

def plot_mechanism_pie(
    events: List[Any],
    title: Optional[str] = None,
) -> Figure:
    """
    Pie chart of mechanism label distribution across a batch of MicroEvents.

    Parameters
    ----------
    events : List of MicroEvent objects (or dicts with a "mechanism" key).
    title  : Plot title.

    Returns
    -------
    fig : matplotlib Figure.
    """
    labels = [_ev_get(e, "mechanism") for e in events]
    counts = Counter(labels)

    fig, ax = plt.subplots(figsize=(6, 6))
    keys = list(counts.keys())
    vals = [counts[k] for k in keys]
    colors = [MECHANISM_COLORS.get(k, "#999999") for k in keys]

    ax.pie(vals, labels=keys, colors=colors, autopct="%1.0f%%", startangle=140)
    ax.set_title(title or "Mechanism Distribution")
    return fig


def plot_collapse_scatter(
    events: List[Any],
    title: Optional[str] = None,
) -> Figure:
    """
    Scatter plot: chi_emp vs diameter, colored by mechanism.

    Parameters
    ----------
    events : List of MicroEvent objects.
    title  : Plot title.

    Returns
    -------
    fig : matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    for e in events:
        mech = _ev_get(e, "mechanism")
        chi = _ev_get(e, "chi_emp")
        d = _ev_get(e, "diameter")
        color = MECHANISM_COLORS.get(mech, "#999")
        ax.scatter(d, chi, c=color, s=50, edgecolors="black", lw=0.5, zorder=3)

    # Legend
    patches = [
        mpatches.Patch(color=c, label=m) for m, c in MECHANISM_COLORS.items()
    ]
    ax.legend(handles=patches, fontsize=8, loc="upper right")
    ax.set_xlabel("Diameter  d")
    ax.set_ylabel(r"$\chi_{emp}$  (collapse time)")
    ax.set_title(title or r"$\chi_{emp}$ vs Diameter by Mechanism")
    ax.grid(True, alpha=0.3, ls="--")
    return fig


def plot_regime_map(
    events: List[Any],
    x_key: str = "diameter",
    y_key: str = "gamma_gate",
    title: Optional[str] = None,
) -> Figure:
    """
    Regime map: colored scatter/grid of mechanism labels over two parameters.

    Parameters
    ----------
    events : List of MicroEvent objects.
    x_key  : Attribute for x-axis (default "diameter").
    y_key  : Attribute for y-axis (default "gamma_gate").
    title  : Plot title.

    Returns
    -------
    fig : matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    gamma_order = {"inward": 0, "tangent": 1, "outward": 2}

    for e in events:
        mech = _ev_get(e, "mechanism")
        x = _ev_get(e, x_key)
        y_raw = _ev_get(e, y_key)

        # Convert gamma_gate to numeric if needed
        if y_key == "gamma_gate":
            y = gamma_order.get(y_raw, 0)
        else:
            y = float(y_raw)

        color = MECHANISM_COLORS.get(mech, "#999")
        ax.scatter(x, y, c=color, s=100, edgecolors="black", lw=0.5, marker="s", zorder=3)

    if y_key == "gamma_gate":
        ax.set_yticks([0, 1, 2])
        ax.set_yticklabels(["inward", "tangent", "outward"])

    patches = [mpatches.Patch(color=c, label=m) for m, c in MECHANISM_COLORS.items()]
    ax.legend(handles=patches, fontsize=8, loc="upper right")
    ax.set_xlabel(x_key.replace("_", " ").title())
    ax.set_ylabel(y_key.replace("_", " ").title())
    ax.set_title(title or f"Regime Map: {y_key} vs {x_key}")
    ax.grid(True, alpha=0.3, ls="--")
    return fig


# ---------------------------------------------------------------------------
# DECISION-TREE SNAPSHOT
# ---------------------------------------------------------------------------

def plot_decision_snapshot(
    event: Any,
    title: Optional[str] = None,
) -> Figure:
    """
    Annotated diagram of the geometric quantities used by the decision tree
    and which classification branch was taken.

    Parameters
    ----------
    event : A MicroEvent object.
    title : Plot title.

    Returns
    -------
    fig : matplotlib Figure.
    """
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.axis("off")

    def _get(attr, default="?"):
        if hasattr(event, attr):
            return getattr(event, attr)
        if isinstance(event, dict):
            return event.get(attr, default)
        return default

    mech = _get("mechanism")
    color = MECHANISM_COLORS.get(mech, "#999")

    lines = [
        f"Mechanism:   {mech}",
        "",
        f"collapse_step:   {_get('collapse_step')}",
        f"chi_emp:         {_get('chi_emp'):.4f}",
        f"d_min:           {_get('d_min_at_collapse'):.6f}",
        f"R:               {_get('R_at_collapse'):.6f}",
        f"delta_R:         {_get('delta_R_at_collapse'):.6f}",
        f"pbc_edge_dist:   {_get('pbc_edge_dist'):.6f}",
        f"pbc_corner_dist: {_get('pbc_corner_dist'):.6f}",
        "",
        f"N={_get('N')}  d={_get('diameter')}  "
        f"K={_get('K'):.2f}  gamma={_get('gamma_gate')}",
    ]

    # Decision path
    lines.append("")
    lines.append("Decision path:")
    collapsed = _get("collapse_step", -1) >= 0
    dR = _get("delta_R_at_collapse", 0.0)

    if not collapsed:
        lines.append("  1. No collapse → DECAY")
    elif dR < 0:
        lines.append(f"  1. Collapsed  →  2. dR={dR:.4f} < 0 → INWARD_COLLAPSE")
    elif dR > 0:
        edge = getattr(event, "pbc_edge_dist", 999)
        corner = getattr(event, "pbc_corner_dist", 999)
        lines.append(f"  1. Collapsed  →  2. dR={dR:.4f} > 0  →  3. PBC check")
        lines.append(f"     edge={edge:.4f}  corner={corner:.4f}")
        if mech == "PBC-corner":
            lines.append("     → near corner → PBC_CORNER")
        elif mech == "outward-PBC":
            lines.append("     → near edge → OUTWARD_PBC")
        else:
            lines.append("     → fallback → OTHER_LATE")
    else:
        lines.append(f"  1. Collapsed  →  2. dR≈0  →  4. OTHER_LATE")

    text = "\n".join(lines)
    ax.text(
        0.05, 0.95, text,
        transform=ax.transAxes, fontsize=10, va="top", fontfamily="monospace",
        bbox=dict(boxstyle="round,pad=0.5", fc=color, alpha=0.15),
    )
    ax.set_title(title or f"Decision Snapshot: {mech}", fontweight="bold")
    return fig


# ---------------------------------------------------------------------------
# Utility
# ---------------------------------------------------------------------------

def save_figure(fig: Figure, path: str, dpi: int = 150, **kwargs) -> None:
    """Save a matplotlib Figure to a file."""
    fig.savefig(path, dpi=dpi, bbox_inches="tight", **kwargs)
    print(f"Saved: {path}")
