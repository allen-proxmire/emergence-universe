"""
run_arch_harness.py
===================
CLI runner for the ED-Arch validation harness.  Orchestrates law-surface
validation (Arch-20) and architectural invariant tests (Arch-21), produces
summary reports and diagnostic plots.

SUITES
------
    law-surfaces : Sample analytic χ predictions across the (d, ncl, γ)
                   parameter space using Arch-20 law surfaces.
    invariants   : Run the seven INV-21-* invariant tests from Arch-21.

USAGE
-----
    python run_arch_harness.py --suite law-surfaces
    python run_arch_harness.py --suite invariants
    python run_arch_harness.py --suite invariants --invariant INV-21-3
    python run_arch_harness.py --suite all --outdir results/ --no-show
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Callable, Any

import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Allow imports from the parent ED Simulations directory
# ---------------------------------------------------------------------------
_PARENT = str(Path(__file__).resolve().parent.parent)
if _PARENT not in sys.path:
    sys.path.insert(0, _PARENT)

from event_lattice import RingParams, RingState, make_triad  # noqa: F401
from event_update import apply_pbc, circumradius, d_min, init_ring  # noqa: F401
from micro_event_operator import MicroEvent, detect_micro_event  # noqa: F401

from arch20_law_surfaces import (
    chi_inward,
    chi_outward,
    closing_rate,
    decay_window,
    deviation,
    pythagorean_drift,
    sample_law_surface,
)

from arch21_invariants import (
    InvariantResult,
    test_inv21_1_compositionality,
    test_inv21_2_memoryless_switching,
    test_inv21_3_pythagorean_drift,
    test_inv21_4_programmed_collapse,
    test_inv21_5_perturbation_hardness,
    test_inv21_6_decay_angular_subregime,
    test_inv21_7_ghost_compositionality,
)


# ═══════════════════════════════════════════════════════════════════════════
# Law-surface validation suite (analytic sampling)
# ═══════════════════════════════════════════════════════════════════════════

def run_law_surface_validation(
    ncl_range: Optional[List[int]] = None,
    d_range: Optional[List[float]] = None,
    outdir: Optional[str] = None,
) -> Dict[str, object]:
    """
    Sample the Arch-20 law surfaces across γ = -1, 0, +1.

    For each (d, ncl, γ) triple, record analytic χ_pred and mechanism.
    This version does not yet compare to engine χ_emp; deviations are
    only computed if 'chi_emp' is present in records.

    Parameters
    ----------
    ncl_range : List of cluster sizes to test.
    d_range   : List of initial separations (pixels) to test.
    outdir    : Directory for output files.  None → no file output.

    Returns
    -------
    Dict with keys: 'records', 'max_deviation', 'all_passed'.
    """
    if ncl_range is None:
        ncl_range = [3, 4, 5, 6]
    if d_range is None:
        d_range = [20.0, 30.0, 40.0, 50.0, 60.0]

    d_array = np.array(d_range, dtype=float)

    # Arch-20 uses gamma sheets: -1, 0, +1
    gamma_values = [-1, 0, 1]

    records: List[Dict[str, Any]] = []

    for g in gamma_values:
        if g == -1:
            # Outward sheet requires d_far; choose a reasonable far boundary
            d_far = float(max(d_range) * 1.8)
            points = sample_law_surface(
                gamma=g,
                ncl_range=ncl_range,
                d_range=d_array,
                d_far=d_far,
            )
        else:
            points = sample_law_surface(
                gamma=g,
                ncl_range=ncl_range,
                d_range=d_array,
            )

        # Convert LawSurfacePoint → dict
        for p in points:
            rec = {
                "d": float(p.d),
                "ncl": int(p.ncl),
                "gamma": int(p.gamma),
                "chi_pred": float(p.chi_predicted),
                "mechanism_pred": p.mechanism_predicted,
            }
            records.append(rec)

    # Compute deviations only if empirical χ is present
    devs: List[float] = []
    for rec in records:
        chi_pred = rec.get("chi_pred", None)
        chi_emp = rec.get("chi_emp", None)
        if chi_pred is None or chi_emp is None:
            continue
        devs.append(abs(float(chi_emp) - float(chi_pred)))

    max_dev = float(max(devs)) if devs else 0.0
    all_passed = max_dev <= 0.015 if devs else True

    summary = {
        "records": records,
        "max_deviation": max_dev,
        "all_passed": all_passed,
        "ncl_range": ncl_range,
        "d_range": d_range,
    }

    if outdir is not None:
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)
        with (outdir_path / "law_surface_results.json").open("w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)

    return summary


def run_angular_rate_check(
    ncl_range: Optional[List[int]] = None,
    outdir: Optional[str] = None,
) -> Dict[str, object]:
    """
    Verify the angular rate law K = 2·sin(π/ncl) in analytic form.

    Parameters
    ----------
    ncl_range : List of cluster sizes to test.
    outdir    : Directory for output files.

    Returns
    -------
    Dict with keys: 'records', 'max_deviation', 'all_passed'.
    """
    if ncl_range is None:
        ncl_range = [3, 4, 5, 6]

    records: List[Dict[str, Any]] = []

    for ncl in ncl_range:
        # Analytic K
        k_pred = 2.0 * np.sin(np.pi / float(ncl))
        # Use closing_rate helper (should match exactly)
        k_emp = closing_rate(ncl)

        dev = abs(float(k_emp) - float(k_pred))
        records.append(
            {
                "ncl": ncl,
                "k_pred": float(k_pred),
                "k_emp": float(k_emp),
                "deviation": float(dev),
            }
        )

    max_dev = float(max(r["deviation"] for r in records)) if records else 0.0
    all_passed = max_dev <= 1e-12  # should be analytically exact

    summary = {
        "records": records,
        "max_deviation": max_dev,
        "all_passed": all_passed,
        "ncl_range": ncl_range,
    }

    if outdir is not None:
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)
        with (outdir_path / "angular_rate_results.json").open("w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)

    return summary


def run_decay_window_check(
    ncl_range: Optional[List[int]] = None,
    outdir: Optional[str] = None,
) -> Dict[str, object]:
    """
    Inspect DECAY window boundaries for γ = 0 across cluster sizes.

    This version uses the analytic decay_window() oracle only; it does not
    run the engine. It reports the windows and a trivial boundary_accuracy
    of 1.0 (since no empirical classification is performed yet).

    Parameters
    ----------
    ncl_range : List of cluster sizes to test.
    outdir    : Directory for output files.

    Returns
    -------
    Dict with keys: 'records', 'boundary_accuracy', 'all_passed'.
    """
    if ncl_range is None:
        ncl_range = [3, 4, 5, 6]

    records: List[Dict[str, Any]] = []

    for ncl in ncl_range:
        d_lo, d_hi = decay_window(ncl=ncl)
        records.append(
            {
                "ncl": ncl,
                "d_lower": float(d_lo),
                "d_upper": float(d_hi),
            }
        )

    boundary_accuracy = 1.0
    all_passed = True

    summary = {
        "records": records,
        "boundary_accuracy": boundary_accuracy,
        "all_passed": all_passed,
        "ncl_range": ncl_range,
    }

    if outdir is not None:
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)
        with (outdir_path / "decay_window_results.json").open("w", encoding="utf-8") as f:
            json.dump(summary, f, indent=2)

    return summary


def run_engine_law_surface_validation(
    ncl_range: Optional[List[int]] = None,
    d_range: Optional[List[float]] = None,
    outdir: Optional[str] = None,
) -> Dict[str, object]:
    """
    Engine-integrated law-surface validation.

    For each (d, ncl, gamma) triple:
        - Use Arch-20 to get chi_pred, mechanism_pred
        - Run the ED engine to get chi_emp, mechanism_emp
        - Compute deviation = |chi_emp - chi_pred|

    Writes:
        - law_surface_engine_results.json
        - law_surface_engine_comparison.png
    """
    if ncl_range is None:
        ncl_range = [3, 4, 5, 6]
    if d_range is None:
        d_range = [20.0, 30.0, 40.0, 50.0, 60.0]

    # 1. Analytic sampling (reuse existing function)
    analytic = run_law_surface_validation(
        ncl_range=ncl_range,
        d_range=d_range,
        outdir=None,
    )["records"]

    records: List[Dict[str, Any]] = []

    # 2. Engine sampling — using the real engine runner from arch21_invariants
    from arch21_invariants import _run_single

    for rec in analytic:
        d = rec["d"]
        ncl = rec["ncl"]
        gamma = rec["gamma"]

        ring, event, chi_emp = _run_single(d, ncl, gamma)
        mechanism_emp = event.mechanism

        rec2 = dict(rec)
        rec2["chi_emp"] = float(chi_emp)
        rec2["mechanism_emp"] = mechanism_emp
        # If analytic prediction is DECAY, skip χ deviation entirely
        if rec["mechanism_pred"] == "DECAY":
            rec2["deviation"] = None
        else:
            rec2["deviation"] = abs(float(chi_emp) - float(rec["chi_pred"]))

        records.append(rec2)

    # --- DECAY summary block -----------------------------------------
    from collections import Counter

    decay_predicted = [r for r in records if r["mechanism_pred"] == "DECAY"]
    decay_total = len(decay_predicted)
    decay_collapsed = sum(
        1 for r in decay_predicted if r["mechanism_emp"] == "inward-collapse"
    )
    decay_preserved = decay_total - decay_collapsed
    decay_mech_counts = Counter(r["mechanism_emp"] for r in decay_predicted)

    decay_summary = {
        "decay_total": decay_total,
        "decay_preserved": decay_preserved,
        "decay_collapsed": decay_collapsed,
        "mechanism_counts": dict(decay_mech_counts),
    }
    # -----------------------------------------------------------------

    # Filter None deviations (DECAY points) when computing max
    real_devs = [r["deviation"] for r in records if r["deviation"] is not None]
    max_dev = float(max(real_devs)) if real_devs else 0.0
    all_passed = max_dev <= 0.015

    summary = {
        "records": records,
        "max_deviation": max_dev,
        "all_passed": all_passed,
        "ncl_range": ncl_range,
        "d_range": d_range,
        "decay_summary": decay_summary,
    }

    # 3. JSON output
    if outdir is not None:
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)
        with (outdir_path / "law_surface_engine_results.json").open(
            "w", encoding="utf-8"
        ) as f:
            json.dump(summary, f, indent=2)

    # 4. Scatter plot
    chi_pred_vals = [r["chi_pred"] for r in records]
    chi_emp_vals = [r["chi_emp"] for r in records]
    gamma_vals = [r["gamma"] for r in records]

    if chi_pred_vals and chi_emp_vals:
        chi_pred_arr = np.array(chi_pred_vals, dtype=float)
        chi_emp_arr = np.array(chi_emp_vals, dtype=float)
        gammas_arr = np.array(gamma_vals, dtype=float)

        plt.figure(figsize=(6, 6))
        sc = plt.scatter(
            chi_pred_arr, chi_emp_arr, c=gammas_arr, cmap="viridis", s=20, alpha=0.8
        )
        lims = [
            min(chi_pred_arr.min(), chi_emp_arr.min()),
            max(chi_pred_arr.max(), chi_emp_arr.max()),
        ]
        plt.plot(lims, lims, "k--", linewidth=1.0)
        plt.xlabel(r"$\chi_{\mathrm{pred}}$")
        plt.ylabel(r"$\chi_{\mathrm{emp}}$")
        plt.title("Engine-Integrated Law-Surface Validation")
        cbar = plt.colorbar(sc)
        cbar.set_label("γ-sheet")

        if outdir is not None:
            plt.savefig(
                outdir_path / "law_surface_engine_comparison.png",
                dpi=150,
                bbox_inches="tight",
            )
        plt.close()

    return summary



# ═══════════════════════════════════════════════════════════════════════════
# Invariant test suite
# ═══════════════════════════════════════════════════════════════════════════

_INVARIANT_REGISTRY: Dict[str, Callable[[], InvariantResult]] = {
    "INV-21-1": test_inv21_1_compositionality,
    "INV-21-2": test_inv21_2_memoryless_switching,
    "INV-21-3": test_inv21_3_pythagorean_drift,
    "INV-21-4": test_inv21_4_programmed_collapse,
    "INV-21-5": test_inv21_5_perturbation_hardness,
    "INV-21-6": test_inv21_6_decay_angular_subregime,
    "INV-21-7": test_inv21_7_ghost_compositionality,
}


def run_invariant_suite(
    outdir: Optional[str] = None,
) -> List[InvariantResult]:
    """
    Run all seven Arch-21 invariant tests and collect results.

    Parameters
    ----------
    outdir : Directory for output files.

    Returns
    -------
    List of InvariantResult objects (one per invariant).
    """
    results: List[InvariantResult] = []

    for inv_id, fn in _INVARIANT_REGISTRY.items():
        res = fn()
        results.append(res)

    if outdir is not None:
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)
        serializable = []
        for r in results:
            serializable.append(
                {
                    "id": getattr(r, "name", None),
                    "title": getattr(r, "title", None),
                    "passed": bool(getattr(r, "passed", False)),
                    "max_deviation": float(getattr(r, "deviation", 0.0)),
                    "details": getattr(r, "details", None),
                }
            )
        with (outdir_path / "invariant_results.json").open("w", encoding="utf-8") as f:
            json.dump(serializable, f, indent=2)

    return results


def run_single_invariant(
    invariant_id: str,
    outdir: Optional[str] = None,
) -> InvariantResult:
    """
    Run a single invariant test by its identifier.

    Parameters
    ----------
    invariant_id : One of "INV-21-1" through "INV-21-7".
    outdir       : Directory for output files.

    Returns
    -------
    InvariantResult for the specified invariant.

    Raises
    ------
    ValueError if invariant_id is not recognized.
    """
    if invariant_id not in _INVARIANT_REGISTRY:
        raise ValueError(f"Unknown invariant id: {invariant_id}")

    res = _INVARIANT_REGISTRY[invariant_id]()

    if outdir is not None:
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)
        payload = {
            "id": getattr(res, "id", None),
            "title": getattr(res, "title", None),
            "passed": bool(getattr(res, "passed", False)),
            "max_deviation": float(getattr(res, "max_deviation", 0.0)),
            "details": getattr(res, "details", None),
        }
        fname = f"invariant_{invariant_id}.json"
        with (outdir_path / fname).open("w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)

    return res


# ═══════════════════════════════════════════════════════════════════════════
# Plot helpers
# ═══════════════════════════════════════════════════════════════════════════

def _plot_law_surface_comparison(
    results: Dict[str, object],
    outdir: Optional[str] = None,
    show: bool = True,
) -> None:
    """
    Plot predicted vs empirical χ for the law-surface validation.

    If no empirical χ is present, this plot is skipped.

    Parameters
    ----------
    results : Output from run_law_surface_validation().
    outdir  : Directory for saving the figure.
    show    : Whether to call plt.show().
    """
    records = results.get("records", [])
    if not records:
        return

    chi_pred = []
    chi_emp = []
    gammas = []

    for rec in records:
        if "chi_pred" in rec and "chi_emp" in rec:
            chi_pred.append(float(rec["chi_pred"]))
            chi_emp.append(float(rec["chi_emp"]))
            gammas.append(float(rec.get("gamma", 0.0)))

    if not chi_pred:
        return

    chi_pred = np.array(chi_pred)
    chi_emp = np.array(chi_emp)
    gammas = np.array(gammas)

    plt.figure(figsize=(6, 6))
    sc = plt.scatter(chi_pred, chi_emp, c=gammas, cmap="viridis", s=20, alpha=0.8)
    lims = [
        min(chi_pred.min(), chi_emp.min()),
        max(chi_pred.max(), chi_emp.max()),
    ]
    plt.plot(lims, lims, "k--", linewidth=1.0)
    plt.xlabel(r"$\chi_{\mathrm{pred}}$")
    plt.ylabel(r"$\chi_{\mathrm{emp}}$")
    plt.title("Law-Surface Validation: Predicted vs Empirical χ")
    cbar = plt.colorbar(sc)
    cbar.set_label("γ-sheet")

    if outdir is not None:
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)
        plt.savefig(outdir_path / "law_surface_comparison.png", dpi=150, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()


def _plot_angular_rate_law(
    results: Dict[str, object],
    outdir: Optional[str] = None,
    show: bool = True,
) -> None:
    """
    Plot K_empirical vs K_predicted = 2·sin(π/ncl).

    Parameters
    ----------
    results : Output from run_angular_rate_check().
    outdir  : Directory for saving the figure.
    show    : Whether to call plt.show().
    """
    records = results.get("records", [])
    if not records:
        return

    k_pred = np.array([float(r["k_pred"]) for r in records])
    k_emp = np.array([float(r["k_emp"]) for r in records])
    ncl = np.array([int(r["ncl"]) for r in records])

    plt.figure(figsize=(6, 4))
    plt.plot(ncl, k_pred, "k--", label="K_pred = 2 sin(π/ncl)")
    plt.plot(ncl, k_emp, "o-", label="K_emp (closing_rate)")
    plt.xlabel("ncl")
    plt.ylabel("K")
    plt.title("Angular Rate Law")
    plt.legend()

    if outdir is not None:
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)
        plt.savefig(outdir_path / "angular_rate_law.png", dpi=150, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()


def _plot_invariant_report(
    invariant_results: List[InvariantResult],
    outdir: Optional[str] = None,
    show: bool = True,
) -> None:
    """
    Summary bar chart of invariant pass/fail status with deviation magnitudes.

    Parameters
    ----------
    invariant_results : List of InvariantResult objects.
    outdir            : Directory for saving the figure.
    show              : Whether to call plt.show().
    """
    if not invariant_results:
        return

    ids = [getattr(r, "name", "") for r in invariant_results]
    devs = [float(getattr(r, "deviation", 0.0)) for r in invariant_results]
    passed = [bool(getattr(r, "passed", False)) for r in invariant_results]

    colors = ["tab:green" if p else "tab:red" for p in passed]

    plt.figure(figsize=(8, 4))
    x = np.arange(len(ids))
    plt.bar(x, devs, color=colors)
    plt.xticks(x, ids, rotation=45)
    plt.ylabel("Max deviation")
    plt.title("Arch-21 Invariant Summary")

    if outdir is not None:
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)
        plt.savefig(outdir_path / "invariant_summary.png", dpi=150, bbox_inches="tight")

    if show:
        plt.show()
    else:
        plt.close()


# ═══════════════════════════════════════════════════════════════════════════
# Report generation
# ═══════════════════════════════════════════════════════════════════════════

def _save_harness_report(
    law_surface_results,
    invariant_results,
    engine_law_results=None,
    outdir="arch_harness_output",
) -> None:
    """
    Write a combined text report summarizing all harness results.

    Parameters
    ----------
    law_surface_results : Output from law-surface validation (or None).
    invariant_results   : List of InvariantResult objects (or None).
    engine_law_results  : Output from engine-integrated validation (or None).
    outdir              : Directory where the report file is written.
    """
    outdir_path = Path(outdir)
    outdir_path.mkdir(parents=True, exist_ok=True)
    report_path = outdir_path / "arch_harness_report.txt"

    lines: List[str] = []
    lines.append("ED-Arch Validation Harness Report")
    lines.append("=================================")
    lines.append("")

    # --- Analytic law-surface validation ---
    if law_surface_results is not None:
        lines.append("Law-Surface Validation (Analytic)")
        lines.append("---------------------------------")
        max_dev = law_surface_results.get("max_deviation", None)
        all_passed = law_surface_results.get("all_passed", None)
        n_records = len(law_surface_results.get("records", []))
        lines.append(f"  Points sampled: {n_records}")
        lines.append(f"  Max deviation:  {max_dev}")
        lines.append(f"  All passed:     {all_passed}")
        lines.append("")

    # --- Engine-integrated law-surface validation ---
    if engine_law_results is not None:
        lines.append("Law-Surface Validation (Engine-Integrated)")
        lines.append("------------------------------------------")
        max_dev = engine_law_results.get("max_deviation", None)
        all_passed = engine_law_results.get("all_passed", None)
        n_records = len(engine_law_results.get("records", []))
        lines.append(f"  Points sampled: {n_records}")
        lines.append(f"  Max deviation:  {max_dev}")
        lines.append(f"  All passed:     {all_passed}")
        ds = engine_law_results.get("decay_summary")
        if ds:
            lines.append(f"  DECAY region:   {ds['decay_total']} points "
                         f"({ds['decay_preserved']} preserved, "
                         f"{ds['decay_collapsed']} collapsed)")
        lines.append("")

    # --- Invariant suite ---
    if invariant_results is not None:
        lines.append("Invariant Suite (Arch-21)")
        lines.append("------------------------")
        for r in invariant_results:
            inv_name = getattr(r, "name", "")
            title = getattr(r, "title", "")
            passed = bool(getattr(r, "passed", False))
            dev = float(getattr(r, "deviation", 0.0))
            lines.append(
                f"  {inv_name:10s}  {'PASS' if passed else 'FAIL'}  "
                f"dev={dev:.4f}  {title}"
            )
        n_pass = sum(1 for r in invariant_results if getattr(r, "passed", False))
        n_total = len(invariant_results)
        lines.append("")
        lines.append(f"  Result: {n_pass}/{n_total} invariants passed.")
        lines.append("")

    with report_path.open("w", encoding="utf-8") as f:
        f.write("\n".join(lines))


# ═══════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════

def _parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments for the arch harness runner.

    Returns
    -------
    argparse.Namespace with fields: suite, invariant, outdir, no_show.
    """
    parser = argparse.ArgumentParser(
        description="ED-Arch Validation Harness — law-surface and invariant tests.",
    )
    parser.add_argument(
        "--suite",
        choices=["law-surfaces", "invariants", "all"],
        default="all",
        help="Which test suite to run (default: all).",
    )
    parser.add_argument(
        "--invariant",
        type=str,
        default=None,
        help="Run a single invariant by ID (e.g. INV-21-3).  Implies --suite invariants.",
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="arch_harness_output",
        help="Output directory for reports, plots, and JSON (default: arch_harness_output).",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Suppress interactive plot display (save only).",
    )
    parser.add_argument(
        "--engine-law-surfaces",
        action="store_true",
        help="Also run engine-integrated law-surface validation.",
    )

    return parser.parse_args()


def main() -> None:
    """
    Entry point for the ED-Arch Validation Harness CLI.

    Dispatches to the selected suite(s), generates plots and reports,
    and prints a summary to stdout.
    """
    args = _parse_args()
    outdir = args.outdir
    show = not args.no_show

    law_surface_results: Optional[Dict[str, object]] = None
    invariant_results: Optional[List[InvariantResult]] = None
    engine_law_results: Optional[Dict[str, object]] = None

    # If a single invariant is requested, force suite="invariants"
    if args.invariant is not None:
        inv_id = args.invariant.strip()
        print(f"Running single invariant: {inv_id}")
        invariant_results = [run_single_invariant(inv_id, outdir=outdir)]
        _plot_invariant_report(invariant_results, outdir=outdir, show=show)
        _save_harness_report(
            law_surface_results=None,
            invariant_results=invariant_results,
            engine_law_results=None,
            outdir=outdir,
        )
        return

    # Otherwise, dispatch based on suite
    if args.suite in ("law-surfaces", "all"):
        print("Running law-surface validation...")
        law_surface_results = run_law_surface_validation(outdir=outdir)
        _plot_law_surface_comparison(law_surface_results, outdir=outdir, show=show)

        if args.engine_law_surfaces:
            print("Running engine-integrated law-surface validation...")
            engine_law_results = run_engine_law_surface_validation(outdir=outdir)

        print("Running angular rate check...")
        angular_results = run_angular_rate_check(outdir=outdir)
        _plot_angular_rate_law(angular_results, outdir=outdir, show=show)

        print("Running DECAY window check...")
        decay_results = run_decay_window_check(outdir=outdir)

    if args.suite in ("invariants", "all"):
        print("Running invariant suite...")
        invariant_results = run_invariant_suite(outdir=outdir)
        _plot_invariant_report(invariant_results, outdir=outdir, show=show)

    # Combined report
    _save_harness_report(
        law_surface_results=law_surface_results,
        invariant_results=invariant_results,
        engine_law_results=engine_law_results,
        outdir=outdir,
    )

    print("")
    print(f"ED-Arch harness complete. Results written to: {outdir}")


if __name__ == "__main__":
    main()
