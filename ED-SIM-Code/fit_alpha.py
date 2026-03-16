"""
fit_alpha.py
============
Compute the alpha parameter for the ED-Arch analytic chi law.

Loads enriched MicroEvent JSON files, extracts:
    - chi_emp
    - chi_pred

Then fits alpha in the relation:
    chi_emp ≈ alpha * chi_pred

This version fits alpha using ONLY inward-collapse events.
"""

import json
import numpy as np
from pathlib import Path


def load_events(path: str):
    with open(path, "r") as f:
        return json.load(f)


def fit_alpha(events):
    chi_emp = np.array([e["chi_emp"] for e in events])
    chi_pred = np.array([e["chi_pred"] for e in events])

    # Avoid division by zero
    mask = chi_pred > 1e-12
    chi_emp = chi_emp[mask]
    chi_pred = chi_pred[mask]

    # Fit alpha via least squares:
    #   chi_emp ≈ alpha * chi_pred
    alpha = np.sum(chi_emp * chi_pred) / np.sum(chi_pred ** 2)

    residuals = chi_emp - alpha * chi_pred
    rms = np.sqrt(np.mean(residuals ** 2))

    return alpha, rms, residuals


def main():
    # Load all events from the matrix scenario
    path = Path("me_output/sweep_events.json")
    events = load_events(path)

    # Filter to inward-collapse only
    inward_events = [e for e in events if e["gamma_gate"] == "inward"]

    if not inward_events:
        print("No inward-collapse events found.")
        return

    alpha, rms, residuals = fit_alpha(inward_events)

    print("\n=== Alpha Fit Results (Inward Only) ===")
    print(f"alpha = {alpha:.6f}")
    print(f"RMS residual = {rms:.6f}")
    print("Residuals:", residuals)


if __name__ == "__main__":
    main()
