# run_law_II.py
#
# Reproduce the results of:
# ED‑Arch‑23 — Law II: Islands, Ridges, and Resonance

import subprocess
from pathlib import Path

def main():
    repo_root = Path(__file__).resolve().parents[2]

    cmd = [
        "python",
        str(repo_root / "ED Research" / "ED Simulations" / "ED_Arch_Harness" / "run_arch_harness.py"),
        "--suite", "law-surfaces"
    ]

    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
