# run_law_21.py
#
# Reproduce the results of:
# ED‑Arch‑21 — Compositionality, y‑Switching, and Architectural Robustness

import subprocess
from pathlib import Path

def main():
    repo_root = Path(__file__).resolve().parents[2]

    cmd = [
        "python",
        str(repo_root / "ED Research" / "ED Simulations" / "ED_Arch_Harness" / "run_arch_harness.py"),
        "--suite", "invariants"
    ]

    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
