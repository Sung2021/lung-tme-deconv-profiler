#!/usr/bin/env python3
"""
run_pipeline.py
---------------
Orchestrate the TME-lung-cancer analysis pipeline by executing R scripts
in the required order (01 → 05). Stops on first failure with a clear message.

Usage:
    python python/run_pipeline.py              # run all steps
    python python/run_pipeline.py --steps 1 3  # run only steps 1 and 3
    python python/run_pipeline.py --rscript /path/to/Rscript

Requirements:
    pip install pyyaml
"""

import argparse
import subprocess
import sys
import time
import yaml
from pathlib import Path

# ---------------------------------------------------------------------------
# Step 1: Resolve project root and load config
# ---------------------------------------------------------------------------
ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = ROOT / "config.yaml"

with open(CONFIG_PATH) as f:
    config = yaml.safe_load(f)

R_SCRIPTS = [
    ROOT / "R" / "01_preprocessing.R",
    ROOT / "R" / "02a_signature_matrix.R",
    ROOT / "R" / "02b_deconvolution.R",
    ROOT / "R" / "03_stage_profiling.R",
    ROOT / "R" / "04_outlier_detection.R",
    ROOT / "R" / "05_scoring.R",
]

# ---------------------------------------------------------------------------
# Step 2: Parse arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description="Run the TME lung cancer R pipeline.")
parser.add_argument(
    "--steps",
    nargs="+",
    type=int,
    metavar="N",
    help=(
        "Step numbers to run (1-6). Default: all steps.\n"
        "  1: 01_preprocessing.R\n"
        "  2: 02a_signature_matrix.R\n"
        "  3: 02b_deconvolution.R\n"
        "  4: 03_stage_profiling.R\n"
        "  5: 04_outlier_detection.R\n"
        "  6: 05_scoring.R"
    ),
)
parser.add_argument(
    "--rscript",
    default="Rscript",
    metavar="PATH",
    help="Path to the Rscript executable. Default: 'Rscript' (from PATH).",
)
args = parser.parse_args()

selected = args.steps if args.steps else list(range(1, len(R_SCRIPTS) + 1))
scripts_to_run = [(i, R_SCRIPTS[i - 1]) for i in selected]


# ---------------------------------------------------------------------------
# Step 3: Helper — run a single R script with live output streaming
# ---------------------------------------------------------------------------
def run_r_script(step: int, script_path: Path, rscript: str) -> int:
    """Execute one R script via subprocess and return its exit code."""
    print(f"\n{'=' * 60}")
    print(f"[STEP {step}] {script_path.name}")
    print(f"{'=' * 60}")
    start = time.time()

    if not script_path.exists():
        print(f"[ERROR] Script not found: {script_path}")
        return 1

    result = subprocess.run(
        [rscript, str(script_path)],
        cwd=str(ROOT),
        text=True,
        # inherit stdout/stderr so output streams live to console
    )

    elapsed = time.time() - start
    status = "OK" if result.returncode == 0 else "FAILED"
    print(f"\n[{status}] {script_path.name} finished in {elapsed:.1f}s")
    return result.returncode


# ---------------------------------------------------------------------------
# Step 4: Run selected scripts in sequence
# ---------------------------------------------------------------------------
print("=" * 60)
print("TME Lung Cancer Analysis Pipeline")
print(f"Project root : {ROOT}")
print(f"Config       : {CONFIG_PATH}")
print(f"Steps to run : {[s for s, _ in scripts_to_run]}")
print("=" * 60)

start_total = time.time()

for step, script in scripts_to_run:
    code = run_r_script(step, script, args.rscript)
    if code != 0:
        print(f"\n[PIPELINE STOPPED] Step {step} ({script.name}) failed (exit code {code}).")
        print("Fix the error above, then re-run from this step with:")
        print(f"  python python/run_pipeline.py --steps {step}")
        sys.exit(code)

total = time.time() - start_total
print(f"\n{'=' * 60}")
print(f"Pipeline complete. {len(scripts_to_run)} step(s) finished in {total:.1f}s.")
print("=" * 60)
