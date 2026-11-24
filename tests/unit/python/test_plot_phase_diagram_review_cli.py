from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
PLOTTER = ROOT / "scripts" / "pnjl" / "plot_phase_diagram.py"


def run(*arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(PLOTTER), *arguments],
        cwd=ROOT,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )


@pytest.mark.parametrize(
    "arguments, expected",
    [
        (("--review-c2", "--format", "pdf"), "does not accept --format"),
        (("--review-c2", "--xi", "0.3"), "does not accept --xi"),
        (("--review-c2", "--dpi", "600"), "does not accept --dpi"),
        (("--review-c2", "--type", "combined"), "does not accept --type"),
        (("--review-c2", "--phase-reference-layer", "strict"), "does not accept --phase-reference-layer"),
        (("--review-c2", "--no-spinodal"), "does not accept --no-spinodal"),
        (("--review-c2", "--output-dir", "data/outputs/figures/guard"), "requires --output-dir and --results-dir together"),
        (("--review-c2", "--results-dir", "data/outputs/results/guard"), "requires --output-dir and --results-dir together"),
    ],
)
def test_review_mode_rejects_quick_plot_flags(arguments: tuple[str, ...], expected: str) -> None:
    completed = run(*arguments)
    assert completed.returncode == 2, completed.stderr
    assert expected in completed.stderr


def test_quick_plot_keeps_its_default_output_directory() -> None:
    """The review mode must not move the quick-plot default output directory."""
    completed = run("--help")
    assert completed.returncode == 0
    assert "data/outputs/figures/pnjl" in completed.stdout
