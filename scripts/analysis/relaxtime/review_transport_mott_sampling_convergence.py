"""Audit the retained local Mott convergence study; run on quote-remote."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

POINTS = ("p004", "p005", "p006", "p014", "p015", "p016")
CHANNELS = ("udbar_to_udbar", "dubar_to_dubar", "uubar_to_uubar", "uubar_to_ddbar")
OUTER_CHANNELS = CHANNELS[:2]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def finite(row: dict, key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"Nonfinite {key}: {row}")
    return value


def verify_job(directory: Path) -> None:
    completion = read_json(directory / "completed.json")
    for name, expected in completion["files"].items():
        if digest(directory / name) != expected:
            raise ValueError(f"Retained job hash differs: {directory / name}")


def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit("CASE_ROOT PROPAGATION_RUN OUTPUT_DIR CHANNEL_RUNS")
    case = Path(sys.argv[1]).resolve()
    propagation = case / sys.argv[2]
    output = Path(sys.argv[3]).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("Refusing to overwrite retained review")
    output.mkdir(parents=True, exist_ok=True)
    runs = sys.argv[4].split(",")
    if any(not name.replace("_", "").replace("-", "").isalnum() for name in runs):
        raise ValueError("Invalid run name")
    s_rows, outer_rows = [], []
    origins = {}
    for run in runs:
        summary = read_json(case / run / "summary.json")
        if digest(case / run / "channels.csv") != summary["channels_sha256"]:
            raise ValueError(f"Channel aggregate hash differs: {run}")
        for row in read_csv(case / run / "channels.csv"):
            verify_job(case / run / row["job_id"])
            key = row["point_id"], row["channel"]
            category = "s" if "reference_rate" in row else "outer"
            if (category, *key) in origins:
                raise ValueError(f"Duplicate channel: {key}")
            origins[category, *key] = run
            (s_rows if category == "s" else outer_rows).append(row)
    expected = {(point, channel) for point in POINTS for channel in CHANNELS}
    if {(r["point_id"], r["channel"]) for r in s_rows} != expected or len(s_rows) != 24:
        raise ValueError("Incomplete fixed-outer s convergence study")
    if outer_rows:
        outer_expected = {(point, channel) for point in POINTS for channel in OUTER_CHANNELS}
        if {(r["point_id"], r["channel"]) for r in outer_rows} != outer_expected or len(outer_rows) != 12:
            raise ValueError("Incomplete outer convergence study")
    derived = []
    for row in s_rows:
        rate = finite(row, "reference_rate")
        if rate <= 0:
            raise ValueError("Reference rate must be positive")
        previous = finite(row, "reference_previous_l1")
        shifted = finite(row, "reference_shifted_l1")
        reproduction = finite(row, "reproduction_relative")
        if max(previous, shifted) > 2e-5 or abs(reproduction) > 1e-6:
            raise ValueError("Required reference/reproduction gate failed")
        derived.append({
            "point_id": row["point_id"], "window": row["window"], "xi": finite(row, "xi"),
            "channel": row["channel"], "original_rate_fm2": finite(row, "original_rate"),
            "reference_rate_fm2": rate,
            "original_minus_reference_percent": 100 * (finite(row, "original_rate") / rate - 1),
            "original_reference_l1_percent": 100 * finite(row, "base_reference_l1"),
            "tight2_minus_reference_percent": 100 * (finite(row, "tight2_rate") / rate - 1),
            "tight2_reference_l1_percent": 100 * finite(row, "tight2_reference_l1"),
            "reference_order": int(row["reference_order"]),
            "reference_previous_l1_percent": 100 * previous,
            "reference_shifted_l1_percent": 100 * shifted,
            "direct_events_checked": int(row["direct_events_checked"]),
            "direct_subset_rate_percent": 100 * finite(row, "direct_subset_rate_fraction"),
            "direct_reference_l1_over_total_percent": 100 * finite(row, "direct_reference_l1_over_total"),
        })
    ps = read_json(propagation / "summary.json")
    for name, expected_hash in ps["files"].items():
        if digest(propagation / name) != expected_hash:
            raise ValueError(f"Propagation output hash differs: {name}")
    for run, expected_hash in ps["input_hashes"].items():
        if digest(case / run / "channels.csv") != expected_hash:
            raise ValueError(f"Propagation input hash differs: {run}")
    replay = read_csv(propagation / "formal_replay.csv")
    if len(replay) != 6 or any(r["passed"] != "true" for r in replay):
        raise ValueError("Formal transport replay did not pass")
    observables = read_csv(propagation / "observables.csv")
    shapes = read_csv(propagation / "local_shapes.csv")
    peak_rows = []
    for window in ("alpha1.0", "alpha1.1"):
        original = next(r for r in shapes if r["window"] == window
                        and r["configuration"] == "original" and r["observable"] == "sigma_over_T")
        for name in ("reference", "outer_joint1", "outer_joint2"):
            selected = [r for r in shapes if r["window"] == window
                        and r["configuration"] == name and r["observable"] == "sigma_over_T"]
            if not selected:
                if name == "reference" or outer_rows:
                    raise ValueError(f"Missing local shape: {window}/{name}")
                continue
            current, = selected
            original_bump = finite(original, "center_above_neighbor_chord")
            current_bump = finite(current, "center_above_neighbor_chord")
            peak_rows.append({
                "window": window, "configuration": name,
                "center_xi": finite(current, "center_xi"),
                "formal_center_sigma_over_T": finite(original, "center_value"),
                "current_center_sigma_over_T": finite(current, "center_value"),
                "center_change_percent": 100 * (finite(current, "center_value") / finite(original, "center_value") - 1),
                "formal_height_above_neighbor_chord": original_bump,
                "current_height_above_neighbor_chord": current_bump,
                "height_change_percent": 100 * (current_bump / original_bump - 1),
                "formal_center_minus_right": finite(original, "center_value") - finite(original, "right_value"),
                "current_center_minus_right": finite(current, "center_value") - finite(current, "right_value"),
                "center_is_local_max": current["center_is_local_max"] == "true",
                "sampled_max_xi": finite(current, "sampled_max_xi"),
            })
    outer_derived = []
    for row in outer_rows:
        outer_derived.append({
            "point_id": row["point_id"], "channel": row["channel"],
            "joint1_change_percent": 100 * finite(row, "joint1_relative_to_base"),
            "joint2_change_percent": 100 * finite(row, "joint2_relative_to_base"),
            "joint2_vs_joint1_percent": 100 * finite(row, "joint2_relative_to_joint1"),
            "max_s_l1_percent": 100 * finite(row, "max_s_l1"),
            "outer_converged": row["outer_converged"] == "true",
        })
    observable_changes = []
    for row in observables:
        if row["configuration"] not in ("reference", "outer_joint1", "outer_joint2"):
            continue
        observable_changes.append({
            "point_id": row["point_id"], "window": row["window"], "xi": finite(row, "xi"),
            "configuration": row["configuration"], "sigma_over_T": finite(row, "sigma_over_T"),
            "eta_over_s": finite(row, "eta_over_s"),
            "sigma_change_percent": 100 * finite(row, "sigma_relative_to_formal"),
            "eta_change_percent": 100 * finite(row, "eta_relative_to_formal"),
        })
    write_csv(output / "channel_comparison.csv", derived)
    write_csv(output / "peak_comparison.csv", peak_rows)
    write_csv(output / "observable_changes.csv", observable_changes)
    write_csv(output / "outer_comparison.csv", outer_derived)
    summary = {
        "status": "completed_diagnostic_review", "finished_utc": datetime.now(timezone.utc).isoformat(),
        "source_sha": ps["source_sha"], "review_script_sha256": digest(Path(__file__)),
        "s_jobs": len(s_rows), "outer_jobs": len(outer_rows), "formal_replay_passed": True,
        "all_reference_checks_passed": True,
        "all_base_l1_within_production_rtol": all(finite(r, "base_reference_l1") <= 1e-3 for r in s_rows),
        "max_channel_signed_percent": max(abs(r["original_minus_reference_percent"]) for r in derived),
        "max_channel_l1_percent": max(r["original_reference_l1_percent"] for r in derived),
        "max_tight2_channel_signed_percent": max(abs(r["tight2_minus_reference_percent"]) for r in derived),
        "max_tight2_channel_l1_percent": max(r["tight2_reference_l1_percent"] for r in derived),
        "all_outer_converged": all(r["outer_converged"] for r in outer_derived) if outer_derived else None,
        "max_outer_change_percent": max(abs(r["joint2_change_percent"]) for r in outer_derived) if outer_derived else None,
        "max_outer_last_step_percent": max(abs(r["joint2_vs_joint1_percent"]) for r in outer_derived) if outer_derived else None,
        "observable_maxima": {
            name: {
                "sigma_change_percent": max(abs(r["sigma_change_percent"]) for r in observable_changes if r["configuration"] == name),
                "eta_change_percent": max(abs(r["eta_change_percent"]) for r in observable_changes if r["configuration"] == name),
            }
            for name in sorted({r["configuration"] for r in observable_changes})
        },
        "peaks": peak_rows,
        "limits": [
            "Only six original backgrounds and four light channels were studied for s interpolation.",
            "Outer study, if present, refines two charged-pion-sensitive channels and fixes incoming density normalization.",
            "Direct event audits cover only the recorded subsets, not every thermal event.",
            "No new t/B0/background/terminal convergence claim and no meson-amplitude removal experiment.",
            "No formal CSV, PNG, production parameter or paper artifact was changed.",
        ],
        "input_hashes": {str(case / run / "channels.csv"): digest(case / run / "channels.csv") for run in runs},
        "propagation_summary_sha256": digest(propagation / "summary.json"),
        "files": {path.name: digest(path) for path in sorted(output.glob("*.csv"))},
    }
    (output / "review.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False))


if __name__ == "__main__":
    main()
