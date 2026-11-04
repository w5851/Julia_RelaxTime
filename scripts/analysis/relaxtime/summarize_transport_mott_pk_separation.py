"""Compare existing pi/K crossing brackets with adjacent transport responses.

Reads the fixed-background Mott probe; performs no solver or integration calls.
Process groups identify allowed s-channel exchanges, not meson-only amplitudes.
"""

import argparse
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


SPECIES = ("u", "d", "s", "ubar", "dbar", "sbar")
CONTRIBUTION_SPECIES = ("u", "s", "ubar", "sbar")
PI_S_PROCESSES = (
    "udbar_to_udbar", "dubar_to_dubar", "uubar_to_uubar", "uubar_to_ddbar"
)
K_S_PROCESSES = ("usbar_to_usbar", "subar_to_subar")


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("probe_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    assert not args.output.exists(), "refusing to overwrite an existing analysis"
    inputs = {
        "points": args.probe_root / "inputs/source_rows.csv",
        "contributions": args.probe_root / "full_v1/channel_contributions.csv",
        "crossings": args.probe_root / "final_review_v2/root_threshold_crossings.csv",
    }
    points = read_csv(inputs["points"])
    contributions = read_csv(inputs["contributions"])
    crossings = read_csv(inputs["crossings"])
    index = {(r["point_id"], r["species"], r["channel"]): r for r in contributions}
    assert len(index) == len(contributions), "duplicate process contribution"
    processes = sorted({r["channel"] for r in contributions})
    assert len(processes) == 17
    assert len(points) == 21 and len({r["point_id"] for r in points}) == 21
    assert {r["species"] for r in contributions} == set(CONTRIBUTION_SPECIES)
    assert len(contributions) == len(points) * len(CONTRIBUTION_SPECIES) * len(processes)
    for row in points:
        assert row["phase_curr"] == "quark"
        for species in CONTRIBUTION_SPECIES:
            total = sum(float(index[row["point_id"], species, p]["inverse_tau_contribution_inv_fm"])
                        for p in processes)
            assert math.isclose(total, float(row[f"tauinv_{species}"]), rel_tol=1e-10)
    brackets = {(r["window"], float(r["xi_left"]), float(r["xi_right"])): r["meson"]
                for r in crossings if r["meson"] in ("pi", "K") and r["cut"] == "pair"}
    assert len(brackets) == 4
    segments, groups, process_changes = [], [], []
    point_fields = ("point_id", "window", "T_MeV", "muB_MeV", "xi", "sigma_over_T",
                    "eta_over_s", *[f"tau_{sp}" for sp in SPECIES])
    for window in sorted({r["window"] for r in points}):
        ordered = sorted((r for r in points if r["window"] == window), key=lambda r: float(r["xi"]))
        assert len({r["T_MeV"] for r in ordered}) == 1
        for left, right in zip(ordered, ordered[1:]):
            xl, xr = float(left["xi"]), float(right["xi"])
            crossing = brackets.get((window, xl, xr), "")
            segment = dict(window=window, left_id=left["point_id"], right_id=right["point_id"],
                           xi_left=xl, xi_right=xr, crossing=crossing)
            for key in ("sigma_over_T", "eta_over_s", *[f"tau_{sp}" for sp in SPECIES]):
                a, b = float(left[key]), float(right[key])
                assert math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0
                segment[f"{key}_left"] = a
                segment[f"{key}_right"] = b
                segment[f"{key}_change_percent"] = 100 * (b / a - 1)
            segment["all_tau_increase"] = all(float(right[f"tau_{sp}"]) > float(left[f"tau_{sp}"])
                                              for sp in SPECIES)
            segment["sigma_decreases_while_all_tau_increase"] = (
                segment["sigma_over_T_change_percent"] < 0 and segment["all_tau_increase"])
            segments.append(segment)
            for species in CONTRIBUTION_SPECIES:
                total_delta = float(right[f"tauinv_{species}"]) - float(left[f"tauinv_{species}"])
                delta_by_process = {}
                for process in processes:
                    a = index[left["point_id"], species, process]
                    b = index[right["point_id"], species, process]
                    delta = float(b["inverse_tau_contribution_inv_fm"]) - float(a["inverse_tau_contribution_inv_fm"])
                    delta_by_process[process] = delta
                    if crossing:
                        process_changes.append(dict(window=window, crossing=crossing, xi_left=xl, xi_right=xr,
                            species=species, process=process, rate_left_fm2=float(a["rate_fm2"]),
                            rate_right_fm2=float(b["rate_fm2"]), contribution_change_inv_fm=delta,
                            total_change_inv_fm=total_delta))
                assert math.isclose(sum(delta_by_process.values()), total_delta, rel_tol=1e-9, abs_tol=1e-12)
                item = dict(window=window, crossing=crossing, xi_left=xl, xi_right=xr, species=species,
                            total_change_inv_fm=total_delta)
                for name, selected in (("pi_s_processes", PI_S_PROCESSES), ("K_s_processes", K_S_PROCESSES)):
                    delta = sum(delta_by_process[p] for p in selected)
                    item[f"{name}_change_inv_fm"] = delta
                    item[f"{name}_fraction_of_total_change"] = delta / total_delta if total_delta else None
                groups.append(item)
    crossing_segments = [r for r in segments if r["crossing"]]
    assert len(crossing_segments) == 4
    args.output.mkdir(parents=True)
    for name, rows in (("point_observations.csv", [{k: r[k] for k in point_fields} for r in points]),
                       ("adjacent_changes.csv", segments), ("crossing_responses.csv", crossing_segments),
                       ("channel_group_changes.csv", groups), ("crossing_process_changes.csv", process_changes)):
        write_csv(args.output / name, rows)
    summary = dict(status="completed_existing_data_analysis", finished_utc=datetime.now(timezone.utc).isoformat(),
        script_sha256=digest(Path(__file__)), source_sha="2e5c411c716d5cbc524c3fc6cfb335fe05b39444",
        points=len(points), adjacent_pairs=len(segments),
        contribution_species=CONTRIBUTION_SPECIES,
        reconstructed_species_rates=len(points) * len(CONTRIBUTION_SPECIES),
        input_hashes={str(p): digest(p) for p in inputs.values()}, crossing_responses=crossing_segments,
        crossing_groups=[r for r in groups if r["crossing"] and r["species"] in ("u", "ubar", "s", "sbar")],
        group_definitions={"pi_s_processes": PI_S_PROCESSES, "K_s_processes": K_S_PROCESSES},
        limits=["Original production rates only; no new integration or solver calls.",
                "The input contribution table contains u/s/ubar/sbar; d/dbar are not filled or reconstructed.",
                "Process group changes include other exchanges, interference and density changes.",
                "K crossing neighborhoods did not receive the preceding six-point local convergence study.",
                "No meson-amplitude ablation or species-resolved terminal-kernel decomposition.",
                "The two windows have different fixed temperatures and are not two peaks on one curve."],
        files={p.name: digest(p) for p in args.output.iterdir() if p.is_file()})
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(dict(status=summary["status"], points=len(points), adjacent_pairs=len(segments),
                          crossing_responses=crossing_segments,
                          crossing_groups=summary["crossing_groups"]), allow_nan=False))


if __name__ == "__main__":
    main()
