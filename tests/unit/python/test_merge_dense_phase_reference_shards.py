from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT = REPO_ROOT / "scripts" / "pnjl" / "merge_dense_phase_reference_shards.py"
VALIDATOR = REPO_ROOT / "scripts" / "pnjl" / "validate_dense_reference_artifact.py"
TAG = "test"


HEADERS = {
    "boundary": ["xi", "T_MeV", "mu_transition_MeV", "rho_hadron", "rho_quark", "area_residual", "converged", "curve_parameter", "plot_order_key"],
    "cep": ["xi", "T_CEP_MeV", "muq_CEP_MeV", "muB_CEP_MeV", "uncertainty_T_MeV", "T_bracket_low_MeV", "T_bracket_high_MeV", "bracket_width_T_MeV"],
    "spinodals": ["xi", "T_MeV", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV", "rho_spinodal_hadron", "rho_spinodal_quark", "curve_parameter", "plot_order_key"],
    "crossover": ["xi", "mu_MeV", "T_crossover_MeV", "rho", "method", "converged", "derivative", "variable", "curve_parameter", "plot_order_key"],
    "phase_grid_convergence": ["axis", "xi", "T_MeV", "level", "left", "right", "midpoint", "position_error_MeV", "density_error", "maxwell_area", "response_rtol", "converged", "reason"],
}


def write_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_shard(root: Path, name: str, xis: list[float], conflict_at_zero: bool = False) -> None:
    shard = root / name
    shard.mkdir(parents=True)
    boundary_rows: list[list[str]] = []
    cep_rows: list[list[str]] = []
    spinodal_rows: list[list[str]] = []
    crossover_rows: list[list[str]] = []
    grid_rows: list[list[str]] = []
    runs = []
    for xi in xis:
        offset = 99.0 if conflict_at_zero and xi == 0.0 else 2.0 * xi
        boundary_rows.append([str(xi), "100", str(300 + offset), "1", "2", "0.00005", "true", "100", "100"])
        cep_rows.append([str(xi), str(130 + offset), str(295 + offset), str(3 * (295 + offset)), "0.05", "129.95", "130.05", "0.1"])
        spinodal_rows.append([str(xi), "100", str(310 + offset), str(290 + offset), "0.8", "2.2", "100", "100"])
        crossover_rows.append([str(xi), "0", str(180 + offset), "0.3", "peak", "true", "4", "phi_u", "0", "0"])
        grid_rows.append(["rho", str(xi), "100", "1", "0", "1", "1", "0.01", "0.001", "0.00005", "0", "true", "converged"])
        runs.append({
            "xi": xi,
            "run_id": f"xi_{xi}",
            "run_dir": f"processed/{name}/{xi}",
            "boundary_count": 1,
            "spinodal_count": 1,
            "crossover_count": 1,
            "cep_found": True,
        })

    write_csv(shard / f"boundary_{TAG}.csv", HEADERS["boundary"], boundary_rows)
    write_csv(shard / f"cep_{TAG}.csv", HEADERS["cep"], cep_rows)
    write_csv(shard / f"spinodals_{TAG}.csv", HEADERS["spinodals"], spinodal_rows)
    write_csv(shard / f"crossover_{TAG}.csv", HEADERS["crossover"], crossover_rows)
    write_csv(shard / f"phase_grid_convergence_{TAG}.csv", HEADERS["phase_grid_convergence"], grid_rows)

    meta = {
        "schema_version": "v1",
        "generator": {"git_commit": "abc123", "generated_at": "2026-01-01T00:00:00Z"},
        "column_definitions": [{"name": name} for name in HEADERS["crossover"]],
        "dense_meaning": {"xi_sampling": {"values": xis}},
    }
    (shard / f"crossover_{TAG}.meta.json").write_text(json.dumps(meta), encoding="utf-8")
    manifest = {
        "schema_version": "v1",
        "generator": {"git_commit": "abc123", "generated_at": "2026-01-01T00:00:00Z"},
        "config": {"tag": TAG, "xi_values": xis, "requested_xi_values": xis, "crossover_only": False},
        "output_root": "processed/test",
        "runs": runs,
    }
    (shard / f"phase_reference_{TAG}_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


def run_merge(shards: Path, output: Path, xi_convergence_root: Path | None = None) -> subprocess.CompletedProcess[str]:
    command = [
            sys.executable,
            str(SCRIPT),
            "--shards-root",
            str(shards),
            "--reference-root",
            str(output),
            "--tag",
            TAG,
            "--expected-xi-list=-0.1,0.0,0.1",
            "--overwrite",
        ]
    if xi_convergence_root is not None:
        command.extend(["--xi-convergence-root", str(xi_convergence_root)])
    return subprocess.run(
        command,
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )


def test_merge_is_deterministic_and_deduplicates_interval_endpoints(tmp_path: Path) -> None:
    shards = tmp_path / "shards"
    write_shard(shards, "left", [-0.1, 0.0])
    write_shard(shards, "right", [0.0, 0.1])
    output = tmp_path / "merged"

    first = run_merge(shards, output)
    assert first.returncode == 0, first.stderr
    manifest_path = output / f"phase_reference_{TAG}_manifest.json"
    first_hash = hashlib.sha256(manifest_path.read_bytes()).hexdigest()

    second = run_merge(shards, output)
    assert second.returncode == 0, second.stderr
    second_hash = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    assert second_hash == first_hash

    with (output / f"boundary_{TAG}.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert [float(row["xi"]) for row in rows] == [-0.1, 0.0, 0.1]
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["config"]["xi_values"] == [-0.1, 0.0, 0.1]
    assert len(manifest["shards"]) == 2
    validation = subprocess.run(
        [
            sys.executable,
            str(VALIDATOR),
            "--reference-root",
            str(output),
            "--tag",
            TAG,
            "--expect-full-reference",
        ],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )
    assert validation.returncode == 0, validation.stderr


def test_merge_rejects_conflicting_duplicate_endpoint(tmp_path: Path) -> None:
    shards = tmp_path / "shards"
    write_shard(shards, "left", [-0.1, 0.0])
    write_shard(shards, "right", [0.0, 0.1], conflict_at_zero=True)
    result = run_merge(shards, tmp_path / "merged")
    assert result.returncode != 0
    assert "conflicting duplicate boundary row" in (result.stdout + result.stderr)


def test_merge_includes_staged_xi_convergence_records(tmp_path: Path) -> None:
    shards = tmp_path / "shards"
    write_shard(shards, "left", [-0.1, 0.0])
    write_shard(shards, "right", [0.1])
    xi_plan = tmp_path / "xi_plan"
    xi_plan.mkdir()
    write_csv(
        xi_plan / f"xi_grid_convergence_{TAG}_level1.csv",
        HEADERS["phase_grid_convergence"],
        [["xi", "0", "", "1", "-0.1", "0.1", "0", "0.02", "0.001", "0.00005", "0.01", "true", "converged"]],
    )

    output = tmp_path / "merged"
    result = run_merge(shards, output, xi_plan)
    assert result.returncode == 0, result.stderr
    with (output / f"phase_grid_convergence_{TAG}.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert any(row["axis"] == "xi" and row["level"] == "1" for row in rows)
    manifest = json.loads((output / f"phase_reference_{TAG}_manifest.json").read_text(encoding="utf-8"))
    assert len(manifest["xi_refinement_records"]) == 1
