from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).parents[3]
SCRIPT = ROOT / "scripts" / "analysis" / "relaxtime" / "formalize_issue130_rs_candidate_figures.py"
SPEC = importlib.util.spec_from_file_location("issue130_rs_candidate_figure_formalization", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _write_fixture(root: Path, *, mode: str) -> tuple[Path, dict[str, str]]:
    case = root / mode / "candidate_runtime" / MODULE.CASE_NAME
    hashes: dict[str, str] = {}
    for index in range(MODULE.EXPECTED_FIGURE_COUNT):
        relative = f"plot_panel=panel{index % 3}/observable_{index}.png"
        path = case / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(f"png-{index}".encode("ascii"))
        hashes[relative] = MODULE.sha256(path)
    plot_manifest = {
        "schema_version": "v2",
        "count": MODULE.EXPECTED_FIGURE_COUNT,
        "figure_hashes": hashes,
        "source_csv": "external.csv",
        "source_csv_sha256": "a" * 64,
        "plot_entrypoint": "scripts/relaxtime/run_phase_guided_transport_plots.jl",
        "x": "xi",
        "y_columns": ["tau_u"],
        "format": "png",
        "dpi": 600,
    }
    (case / "plot_manifest.json").write_text(
        json.dumps(plot_manifest, sort_keys=True) + "\n", encoding="utf-8"
    )
    return case, hashes


def test_validate_candidate_tree_rejects_hash_drift(tmp_path: Path) -> None:
    case, hashes = _write_fixture(tmp_path, mode=MODULE.MODES[0])
    rendered = {
        "mode": MODULE.MODES[0],
        "reference": "candidate_runtime",
        "figure_root": str(case),
        "figure_count": MODULE.EXPECTED_FIGURE_COUNT,
        "figure_hashes": hashes,
    }
    (case / "observable_drift.png").write_bytes(b"not-in-manifest")
    try:
        MODULE.validate_candidate_tree(source_root=tmp_path, mode=MODULE.MODES[0], rendered=rendered)
    except ValueError as exc:
        assert "hash map" in str(exc)
    else:
        raise AssertionError("unlisted figure must be rejected")


def test_import_candidate_tree_preserves_bytes_and_marks_diagnostic(tmp_path: Path) -> None:
    case, hashes = _write_fixture(tmp_path / "source", mode=MODULE.MODES[0])
    source = {
        "mode": MODULE.MODES[0],
        "source_case_root": str(case),
        "source_plot_manifest": str(case / "plot_manifest.json"),
        "source_plot_manifest_sha256": MODULE.sha256(case / "plot_manifest.json"),
        "source_csv": "external.csv",
        "source_csv_sha256": "a" * 64,
        "plot_contract": {"entrypoint": "plot.jl", "x": "xi", "y_columns": ["tau_u"], "format": "png", "dpi": 600},
        "figure_hashes": hashes,
    }
    result = MODULE.import_candidate_tree(
        repo_root=tmp_path,
        target_root=tmp_path / "data" / "outputs" / "figures",
        source=source,
        source_manifest_sha256="b" * 64,
        workflow_head_sha="c" * 40,
        author_review_note="accepted",
    )
    manifest = json.loads(Path(result["plot_manifest"]).read_text(encoding="utf-8"))
    assert manifest["figure_status"] == "author_accepted_formal_layout"
    assert manifest["numerical_status"] == "diagnostic_only"
    assert manifest["solver_called"] is False
    assert manifest["figure_count"] == MODULE.EXPECTED_FIGURE_COUNT
    assert MODULE.sha256(Path(result["target_case_root"]) / "plot_panel=panel0" / "observable_0.png") == hashes["plot_panel=panel0/observable_0.png"]
