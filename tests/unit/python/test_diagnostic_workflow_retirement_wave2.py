from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[3]
PACKAGE = ROOT / "docs" / "analysis" / "governance" / "diagnostic_workflow_retirement_wave2_v1"
DEFINITIONS = PACKAGE / "definitions"
ACTIVE_WORKFLOW_DIR = ROOT / ".github" / "workflows"

EXPECTED_RETIRED = {
    "pnjl-c2-cep-limited-feasibility.yml",
    "pnjl-c2-cep-three-midpoint.yml",
    "pnjl-c2-cep-xi05-high-side-extension.yml",
    "pnjl-c2-limited-feasibility.yml",
    "pnjl-c2-targeted-closure-v1.yml",
    "pnjl-cep-deep-oracle.yml",
    "pnjl-cep-hybrid-production-shadow.yml",
    "pnjl-cep-narrow-pilot.yml",
    "pnjl-cep-production-shadow.yml",
    "pnjl-issue130-crossover-mu-endpoint-pilot.yml",
    "pnjl-issue130-maxwell-cep-local-pilot.yml",
    "pnjl-maxwell-endpoint-local-production-shadow-v4.yml",
    "pnjl-maxwell-endpoint-production-shadow.yml",
    "pnjl-maxwell-endpoint-refinement.yml",
    "pnjl-phase-diagram.yml",
    "pnjl-stagec-density-certificate-feasibility-v2.yml",
}


def test_wave2_manifest_and_definition_hashes_are_consistent() -> None:
    manifest = json.loads((PACKAGE / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "diagnostic_workflow_retirement_wave2_v1"
    assert manifest["solver_called"] is False
    assert manifest["github_runs_deleted"] is False
    assert manifest["github_artifacts_deleted"] is False
    assert manifest["numerical_results_changed"] is False
    assert manifest["counts"]["retired_workflows"] == len(EXPECTED_RETIRED)
    assert manifest["counts"]["retained_parameterized_workflows"] == 1

    with (PACKAGE / "retired_workflows.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert {row["workflow_file"] for row in rows} == EXPECTED_RETIRED
    for row in rows:
        path = DEFINITIONS / row["workflow_file"]
        assert path.is_file()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["source_sha256"]
        assert path.stat().st_size == int(row["bytes"])
        assert not (ACTIVE_WORKFLOW_DIR / row["workflow_file"]).exists()


def test_only_maxwell_expansion_remains_active_in_retired_family() -> None:
    retained = list(csv.DictReader((PACKAGE / "retained_parameterized_workflows.csv").open(newline="", encoding="utf-8")))
    assert [row["workflow_file"] for row in retained] == ["pnjl-issue130-maxwell-cep-local-expansion.yml"]
    workflow = ACTIVE_WORKFLOW_DIR / retained[0]["workflow_file"]
    assert workflow.is_file()
    payload = yaml.load(workflow.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)
    assert "workflow_dispatch" in payload["on"]
    text = workflow.read_text(encoding="utf-8")
    assert "target" in text.lower()
    assert "failed" in text.lower()


def test_wave2_reference_migrations_do_not_leave_active_retired_paths() -> None:
    migrations = list(csv.DictReader((PACKAGE / "reference_migrations.csv").open(newline="", encoding="utf-8")))
    assert len(migrations) >= 16
    for row in migrations:
        new_reference = row["new_reference"]
        if "<same>" in new_reference:
            continue
        assert new_reference.startswith("docs/analysis/governance/diagnostic_workflow_retirement_wave2_v1/definitions/")
        assert (ROOT / new_reference).is_file()
