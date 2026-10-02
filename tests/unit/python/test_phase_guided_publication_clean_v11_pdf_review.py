from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import tomllib

from scripts.plotting.plot_provenance import validate_hash_record, validate_snapshot

ROOT = Path(__file__).resolve().parents[3]
CODE_REF = tomllib.loads((ROOT / "config/plotting/historical_snapshots.toml").read_text(encoding="utf-8"))["publication_clean_v11_stage"]["code_commit"]

import pytest

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts/analysis/relaxtime/export_phase_guided_publication_clean_v11_pdf_review.py"


@pytest.fixture(scope="module")
def module():
    spec = importlib.util.spec_from_file_location("v11_pdf_review_test", SCRIPT)
    assert spec is not None and spec.loader is not None
    loaded = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(loaded)
    return loaded


def test_v11_pdf_review_uses_the_frozen_renderer_and_preserves_source_hashes(module):
    package, index = module.read_json(module.PNG_PACKAGE), module.read_json(module.PNG_INDEX)
    assert package["delivery_stage"] == "png_review"
    assert len(index["charts"]) == 74
    assert validate_snapshot(module.PNG_PACKAGE, root=ROOT, code_ref=CODE_REF) == []
    assert module.V11_SCRIPT.name == "build_phase_guided_publication_clean_v11.py"
    pointer = module.read_json(module.V11.CURRENT_POINTER)
    assert pointer["current_analysis_package"].endswith("publication_clean_v5")


def test_v11_pdf_review_artifact_counts_technical_checks_and_retained_failures(module):
    index = module.read_json(module.FIGURE_ROOT / "pdf_review_index.json")
    assert index["manuscript_eligible"] is False
    assert index["formal_vector_delivery_complete"] is False
    assert index["single_figure_count"] == 72
    assert index["composite_figure_count"] == 2
    assert index["mode_counts"] == {"mode_a": 36, "mode_b": 36}
    assert len(index["charts"]) == len(list(module.FIGURE_ROOT.rglob("*.pdf"))) == 74
    modes = {"mode_a": 0, "mode_b": 0}
    composites = 0
    for chart in index["charts"]:
        path = ROOT / chart["manifest"]
        assert module.sha256_file(path) == chart["manifest_sha256"]
        assert module.verify_companion(path) == []
        manifest = module.read_json(path)
        assert module.sha256_file(ROOT / manifest["generator"]["path"]) == manifest["generator"]["sha256"]
        source = module.read_json(ROOT / manifest["source_png_manifest"]["path"])
        for field in ("axes", "series", "selection_rule", "interpolation_policy", "connector_policy", "missing_value_policy"):
            assert manifest[field] == source[field]
        assert manifest["submission_preflight"]["typography_exception_applied"] is False
        for field in ("manuscript_eligible", "current_publication_layer", "formal_vector_delivery_complete",
                      "solver_called", "canonical_data_modified", "new_display_values"):
            assert manifest[field] is False
        if chart["kind"] == "composite":
            composites += 1
            assert manifest["source_png_pixel_comparison"]["pixel_identical"] is True
            assert manifest["submission_preflight"]["violations"] == [module.GLYPH_BLOCKER]
            assert manifest["submission_preflight"]["passed"] is False
        else:
            modes[chart["mode_key"]] += 1
            assert manifest["submission_preflight"]["violations"] == []
            assert manifest["submission_preflight"]["passed"] is True
    assert modes == {"mode_a": 36, "mode_b": 36}
    assert composites == 2


def test_v11_pdf_package_retains_full_png_contract_and_prevents_overwrite(module):
    package = module.read_json(module.ANALYSIS_ROOT / "manifest.json")
    assert package["native_size_preflight_failure_count"] == 2
    for record in [*package["inputs"], *package["outputs"]]:
        assert validate_hash_record(record, root=ROOT, label="package", code_ref=CODE_REF) == []
    assert package["source_package_sha256"] == module.sha256_file(module.PNG_PACKAGE)
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        module.main()
