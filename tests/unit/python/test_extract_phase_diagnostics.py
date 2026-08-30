from __future__ import annotations

import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "scripts" / "pnjl" / "extract_phase_diagnostics.py"


def load_module():
    spec = importlib.util.spec_from_file_location("extract_phase_diagnostics", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def args_for(module, process_root: Path, output: Path):
    return module.argparse.Namespace(
        process_root=process_root,
        output=output,
        tag="c0",
        stage="initial",
        shard_id="s000",
        xi=0.0,
        calculation_sha="a" * 40,
        postprocess_sha="b" * 40,
        source_workflow_run_id="123",
    )


def test_missing_summary_is_explicitly_unavailable(tmp_path: Path):
    module = load_module()
    payload = module.build_payload(args_for(module, tmp_path / "process", tmp_path / "diag.json"))
    assert payload["availability"] == "telemetry_unavailable"
    assert payload["unavailable_reason"] == "phase_summary_missing"
    assert payload["diagnostics"] is None
    assert payload["missing_fields"] == list(module.COUNTER_FIELDS)


def test_projection_preserves_counters_and_classification_counts(tmp_path: Path):
    module = load_module()
    summary_dir = tmp_path / "process" / "xi_0p0"
    summary_dir.mkdir(parents=True)
    (summary_dir / "phase_summary.json").write_text(
        json.dumps(
            {
                "schema_version": "phase-v2",
                "config_hash": "config-hash",
                "stats": {
                    "scan_total": 12,
                    "scan_success": 11,
                    "scan_failure": 1,
                    "rho_support_cache": {
                        "point_requests": 20,
                        "cache_hits": 4,
                        "unique_solves": 16,
                        "targeted_additions": 3,
                        "failed_points": 1,
                    },
                    "sweep_records": [
                        {
                            "status": "confirmed_first_order",
                            "stage_used": "stage_b_dense",
                            "hybrid_certificate_type": "endpoint_limited_first_order",
                            "geometry_missing_fields": [],
                            "geometry_normalization_version": "nothing_to_inf_v1",
                        },
                        {
                            "status": "ambiguous_near_critical",
                            "stage_used": "stage_c_local_oracle",
                            "hybrid_certificate_type": "none",
                            "geometry_missing_fields": ["density_error"],
                            "geometry_normalization_version": "nothing_to_inf_v1",
                        },
                    ],
                    "rho_unconverged_count": 2,
                    "temperature_unconverged_count": 1,
                },
            }
        ),
        encoding="utf-8",
    )
    payload = module.build_payload(args_for(module, tmp_path / "process", tmp_path / "diag.json"))
    assert payload["availability"] == "available"
    assert payload["config_hash"] == "config-hash"
    assert payload["diagnostics"]["counters"]["unique_solves"] == 16
    assert payload["diagnostics"]["status_counts"] == {
        "ambiguous_near_critical": 1,
        "confirmed_first_order": 1,
    }
    assert payload["diagnostics"]["geometry_missing_record_count"] == 1
    assert payload["missing_fields"] == []
