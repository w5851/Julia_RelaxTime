#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Any


MAX_ACTION_XI_REFINE_LEVEL = 3
IGNORED_CONFIG_FIELDS = {"xi_values", "requested_xi_values"}


def fail(message: str) -> None:
    raise ValueError(message)


def parse_xi_list(raw: str) -> list[float]:
    values = sorted({_xi_key(float(token.strip())) for token in raw.split(",") if token.strip()})
    if not values:
        fail("expected xi list is empty")
    if len(values) < 2:
        fail("expected xi list must contain at least two anchors")
    if values[0] <= -1 or any(not math.isfinite(value) for value in values):
        fail("expected xi values must be finite and satisfy xi > -1")
    return values


def _xi_key(value: float) -> float:
    return round(value, 12)


def _xi_text(value: float) -> str:
    rendered = format(_xi_key(value), ".12g")
    return "0" if rendered == "-0" else rendered


def _number(config: dict[str, Any], key: str) -> int | float:
    value = config.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        fail(f"source config field {key} must be numeric")
    number = float(value)
    if not math.isfinite(number):
        fail(f"source config field {key} must be finite")
    return value


def _integer(config: dict[str, Any], key: str) -> int:
    value = _number(config, key)
    if int(value) != value:
        fail(f"source config field {key} must be an integer")
    return int(value)


def _boolean(config: dict[str, Any], key: str) -> bool:
    value = config.get(key)
    if not isinstance(value, bool):
        fail(f"source config field {key} must be boolean")
    return value


def _text(value: Any) -> str:
    if isinstance(value, bool):
        fail("boolean cannot be rendered as a numeric CLI value")
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            fail("non-finite CLI value")
        return format(value, ".17g")
    return str(value)


def normalize_config(config: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in config.items() if key not in IGNORED_CONFIG_FIELDS}


def common_cli_args(config: dict[str, Any]) -> list[str]:
    if config.get("mode") != "production":
        fail("source shard mode must be production")
    if not _boolean(config, "compute_crossover"):
        fail("source shard must include crossover results")
    if config.get("crossover_method") != "peak" or config.get("crossover_variable") != "phi_u":
        fail("source crossover method/variable is not supported by staged resume")

    args = [
        "--T-min", _text(_number(config, "T_min_MeV")),
        "--T-max", _text(_number(config, "T_max_MeV")),
        "--T-step", _text(_number(config, "T_step_MeV")),
        "--T-refine-levels", _text(_integer(config, "temperature_max_refine_level")),
        "--T-position-tol", _text(_number(config, "temperature_position_tol_MeV")),
        "--T-density-tol", _text(_number(config, "temperature_density_tol")),
        "--T-maxwell-area-tol", _text(_number(config, "temperature_maxwell_area_tol")),
        "--rho-min", _text(_number(config, "rho_min")),
        "--rho-max", _text(_number(config, "rho_max")),
        "--rho-step", _text(_number(config, "rho_step")),
        "--rho-position-tol", _text(_number(config, "rho_position_tol_MeV")),
        "--rho-density-tol", _text(_number(config, "rho_density_tol")),
        "--rho-maxwell-area-tol", _text(_number(config, "rho_maxwell_area_tol")),
        "--cep-tol", _text(_number(config, "cep_tol_MeV")),
        "--p-num", _text(_integer(config, "p_num")),
        "--t-num", _text(_integer(config, "t_num")),
        "--thermo-quadrature-policy", str(config["thermo_quadrature_policy"]),
        "--thermo-quadrature-rtol", _text(_number(config, "thermo_quadrature_rtol")),
        "--thermo-quadrature-atol", _text(_number(config, "thermo_quadrature_atol")),
        "--thermo-quadrature-maxevals", _text(_integer(config, "thermo_quadrature_maxevals")),
        "--iterations", _text(_integer(config, "iterations")),
        "--crossover-n-mu", _text(_integer(config, "crossover_n_mu")),
        "--crossover-mu-max", _text(_number(config, "crossover_mu_max_MeV")),
    ]
    # Preserve rho-policy provenance when a staged resume is reconstructed
    # from a source manifest.  Older manifests predate these fields, so they
    # intentionally retain the historical uniform defaults.
    if "rho_refinement_policy" in config:
        rho_policy = str(config["rho_refinement_policy"])
        if rho_policy not in {"uniform_nested", "rho_support_cascade", "rho_support_hybrid"}:
            fail(f"unsupported source rho_refinement_policy: {rho_policy}")
        args.extend(["--rho-refinement-policy", rho_policy])
        if "rho_refine_levels" in config:
            args.extend(["--rho-refine-levels", _text(_integer(config, "rho_refine_levels"))])
        if "rho_support_fine_step" in config:
            args.extend(["--rho-support-fine-step", _text(_number(config, "rho_support_fine_step"))])
        if "rho_support_target_point_count" in config:
            args.extend(["--rho-support-target-point-count", _text(_integer(config, "rho_support_target_point_count"))])
        if "rho_support_targeted_cap" in config:
            args.extend(["--rho-support-targeted-cap", _text(_integer(config, "rho_support_targeted_cap"))])
        endpoint_policy = config.get("rho_hybrid_endpoint_policy")
        if endpoint_policy is not None:
            if endpoint_policy not in {"bounded_zero_density_v1", "three_crossing_endpoint_local_v2"}:
                fail(f"unsupported source rho_hybrid_endpoint_policy: {endpoint_policy}")
            args.extend(["--rho-hybrid-endpoint-policy", str(endpoint_policy)])
    if _boolean(config, "adaptive_xi"):
        args.extend([
            "--adaptive-xi",
            "--xi-refine-levels", _text(_integer(config, "xi_max_refine_level")),
            "--xi-position-tol", _text(_number(config, "xi_position_tol_MeV")),
            "--xi-density-tol", _text(_number(config, "xi_density_tol")),
            "--xi-maxwell-area-tol", _text(_number(config, "xi_maxwell_area_tol")),
            "--xi-response-rtol", _text(_number(config, "xi_response_rtol")),
        ])
    if not _boolean(config, "adaptive_temperature"):
        args.append("--no-adaptive-T")
    if not _boolean(config, "rho_geometry_convergence"):
        args.append("--no-rho-geometry-convergence")
    if _boolean(config, "crossover_only"):
        args.append("--crossover-only")
    if _boolean(config, "crossover_mu0_only"):
        args.append("--crossover-mu0-only")
    crossover_T_max = config.get("crossover_T_max_MeV")
    if crossover_T_max is not None:
        args.extend(["--crossover-T-max", _text(_number(config, "crossover_T_max_MeV"))])
    args.append("--overwrite")
    return args


def build_resume_plan(
    shards_root: Path,
    tag: str,
    source_calculation_sha: str,
    expected_xi_list: str,
) -> dict[str, Any]:
    manifests = sorted(shards_root.rglob(f"phase_reference_{tag}_manifest.json"))
    if not manifests:
        fail(f"no source shard manifests found below {shards_root}")

    payloads = [json.loads(path.read_text(encoding="utf-8")) for path in manifests]
    commits = {payload.get("generator", {}).get("git_commit") for payload in payloads}
    if commits != {source_calculation_sha}:
        fail(f"source shard calculation commits differ from {source_calculation_sha}: {sorted(map(str, commits))}")

    configs = [payload.get("config") for payload in payloads]
    if any(not isinstance(config, dict) for config in configs):
        fail("source shard manifest lacks a config object")
    normalized = [normalize_config(config) for config in configs]
    if any(config != normalized[0] for config in normalized[1:]):
        fail("source shard configuration mismatch outside xi sampling fields")
    config = normalized[0]
    if config.get("tag") != tag:
        fail(f"source shard tag mismatch: expected {tag}, got {config.get('tag')}")
    if not _boolean(config, "adaptive_xi"):
        fail("source shard config must enable adaptive_xi")
    xi_refine_levels = _integer(config, "xi_max_refine_level")
    if not 1 <= xi_refine_levels <= MAX_ACTION_XI_REFINE_LEVEL:
        fail(f"source xi refinement level must be in [1, {MAX_ACTION_XI_REFINE_LEVEL}]")

    anchors = parse_xi_list(expected_xi_list)
    intervals = list(zip(anchors, anchors[1:]))
    level1 = [_xi_key(0.5 * (left + right)) for left, right in intervals]
    expected_initial = sorted(set(anchors + level1))
    resolved: list[float] = []
    for source_config in configs:
        values = source_config.get("xi_values")
        if not isinstance(values, list) or len(values) != 1:
            fail("each source shard must contain exactly one xi value")
        resolved.append(_xi_key(float(values[0])))
    if len(resolved) != len(set(resolved)):
        fail("source run contains duplicate xi shard manifests")
    missing = sorted(set(expected_initial) - set(resolved))
    unexpected = sorted(set(resolved) - set(expected_initial))
    if missing or unexpected:
        fail(f"source run is not an initial staged grid: missing={missing}, unexpected={unexpected}")

    config_digest = hashlib.sha256(
        (json.dumps(config, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")
    ).hexdigest()
    return {
        "expected_xi_csv": ",".join(_xi_text(value) for value in anchors),
        "initial_intervals": [[_xi_text(left), _xi_text(right)] for left, right in intervals],
        "source_shard_count": len(manifests),
        "source_xi_values": sorted(resolved),
        "source_config_sha256": config_digest,
        "xi_refine_levels": xi_refine_levels,
        "xi_position_tol": _text(_number(config, "xi_position_tol_MeV")),
        "xi_density_tol": _text(_number(config, "xi_density_tol")),
        "xi_maxwell_area_tol": _text(_number(config, "xi_maxwell_area_tol")),
        "xi_response_rtol": _text(_number(config, "xi_response_rtol")),
        "common_args": common_cli_args(config),
        "crossover_only": _boolean(config, "crossover_only"),
        "crossover_mu0_only": _boolean(config, "crossover_mu0_only"),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plan a PNJL dense-reference staged source-run resume.")
    parser.add_argument("--shards-root", type=Path, required=True)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--source-calculation-sha", required=True)
    parser.add_argument("--expected-xi-list", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        plan = build_resume_plan(
            args.shards_root,
            args.tag,
            args.source_calculation_sha,
            args.expected_xi_list,
        )
    except (ValueError, json.JSONDecodeError) as exc:
        raise SystemExit(f"[dense-reference-resume-plan] {exc}") from exc
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
