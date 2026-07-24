#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import sys
from decimal import Decimal
from pathlib import Path
from typing import Any


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from resolve_dense_reference_action_config import resolve_action_config  # noqa: E402


MAX_ACTION_XI_REFINE_LEVEL = 3
FULL_GIT_SHA_PATTERN = re.compile(r"[0-9a-f]{40}")


def fail(message: str) -> None:
    raise ValueError(message)


def _value(payload: dict[str, Any], key: str, default: Any) -> Any:
    value = payload.get(key, default)
    return default if value is None or value == "" else value


def _boolean(payload: dict[str, Any], key: str, default: bool) -> bool:
    value = _value(payload, key, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.lower() in {"true", "false"}:
        return value.lower() == "true"
    fail(f"{key} must be boolean")


def _validate_calculation_ref(payload: dict[str, Any]) -> None:
    calculation_ref = str(payload.get("calculation_ref", "") or "").strip()
    if calculation_ref and FULL_GIT_SHA_PATTERN.fullmatch(calculation_ref) is None:
        fail("calculation_ref must be an immutable lowercase 40-character Git SHA")


def _decimal_text(value: Decimal) -> str:
    rendered = format(value, "f")
    return rendered.rstrip("0").rstrip(".") if "." in rendered else rendered


def _xi_anchors(payload: dict[str, Any]) -> list[Decimal]:
    raw = str(_value(payload, "xi_values", "")).strip()
    if raw:
        values = sorted({Decimal(token.strip()) for token in raw.split(",") if token.strip()})
    else:
        start = Decimal(str(_value(payload, "xi_min", "-0.5")))
        stop = Decimal(str(_value(payload, "xi_max", "0.5")))
        step = Decimal(str(_value(payload, "xi_step", "0.05")))
        if step <= 0 or stop < start:
            fail("invalid xi range")
        values = []
        current = start
        while current <= stop:
            values.append(current)
            current += step
        if values[-1] != stop:
            values.append(stop)
    if not values:
        fail("empty xi plan")
    if values[0] <= Decimal("-1"):
        fail("xi values must satisfy xi > -1")
    return values


def _advanced_payload(payload: dict[str, Any]) -> dict[str, Any]:
    raw = str(_value(payload, "advanced_config_json", "{}"))
    # Resolve first so this planner shares the workflow whitelist and numeric validation.
    resolve_action_config(raw)
    advanced = json.loads(raw or "{}")
    if not isinstance(advanced, dict):
        fail("advanced_config_json must be a JSON object")
    return advanced


def _numeric_text(value: Any) -> str:
    if isinstance(value, bool):
        fail("numeric workflow value cannot be boolean")
    return str(value)


def _common_cli_args(
    payload: dict[str, Any],
    advanced_args: list[str],
    *,
    adaptive_xi: bool,
) -> list[str]:
    args = [
        "--T-min", _numeric_text(_value(payload, "T_min", "60")),
        "--T-max", _numeric_text(_value(payload, "T_max", "240")),
        "--T-step", _numeric_text(_value(payload, "T_step", "5")),
        "--T-refine-levels", "2",
        "--T-position-tol", "0.10",
        "--T-density-tol", "0.01",
        "--T-maxwell-area-tol", "1e-4",
        "--rho-min", _numeric_text(_value(payload, "rho_min", "0.0")),
        "--rho-max", _numeric_text(_value(payload, "rho_max", "4.0")),
        "--rho-step", _numeric_text(_value(payload, "rho_step", "0.05")),
        "--rho-position-tol", "0.05",
        "--rho-density-tol", "0.005",
        "--rho-maxwell-area-tol", "1e-4",
        "--cep-tol", "0.1",
        "--p-num", _numeric_text(_value(payload, "p_num", "24")),
        "--t-num", _numeric_text(_value(payload, "t_num", "8")),
        "--thermo-quadrature-policy", str(_value(payload, "thermo_quadrature_policy", "tensor_gauss")),
        "--thermo-quadrature-rtol", "1e-8",
        "--thermo-quadrature-atol", "1e-10",
        "--thermo-quadrature-maxevals", "10000000",
        "--iterations", _numeric_text(_value(payload, "iterations", "80")),
        "--crossover-n-mu", _numeric_text(_value(payload, "crossover_n_mu", "16")),
        "--crossover-mu-max", _numeric_text(_value(payload, "crossover_mu_max", "450")),
    ]
    if adaptive_xi:
        args.extend([
            "--adaptive-xi",
            "--xi-refine-levels", "2",
            "--xi-position-tol", "0.10",
            "--xi-density-tol", "0.01",
            "--xi-maxwell-area-tol", "1e-4",
            "--xi-response-rtol", "0.05",
        ])
    if not _boolean(payload, "adaptive_T", True):
        args.append("--no-adaptive-T")
    if not _boolean(payload, "rho_geometry_convergence", True):
        args.append("--no-rho-geometry-convergence")
    if _boolean(payload, "crossover_only", True):
        args.append("--crossover-only")
    if _boolean(payload, "crossover_mu0_only", False):
        args.append("--crossover-mu0-only")
    if _boolean(payload, "overwrite", True):
        args.append("--overwrite")
    args.extend(advanced_args)
    return args


def build_action_plan(payload: dict[str, Any]) -> dict[str, Any]:
    _validate_calculation_ref(payload)
    anchors = _xi_anchors(payload)
    adaptive_xi = _boolean(payload, "adaptive_xi", False)
    crossover_only = _boolean(payload, "crossover_only", True)
    if adaptive_xi and crossover_only:
        fail("adaptive_xi requires crossover_only=false")

    advanced = _advanced_payload(payload)
    xi_refine_levels = int(advanced.get("xi_refine_levels", 2)) if adaptive_xi else 0
    if xi_refine_levels > MAX_ACTION_XI_REFINE_LEVEL:
        fail(
            f"GitHub Actions staged xi refinement supports at most "
            f"{MAX_ACTION_XI_REFINE_LEVEL} levels, got {xi_refine_levels}"
        )

    intervals = list(zip(anchors, anchors[1:]))
    initial_entries = [
        {"stage": "anchor", "shard_id": f"a{index:03d}", "xi": _decimal_text(value)}
        for index, value in enumerate(anchors)
    ]
    if xi_refine_levels >= 1:
        midpoints = sorted({(left + right) / Decimal(2) for left, right in intervals})
        initial_entries.extend(
            {"stage": "level1", "shard_id": f"l1_{index:03d}", "xi": _decimal_text(value)}
            for index, value in enumerate(midpoints)
        )

    advanced_raw = str(_value(payload, "advanced_config_json", "{}"))
    advanced_args = resolve_action_config(advanced_raw)
    common_args = _common_cli_args(payload, advanced_args, adaptive_xi=adaptive_xi)

    def effective(name: str, default: Any) -> str:
        return _numeric_text(advanced.get(name, default))

    return {
        "initial_matrix": {"include": initial_entries},
        "expected_xi_csv": ",".join(_decimal_text(value) for value in anchors),
        "initial_intervals": [[_decimal_text(left), _decimal_text(right)] for left, right in intervals],
        "initial_shard_count": len(initial_entries),
        "xi_refine_levels": xi_refine_levels,
        "xi_position_tol": effective("xi_position_tol", 0.10),
        "xi_density_tol": effective("xi_density_tol", 0.01),
        "xi_maxwell_area_tol": effective("xi_maxwell_area_tol", 1e-4),
        "xi_response_rtol": effective("xi_response_rtol", 0.05),
        "common_args": common_args,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plan staged one-xi PNJL dense-reference Action jobs.")
    parser.add_argument("--inputs-json", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        payload = json.loads(args.inputs_json)
        if not isinstance(payload, dict):
            fail("workflow inputs must be a JSON object")
        plan = build_action_plan(payload)
    except (ValueError, json.JSONDecodeError) as exc:
        raise SystemExit(f"[dense-reference-action-plan] {exc}") from exc
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
