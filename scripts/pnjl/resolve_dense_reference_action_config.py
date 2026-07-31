#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


FIELD_SPECS: dict[str, tuple[str, str, bool]] = {
    "xi_refine_levels": ("--xi-refine-levels", "int", True),
    "xi_position_tol": ("--xi-position-tol", "float", False),
    "xi_density_tol": ("--xi-density-tol", "float", False),
    "xi_maxwell_area_tol": ("--xi-maxwell-area-tol", "float", False),
    "xi_response_rtol": ("--xi-response-rtol", "float", False),
    "T_refine_levels": ("--T-refine-levels", "int", True),
    "T_position_tol": ("--T-position-tol", "float", False),
    "T_density_tol": ("--T-density-tol", "float", False),
    "T_maxwell_area_tol": ("--T-maxwell-area-tol", "float", False),
    "rho_position_tol": ("--rho-position-tol", "float", False),
    "rho_density_tol": ("--rho-density-tol", "float", False),
    "rho_maxwell_area_tol": ("--rho-maxwell-area-tol", "float", False),
    "cep_tol": ("--cep-tol", "float", False),
    "thermo_quadrature_rtol": ("--thermo-quadrature-rtol", "float", False),
    "thermo_quadrature_atol": ("--thermo-quadrature-atol", "float", False),
    "thermo_quadrature_maxevals": ("--thermo-quadrature-maxevals", "int", False),
    "crossover_T_max": ("--crossover-T-max", "float", False),
    "rho_refinement_policy": ("--rho-refinement-policy", "enum", False),
    "rho_refine_levels": ("--rho-refine-levels", "int", True),
    "rho_support_fine_step": ("--rho-support-fine-step", "float", False),
    "rho_support_target_point_count": ("--rho-support-target-point-count", "int", False),
    "rho_support_targeted_cap": ("--rho-support-targeted-cap", "int", False),
}


def fail(message: str) -> None:
    raise ValueError(message)


def _normalized_value(key: str, value: Any, kind: str, allow_zero: bool) -> str:
    if kind == "enum":
        if not isinstance(value, str) or value not in {
            "uniform_nested",
            "rho_support_cascade",
            "rho_support_hybrid",
        }:
            fail(f"{key} must be uniform_nested, rho_support_cascade, or rho_support_hybrid")
        return value
    if isinstance(value, bool):
        fail(f"{key} must be a numeric value, not boolean")
    if kind == "int":
        if not isinstance(value, int):
            fail(f"{key} must be an integer")
        if value < 0 or (value == 0 and not allow_zero):
            qualifier = "nonnegative" if allow_zero else "positive"
            fail(f"{key} must be {qualifier}")
        return str(value)
    if not isinstance(value, (int, float)):
        fail(f"{key} must be numeric")
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        fail(f"{key} must be finite and positive")
    return format(number, ".17g")


def resolve_action_config(raw: str) -> list[str]:
    try:
        payload = json.loads(raw or "{}")
    except json.JSONDecodeError as exc:
        fail(f"advanced_config_json is not valid JSON: {exc}")
    if not isinstance(payload, dict):
        fail("advanced_config_json must be a JSON object")
    unknown = sorted(set(payload) - set(FIELD_SPECS))
    if unknown:
        fail(f"unsupported advanced_config_json keys: {unknown}")

    args: list[str] = []
    for key in FIELD_SPECS:
        if key not in payload:
            continue
        flag, kind, allow_zero = FIELD_SPECS[key]
        args.extend([flag, _normalized_value(key, payload[key], kind, allow_zero)])
    return args


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Resolve whitelisted PNJL dense-reference Action controls.")
    parser.add_argument("--config-json", default="{}")
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        resolved = resolve_action_config(args.config_json)
    except ValueError as exc:
        raise SystemExit(f"[dense-reference-action-config] {exc}") from exc
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("".join(f"{token}\n" for token in resolved), encoding="utf-8")


if __name__ == "__main__":
    main()
