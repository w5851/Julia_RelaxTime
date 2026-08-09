#!/usr/bin/env python3

import argparse
import csv
import json
import math
from pathlib import Path


ARTIFACT_FILES = {
    "boundary": "boundary_{tag}.csv",
    "cep": "cep_{tag}.csv",
    "spinodals": "spinodals_{tag}.csv",
    "crossover": "crossover_{tag}.csv",
}

OPTIONAL_ARTIFACT_FILES = {
    "grid_convergence": "phase_grid_convergence_{tag}.csv",
}

# The public comparison grid is intentionally independent of whichever
# adaptive xi plan a particular C1/C2 run used.  ``1:5:236`` is Julia's
# inclusive range (1, 6, ..., 236), with the explicit 240 MeV endpoint.
PUBLIC_XI_VALUES = tuple(round(-0.5 + 0.05 * index, 10) for index in range(21))
PUBLIC_T_VALUES = tuple(list(range(1, 237, 5)) + [240])

GEOMETRY_TOLERANCES = {
    "mu_transition_MeV": 0.025,
    "rho_hadron": 0.0025,
    "rho_quark": 0.0025,
    "area_residual": 5.0e-5,
    "rho_spinodal_hadron": 0.0025,
    "rho_spinodal_quark": 0.0025,
}
CEP_BRACKET_TOLERANCE_MEV = 0.1
CROSSOVER_TOLERANCES = {
    "T_crossover_MeV": 0.05,
    "rho": 0.005,
    # The derivative field is a response-like quantity.  Its gate is
    # relative and is applied after interpolation on physical mu points.
    "derivative_rel": 0.025,
}

STATE_FIRST_ORDER = "confirmed_first_order"
STATE_MONOTONE = "confirmed_monotone"
STATE_AMBIGUOUS = "ambiguous_near_critical"
STATE_NOT_EVALUATED = "not_evaluated"

COMPARE_COLUMNS = {
    "boundary": ["mu_transition_MeV", "rho_hadron", "rho_quark", "area_residual"],
    "cep": ["T_CEP_MeV", "muq_CEP_MeV", "muB_CEP_MeV"],
    "spinodals": [
        "mu_spinodal_hadron_MeV",
        "mu_spinodal_quark_MeV",
        "rho_spinodal_hadron",
        "rho_spinodal_quark",
    ],
    "crossover": ["mu_MeV", "T_crossover_MeV", "rho", "derivative"],
}

CEP_NUMERIC_COLUMNS = [
    "T_CEP_MeV",
    "muq_CEP_MeV",
    "muB_CEP_MeV",
    "T_bracket_low_MeV",
    "T_bracket_high_MeV",
    "bracket_width_T_MeV",
    "T_last_first_order_MeV",
    "muq_last_first_order_MeV",
    "muB_last_first_order_MeV",
    "T_first_monotone_MeV",
    "ambiguity_width_T_MeV",
    "temperature_resolution_target_MeV",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare PNJL dense phase reference artifacts for convergence audits."
    )
    parser.add_argument("--candidate-root", required=True, type=Path)
    parser.add_argument("--candidate-tag", required=True)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--reference-tag", required=True)
    parser.add_argument("--candidate-label", default="candidate")
    parser.add_argument("--reference-label", default="reference")
    parser.add_argument(
        "--xi-values",
        help="Optional comma-separated xi anchors to compare; rows outside this set are ignored.",
    )
    parser.add_argument(
        "--public-grid",
        action="store_true",
        help="Also emit the fixed public xi/T anchor state and geometry tables.",
    )
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


def fail(message):
    raise SystemExit(f"[phase-reference-convergence] {message}")


def parse_float_set(raw):
    if raw is None or not raw.strip():
        return None
    values = set()
    for token in raw.split(","):
        token = token.strip()
        if token:
            values.add(round(float(token), 10))
    if not values:
        return None
    return values


def read_csv(path):
    if not path.is_file():
        fail(f"missing csv: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames or [], list(reader)


def try_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def finite_float(row, column):
    value = try_float(row.get(column))
    if value is None or not math.isfinite(value):
        return None
    return value


def key_float(row, column, ndigits=10):
    value = finite_float(row, column)
    if value is None:
        return None
    return round(value, ndigits)


def crossover_sort_key(row):
    plot_order_key = finite_float(row, "plot_order_key")
    if plot_order_key is not None:
        return plot_order_key
    mu_value = finite_float(row, "mu_MeV")
    if mu_value is not None:
        return mu_value
    return math.inf


def summary_sort_key(item):
    return (item["artifact"], item["metric"])


def artifact_path(root, tag, artifact):
    return root / ARTIFACT_FILES[artifact].format(tag=tag)


def filter_rows_by_xi(rows, xi_filter):
    if xi_filter is None:
        return rows
    filtered = []
    for row in rows:
        xi_value = key_float(row, "xi")
        if xi_value in xi_filter:
            filtered.append(row)
    return filtered


def load_artifacts(root, tag, xi_filter=None):
    loaded = {}
    for artifact in ARTIFACT_FILES:
        fields, rows = read_csv(artifact_path(root, tag, artifact))
        rows = filter_rows_by_xi(rows, xi_filter)
        loaded[artifact] = {"fields": fields, "rows": rows}
    return loaded


def load_optional_artifacts(root, tag, xi_filter=None):
    loaded = {}
    for artifact, template in OPTIONAL_ARTIFACT_FILES.items():
        path = root / template.format(tag=tag)
        if not path.is_file():
            loaded[artifact] = {"fields": [], "rows": [], "available": False, "path": path.as_posix()}
            continue
        fields, rows = read_csv(path)
        loaded[artifact] = {
            "fields": fields,
            "rows": filter_rows_by_xi(rows, xi_filter),
            "available": True,
            "path": path.as_posix(),
        }
    return loaded


def artifact_key(artifact, row):
    if artifact == "cep":
        return (key_float(row, "xi"),)
    if artifact in ("boundary", "spinodals", "grid_convergence"):
        if artifact == "grid_convergence":
            return (
                str(row.get("axis", "")),
                key_float(row, "xi"),
                key_float(row, "T_MeV"),
                key_float(row, "level"),
                key_float(row, "left"),
                key_float(row, "right"),
                key_float(row, "midpoint"),
            )
        return (key_float(row, "xi"), key_float(row, "T_MeV"))
    if artifact == "crossover":
        return (key_float(row, "xi"), key_float(row, "mu_MeV"))
    return tuple()


def summarize_artifact(artifact, fields, rows):
    xi_values = sorted({finite_float(row, "xi") for row in rows if finite_float(row, "xi") is not None})
    naninf_count = 0
    negative_by_column = {}
    for row in rows:
        for column, raw in row.items():
            value = try_float(raw)
            if value is None:
                continue
            if not math.isfinite(value):
                naninf_count += 1
            if value < 0 and column != "xi":
                negative_by_column[column] = negative_by_column.get(column, 0) + 1
    converged_false = None
    if "converged" in fields:
        converged_false = sum(1 for row in rows if str(row.get("converged", "")).lower() != "true")
    keys = [artifact_key(artifact, row) for row in rows]
    duplicate_key_count = len(keys) - len(set(keys))
    required_columns = {
        "boundary": {"xi", "T_MeV", "mu_transition_MeV", "rho_hadron", "rho_quark", "area_residual"},
        "cep": {"xi", "result_status", "T_bracket_low_MeV", "T_bracket_high_MeV", "bracket_width_T_MeV", "ambiguity_width_T_MeV"},
        "spinodals": {"xi", "T_MeV", "mu_spinodal_hadron_MeV", "mu_spinodal_quark_MeV", "rho_spinodal_hadron", "rho_spinodal_quark"},
        "crossover": {"xi", "mu_MeV", "T_crossover_MeV", "rho", "derivative"},
        "grid_convergence": {"axis", "xi", "T_MeV", "converged", "reason"},
    }
    missing_columns = sorted(required_columns.get(artifact, set()) - set(fields))
    return {
        "artifact": artifact,
        "row_count": len(rows),
        "columns": fields,
        "xi_count": len(xi_values),
        "xi_values": xi_values,
        "naninf_count": naninf_count,
        "negative_by_column": negative_by_column,
        "converged_false_count": converged_false,
        "duplicate_key_count": duplicate_key_count,
        "missing_columns": missing_columns,
    }


def keyed_rows(artifact, rows):
    if artifact == "cep":
        return {(key_float(row, "xi"),): row for row in rows}
    if artifact in ("boundary", "spinodals"):
        return {(key_float(row, "xi"), key_float(row, "T_MeV")): row for row in rows}
    if artifact == "crossover":
        return {(key_float(row, "xi"), key_float(row, "mu_MeV")): row for row in rows}
    fail(f"unknown artifact: {artifact}")


def rel_diff(candidate, reference):
    denom = max(abs(reference), 1.0e-12)
    return abs(candidate - reference) / denom


def append_missing_rows(comparison_rows, artifact, keys, side, *, adaptive_xi=False, grid_change_kind=None):
    for key in sorted(keys):
        comparison_rows.append(
            {
                "artifact": artifact,
                "match_status": "adaptive_xi_added" if adaptive_xi else f"missing_in_{side}",
                "grid_change_kind": grid_change_kind or ("new_xi" if adaptive_xi else "missing_key"),
                "xi": key[0] if key else "",
                "match_key": "|".join(str(part) for part in key),
                "metric": "",
                "candidate_value": "",
                "reference_value": "",
                "abs_diff": "",
                "rel_diff": "",
            }
        )


def _append_key_gaps(comparison_rows, artifact, candidate_keys, reference_keys):
    # Keep xi-plan changes separate from a changed T grid on a shared xi.  A
    # candidate-only xi is an adaptive addition; a candidate-only T is a
    # denser/replaced temperature grid and must not be counted as an xi gap.
    candidate_only = candidate_keys - reference_keys
    reference_only = reference_keys - candidate_keys
    candidate_xi = {key[0] for key in candidate_only}
    reference_xi = {key[0] for key in reference_only}
    shared_xi = candidate_xi & reference_xi
    candidate_xi_only = {key for key in candidate_only if key[0] not in reference_xi}
    reference_xi_only = {key for key in reference_only if key[0] not in candidate_xi}
    candidate_t_added = {key for key in candidate_only if key[0] in shared_xi}
    reference_t_missing = {key for key in reference_only if key[0] in shared_xi}
    append_missing_rows(
        comparison_rows,
        artifact,
        candidate_xi_only,
        "reference",
        adaptive_xi=True,
        grid_change_kind="new_xi",
    )
    append_missing_rows(
        comparison_rows,
        artifact,
        reference_xi_only,
        "candidate",
        grid_change_kind="reference_xi_missing",
    )
    append_missing_rows(
        comparison_rows,
        artifact,
        candidate_t_added,
        "reference",
        grid_change_kind="shared_xi_new_T",
    )
    append_missing_rows(
        comparison_rows,
        artifact,
        reference_t_missing,
        "candidate",
        grid_change_kind="shared_xi_missing_T",
    )


def _value_text(row, column):
    value = row.get(column, "")
    return "" if value is None else str(value).strip()


def compare_cep_artifact(candidate_rows, reference_rows):
    candidate_keyed = keyed_rows("cep", candidate_rows)
    reference_keyed = keyed_rows("cep", reference_rows)
    candidate_keys = set(candidate_keyed)
    reference_keys = set(reference_keyed)
    comparison_rows = []
    _append_key_gaps(comparison_rows, "cep", candidate_keys, reference_keys)

    for key in sorted(candidate_keys & reference_keys):
        candidate = candidate_keyed[key]
        reference = reference_keyed[key]
        candidate_status = _value_text(candidate, "result_status")
        reference_status = _value_text(reference, "result_status")
        comparison_rows.append(
            {
                "artifact": "cep",
                "match_status": "matched" if candidate_status == reference_status else "status_changed",
                "xi": key[0],
                "match_key": "|".join(str(part) for part in key),
                "metric": "result_status",
                "candidate_value": candidate_status,
                "reference_value": reference_status,
                "abs_diff": "",
                "rel_diff": "",
            }
        )
        for column in CEP_NUMERIC_COLUMNS:
            candidate_text = _value_text(candidate, column)
            reference_text = _value_text(reference, column)
            if not candidate_text and not reference_text:
                comparison_rows.append(
                    {
                        "artifact": "cep",
                        "match_status": "not_applicable",
                        "xi": key[0],
                        "match_key": "|".join(str(part) for part in key),
                        "metric": column,
                        "candidate_value": "",
                        "reference_value": "",
                        "abs_diff": "",
                        "rel_diff": "",
                    }
                )
                continue
            candidate_value = try_float(candidate_text)
            reference_value = try_float(reference_text)
            if candidate_value is None or reference_value is None or not math.isfinite(candidate_value) or not math.isfinite(reference_value):
                status = "non_numeric"
                abs_diff = ""
                rel_value = ""
            else:
                status = "matched"
                abs_diff = abs(candidate_value - reference_value)
                rel_value = rel_diff(candidate_value, reference_value)
            comparison_rows.append(
                {
                    "artifact": "cep",
                    "match_status": status,
                    "xi": key[0],
                    "match_key": "|".join(str(part) for part in key),
                    "metric": column,
                    "candidate_value": candidate_text,
                    "reference_value": reference_text,
                    "abs_diff": abs_diff,
                    "rel_diff": rel_value,
                }
            )
    return comparison_rows


def _is_physical_crossover_row(row):
    """Return whether a crossover sample belongs to a comparable branch."""
    status = _value_text(row, "phase_status") or _value_text(row, "result_status")
    if status:
        normalized = status.lower()
        if normalized not in {
            "confirmed_monotone",
            "monotone",
            "crossover",
            "confirmed_crossover",
        }:
            return False
    return (
        finite_float(row, "mu_MeV") is not None
        and finite_float(row, "T_crossover_MeV") is not None
        and finite_float(row, "rho") is not None
        and finite_float(row, "derivative") is not None
        and str(row.get("converged", "true")).lower() == "true"
    )


def _group_crossover_rows(rows, min_temperature_by_xi=None):
    grouped = {}
    for row in rows:
        if not _is_physical_crossover_row(row):
            continue
        xi = key_float(row, "xi")
        mu = finite_float(row, "mu_MeV")
        if xi is None or mu is None:
            continue
        if min_temperature_by_xi is not None:
            threshold = min_temperature_by_xi.get(xi)
            temperature = finite_float(row, "T_crossover_MeV")
            if threshold is not None and (temperature is None or temperature < threshold - 1.0e-9):
                continue
        grouped.setdefault(xi, {})[mu] = row
    return grouped


def _interpolate_row(rows_by_mu, mu):
    points = sorted(rows_by_mu)
    if not points or mu < points[0] - 1.0e-10 or mu > points[-1] + 1.0e-10:
        return None, False
    if mu in rows_by_mu:
        row = rows_by_mu[mu]
        return {
            "mu_MeV": mu,
            "T_crossover_MeV": finite_float(row, "T_crossover_MeV"),
            "rho": finite_float(row, "rho"),
            "derivative": finite_float(row, "derivative"),
        }, True
    right_index = next((index for index, value in enumerate(points) if value > mu), len(points))
    if right_index == 0 or right_index == len(points):
        return None, False
    left_mu, right_mu = points[right_index - 1], points[right_index]
    weight = (mu - left_mu) / (right_mu - left_mu)
    left = rows_by_mu[left_mu]
    right = rows_by_mu[right_mu]
    values = {}
    for column in ("T_crossover_MeV", "rho", "derivative"):
        left_value = finite_float(left, column)
        right_value = finite_float(right, column)
        if left_value is None or right_value is None:
            return None, False
        values[column] = left_value + weight * (right_value - left_value)
    values["mu_MeV"] = mu
    return values, False


def _crossover_thresholds(cep_rows):
    thresholds = {}
    for row in cep_rows:
        xi = key_float(row, "xi")
        if xi is None:
            continue
        threshold = finite_float(row, "T_first_monotone_MeV")
        if threshold is None:
            threshold = finite_float(row, "T_bracket_high_MeV")
        if threshold is not None:
            thresholds[xi] = threshold
    return thresholds


def compare_crossover_artifact(candidate_rows, reference_rows, candidate_cep_rows=None, reference_cep_rows=None):
    """Compare crossover curves on the union of native physical-mu points.

    The old implementation sorted rows and paired ordinal N with ordinal N.
    Adaptive runs do not share those nodes, so this function interpolates each
    side on the union inside the common physical interval and performs the
    comparison in both directions without extrapolation.
    """
    candidate_thresholds = _crossover_thresholds(candidate_cep_rows or [])
    reference_thresholds = _crossover_thresholds(reference_cep_rows or [])
    candidate_grouped = _group_crossover_rows(candidate_rows, candidate_thresholds or None)
    reference_grouped = _group_crossover_rows(reference_rows, reference_thresholds or None)
    comparison_rows = []
    common_xi = sorted(set(candidate_grouped) & set(reference_grouped))
    for xi in common_xi:
        candidate_mu = sorted(candidate_grouped[xi])
        reference_mu = sorted(reference_grouped[xi])
        lower = max(candidate_mu[0], reference_mu[0])
        upper = min(candidate_mu[-1], reference_mu[-1])
        union_mu = sorted(
            {
                mu
                for mu in candidate_mu + reference_mu
                if lower - 1.0e-10 <= mu <= upper + 1.0e-10
            }
        )
        for mu in union_mu:
            candidate, candidate_native = _interpolate_row(candidate_grouped[xi], mu)
            reference, reference_native = _interpolate_row(reference_grouped[xi], mu)
            if candidate is None or reference is None:
                continue
            status = "matched" if candidate_native and reference_native else "interpolated"
            for column in ("T_crossover_MeV", "rho", "derivative"):
                candidate_value = candidate[column]
                reference_value = reference[column]
                if column == "derivative":
                    abs_diff = abs(candidate_value - reference_value)
                    rel_value = rel_diff(candidate_value, reference_value)
                else:
                    abs_diff = abs(candidate_value - reference_value)
                    rel_value = rel_diff(candidate_value, reference_value)
                comparison_rows.append(
                    {
                        "artifact": "crossover",
                        "match_status": status,
                        "grid_change_kind": "physical_mu_union",
                        "xi": xi,
                        "match_key": f"{xi}|{mu:.12g}",
                        "metric": column,
                        "candidate_value": candidate_value,
                        "reference_value": reference_value,
                        "abs_diff": abs_diff,
                        "rel_diff": rel_value,
                    }
                )
    return comparison_rows


def compare_artifact(artifact, candidate_rows, reference_rows, *, candidate_cep_rows=None, reference_cep_rows=None):
    if artifact == "cep":
        return compare_cep_artifact(candidate_rows, reference_rows)
    if artifact == "crossover":
        return compare_crossover_artifact(
            candidate_rows,
            reference_rows,
            candidate_cep_rows=candidate_cep_rows,
            reference_cep_rows=reference_cep_rows,
        )
    candidate_keyed = keyed_rows(artifact, candidate_rows)
    reference_keyed = keyed_rows(artifact, reference_rows)
    candidate_keys = set(candidate_keyed)
    reference_keys = set(reference_keyed)
    shared_keys = sorted(candidate_keys & reference_keys)
    comparison_rows = []

    _append_key_gaps(comparison_rows, artifact, candidate_keys, reference_keys)

    for key in shared_keys:
        candidate = candidate_keyed[key]
        reference = reference_keyed[key]
        for column in COMPARE_COLUMNS[artifact]:
            candidate_value = finite_float(candidate, column)
            reference_value = finite_float(reference, column)
            if candidate_value is None or reference_value is None:
                comparison_rows.append(
                    {
                        "artifact": artifact,
                        "match_status": "non_numeric",
                        "xi": key[0] if key else "",
                        "match_key": "|".join(str(part) for part in key),
                        "metric": column,
                        "candidate_value": candidate.get(column, ""),
                        "reference_value": reference.get(column, ""),
                        "abs_diff": "",
                        "rel_diff": "",
                    }
                )
                continue
            comparison_rows.append(
                {
                    "artifact": artifact,
                    "match_status": "matched",
                    "xi": key[0] if key else "",
                    "match_key": "|".join(str(part) for part in key),
                    "metric": column,
                    "candidate_value": candidate_value,
                    "reference_value": reference_value,
                    "abs_diff": abs(candidate_value - reference_value),
                    "rel_diff": rel_diff(candidate_value, reference_value),
                }
            )
    return comparison_rows


def summarize_comparison(rows):
    summary = {}
    for row in rows:
        artifact = row["artifact"]
        metric = row["metric"] or "__missing__"
        entry = summary.setdefault(
            (artifact, metric),
            {
                "artifact": artifact,
                "metric": metric,
                "matched_count": 0,
                "missing_count": 0,
                "non_numeric_count": 0,
                "max_abs_diff": None,
                "max_rel_diff": None,
                "mean_abs_diff": None,
                "mean_rel_diff": None,
                "max_abs_match_key": None,
                "max_rel_match_key": None,
                "status_changed_count": 0,
                "not_applicable_count": 0,
                "adaptive_xi_count": 0,
            },
        )
        status = row["match_status"]
        if status == "adaptive_xi_added":
            entry["adaptive_xi_count"] += 1
            continue
        if status.startswith("missing"):
            entry["missing_count"] += 1
            continue
        if status == "status_changed":
            entry["status_changed_count"] += 1
            continue
        if status == "not_applicable":
            entry["not_applicable_count"] += 1
            continue
        if status == "non_numeric":
            entry["non_numeric_count"] += 1
            continue
        entry["matched_count"] += 1
        if row["abs_diff"] in ("", None):
            # String/status matches have no numeric difference to aggregate.
            continue
        abs_value = float(row["abs_diff"])
        rel_value = float(row["rel_diff"])
        entry.setdefault("_abs_values", []).append(abs_value)
        entry.setdefault("_rel_values", []).append(rel_value)
        if entry["max_abs_diff"] is None or abs_value > entry["max_abs_diff"]:
            entry["max_abs_diff"] = abs_value
            entry["max_abs_match_key"] = row["match_key"]
        if entry["max_rel_diff"] is None or rel_value > entry["max_rel_diff"]:
            entry["max_rel_diff"] = rel_value
            entry["max_rel_match_key"] = row["match_key"]

    values = []
    for entry in summary.values():
        abs_values = entry.pop("_abs_values", [])
        rel_values = entry.pop("_rel_values", [])
        if abs_values:
            entry["mean_abs_diff"] = sum(abs_values) / len(abs_values)
            entry["mean_rel_diff"] = sum(rel_values) / len(rel_values)
        values.append(entry)
    decorated_values = [(summary_sort_key(entry), index, entry) for index, entry in enumerate(values)]
    decorated_values.sort()
    return [entry for _, _, entry in decorated_values]


def _row_key_set(artifact_rows, artifact):
    return {artifact_key(artifact, row) for row in artifact_rows}


def _state_for_anchor(xi, temperature, boundary_rows, spinodal_rows, grid_rows=None):
    key = (round(float(xi), 10), round(float(temperature), 10))
    boundary_keys = _row_key_set(boundary_rows, "boundary")
    spinodal_keys = _row_key_set(spinodal_rows, "spinodals")
    if key in boundary_keys:
        return STATE_FIRST_ORDER
    if key in spinodal_keys:
        return STATE_AMBIGUOUS
    if grid_rows:
        if any(
            key_float(row, "xi") == key[0] and key_float(row, "T_MeV") == key[1]
            for row in grid_rows
        ):
            return STATE_AMBIGUOUS
    # The dense reference emits a boundary/spinodal row only for a first-order
    # slice.  Absence on the public scan grid therefore represents the
    # monotone branch, while a first-order row disappearing between stages is
    # recorded separately below as a classification regression.
    return STATE_MONOTONE


def build_public_anchor_state_table(candidate, reference):
    rows = []
    for xi in PUBLIC_XI_VALUES:
        for temperature in PUBLIC_T_VALUES:
            candidate_state = _state_for_anchor(
                xi,
                temperature,
                candidate["boundary"]["rows"],
                candidate["spinodals"]["rows"],
                candidate.get("grid_convergence", {}).get("rows", []),
            )
            reference_state = _state_for_anchor(
                xi,
                temperature,
                reference["boundary"]["rows"],
                reference["spinodals"]["rows"],
                reference.get("grid_convergence", {}).get("rows", []),
            )
            regression = (
                reference_state == STATE_FIRST_ORDER
                and candidate_state in {STATE_AMBIGUOUS, STATE_MONOTONE, STATE_NOT_EVALUATED}
            )
            rows.append(
                {
                    "xi": xi,
                    "T_MeV": temperature,
                    "candidate_state": candidate_state,
                    "reference_state": reference_state,
                    "classification_regression": str(regression).lower(),
                }
            )
    return rows


def build_geometry_gate_table(candidate, reference):
    rows = []
    candidate_boundary = keyed_rows("boundary", candidate["boundary"]["rows"])
    reference_boundary = keyed_rows("boundary", reference["boundary"]["rows"])
    candidate_spinodal = keyed_rows("spinodals", candidate["spinodals"]["rows"])
    reference_spinodal = keyed_rows("spinodals", reference["spinodals"]["rows"])
    for key in sorted(set(candidate_boundary) & set(reference_boundary)):
        c_boundary = candidate_boundary[key]
        r_boundary = reference_boundary[key]
        c_spin = candidate_spinodal.get(key)
        r_spin = reference_spinodal.get(key)
        checks = {}
        for column, tolerance in GEOMETRY_TOLERANCES.items():
            if column in c_boundary and column in r_boundary:
                c_value = finite_float(c_boundary, column)
                r_value = finite_float(r_boundary, column)
            else:
                spin_column = column
                c_value = finite_float(c_spin, spin_column) if c_spin else None
                r_value = finite_float(r_spin, spin_column) if r_spin else None
            checks[column] = {
                "candidate": c_value,
                "reference": r_value,
                "abs_diff": abs(c_value - r_value) if c_value is not None and r_value is not None else None,
                "tolerance": tolerance,
                "pass": c_value is not None and r_value is not None and abs(c_value - r_value) <= tolerance,
            }
        rows.append(
            {
                "xi": key[0],
                "T_MeV": key[1],
                "mu_pass": checks["mu_transition_MeV"]["pass"],
                "density_pass": checks["rho_hadron"]["pass"] and checks["rho_quark"]["pass"],
                "area_pass": checks["area_residual"]["pass"],
                "spinodal_pass": checks["rho_spinodal_hadron"]["pass"] and checks["rho_spinodal_quark"]["pass"],
                "all_pass": all(check["pass"] for check in checks.values()),
                "mu_abs_diff": checks["mu_transition_MeV"]["abs_diff"],
                "rho_hadron_abs_diff": checks["rho_hadron"]["abs_diff"],
                "rho_quark_abs_diff": checks["rho_quark"]["abs_diff"],
                "area_abs_diff": checks["area_residual"]["abs_diff"],
                "rho_spinodal_hadron_abs_diff": checks["rho_spinodal_hadron"]["abs_diff"],
                "rho_spinodal_quark_abs_diff": checks["rho_spinodal_quark"]["abs_diff"],
            }
        )
    return rows


def build_cep_gate_table(candidate, reference):
    rows = []
    candidate_rows = keyed_rows("cep", candidate["cep"]["rows"])
    reference_rows = keyed_rows("cep", reference["cep"]["rows"])
    for key in sorted(set(candidate_rows) & set(reference_rows)):
        c_row = candidate_rows[key]
        r_row = reference_rows[key]
        c_bracket = finite_float(c_row, "bracket_width_T_MeV")
        r_bracket = finite_float(r_row, "bracket_width_T_MeV")
        c_ambiguity = finite_float(c_row, "ambiguity_width_T_MeV")
        r_ambiguity = finite_float(r_row, "ambiguity_width_T_MeV")
        c_resolution = finite_float(c_row, "temperature_resolution_target_MeV")
        r_resolution = finite_float(r_row, "temperature_resolution_target_MeV")
        rows.append(
            {
                "xi": key[0],
                "candidate_result_status": _value_text(c_row, "result_status"),
                "reference_result_status": _value_text(r_row, "result_status"),
                "candidate_bracket_width_T_MeV": c_bracket,
                "reference_bracket_width_T_MeV": r_bracket,
                "candidate_ambiguity_width_T_MeV": c_ambiguity,
                "reference_ambiguity_width_T_MeV": r_ambiguity,
                "candidate_endpoint_resolution_MeV": c_resolution,
                "reference_endpoint_resolution_MeV": r_resolution,
                "candidate_bracket_pass": c_bracket is not None and c_bracket <= CEP_BRACKET_TOLERANCE_MEV,
                "reference_bracket_pass": r_bracket is not None and r_bracket <= CEP_BRACKET_TOLERANCE_MEV,
            }
        )
    return rows


def classification_regressions(anchor_rows):
    return [row for row in anchor_rows if row["classification_regression"] == "true"]


def cep_bracket_failures(cep_rows):
    return [
        row
        for row in cep_rows
        if row["candidate_bracket_pass"] is False
        or row["reference_bracket_pass"] is False
    ]


def geometry_failures(geometry_rows):
    return [row for row in geometry_rows if not row["all_pass"]]


def crossover_failures(comparison_rows):
    failures = []
    for row in comparison_rows:
        if row["artifact"] != "crossover" or row["metric"] not in {"T_crossover_MeV", "rho", "derivative"}:
            continue
        if row["metric"] == "derivative":
            bad = float(row["rel_diff"]) > CROSSOVER_TOLERANCES["derivative_rel"]
        else:
            bad = float(row["abs_diff"]) > CROSSOVER_TOLERANCES[row["metric"]]
        if bad:
            failures.append(row)
    return failures


def determine_verdict(*, bad_inventory, anchor_rows, cep_rows, geometry_rows, comparison_rows):
    secondary = []
    if bad_inventory:
        secondary.append("artifact_invalid")
    if classification_regressions(anchor_rows):
        secondary.append("classification_regression")
    if cep_bracket_failures(cep_rows):
        secondary.append("cep_bracket_unresolved")
    if geometry_failures(geometry_rows):
        secondary.append("first_order_geometry_unstable")
    if crossover_failures(comparison_rows):
        secondary.append("crossover_unstable")
    precedence = [
        "artifact_invalid",
        "classification_regression",
        "cep_bracket_unresolved",
        "first_order_geometry_unstable",
        "crossover_unstable",
        "convergence_candidate",
    ]
    verdict = next((item for item in precedence if item in secondary), "convergence_candidate")
    return verdict, secondary


def write_csv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_markdown(path, payload):
    lines = [
        "# PNJL Phase Reference Convergence Summary",
        "",
        f"- candidate: `{payload['candidate_label']}`",
        f"- reference: `{payload['reference_label']}`",
        f"- verdict: `{payload['verdict']}`",
        f"- secondary risks: `{', '.join(payload.get('secondary_failures', [])) or 'none'}`",
        "",
        "## Artifact Inventory",
        "",
        "| side | artifact | rows | xi_count | naninf | converged_false |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for side in ("candidate", "reference"):
        for artifact, info in payload[f"{side}_inventory"].items():
            lines.append(
                f"| {side} | {artifact} | {info['row_count']} | {info['xi_count']} | "
                f"{info['naninf_count']} | {info['converged_false_count']} |"
            )

    lines.extend(["", "## Max Differences", "", "| artifact | metric | n | max_abs | max_rel | max_abs_key | max_rel_key |", "|---|---|---:|---:|---:|---|---|"])
    for row in payload["comparison_summary"]:
        if row["metric"] == "__missing__":
            continue
        lines.append(
            f"| {row['artifact']} | {row['metric']} | {row['matched_count']} | "
            f"{row['max_abs_diff']} | {row['max_rel_diff']} | "
            f"{row['max_abs_match_key']} | {row['max_rel_match_key']} |"
        )

    if payload["missing_or_non_numeric"]:
        lines.extend(["", "## Missing Or Non-Numeric", ""])
        for row in payload["missing_or_non_numeric"]:
            lines.append(
                f"- {row['artifact']} {row['metric']}: missing={row['missing_count']}, "
                f"non_numeric={row['non_numeric_count']}"
            )

    if payload.get("adaptive_xi_rows"):
        lines.extend(["", "## Adaptive Xi Additions", ""])
        for row in payload["adaptive_xi_rows"]:
            lines.append(
                f"- {row['artifact']}: {row['adaptive_xi_count']} candidate-only xi rows; "
                "reported separately from missing matches"
            )

    if payload.get("classification_regressions"):
        lines.extend(["", "## Classification Regressions", ""])
        for row in payload["classification_regressions"]:
            lines.append(
                f"- xi={row['xi']}, T={row['T_MeV']}: "
                f"{row['reference_state']} -> {row['candidate_state']}"
            )
    if payload.get("cep_bracket_failures"):
        lines.extend(["", "## CEP Bracket Failures", ""])
        for row in payload["cep_bracket_failures"]:
            lines.append(
                f"- xi={row['xi']}: candidate bracket={row['candidate_bracket_width_T_MeV']} MeV, "
                f"ambiguity={row['candidate_ambiguity_width_T_MeV']} MeV"
            )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    xi_filter = parse_float_set(args.xi_values)
    candidate = load_artifacts(args.candidate_root, args.candidate_tag, xi_filter)
    reference = load_artifacts(args.reference_root, args.reference_tag, xi_filter)
    candidate_optional = load_optional_artifacts(args.candidate_root, args.candidate_tag, xi_filter)
    reference_optional = load_optional_artifacts(args.reference_root, args.reference_tag, xi_filter)
    if candidate_optional["grid_convergence"]["available"]:
        candidate["grid_convergence"] = candidate_optional["grid_convergence"]
    if reference_optional["grid_convergence"]["available"]:
        reference["grid_convergence"] = reference_optional["grid_convergence"]

    candidate_inventory = {
        artifact: summarize_artifact(artifact, data["fields"], data["rows"])
        for artifact, data in candidate.items()
    }
    reference_inventory = {
        artifact: summarize_artifact(artifact, data["fields"], data["rows"])
        for artifact, data in reference.items()
    }
    for side, optional, inventory in (
        ("candidate", candidate_optional, candidate_inventory),
        ("reference", reference_optional, reference_inventory),
    ):
        for artifact, data in optional.items():
            if data["available"]:
                inventory[artifact] = summarize_artifact(artifact, data["fields"], data["rows"])

    comparison_rows = []
    for artifact in ARTIFACT_FILES:
        comparison_rows.extend(
            compare_artifact(
                artifact,
                candidate[artifact]["rows"],
                reference[artifact]["rows"],
                candidate_cep_rows=candidate["cep"]["rows"],
                reference_cep_rows=reference["cep"]["rows"],
            )
        )

    comparison_summary = summarize_comparison(comparison_rows)
    missing_or_non_numeric = [
        row
        for row in comparison_summary
        if row["missing_count"] > 0 or row["non_numeric_count"] > 0
    ]
    adaptive_xi_rows = [
        row
        for row in comparison_summary
        if row.get("adaptive_xi_count", 0) > 0
    ]
    bad_inventory = []
    for side, inventory in (("candidate", candidate_inventory), ("reference", reference_inventory)):
        for artifact, info in inventory.items():
            if (
                info["naninf_count"] > 0
                or info["duplicate_key_count"] > 0
                or info["missing_columns"]
                or (artifact != "grid_convergence" and (info["converged_false_count"] or 0) > 0)
            ):
                bad_inventory.append({"side": side, "artifact": artifact, **info})
    non_numeric_comparisons = [row for row in comparison_rows if row["match_status"] == "non_numeric"]
    if non_numeric_comparisons:
        bad_inventory.append(
            {
                "side": "comparison",
                "artifact": "matched_artifact_values",
                "non_numeric_count": len(non_numeric_comparisons),
                "examples": non_numeric_comparisons[:5],
            }
        )

    anchor_rows = build_public_anchor_state_table(candidate, reference)
    geometry_rows = build_geometry_gate_table(candidate, reference)
    cep_rows = build_cep_gate_table(candidate, reference)
    verdict, secondary_failures = determine_verdict(
        bad_inventory=bad_inventory,
        anchor_rows=anchor_rows,
        cep_rows=cep_rows,
        geometry_rows=geometry_rows,
        comparison_rows=comparison_rows,
    )
    classification_rows = classification_regressions(anchor_rows)
    cep_failures = cep_bracket_failures(cep_rows)

    payload = {
        "candidate_label": args.candidate_label,
        "reference_label": args.reference_label,
        "candidate_root": args.candidate_root.as_posix(),
        "reference_root": args.reference_root.as_posix(),
        "candidate_tag": args.candidate_tag,
        "reference_tag": args.reference_tag,
        "xi_filter": sorted(xi_filter) if xi_filter is not None else None,
        "verdict": verdict,
        "secondary_failures": secondary_failures,
        "public_grid": {
            "xi_values": list(PUBLIC_XI_VALUES),
            "T_values_MeV": list(PUBLIC_T_VALUES),
            "anchor_count": len(PUBLIC_XI_VALUES) * len(PUBLIC_T_VALUES),
        },
        "candidate_inventory": candidate_inventory,
        "reference_inventory": reference_inventory,
        "comparison_summary": comparison_summary,
        "missing_or_non_numeric": missing_or_non_numeric,
        "adaptive_xi_rows": adaptive_xi_rows,
        "bad_inventory": bad_inventory,
        "public_anchor_states": anchor_rows,
        "classification_regressions": classification_rows,
        "geometry_gate_rows": geometry_rows,
        "cep_gate_rows": cep_rows,
        "cep_bracket_failures": cep_failures,
        "crossover_failure_count": len(crossover_failures(comparison_rows)),
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        args.out_dir / "phase_reference_convergence_comparison.csv",
        comparison_rows,
        [
            "artifact",
            "match_status",
            "grid_change_kind",
            "xi",
            "match_key",
            "metric",
            "candidate_value",
            "reference_value",
            "abs_diff",
            "rel_diff",
        ],
    )
    inventory_rows = []
    for side, inventory in (("candidate", candidate_inventory), ("reference", reference_inventory)):
        for artifact, info in inventory.items():
            inventory_rows.append(
                {
                    "side": side,
                    "artifact": artifact,
                    "row_count": info["row_count"],
                    "xi_count": info["xi_count"],
                    "xi_values": json.dumps(info["xi_values"]),
                    "naninf_count": info["naninf_count"],
                    "negative_by_column": json.dumps(info["negative_by_column"], sort_keys=True),
                    "converged_false_count": info["converged_false_count"],
                    "duplicate_key_count": info["duplicate_key_count"],
                    "missing_columns": json.dumps(info["missing_columns"]),
                }
            )
    write_csv(
        args.out_dir / "phase_reference_artifact_inventory.csv",
        inventory_rows,
        [
            "side",
            "artifact",
            "row_count",
            "xi_count",
            "xi_values",
            "naninf_count",
            "negative_by_column",
            "converged_false_count",
            "duplicate_key_count",
            "missing_columns",
        ],
    )
    write_csv(
        args.out_dir / "phase_reference_anchor_states.csv",
        anchor_rows,
        ["xi", "T_MeV", "candidate_state", "reference_state", "classification_regression"],
    )
    write_csv(
        args.out_dir / "phase_reference_geometry_gates.csv",
        geometry_rows,
        [
            "xi", "T_MeV", "mu_pass", "density_pass", "area_pass", "spinodal_pass", "all_pass",
            "mu_abs_diff", "rho_hadron_abs_diff", "rho_quark_abs_diff", "area_abs_diff",
            "rho_spinodal_hadron_abs_diff", "rho_spinodal_quark_abs_diff",
        ],
    )
    write_csv(
        args.out_dir / "phase_reference_cep_gates.csv",
        cep_rows,
        [
            "xi", "candidate_result_status", "reference_result_status",
            "candidate_bracket_width_T_MeV", "reference_bracket_width_T_MeV",
            "candidate_ambiguity_width_T_MeV", "reference_ambiguity_width_T_MeV",
            "candidate_endpoint_resolution_MeV", "reference_endpoint_resolution_MeV",
            "candidate_bracket_pass", "reference_bracket_pass",
        ],
    )
    (args.out_dir / "phase_reference_convergence_summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    write_markdown(args.out_dir / "phase_reference_convergence_summary.md", payload)
    print(json.dumps({"verdict": verdict, "out_dir": args.out_dir.as_posix()}, ensure_ascii=False))


if __name__ == "__main__":
    main()
