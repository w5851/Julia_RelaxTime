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

COMPARE_COLUMNS = {
    "boundary": ["mu_transition_MeV", "rho_hadron", "rho_quark"],
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


def summarize_artifact(fields, rows):
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
    return {
        "row_count": len(rows),
        "columns": fields,
        "xi_count": len(xi_values),
        "xi_values": xi_values,
        "naninf_count": naninf_count,
        "negative_by_column": negative_by_column,
        "converged_false_count": converged_false,
    }


def keyed_rows(artifact, rows):
    if artifact == "cep":
        return {(key_float(row, "xi"),): row for row in rows}
    if artifact in ("boundary", "spinodals"):
        return {(key_float(row, "xi"), key_float(row, "T_MeV")): row for row in rows}
    if artifact == "crossover":
        by_xi = {}
        for row in rows:
            xi = key_float(row, "xi")
            by_xi.setdefault(xi, []).append(row)
        keyed = {}
        for xi, xi_rows in by_xi.items():
            decorated_rows = [(crossover_sort_key(row), index, row) for index, row in enumerate(xi_rows)]
            decorated_rows.sort()
            for index, (_, _, row) in enumerate(decorated_rows, start=1):
                keyed[(xi, index)] = row
        return keyed
    fail(f"unknown artifact: {artifact}")


def rel_diff(candidate, reference):
    denom = max(abs(reference), 1.0e-12)
    return abs(candidate - reference) / denom


def append_missing_rows(comparison_rows, artifact, keys, side, *, adaptive_xi=False):
    for key in sorted(keys):
        comparison_rows.append(
            {
                "artifact": artifact,
                "match_status": "adaptive_xi_added" if adaptive_xi else f"missing_in_{side}",
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
    # A denser C1/C2 grid legitimately adds xi midpoints.  Keep these rows in
    # the audit output without treating them as a failed one-to-one match.
    append_missing_rows(
        comparison_rows,
        artifact,
        candidate_keys - reference_keys,
        "reference",
        adaptive_xi=True,
    )
    append_missing_rows(comparison_rows, artifact, reference_keys - candidate_keys, "candidate")


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


def compare_artifact(artifact, candidate_rows, reference_rows):
    if artifact == "cep":
        return compare_cep_artifact(candidate_rows, reference_rows)
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

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    xi_filter = parse_float_set(args.xi_values)
    candidate = load_artifacts(args.candidate_root, args.candidate_tag, xi_filter)
    reference = load_artifacts(args.reference_root, args.reference_tag, xi_filter)

    candidate_inventory = {
        artifact: summarize_artifact(data["fields"], data["rows"])
        for artifact, data in candidate.items()
    }
    reference_inventory = {
        artifact: summarize_artifact(data["fields"], data["rows"])
        for artifact, data in reference.items()
    }

    comparison_rows = []
    for artifact in ARTIFACT_FILES:
        comparison_rows.extend(
            compare_artifact(artifact, candidate[artifact]["rows"], reference[artifact]["rows"])
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
            if info["naninf_count"] > 0 or (info["converged_false_count"] or 0) > 0:
                bad_inventory.append({"side": side, "artifact": artifact, **info})

    verdict = "needs_manual_review"
    if missing_or_non_numeric:
        verdict = "blocked_missing_matches"
    elif bad_inventory:
        verdict = "blocked_invalid_values"

    payload = {
        "candidate_label": args.candidate_label,
        "reference_label": args.reference_label,
        "candidate_root": args.candidate_root.as_posix(),
        "reference_root": args.reference_root.as_posix(),
        "candidate_tag": args.candidate_tag,
        "reference_tag": args.reference_tag,
        "xi_filter": sorted(xi_filter) if xi_filter is not None else None,
        "verdict": verdict,
        "candidate_inventory": candidate_inventory,
        "reference_inventory": reference_inventory,
        "comparison_summary": comparison_summary,
        "missing_or_non_numeric": missing_or_non_numeric,
        "adaptive_xi_rows": adaptive_xi_rows,
        "bad_inventory": bad_inventory,
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(
        args.out_dir / "phase_reference_convergence_comparison.csv",
        comparison_rows,
        [
            "artifact",
            "match_status",
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
