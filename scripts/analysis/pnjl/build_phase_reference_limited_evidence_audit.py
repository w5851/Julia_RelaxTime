#!/usr/bin/env python3
"""Audit the 763 C2 unresolved rho rows against the external raw oracle archive.

This is deliberately a solver-free, diagnostic-only audit.  It reads the C2
grid-convergence CSV, retrieves only the matching raw fixed-rho curves from the
Zenodo ZIP by HTTP range requests, and checks coverage, finite/converged status,
and the presence of a discrete ``+ -> - -> +`` S-shape in ``mu_B(rho)``.
It never calls Julia, rewrites C2 labels, or writes reference/production data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import struct
import time
import urllib.request
import zlib
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


CALCULATION_SHA = "3c5f6b3c9bd535cff7657364dadb2efc31f2ea48"
ARCHIVE_RECORD = "21980679"
ARCHIVE_SHA256 = "467be7fb1075d1a5f0de3dd0d8afe29d9206a156c0ca7135a1e50967a4f18ccc"
ARCHIVE_URL = (
    "https://zenodo.org/api/records/21980679/files/"
    "pnjl-raw-curve-archive-issue130-c2-oracle-raw-v1.zip/content"
)
EXPECTED_TARGET_ROWS = 763
EXPECTED_TARGET_COORDINATES = 761
EXPECTED_CURVE_ROWS = 1281
RHO_STEP = 0.003125
SLOPE_EPS = 1.0e-8
NEGATIVE_SLOPE_REQUIRED = -1.0e-5
MIN_SIGN_RUN = 2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--c2-root", type=Path, required=True)
    parser.add_argument("--curve-index", type=Path, required=True)
    parser.add_argument("--central-directory", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--cache-root", type=Path, required=True)
    parser.add_argument("--archive-url", default=ARCHIVE_URL)
    parser.add_argument("--workers", type=int, default=8)
    return parser.parse_args()


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"CSV has no rows: {path}")
    return rows


def write_csv(path: Path, fieldnames: Iterable[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        writer.writerows(rows)


def finite(value: object) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def float_key(value: str) -> float:
    return round(float(value), 10)


def reason_class_for_reason(reason: str) -> str:
    if reason == "ok":
        return "unresolved_unclassified"
    if reason == "hybrid_stage_c_not_converged":
        return "hybrid_stage_c_unresolved"
    if reason.startswith("classification_changed_or_unresolved:"):
        return "classification_transition"
    return "other_unresolved"


def xi_token(xi: float) -> str:
    if abs(xi) < 1.0e-12:
        return "0"
    text = f"{xi:g}".replace("-", "m").replace(".", "p")
    return text


def temperature_token(T: float) -> str:
    text = f"{T:g}".replace(".", "p")
    return text


def parse_central_directory(path: Path) -> dict[str, dict[str, int]]:
    data = path.read_bytes()
    entries: dict[str, dict[str, int]] = {}
    position = 0
    while position < len(data):
        if data[position : position + 4] != b"PK\x01\x02":
            raise ValueError(f"invalid central-directory signature at {position}")
        values = struct.unpack_from("<4s6H3I5H2I", data, position)
        filename_len, extra_len, comment_len = values[10:13]
        name = data[position + 46 : position + 46 + filename_len].decode("utf-8")
        entries[name] = {
            "method": values[4],
            "compressed_size": values[8],
            "uncompressed_size": values[9],
            "local_offset": values[16],
        }
        position += 46 + filename_len + extra_len + comment_len
    if position != len(data):
        raise ValueError("central-directory parse did not consume the input")
    return entries


def range_get(url: str, start: int, end: int, attempts: int = 4) -> bytes:
    expected = end - start + 1
    last_error: Exception | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(request, timeout=180) as response:
                data = response.read()
            if len(data) == expected:
                return data
            # A few HTTP proxies ignore Range and return the whole object.  Do
            # not accept a truncated response, but safely slice a full response.
            if len(data) > end:
                return data[start : end + 1]
            raise IOError(f"range length {len(data)} != {expected}")
        except Exception as error:  # pragma: no cover - network failures vary
            last_error = error
            time.sleep(2**attempt)
    raise RuntimeError(f"range request failed {start}-{end}: {last_error}")


def extract_entry(name: str, entry: dict[str, int], url: str) -> bytes:
    offset = entry["local_offset"]
    header = range_get(url, offset, offset + 99)
    if header[:4] != b"PK\x03\x04":
        raise ValueError(f"invalid local header for {name}")
    values = struct.unpack_from("<4s5H3I2H", header, 0)
    _, _, _, method, _, _, _, _, _, filename_len, extra_len = values
    data_start = offset + 30 + filename_len + extra_len
    compressed = range_get(url, data_start, data_start + entry["compressed_size"] - 1)
    if method == 8:
        return zlib.decompress(compressed, -15)
    if method == 0:
        return compressed
    raise ValueError(f"unsupported ZIP method {method} for {name}")


def cache_path(cache_root: Path, raw_curve_file: str) -> Path:
    relative = Path(raw_curve_file)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"unsafe raw curve path: {raw_curve_file}")
    return cache_root / relative


def load_raw_curve(
    index_row: dict[str, str],
    central: dict[str, dict[str, int]],
    cache_root: Path,
    url: str,
) -> tuple[bytes, Path, str]:
    archive_name = "archive/" + index_row["raw_curve_file"].replace("\\", "/")
    if archive_name not in central:
        raise FileNotFoundError(archive_name)
    output = cache_path(cache_root, index_row["raw_curve_file"])
    output.parent.mkdir(parents=True, exist_ok=True)
    expected_sha = index_row["raw_curve_sha256"].lower()
    if output.is_file():
        raw = output.read_bytes()
        if sha256_bytes(raw) == expected_sha:
            return raw, output, "cache"
    raw = extract_entry(archive_name, central[archive_name], url)
    if sha256_bytes(raw) != expected_sha:
        raise ValueError(f"raw curve hash mismatch: {archive_name}")
    output.write_bytes(raw)
    return raw, output, "range_fetch"


def moving_average(values: list[float], width: int = 5) -> list[float]:
    half = width // 2
    result = []
    for index in range(len(values)):
        left = max(0, index - half)
        right = min(len(values), index + half + 1)
        result.append(sum(values[left:right]) / (right - left))
    return result


def sign_runs(slopes: list[float]) -> list[tuple[int, int, int]]:
    signs: list[int] = []
    for slope in slopes:
        if slope > SLOPE_EPS:
            signs.append(1)
        elif slope < -SLOPE_EPS:
            signs.append(-1)
        else:
            signs.append(0)
    # Zero slopes are not evidence for a turning point.  Keep the original
    # interval index while collapsing them for the topology observation.
    nonzero = [(index, sign) for index, sign in enumerate(signs) if sign]
    runs: list[tuple[int, int, int]] = []
    for index, sign in nonzero:
        if not runs or runs[-1][2] != sign:
            runs.append((index, index, sign))
        else:
            runs[-1] = (runs[-1][0], index, sign)
    return runs


def classify_curve(raw: bytes) -> tuple[dict[str, object], list[tuple[float, float]]]:
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8")))
    required = {"T_MeV", "rho", "xi", "mu_B_MeV", "residual_norm", "converged"}
    if not required.issubset(reader.fieldnames or set()):
        raise ValueError(f"raw curve missing required columns: {sorted(required)}")
    rows = list(reader)
    errors: list[str] = []
    if len(rows) != EXPECTED_CURVE_ROWS:
        errors.append(f"row_count={len(rows)}")
    values: list[tuple[float, float]] = []
    seen_rho: set[float] = set()
    for row in rows:
        if not all(finite(row.get(column)) for column in ("T_MeV", "rho", "xi", "mu_B_MeV", "residual_norm")):
            errors.append("nonfinite_required_field")
            continue
        rho = float(row["rho"])
        mu = float(row["mu_B_MeV"])
        if rho in seen_rho:
            errors.append("duplicate_rho")
        seen_rho.add(rho)
        if row["converged"].strip().lower() != "true":
            errors.append("nonconverged_row")
        values.append((rho, mu))
    values.sort()
    if values:
        if abs(values[0][0]) > 1.0e-10 or abs(values[-1][0] - 4.0) > 1.0e-10:
            errors.append("rho_range_not_0_to_4")
        gaps = [b[0] - a[0] for a, b in zip(values, values[1:])]
        if gaps and max(abs(gap - RHO_STEP) for gap in gaps) > 1.0e-10:
            errors.append("rho_step_not_0.003125")
    if errors:
        return {
            "finite_and_converged": False,
            "raw_curve_status": "invalid_raw_curve",
            "validation_errors": ";".join(sorted(set(errors))),
        }, values

    mu_values = [mu for _, mu in values]
    rho_values = [rho for rho, _ in values]
    raw_slopes = [
        (right - left) / RHO_STEP for left, right in zip(mu_values, mu_values[1:])
    ]
    smooth_mu = moving_average(mu_values)
    smooth_slopes = [
        (right - left) / RHO_STEP for left, right in zip(smooth_mu, smooth_mu[1:])
    ]
    runs = sign_runs(smooth_slopes)
    patterns: list[tuple[tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]]] = []
    for first in range(len(runs) - 2):
        candidate = runs[first : first + 3]
        if [run[2] for run in candidate] == [1, -1, 1]:
            negative_slopes = smooth_slopes[candidate[1][0] : candidate[1][1] + 1]
            if len(negative_slopes) >= MIN_SIGN_RUN and min(negative_slopes) <= NEGATIVE_SLOPE_REQUIRED:
                patterns.append((candidate[0], candidate[1], candidate[2]))
    if patterns:
        first = patterns[0]
        extrema_rho = (rho_values[first[0][1]], rho_values[first[1][1]], rho_values[first[2][1]])
        shape = "s_shape_present"
        interpretation = "raw oracle curve contains a strict +→−→+ slope topology; this is diagnostic evidence only"
    else:
        extrema_rho = ("", "", "")
        shape = "monotone_or_no_strict_s_shape"
        interpretation = "raw oracle curve has no strict +→−→+ topology under the stated discrete test"
    return {
        "finite_and_converged": True,
        "raw_curve_status": "valid",
        "validation_errors": "",
        "oracle_shape_observation": shape,
        "s_shape_candidate_count": len(patterns),
        "turning_run_count": len(runs),
        "min_smoothed_slope": min(smooth_slopes),
        "max_smoothed_slope": max(smooth_slopes),
        "min_raw_slope": min(raw_slopes),
        "max_raw_slope": max(raw_slopes),
        "rho_extremum_1": extrema_rho[0],
        "rho_extremum_2": extrema_rho[1],
        "rho_extremum_3": extrema_rho[2],
        "interpretation": interpretation,
    }, values


def choose_figure_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[(str(row["reason_class"]), str(row.get("oracle_shape_observation", "unknown")))].append(row)
    selected: list[dict[str, object]] = []
    for key in sorted(groups):
        candidates = sorted(groups[key], key=lambda row: (-abs(float(row.get("min_smoothed_slope") or 0.0)), float(row["xi"]), float(row["T_MeV"])))
        selected.append(candidates[0])
    return selected[:12]


def make_figures(output_root: Path, rows: list[dict[str, object]], curves: dict[tuple[float, float], list[tuple[float, float]]], input_sha: str) -> list[dict[str, object]]:
    try:
        import matplotlib.pyplot as plt
    except Exception as error:  # pragma: no cover - environment dependent
        return [{"status": "skipped", "reason": f"matplotlib unavailable: {error}"}]
    figure_root = output_root / "figures"
    figure_root.mkdir(parents=True, exist_ok=True)
    manifest: list[dict[str, object]] = []
    for row in choose_figure_rows(rows):
        key = (float(row["xi"]), float(row["T_MeV"]))
        if key not in curves:
            continue
        values = curves[key]
        rho = [item[0] for item in values]
        mu = [item[1] for item in values]
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
        axes[0].plot(rho, mu, color="#1f4e79", linewidth=1.0)
        axes[0].set_xlabel(r"$\rho$ [fm$^{-3}$]")
        axes[0].set_ylabel(r"$\mu_B$ [MeV]")
        axes[0].set_title("full raw curve")
        if row.get("oracle_shape_observation") == "s_shape_present":
            center = [float(row["rho_extremum_1"]), float(row["rho_extremum_2"]), float(row["rho_extremum_3"])]
            left = max(0.0, min(center) - 0.15)
            right = min(4.0, max(center) + 0.15)
            selected = [index for index, value in enumerate(rho) if left <= value <= right]
            if selected:
                y_values = [mu[index] for index in selected]
                y_margin = max(0.002, 0.12 * (max(y_values) - min(y_values)))
                axes[1].set_xlim(left, right)
                axes[1].set_ylim(min(y_values) - y_margin, max(y_values) + y_margin)
        else:
            center_index = max(range(1, len(mu) - 1), key=lambda index: abs(mu[index + 1] - 2 * mu[index] + mu[index - 1]))
            left_index = max(0, center_index - 40)
            right_index = min(len(rho) - 1, center_index + 40)
            y_values = mu[left_index : right_index + 1]
            y_margin = max(0.002, 0.12 * (max(y_values) - min(y_values)))
            axes[1].set_xlim(rho[left_index], rho[right_index])
            axes[1].set_ylim(min(y_values) - y_margin, max(y_values) + y_margin)
        axes[1].plot(rho, mu, color="#b23a48", linewidth=1.2)
        axes[1].set_xlabel(r"$\rho$ [fm$^{-3}$]")
        axes[1].set_ylabel(r"$\mu_B$ [MeV]")
        axes[1].set_title("local diagnostic zoom")
        filename = f"rho_mu_xi_{xi_token(key[0])}_T_{temperature_token(key[1])}_raw_audit.png"
        path = figure_root / filename
        fig.savefig(path, dpi=180)
        plt.close(fig)
        manifest.append({"path": f"figures/{filename}", "xi": key[0], "T_MeV": key[1], "input_sha256": input_sha, "sha256": sha256_file(path)})
    (figure_root / "plot_manifest.json").write_text(json.dumps({"schema_version": "pnjl_phase_reference_limited_evidence_audit_plots_v1", "figures": manifest}, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    args = parse_args()
    c2_root = args.c2_root.resolve()
    curve_index_path = args.curve_index.resolve()
    central_path = args.central_directory.resolve()
    output_root = args.output_root.resolve()
    cache_root = args.cache_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    tables_root = output_root / "tables"
    tables_root.mkdir(parents=True, exist_ok=True)

    grid_path = next((c2_root / "reference").glob("phase_grid_convergence_*.csv"))
    c2_manifest_path = next((c2_root / "reference").glob("phase_reference_*_manifest.json"))
    validation_path = c2_root / "validation_report.json"
    grid_rows = read_csv(grid_path)
    c2_manifest = read_json(c2_manifest_path)
    validation_payload = read_json(validation_path) if validation_path.is_file() else {}
    c2_sha = (
        c2_manifest.get("calculation_git_commit")
        or c2_manifest.get("git_commit")
        or c2_manifest.get("calculation_sha")
        or validation_payload.get("calculation_git_commit")
        or validation_payload.get("git_commit")
        or validation_payload.get("calculation_sha")
        or ""
    )
    if str(c2_sha).lower() != CALCULATION_SHA:
        raise ValueError("C2 calculation SHA does not match the raw archive contract")
    unresolved = [row for row in grid_rows if row.get("converged", "").strip().lower() != "true"]
    target_grid_rows = [
        row for row in unresolved
        if row.get("axis") == "rho"
        and (row.get("reason") in {"ok", "hybrid_stage_c_not_converged"} or row.get("reason", "").startswith("classification_changed_or_unresolved:"))
    ]
    target_coordinates = sorted({(float_key(row["xi"]), float_key(row["T_MeV"])) for row in target_grid_rows})
    if len(target_grid_rows) != EXPECTED_TARGET_ROWS or len(target_coordinates) != EXPECTED_TARGET_COORDINATES:
        raise ValueError(f"target scope changed: rows={len(target_grid_rows)}, coordinates={len(target_coordinates)}")

    index_rows = read_csv(curve_index_path)
    index_by_coordinate = {
        (float_key(row["xi"]), float_key(row["T_MeV"])): row
        for row in index_rows if row.get("method") == "independent_oracle"
    }
    missing = [coordinate for coordinate in target_coordinates if coordinate not in index_by_coordinate]
    if missing:
        raise ValueError(f"raw archive lacks {len(missing)} target coordinates")
    central = parse_central_directory(central_path)

    coordinate_reasons: dict[tuple[float, float], list[dict[str, str]]] = defaultdict(list)
    for row in target_grid_rows:
        coordinate_reasons[(float_key(row["xi"]), float_key(row["T_MeV"]))].append(row)

    cache_root.mkdir(parents=True, exist_ok=True)
    payloads: dict[tuple[float, float], tuple[bytes, Path, str]] = {}
    errors: dict[tuple[float, float], str] = {}
    with ThreadPoolExecutor(max_workers=max(1, min(args.workers, 16))) as executor:
        futures = {
            executor.submit(load_raw_curve, index_by_coordinate[coordinate], central, cache_root, args.archive_url): coordinate
            for coordinate in target_coordinates
        }
        for future in as_completed(futures):
            coordinate = futures[future]
            try:
                payloads[coordinate] = future.result()
            except Exception as error:  # pragma: no cover - network failures vary
                errors[coordinate] = repr(error)

    audit_rows: list[dict[str, object]] = []
    curves: dict[tuple[float, float], list[tuple[float, float]]] = {}
    for coordinate in target_coordinates:
        index_row = index_by_coordinate[coordinate]
        reason_rows = coordinate_reasons[coordinate]
        reasons = sorted({row.get("reason", "") for row in reason_rows})
        reason_classes = sorted({reason_class_for_reason(reason) for reason in reasons})
        reason_class = reason_classes[0] if len(reason_classes) == 1 else "mixed_reasons"
        base = {
            "xi": coordinate[0],
            "T_MeV": coordinate[1],
            "c2_unresolved_row_count": len(reason_rows),
            "c2_reasons": "|".join(reasons),
            "c2_levels": "|".join(sorted({row.get("level", "") for row in reason_rows})),
            "reason_class": reason_class,
            "raw_curve_file": index_row["raw_curve_file"],
            "raw_curve_sha256": index_row["raw_curve_sha256"],
            "raw_curve_rows_expected": index_row["raw_curve_rows"],
            "raw_source_run_id": index_row["source_artifact_run_id"],
            "raw_archive_record": ARCHIVE_RECORD,
            "oracle_labels_used_for_routing": False,
            "production_override": False,
        }
        if coordinate in errors:
            base.update({"finite_and_converged": False, "raw_curve_status": "range_fetch_failed", "validation_errors": errors[coordinate], "oracle_shape_observation": "inconclusive"})
            audit_rows.append(base)
            continue
        raw, cache_path_value, fetch_mode = payloads[coordinate]
        result, values = classify_curve(raw)
        curves[coordinate] = values
        base.update(result)
        base.update({"cache_path": str(cache_path_value), "fetch_mode": fetch_mode})
        audit_rows.append(base)

    audit_rows.sort(key=lambda row: (float(row["xi"]), float(row["T_MeV"])))
    write_csv(tables_root / "target_coordinates.csv", ("xi", "T_MeV", "c2_unresolved_row_count", "c2_reasons", "c2_levels", "reason_class"), [
        {key: row[key] for key in ("xi", "T_MeV", "c2_unresolved_row_count", "c2_reasons", "c2_levels", "reason_class")} for row in audit_rows
    ])
    fields = sorted({key for row in audit_rows for key in row})
    write_csv(tables_root / "raw_curve_audit.csv", fields, audit_rows)

    valid_rows = [row for row in audit_rows if row.get("finite_and_converged") is True]
    shape_counts = Counter(str(row.get("oracle_shape_observation", "inconclusive")) for row in audit_rows)
    reason_shape_counts = Counter((str(row["reason_class"]), str(row.get("oracle_shape_observation", "inconclusive"))) for row in audit_rows)
    audit_by_coordinate = {(float(row["xi"]), float(row["T_MeV"])): row for row in audit_rows}
    reason_shape_row_counts: Counter[tuple[str, str]] = Counter()
    reason_shape_coordinate_sets: defaultdict[tuple[str, str], set[tuple[float, float]]] = defaultdict(set)
    for grid_row in target_grid_rows:
        coordinate = (float_key(grid_row["xi"]), float_key(grid_row["T_MeV"]))
        audit_row = audit_by_coordinate[coordinate]
        key = (reason_class_for_reason(grid_row.get("reason", "")), str(audit_row.get("oracle_shape_observation", "inconclusive")))
        reason_shape_row_counts[key] += 1
        reason_shape_coordinate_sets[key].add(coordinate)
    write_csv(tables_root / "coverage_summary.csv", ("metric", "value", "interpretation"), [
        {"metric": "target_unresolved_rows", "value": len(target_grid_rows), "interpretation": "C2 rows in the limited audit scope"},
        {"metric": "target_unique_coordinates", "value": len(target_coordinates), "interpretation": "deduplicated (xi,T) coordinates"},
        {"metric": "duplicate_unresolved_rows_collapsed", "value": len(target_grid_rows) - len(target_coordinates), "interpretation": "multiple unresolved levels for the same (xi,T) coordinate"},
        {"metric": "raw_curve_valid_finite_converged", "value": len(valid_rows), "interpretation": "raw oracle input validity only"},
        {"metric": "range_fetch_or_parse_failures", "value": len(audit_rows) - len(valid_rows), "interpretation": "must remain zero for a usable diagnostic package"},
        {"metric": "raw_archive_curve_rows", "value": EXPECTED_CURVE_ROWS, "interpretation": "expected per-curve rho coverage"},
        {"metric": "solver_called", "value": False, "interpretation": "solver-free audit"},
    ])
    write_csv(tables_root / "shape_summary.csv", ("oracle_shape_observation", "row_count", "note"), [
        {"oracle_shape_observation": key, "row_count": value, "note": "raw mu_B(rho) topology; not a Maxwell certificate"}
        for key, value in sorted(shape_counts.items())
    ])
    write_csv(tables_root / "shape_threshold_sensitivity.csv", ("negative_slope_threshold", "coordinates_with_min_smoothed_slope_below_threshold", "coordinate_total", "note"), [
        {"negative_slope_threshold": threshold, "coordinates_with_min_smoothed_slope_below_threshold": sum(
            1 for row in audit_rows if finite(row.get("min_smoothed_slope")) and float(row["min_smoothed_slope"]) <= threshold
        ), "coordinate_total": len(audit_rows), "note": "sensitivity of the raw topology observation; production tolerances are not changed"}
        for threshold in (-1.0e-5, -1.0e-4, -1.0e-3, -1.0e-2, -1.0e-1)
    ])
    write_csv(tables_root / "reason_shape_summary.csv", ("reason_class", "oracle_shape_observation", "coordinate_count", "row_count", "attribution_status", "note"), [
        {"reason_class": reason[0], "oracle_shape_observation": reason[1], "coordinate_count": len(reason_shape_coordinate_sets[reason]), "row_count": reason_shape_row_counts[reason],
         "attribution_status": "candidate_stage_c_or_certificate_issue" if reason[1] == "s_shape_present" else "not_explained_by_raw_s_shape",
         "note": "diagnostic attribution only; hybrid/C2 labels are unchanged"}
        for reason in sorted(reason_shape_row_counts)
    ])

    input_sha = sha256_file(curve_index_path)
    figure_manifest = make_figures(output_root, audit_rows, curves, input_sha) if curves else []
    coverage_ok = len(valid_rows) == len(target_coordinates) and not errors
    verdict = "raw_curve_coverage_complete_diagnostic_only" if coverage_ok else "raw_curve_audit_incomplete"
    decision = {
        "schema_version": "pnjl_phase_reference_limited_evidence_audit_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "calculation_sha": CALCULATION_SHA,
        "solver_called": False,
        "reference_write": False,
        "production_override": False,
        "verdict": verdict,
        "promotion_effect": "none",
        "scope": {
            "c2_unresolved_rows": len(target_grid_rows),
            "unique_coordinates": len(target_coordinates),
            "raw_curve_rows_per_coordinate": EXPECTED_CURVE_ROWS,
            "reason_classes": sorted({str(row["reason_class"]) for row in audit_rows}),
        },
        "shape_counts": dict(shape_counts),
        "interpretation_boundary": "A raw S-shape is evidence that a first-order topology is present in the independent curve; it is not a Maxwell/geometry certificate and cannot replace the hybrid result.",
    }
    (output_root / "decision.json").write_text(json.dumps(decision, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    write_csv(tables_root / "claim_ledger.csv", ("claim_id", "claim", "status", "evidence", "boundary"), [
        {"claim_id": "raw_coverage", "claim": "All 761 target coordinates have complete finite/converged 1281-point raw curves", "status": "supported" if coverage_ok else "inconclusive", "evidence": "tables/raw_curve_audit.csv", "boundary": "raw input validity only"},
        {"claim_id": "s_shape_observation", "claim": "Some target raw curves exhibit a discrete +→−→+ topology", "status": "supported" if shape_counts.get("s_shape_present", 0) else "inconclusive", "evidence": "tables/shape_summary.csv and figures/", "boundary": "not a Maxwell certificate"},
        {"claim_id": "stage_c_attribution", "claim": "S-shaped raw curves make Stage-C/certificate insufficiency a plausible explanation for selected unresolved rows", "status": "candidate", "evidence": "tables/reason_shape_summary.csv", "boundary": "does not prove the exact failing component"},
        {"claim_id": "promotion", "claim": "The raw archive authorizes phase-reference promotion", "status": "blocked", "evidence": "decision.json", "boundary": "raw oracle cannot silently override hybrid/C2"},
    ])
    execution_log = (
        "# Phase-reference limited evidence audit v1\n\n"
        f"- Generated (UTC): `{decision['generated_at_utc']}`\n"
        f"- Calculation SHA: `{CALCULATION_SHA}`\n"
        f"- Zenodo record: `{ARCHIVE_RECORD}`\n"
        f"- Declared ZIP SHA-256: `{ARCHIVE_SHA256}`\n"
        f"- C2 scope: `{len(target_grid_rows)}` rows / `{len(target_coordinates)}` unique coordinates\n"
        "- Solver called: `false`\n"
        "- C2 rows, hybrid labels and reference data were not modified.\n"
        f"- Verdict: `{verdict}`\n"
    )
    (output_root / "execution_log.md").write_text(execution_log, encoding="utf-8")
    (output_root / "README.md").write_text(
        "# Phase-reference limited evidence audit v1\n\n"
        "本包对 C2 中 452 个 `ok`、272 个 `hybrid_stage_c_not_converged` 和 39 个 classification-transition rho 记录进行限定范围的 solver-free raw-oracle 审计；763 行去重后为 761 个坐标，其中 2 个坐标含多个 unresolved level。每个坐标从 Zenodo 归档读取完整 `rho=0:0.003125:4` 曲线，检查 1,281 个点、finite/converged、重复 rho 和 `mu_B(rho)` 的离散 `+→−→+` 拓扑。\n\n"
        f"当前 verdict 为 `{verdict}`。`s_shape_present` 只说明独立 raw 曲线存在 S 形拓扑；它不等于 Maxwell 面积/geometry 证书，也不会覆盖 hybrid/C2 标签。因此本包不触发 reference promotion、不重跑 C0/C1/C2，也不调用 equilibrium solver。\n\n"
        "详细结果见 `tables/raw_curve_audit.csv`、`tables/reason_shape_summary.csv`、`tables/shape_threshold_sensitivity.csv`、`tables/claim_ledger.csv`、`decision.json` 和代表性 `figures/`。\n",
        encoding="utf-8",
    )

    output_files: dict[str, str] = {}
    for path in sorted(output_root.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            output_files[path.relative_to(output_root).as_posix()] = sha256_file(path)
    manifest = {
        "schema_version": "pnjl_phase_reference_limited_evidence_audit_v1",
        "generated_at_utc": decision["generated_at_utc"],
        "calculation_sha": CALCULATION_SHA,
        "solver_called": False,
        "reference_write": False,
        "production_override": False,
        "verdict": verdict,
        "source": {
            "zenodo_record": ARCHIVE_RECORD,
            "archive_url": args.archive_url,
            "declared_zip_sha256": ARCHIVE_SHA256,
            "curve_index_path": str(curve_index_path),
            "curve_index_sha256": sha256_file(curve_index_path),
            "central_directory_path": str(central_path),
            "central_directory_sha256": sha256_file(central_path),
            "c2_grid_path": str(grid_path),
            "c2_grid_sha256": sha256_file(grid_path),
            "c2_manifest_path": str(c2_manifest_path),
            "c2_manifest_sha256": sha256_file(c2_manifest_path),
            "c2_validation_path": str(validation_path),
            "c2_validation_sha256": sha256_file(validation_path) if validation_path.is_file() else "",
        },
        "counts": decision["scope"] | {"valid_curves": len(valid_rows), "figure_count": len(figure_manifest)},
        "output_files": output_files,
        "generator": {"path": str(Path(__file__).resolve()), "sha256": sha256_file(Path(__file__).resolve())},
    }
    (output_root / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"verdict": verdict, "counts": manifest["counts"], "output_root": str(output_root)}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
