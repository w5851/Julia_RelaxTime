#!/usr/bin/env python3
"""Generate deterministic full-curve/support-zoom figures for hybrid shadow."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

METHODS = ("production_hybrid", "memoized_dense", "independent_oracle")
XIS = ("-0.5", "0.0", "0.5")
CURVE_T_MATCH_TOL = 1e-6  # trho CSV serializes T to six decimal places
SMOOTH_SLOPE_QUANTILE = 0.05
LOCAL_X_RELATIVE_PADDING = 0.04
LOCAL_X_MIN_PADDING_RHO = 0.01
LOCAL_Y_RELATIVE_PADDING = 0.08
LOCAL_Y_MIN_PADDING_MEV = 0.0002


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _float(value: str | None) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return math.nan


def _status(rows: list[dict[str, str]], xi: str, temperature: float, method: str) -> str:
    for row in rows:
        if row.get("method") == method and row.get("xi") == xi and math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=CURVE_T_MATCH_TOL, rel_tol=0.0):
            return row.get("result_status", "")
    return "missing"


def _select_anchors(rows: list[dict[str, str]]) -> list[tuple[str, float, str]]:
    selected: dict[tuple[str, float], str] = {}
    for xi in XIS:
        oracle = sorted([row for row in rows if row.get("xi") == xi and row.get("method") == "independent_oracle"], key=lambda row: _float(row.get("T_MeV")))
        if not oracle:
            continue
        selected[(xi, _float(oracle[0].get("T_MeV")))] = "low_temperature"
        first = [row for row in oracle if row.get("result_status") == "confirmed_first_order"]
        monotone = [row for row in oracle if row.get("result_status") == "confirmed_monotone"]
        ambiguous = [row for row in oracle if row.get("result_status") == "ambiguous_near_critical"]
        if first:
            selected[(xi, _float(first[-1].get("T_MeV")))] = "first_order"
        if monotone:
            selected[(xi, _float(monotone[0].get("T_MeV")))] = "first_monotone"
        if ambiguous:
            selected[(xi, _float(ambiguous[-1].get("T_MeV")))] = "ambiguous_near_critical"
        for oracle_row in oracle:
            temperature = _float(oracle_row.get("T_MeV"))
            hybrid = _status(rows, xi, temperature, "production_hybrid")
            oracle_status = oracle_row.get("result_status", "")
            if oracle_status == "confirmed_first_order" and hybrid == "ambiguous_near_critical":
                selected[(xi, temperature)] = "oracle_first_order_hybrid_ambiguous"
            elif oracle_status == "ambiguous_near_critical" and hybrid == "confirmed_first_order":
                selected[(xi, temperature)] = "unsupported_hybrid_confirmation"
    return sorted([(xi, temperature, reason) for (xi, temperature), reason in selected.items()], key=lambda item: (XIS.index(item[0]), item[1], item[2]))


def _bounds(curve_rows: list[dict[str, str]], slice_rows: list[dict[str, str]], xi: str, temperature: float) -> tuple[float, float]:
    values: list[float] = []
    for row in slice_rows:
        if row.get("xi") == xi and row.get("method") == "production_hybrid" and math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=CURVE_T_MATCH_TOL, rel_tol=0.0):
            for field in ("rho_hadron", "rho_quark", "rho_spinodal_hadron", "rho_spinodal_quark", "support_low", "support_high"):
                value = _float(row.get(field))
                if math.isfinite(value):
                    values.append(value)
    if len(values) < 2:
        values = [_float(row.get("rho")) for row in curve_rows if math.isfinite(_float(row.get("rho")))]
    if len(values) < 2:
        return 0.0, 4.0
    low, high = min(values), max(values)
    # Keep the phase/support envelope visible while removing the large
    # display-only margin that previously hid most of a weak S-shape.  This
    # is a plotting rule only; it does not change support selection or gates.
    padding = max((high - low) * LOCAL_X_RELATIVE_PADDING, LOCAL_X_MIN_PADDING_RHO)
    return max(0.0, low - padding), min(4.0, high + padding)


def _phase_bound_values(slice_rows: list[dict[str, str]], xi: str, temperature: float) -> list[float]:
    values: list[float] = []
    for row in slice_rows:
        if row.get("xi") == xi and row.get("method") == "production_hybrid" and math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=CURVE_T_MATCH_TOL, rel_tol=0.0):
            for field in ("rho_hadron", "rho_quark", "rho_spinodal_hadron", "rho_spinodal_quark", "support_low", "support_high"):
                value = _float(row.get(field))
                if math.isfinite(value):
                    values.append(value)
    return values


def _quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * fraction
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def _smooth_window(curve_rows: list[dict[str, str]], xi: str, temperature: float) -> tuple[float, float, dict[str, object]]:
    """Find a display-only smooth window when no physical support exists.

    This route is used for monotone slices only as a visualization aid.  It
    reads the production curve, never an oracle label, and does not certify a
    physical state.  The longest contiguous run of the lowest-slope secants is
    selected so the right panel exposes a representative smooth region.
    """

    points = sorted(
        [
            (_float(row.get("rho")), _float(row.get("muq_MeV")))
            for row in curve_rows
            if row.get("method") == "production_hybrid"
            and row.get("xi") == xi
            and math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=CURVE_T_MATCH_TOL, rel_tol=0.0)
        ],
        key=lambda point: point[0],
    )
    points = [(x, y) for x, y in points if math.isfinite(x) and math.isfinite(y)]
    if len(points) < 4:
        return 0.0, 4.0, {"smooth_window_status": "insufficient_points"}
    slopes: list[float] = []
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x1 <= x0:
            slopes.append(math.inf)
        else:
            slopes.append(abs((y1 - y0) / (x1 - x0)))
    finite_slopes = [slope for slope in slopes if math.isfinite(slope)]
    threshold = _quantile(finite_slopes, SMOOTH_SLOPE_QUANTILE)
    qualifying = [math.isfinite(slope) and slope <= threshold for slope in slopes]
    runs: list[tuple[int, int]] = []
    index = 0
    while index < len(qualifying):
        if not qualifying[index]:
            index += 1
            continue
        end = index
        while end + 1 < len(qualifying) and qualifying[end + 1]:
            end += 1
        runs.append((index, end))
        index = end + 1
    if not runs:
        return 0.0, 4.0, {"smooth_window_status": "no_low_slope_run", "smooth_slope_threshold": threshold}
    start, end = max(runs, key=lambda run: (run[1] - run[0], -run[0]))
    low, high = points[start][0], points[end + 1][0]
    return low, high, {
        "smooth_window_status": "selected",
        "smooth_window_method": "longest_low_slope_run",
        "smooth_slope_threshold_MeV_per_rho": threshold,
        "smooth_slope_quantile": SMOOTH_SLOPE_QUANTILE,
        "smooth_window_rho": [low, high],
    }


def _local_y_bounds(points: list[tuple[float, float]]) -> tuple[float, float]:
    """Return an independently padded y range for the rho-local panel.

    The full panel spans hundreds of MeV while CEP-local S-shape excursions
    can be only a few hundredths of a MeV.  Sharing its y-axis would therefore
    render the local curve as a horizontal line.  Keep a small absolute floor
    so a nearly flat control remains inspectable without inventing a physical
    tolerance.
    """

    values = [y for _, y in points if math.isfinite(y)]
    if len(values) < 2:
        return 0.0, 1.0
    low, high = min(values), max(values)
    padding = max((high - low) * LOCAL_Y_RELATIVE_PADDING, LOCAL_Y_MIN_PADDING_MEV)
    return low - padding, high + padding


def _slice_row(slice_rows: list[dict[str, str]], xi: str, temperature: float) -> dict[str, str] | None:
    for row in slice_rows:
        if (
            row.get("xi") == xi
            and row.get("method") == "production_hybrid"
            and math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=CURVE_T_MATCH_TOL, rel_tol=0.0)
        ):
            return row
    return None


def _plot_anchor(
    output: Path,
    curve_rows: list[dict[str, str]],
    slice_rows: list[dict[str, str]],
    xi: str,
    temperature: float,
    reason: str,
) -> dict[str, object]:
    import matplotlib.pyplot as plt

    subset = [row for row in curve_rows if row.get("xi") == xi and math.isclose(_float(row.get("T_MeV")), temperature, abs_tol=CURVE_T_MATCH_TOL, rel_tol=0.0)]
    phase_values = _phase_bound_values(slice_rows, xi, temperature)
    if len(phase_values) >= 2:
        low, high = _bounds(subset, slice_rows, xi, temperature)
        local_policy = "independent_rho_mu_zoom_with_phase_markers_v3_tight_envelope"
        local_metadata: dict[str, object] = {"smooth_window_status": "not_used"}
    else:
        low, high, local_metadata = _smooth_window(subset, xi, temperature)
        local_policy = "smooth_window_rho_mu_zoom_v1"
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=False)
    styles = {
        "production_hybrid": {"color": "tab:orange", "linestyle": "-", "marker": "o"},
        "memoized_dense": {"color": "tab:blue", "linestyle": "--", "marker": "s"},
        "independent_oracle": {"color": "tab:green", "linestyle": ":", "marker": "^"},
    }
    local_points: list[tuple[float, float]] = []
    for method in METHODS:
        points = sorted([row for row in subset if row.get("method") == method], key=lambda row: _float(row.get("rho")))
        finite = [(_float(row.get("rho")), _float(row.get("muq_MeV"))) for row in points]
        finite = [(x, y) for x, y in finite if math.isfinite(x) and math.isfinite(y)]
        if not finite:
            continue
        x, y = zip(*finite)
        label = f"{method} ({_status(slice_rows, xi, temperature, method)})"
        style = styles[method]
        axes[0].plot(x, y, ms=2, lw=0.7, label=label, **style)
        zoom = [(xx, yy) for xx, yy in finite if low <= xx <= high]
        if zoom:
            zx, zy = zip(*zoom)
            local_points.extend(zoom)
            axes[1].plot(zx, zy, ms=3, lw=1.0, label=method, **style)
    axes[0].set_xlim(0.0, 4.0)
    axes[1].set_xlim(low, high)
    local_y_low, local_y_high = _local_y_bounds(local_points)
    axes[1].set_ylim(local_y_low, local_y_high)
    for axis in axes:
        axis.set_xlabel(r"$\rho$")
        axis.grid(alpha=0.2)
    axes[0].set_ylabel(r"$\mu_q$ [MeV]")
    axes[1].set_ylabel(r"$\mu_q$ [MeV]")
    axes[1].ticklabel_format(axis="y", style="plain", useOffset=False)
    phase_row = _slice_row(slice_rows, xi, temperature)
    if phase_row is not None:
        mu_transition = _float(phase_row.get("mu_transition_MeV"))
        if math.isfinite(mu_transition) and local_y_low <= mu_transition <= local_y_high:
            axes[1].axhline(mu_transition, color="black", lw=0.8, ls="-.", label=r"Maxwell $\mu_M$")
        for field, label in (
            ("rho_hadron", r"$\rho_h$"),
            ("rho_quark", r"$\rho_q$"),
            ("rho_spinodal_hadron", r"$\rho_{s,h}$"),
            ("rho_spinodal_quark", r"$\rho_{s,q}$"),
        ):
            density = _float(phase_row.get(field))
            if math.isfinite(density) and low <= density <= high:
                axes[1].axvline(density, color="0.35", lw=0.7, ls=":" if "spinodal" in field else "--", label=label)
    if axes[0].lines:
        axes[0].legend(fontsize=7, loc="best")
    if axes[1].lines:
        axes[1].legend(fontsize=7, loc="best")
    fig.suptitle(f"CEP hybrid shadow: xi={xi}, T={temperature:g} MeV ({reason})")
    fig.tight_layout()
    filename = f"rho_mu_xi_{xi.replace('.', 'p').replace('-', 'm')}_T_{temperature:g}_{reason}.png"
    fig.savefig(output / filename, dpi=150)
    plt.close(fig)
    return {
        "file": filename,
        "xi": xi,
        "T_MeV": temperature,
        "reason": reason,
        "full_xlim": [0.0, 4.0],
        "local_xlim": [low, high],
        "local_ylim": [local_y_low, local_y_high],
        "local_curve_mu_span_MeV": max((y for _, y in local_points), default=math.nan)
        - min((y for _, y in local_points), default=math.nan),
        "right_panel_policy": local_policy,
        **local_metadata,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    try:
        import matplotlib

        matplotlib.use("Agg")
    except ImportError:
        schema = "cep_cascade_production_shadow_v2"
        manifest_path = args.input_dir / "manifest.json"
        if manifest_path.is_file():
            try:
                schema = json.loads(manifest_path.read_text(encoding="utf-8")).get("schema_version", schema)
            except (OSError, json.JSONDecodeError):
                pass
        (args.output_dir / "plot_manifest.json").write_text(json.dumps({"schema_version": schema, "figures": [], "reason": "matplotlib unavailable"}, indent=2) + "\n", encoding="utf-8")
        return 0
    curve_path = args.input_dir / "curve_points.csv"
    slice_path = args.input_dir / "slice_metrics.csv"
    curve_rows, slice_rows = _rows(curve_path), _rows(slice_path)
    figures = []
    for xi, temperature, reason in _select_anchors(slice_rows):
        figures.append(_plot_anchor(args.output_dir, curve_rows, slice_rows, xi, temperature, reason))
    schema = "cep_cascade_production_shadow_v2"
    aggregate_manifest = args.input_dir / "manifest.json"
    if aggregate_manifest.is_file():
        try:
            schema = json.loads(aggregate_manifest.read_text(encoding="utf-8")).get("schema_version", schema)
        except (OSError, json.JSONDecodeError):
            pass
    manifest = {
        "schema_version": schema,
        "source_sha256": {"curve_points.csv": hashlib.sha256(curve_path.read_bytes()).hexdigest(), "slice_metrics.csv": hashlib.sha256(slice_path.read_bytes()).hexdigest()},
        "plot_policy": {
            "full_panel": "rho_0_to_4_with_shared_method_styles",
            "local_panel": "independent_rho_mu_zoom_with_phase_markers_v3_tight_envelope",
            "local_x_padding_rule": "max(4% of phase/support envelope, 0.01 rho)",
            "local_y_padding_rule": "max(8% of local mu span, 0.0002 MeV)",
            "local_y_axis_is_independent": True,
            "no_support_panel": "smooth_window_rho_mu_zoom_v1",
        },
        "figures": figures,
    }
    for figure in figures:
        figure["sha256"] = hashlib.sha256((args.output_dir / figure["file"]).read_bytes()).hexdigest()
    (args.output_dir / "plot_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
