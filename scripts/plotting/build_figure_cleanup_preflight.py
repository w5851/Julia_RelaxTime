#!/usr/bin/env python3
"""Build a tracked-only preflight manifest for historical figure cleanup.

This command only reads the worktree. It does not delete, move, rename, or
rewrite any result or figure file. The generated manifest is the review gate
for the separate cleanup PR.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[2]
REGISTRY_PATH = Path("docs/analysis/figure_asset_registry_v1/asset_registry.json")
DECISION_PATH = Path("docs/analysis/figure_asset_registry_v1/author_review_decisions.md")
OUTPUT_PATH = Path("docs/analysis/figure_asset_registry_v1/cleanup_preflight_v1.json")

DELETE_FILES = [
    Path("data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram_prod_v1/figure4_phase_diagram_TmuB_Trho.pdf"),
]

DELETE_ROOTS = [
    Path("data/outputs/figures/pnjl/phase_diagrams"),
    Path("data/outputs/figures/relaxtime/fixed_temperature_xi_scan_muB0"),
    Path("data/outputs/figures/relaxtime/gap_transport_by_xi_muB0p0"),
    Path("data/outputs/figures/relaxtime/gap_transport_by_xi_muB800"),
    Path("data/outputs/figures/relaxtime/meson_density/trho_asymmetric_kplus_piplus_scan_v1/convergence/convergence_low"),
    Path("data/outputs/figures/relaxtime/temperature_scan_muB0_xi0"),
    Path("data/outputs/figures/relaxtime/xi_smoothness_sampling"),
]

MOVE_FILES = {
    **{
        Path("data/outputs/figures/relaxtime/literature") / name: Path("docs/analysis/relaxtime/historical/literature_comparison") / name
        for name in (
            "k_mass_literature_julia_fortranmott_comparison.png",
            "meson_mass_julia_vs_literature_comparison.png",
            "sigma_literature_error_by_process.png",
            "sigma_literature_overlay_by_process.png",
            "tau_literature_comparison.png",
        )
    },
    Path("data/outputs/figures/relaxtime/validation/mixed_identity_track_fortran_swap_compare_muB600.png"):
    Path("docs/analysis/relaxtime/historical/validation/mixed_identity_track_fortran_swap_compare_muB600.png"),
    Path("data/outputs/figures/relaxtime/meson_density/plot_review/freezeout_kminus_piminus_mu_pi_100/overlay_kminus_piminus_mu_pi_100.png"):
    Path("docs/analysis/relaxtime/historical/meson_density/freezeout_kminus_piminus_mu_pi_100_analysis/overlay_kminus_piminus_mu_pi_100.png"),
    Path("data/outputs/figures/relaxtime/meson_density/plot_review/freezeout_kminus_piminus_mu_pi_100/residual_kminus_piminus_mu_pi_100.png"):
    Path("docs/analysis/relaxtime/historical/meson_density/freezeout_kminus_piminus_mu_pi_100_analysis/residual_kminus_piminus_mu_pi_100.png"),
    Path("data/outputs/figures/relaxtime/meson_density/plot_review/freezeout_kminus_piminus_mu_pi_100/plot_manifest.json"):
    Path("docs/analysis/relaxtime/historical/meson_density/freezeout_kminus_piminus_mu_pi_100_analysis/plot_manifest.json"),
}

RENAME_ROOTS = {
    Path("data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_validated_anchored_prod_v1"):
    Path("data/outputs/results/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi005_validated_anchored_prod_v1"),
    Path("data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1_p128_validated_anchored_prod_v1"):
    Path("data/outputs/results/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1_p128_xi005_validated_anchored_prod_v1"),
    Path("data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_validated_anchored_prod_v1"):
    Path("data/outputs/figures/relaxtime/transport/phase_guided/mode_a_fixed_muB_phase_scaled/first_canonical_v1_p128_xi005_validated_anchored_prod_v1"),
    Path("data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1_p128_validated_anchored_prod_v1"):
    Path("data/outputs/figures/relaxtime/transport/phase_guided/mode_b_fixed_T_sparse_muB/first_canonical_v1_p128_xi005_validated_anchored_prod_v1"),
}

UNTRACKED_EXCLUSION_ROOTS = [
    Path("data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram__plotv1__audit"),
    Path("data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram__plotv1__audit__line_first_v1"),
    Path("data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram__plotv1__audit__line_first_v2"),
    Path("data/outputs/figures/pnjl/phase_diagram/figure4_phase_diagram__plotv1__audit__pilot_v2"),
    Path("data/outputs/figures/relaxtime/meson_mass/path_scan/freezeout/default_baseline_freezeout_xi0__plotv1__strict"),
    Path("data/outputs/figures/relaxtime/meson_mass/path_scan/freezeout/default_baseline_freezeout_xi0__plotv1__strict__line_first_v1"),
    Path("data/outputs/figures/relaxtime/meson_mass/path_scan/freezeout/default_baseline_freezeout_xi0__plotv1__strict__line_first_v2"),
    Path("data/outputs/figures/relaxtime/meson_mass/path_scan/freezeout/default_baseline_freezeout_xi0__plotv1__strict__pilot_v2"),
    Path("docs/analysis/plotting_pilot_c1_phase_surface__plotv1__estimated_midpoint"),
    Path("docs/analysis/plotting_pilot_c1_phase_surface__plotv1__estimated_midpoint__line_first_v1"),
    Path("docs/analysis/plotting_pilot_c1_phase_surface__plotv1__estimated_midpoint__line_first_v2"),
    Path("docs/analysis/plotting_pilot_c1_phase_surface__plotv1__estimated_midpoint__pilot_v2"),
    Path("docs/analysis/pnjl/c2_blocking_audit_v2"),
    Path("docs/analysis/pnjl/c1_surface_views/pnjl_c1_mu_xi_T_phase_surfaces_diagnostic_v2"),
    Path("docs/analysis/pnjl/c1_surface_views/pnjl_c1_xi_t_mu_phase_surfaces_diagnostic_v1"),
    Path("scripts/analysis/pnjl/build_c2_blocking_audit_v2.py"),
    Path("scripts/analysis/pnjl/plot_c1_cep_bracket.py"),
    Path("scripts/analysis/pnjl/plot_c1_xi_t_mu_phase_surfaces.py"),
    Path("tmp"),
]

HISTORICAL_REFERENCE_PATHS = {
    REGISTRY_PATH,
    Path("docs/analysis/figure_asset_registry_v1/cleanup_candidates.csv"),
    DECISION_PATH,
    OUTPUT_PATH,
    Path("scripts/plotting/build_figure_cleanup_preflight.py"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT_PATH,
    )
    return parser.parse_args()


def rel(path: Path) -> str:
    return path.as_posix()


def absolute(path: Path) -> Path:
    return PROJECT_ROOT / path


def sha256(path: Path) -> str | None:
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_files() -> set[Path]:
    completed = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=PROJECT_ROOT,
        check=True,
        stdout=subprocess.PIPE,
    )
    return {Path(value) for value in completed.stdout.decode("utf-8").split("\0") if value}


def git_head() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=PROJECT_ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return completed.stdout.strip()


def file_record(
    path: Path,
    tracked: set[Path],
    *,
    action: str,
    operation_id: str,
    destination: Path | None = None,
    content_policy: str = "byte_preserve",
) -> dict[str, Any]:
    target = absolute(path)
    record: dict[str, Any] = {
        "operation_id": operation_id,
        "action": action,
        "path": rel(path),
        "path_type": "file" if target.is_file() else "missing",
        "exists": target.is_file(),
        "tracked": path in tracked,
        "size_bytes": target.stat().st_size if target.is_file() else None,
        "sha256": sha256(target),
        "content_policy": content_policy,
    }
    if destination is not None:
        record["destination"] = rel(destination)
        record["destination_exists"] = absolute(destination).exists()
    return record


def files_under(root: Path) -> list[Path]:
    target = absolute(root)
    if not target.is_dir():
        return []
    return sorted(
        (path.relative_to(PROJECT_ROOT) for path in target.rglob("*") if path.is_file()),
        key=rel,
    )


def reference_scan(tracked: set[Path], old_token: str) -> list[dict[str, Any]]:
    matches: list[dict[str, Any]] = []
    for path in sorted(tracked, key=rel):
        if path in HISTORICAL_REFERENCE_PATHS or path.suffix.lower() in {".png", ".pdf", ".svg", ".jpg", ".jpeg"}:
            continue
        target = absolute(path)
        try:
            text = target.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        count = text.count(old_token)
        if count:
            matches.append(
                {
                    "path": rel(path),
                    "sha256": sha256(target),
                    "match_count": count,
                    "classification": "path_update_candidate",
                }
            )
    return matches


def build_manifest() -> dict[str, Any]:
    tracked = git_files()
    entries: list[dict[str, Any]] = []
    root_records: list[dict[str, Any]] = []

    for path in DELETE_FILES:
        entries.append(
            file_record(
                path,
                tracked,
                action="delete",
                operation_id="delete_explicit_file",
            )
        )

    for root in DELETE_ROOTS:
        members = files_under(root)
        root_records.append(
            {
                "action": "delete_tree",
                "root": rel(root),
                "exists": absolute(root).is_dir(),
                "member_count": len(members),
                "all_members_tracked": all(path in tracked for path in members),
            }
        )
        entries.extend(
            file_record(
                path,
                tracked,
                action="delete",
                operation_id=f"delete_tree::{rel(root)}",
            )
            for path in members
        )

    for source, destination in sorted(MOVE_FILES.items(), key=lambda item: rel(item[0])):
        entries.append(
            file_record(
                source,
                tracked,
                action="move",
                operation_id=f"move::{rel(source.parent)}",
                destination=destination,
            )
        )

    rename_records: list[dict[str, Any]] = []
    for source_root, destination_root in sorted(RENAME_ROOTS.items(), key=lambda item: rel(item[0])):
        members = files_under(source_root)
        rename_records.append(
            {
                "action": "rename_tree",
                "source_root": rel(source_root),
                "destination_root": rel(destination_root),
                "source_exists": absolute(source_root).is_dir(),
                "destination_exists": absolute(destination_root).exists(),
                "member_count": len(members),
                "all_members_tracked": all(path in tracked for path in members),
            }
        )
        for path in members:
            destination = destination_root / path.relative_to(source_root)
            policy = "path_only_metadata" if path.suffix.lower() in {".json", ".md"} else "byte_preserve"
            entries.append(
                file_record(
                    path,
                    tracked,
                    action="rename",
                    operation_id=f"rename_tree::{rel(source_root)}",
                    destination=destination,
                    content_policy=policy,
                )
            )

    old_token = "first_canonical_v1_p128_validated_anchored_prod_v1"
    new_token = "first_canonical_v1_p128_xi005_validated_anchored_prod_v1"
    reference_matches = reference_scan(tracked, old_token)

    exclusions: list[dict[str, Any]] = []
    for root in UNTRACKED_EXCLUSION_ROOTS:
        target = absolute(root)
        members = [path for path in target.rglob("*") if path.is_file()] if target.is_dir() else []
        exclusions.append(
            {
                "path": rel(root),
                "exists": target.exists(),
                "tracked_root": root in tracked,
                "file_count": len(members),
                "reason": "untracked C1/C2/pilot or user worktree asset; excluded from PR B",
            }
        )

    registry = json.loads(absolute(REGISTRY_PATH).read_text(encoding="utf-8"))
    return {
        "schema_version": "figure_cleanup_preflight_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "base_git_commit": git_head(),
        "source_registry": {
            "path": rel(REGISTRY_PATH),
            "sha256": sha256(absolute(REGISTRY_PATH)),
            "asset_count": registry["discovery"]["asset_count"],
            "review_group_count": registry["discovery"]["review_group_count"],
        },
        "author_decision_record": {
            "path": rel(DECISION_PATH),
            "sha256": sha256(absolute(DECISION_PATH)),
            "status": "author_review_complete",
        },
        "execution_gate": {
            "manual_review_required": True,
            "delete_performed": False,
            "move_performed": False,
            "rename_performed": False,
            "tracked_only": True,
            "no_solver_or_numeric_data_changes": True,
        },
        "summary": {
            "delete_explicit_file_count": len(DELETE_FILES),
            "delete_tree_count": len(DELETE_ROOTS),
            "delete_member_count": sum(record["member_count"] for record in root_records),
            "move_file_count": len(MOVE_FILES),
            "rename_tree_count": len(RENAME_ROOTS),
            "rename_member_count": sum(record["member_count"] for record in rename_records),
            "reference_update_candidate_count": len(reference_matches),
            "deferred_untracked_root_count": len(exclusions),
        },
        "delete_trees": root_records,
        "rename_trees": rename_records,
        "entries": entries,
        "reference_update": {
            "old_token": old_token,
            "new_token": new_token,
            "historical_snapshot_paths_excluded": sorted(rel(path) for path in HISTORICAL_REFERENCE_PATHS),
            "candidates": reference_matches,
        },
        "untracked_exclusions": exclusions,
    }


def main() -> None:
    args = parse_args()
    output = absolute(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(build_manifest(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
