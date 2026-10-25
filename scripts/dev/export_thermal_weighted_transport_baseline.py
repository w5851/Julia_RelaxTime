"""Promote already qualified values to a versioned regression fixture.

This does not calculate or fit values. Historical v1 mode-B baselines stay intact.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import zipfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def production_values(qualified_rows, comparisons):
    """Keep independent qualification targets while selecting actual CLI values.

    The original fixture fixes a locally seeded background. A production
    fixture must instead describe the stable-background production entrypoint.
    This selection never changes either acceptance or regression tolerances.
    """
    targets = {(r["point_id"], r["observable"]): r for r in qualified_rows}
    if len(targets) != len(qualified_rows):
        raise ValueError("duplicate qualified target")
    selected = {}
    for row in comparisons:
        key = (row["point_id"], row["observable"])
        if key in selected or key not in targets:
            raise ValueError("duplicate or unknown production comparison")
        target = targets[key]
        candidate, reference, value = (float(row[k]) for k in
            ("qualified_candidate", "qualified_reference", "production"))
        if (candidate, reference) != (float(target["candidate"]), float(target["reference"])):
            raise ValueError("production comparison changes a qualified target")
        if not (row["passed"].lower() == "true" and value > 0 and
                abs(value / candidate - 1) <= 1e-4 and abs(value / reference - 1) <= 0.01):
            raise ValueError("production value is outside the qualified acceptance limits")
        selected[key] = {**target, "candidate": row["production"],
                         "qualified_candidate": target["candidate"]}
    if set(selected) != set(targets):
        raise ValueError("production replay does not cover every qualified target")
    return [selected[(r["point_id"], r["observable"])] for r in qualified_rows]


def load_helper(root, relative, name):
    spec = importlib.util.spec_from_file_location(name, root / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def replay_artifact_commit(directory, manifest, source_archive=None):
    """Verify an archived source before accepting the scan's ``unknown`` commit.

    Numerical sources are separately checked against the current checkout by
    ``verified_sources``. The archived package and every retained snapshot must
    also agree; its worktree base must never be used as a calculation commit.
    """
    archive = manifest.get("source_archive")
    if archive is None:
        head = manifest.get("head")
        if source_archive is not None or not isinstance(head, str) or not re.fullmatch(r"[0-9a-f]{40}", head):
            raise ValueError("production replay needs a Git head or a verified source archive")
        return head, None
    if (manifest.get("head") is not None or archive.get("calculation_commit") is not None
            or archive.get("uncommitted_source") is not True):
        raise ValueError("uncommitted source archive cannot declare a calculation commit")
    directory = directory.resolve()
    snapshot = (directory / archive["package_manifest_snapshot"]).resolve()
    if not snapshot.is_relative_to(directory) or digest(snapshot) != archive["package_manifest_sha256"]:
        raise ValueError("source package manifest snapshot changed")
    package = json.loads(snapshot.read_text(encoding="utf-8-sig"))
    for field in ("head_base", "source_branch", "calculation_commit", "uncommitted_source"):
        if package.get(field) != archive.get(field):
            raise ValueError(f"source package identity differs from the replay: {field}")
    files = package.get("source_files", {})
    if not files or len(files) != package.get("files") or len(files) != archive.get("verified_files"):
        raise ValueError("incomplete archived source file map")
    path = Path(source_archive if source_archive is not None else archive["path"])
    if digest(path) != archive["sha256"] or archive["sha256"] != package["zip_sha256"]:
        raise ValueError("frozen source archive SHA256 mismatch")
    with zipfile.ZipFile(path) as frozen:
        entries = [entry for entry in frozen.infolist() if not entry.is_dir()]
        if len(entries) != len(files) or {entry.filename for entry in entries} != set(files):
            raise ValueError("source archive file set differs from its manifest")
        for entry in entries:
            if hashlib.sha256(frozen.read(entry)).hexdigest() != files[entry.filename]:
                raise ValueError(f"archived source file changed: {entry.filename}")
    seen = set()
    for entry in manifest["source_files"]:
        relative = entry["path"].replace("\\", "/")
        frozen = (directory / "source" / relative).resolve()
        if (relative in seen or not frozen.is_relative_to(directory / "source")
                or entry["sha256"] != files.get(relative) or digest(frozen) != entry["sha256"]):
            raise ValueError(f"replay source snapshot differs from its archive: {relative}")
        seen.add(relative)
    if not seen:
        raise ValueError("missing replay source snapshots")
    return "unknown", {
        "sha256": archive["sha256"],
        "package_manifest_sha256": archive["package_manifest_sha256"],
        "source_files_verified": len(files),
        "head_base": archive.get("head_base"),
        "calculation_commit": None,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("qualified_dir", type=Path)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_csv", type=Path)
    parser.add_argument("--production-replay", type=Path,
                        help="select qualified real production values instead of the locally seeded candidate")
    parser.add_argument("--source-archive", type=Path,
                        help="relocated immutable source.zip for an archived production replay")
    args = parser.parse_args()
    if args.source_archive is not None and args.production_replay is None:
        parser.error("--source-archive requires --production-replay")
    root = Path(__file__).resolve().parents[2]
    summary_path = args.qualified_dir / "summary.json"
    comparisons_path = args.qualified_dir / "comparisons.csv"
    input_path = args.input_dir / "review_points.csv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if not (summary["joint_target_passed"] and summary["completed_points"] == 12
            and summary["comparisons"] == 108 and digest(input_path) == summary["input_sha256"]):
        raise ValueError("incomplete qualification or changed input")
    with input_path.open(encoding="utf-8", newline="") as stream:
        points = {row["point_id"]: row for row in csv.DictReader(stream)}
    with comparisons_path.open(encoding="utf-8", newline="") as stream:
        comparisons = list(csv.DictReader(stream))
    observables = {"tau_u", "tau_d", "tau_s", "tau_ubar", "tau_dbar", "tau_sbar",
                   "eta_over_s", "zeta_over_s", "sigma_over_T"}
    grouped = {}
    rows = []
    for row in comparisons:
        point = points[row["point_id"]]
        value, reference = float(row["candidate"]), float(row["reference"])
        if not (row["passed"].lower() == "true" and value > 0 and reference > 0
                and abs(value / reference - 1) <= 0.01):
            raise ValueError("failed or nonfinite qualification row")
        grouped.setdefault(row["point_id"], []).append(row["observable"])
        rows.append({**{key: point[key] for key in ("point_id", "T_MeV", "muB_MeV", "xi")},
                     "observable": row["observable"], "candidate": row["candidate"],
                     "reference": row["reference"]})
    if len(grouped) != 12 or any(len(v) != 9 or set(v) != observables for v in grouped.values()):
        raise ValueError("missing or duplicate point/observable")
    profile = root / "config/models/pnjl/transport_sampling/thermal_weighted_v1.toml"
    metadata = {
        "schema_version": 1,
        "purpose": "native thermal-weighted regression; supersedes active use of legacy mode-B v1 values",
        "profile_sha256": digest(profile),
        "qualified_directory": str(args.qualified_dir.resolve()),
        "sources": [{"path": str(p.resolve()), "sha256": digest(p)}
                    for p in (summary_path, comparisons_path, input_path)],
        "configurations": summary["configurations"],
        "points": 12, "comparisons": 108,
        "regression_rtol": 1e-6, "regression_atol": 1e-8,
        "independent_relative_target": 0.01,
    }
    if args.production_replay is not None:
        replay = args.production_replay.resolve()
        status_path, manifest_path = replay / "status.json", replay / "manifest.json"
        status = json.loads(status_path.read_text(encoding="utf-8"))
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if not (status.get("complete") and status.get("production_replay_passed") and
                status.get("failures") == 0 and status.get("completed_points") == 12 and
                status.get("points") == 12 and status.get("comparisons") == 108):
            raise ValueError("production entrypoint replay is incomplete")
        inputs = load_helper(root, "scripts/analysis/relaxtime/prepare_transport_production_inputs.py",
                             "transport_production_inputs")
        inputs.verified_sources(replay, manifest, root.resolve())
        reference_hash, _ = inputs.verified_reference(replay, manifest)
        if reference_hash != digest(comparisons_path):
            raise ValueError("production replay uses a different independent qualification")
        if manifest["parameters"]["profile"]["source_sha256"] != digest(profile):
            raise ValueError("production replay profile changed")
        replay_comparisons_path = replay / "comparisons.csv"
        with replay_comparisons_path.open(encoding="utf-8", newline="") as stream:
            replay_comparisons = list(csv.DictReader(stream))
        for row in replay_comparisons:
            inputs.verify_comparison(row)
        rows = production_values(rows, replay_comparisons)
        artifact_commit, archive_identity = replay_artifact_commit(replay, manifest, args.source_archive)
        artifact = load_helper(root, "scripts/relaxtime/check_transport_production_artifact.py",
                               "transport_production_artifact")
        for point_id in grouped:
            folder = replay / "points" / point_id
            artifact.check_artifact(folder, profile, artifact_commit, "production", root)
            actual = artifact.read_rows(folder / "phase_guided_transport_scan.csv")
            if len(actual) != 1:
                raise ValueError("production baseline requires one actual row per point")
            actual = actual[0]
            for row in (r for r in rows if r["point_id"] == point_id):
                if (tuple(float(actual[k]) for k in ("T_MeV", "muB_MeV", "xi")) !=
                        tuple(float(row[k]) for k in ("T_MeV", "muB_MeV", "xi")) or
                        float(actual[row["observable"]]) != float(row["candidate"])):
                    raise ValueError("comparison differs from the actual production output")
        metadata.update({
            "purpose": "stable-background production entrypoint regression; independently qualified against retained joint references",
            "baseline_semantics": "production_entrypoint",
            "production_replay_directory": str(replay),
            "qualification_replay_relative_limit": 1e-4,
            "maximum_qualified_candidate_change": max(
                abs(float(r["candidate"]) / float(r["qualified_candidate"]) - 1) for r in rows),
        })
        metadata["sources"].extend({"path": str(p.resolve()), "sha256": digest(p)}
                                   for p in (status_path, manifest_path, replay_comparisons_path))
        if archive_identity is not None:
            metadata["production_source_archive"] = archive_identity
    if args.output_csv.exists() or args.output_csv.with_suffix(".json").exists():
        raise FileExistsError("choose a new version; never overwrite an existing baseline")
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata["baseline_sha256"] = digest(args.output_csv)
    args.output_csv.with_suffix(".json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Exported {len(rows)} qualified comparisons to {args.output_csv}")


if __name__ == "__main__":
    main()
