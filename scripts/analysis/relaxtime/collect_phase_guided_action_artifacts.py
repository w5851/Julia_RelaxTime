#!/usr/bin/env python3
"""Collect one phase-guided transport Action case and write its run manifest.

The collector is fail-closed: it downloads nothing until the expected unique
shard count is present and every selected run completed successfully at the
requested source commit.  Repeated shard attempts are deduplicated by display
title, preferring the newest successful attempt.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from collections import defaultdict
from pathlib import Path
from typing import Any


WORKFLOW = "relaxtime-phase-guided-transport-production.yml"
REPOSITORY = "w5851/Julia_RelaxTime"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-name", required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--expected-count", required=True, type=int)
    parser.add_argument("--download-root", required=True, type=Path)
    parser.add_argument("--manifest-out", required=True, type=Path)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--retry-count", type=int, default=3)
    return parser.parse_args()


def run_gh(args: list[str], retry_count: int) -> str:
    last_error: subprocess.CalledProcessError | None = None
    for attempt in range(1, retry_count + 1):
        try:
            result = subprocess.run(
                ["gh", *args],
                check=True,
                capture_output=True,
                text=True,
                encoding="utf-8",
            )
            return result.stdout
        except subprocess.CalledProcessError as error:
            last_error = error
            if attempt < retry_count:
                time.sleep(float(attempt))
    assert last_error is not None
    detail = (last_error.stderr or last_error.stdout or "").strip()
    raise RuntimeError(f"gh command failed after {retry_count} attempts: {detail}")


def list_runs(case_name: str, retry_count: int) -> list[dict[str, Any]]:
    raw = run_gh(
        [
            "run",
            "list",
            "--workflow",
            WORKFLOW,
            "--limit",
            "100",
            "--json",
            "databaseId,displayTitle,status,conclusion,headSha,url,createdAt",
        ],
        retry_count,
    )
    prefix = f"phase-guided transport {case_name} / "
    return [row for row in json.loads(raw) if str(row["displayTitle"]).startswith(prefix)]


def shard_label(display_title: str) -> str:
    parts = display_title.split(" / ")
    if len(parts) != 3:
        raise ValueError(f"unexpected phase-guided run title: {display_title}")
    return f"{parts[1]} / {parts[2]}"


def select_unique_runs(runs: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, list[int]]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for run in runs:
        groups[shard_label(str(run["displayTitle"]))].append(run)

    selected: list[dict[str, Any]] = []
    duplicates: dict[str, list[int]] = {}
    for label, attempts in groups.items():
        attempts.sort(key=lambda row: str(row["createdAt"]), reverse=True)
        successful = [
            row
            for row in attempts
            if row["status"] == "completed" and row["conclusion"] == "success"
        ]
        selected.append(successful[0] if successful else attempts[0])
        if len(attempts) > 1:
            duplicates[label] = [int(row["databaseId"]) for row in attempts]
    selected.sort(key=lambda row: shard_label(str(row["displayTitle"])))
    return selected, duplicates


def result_artifact(run_id: int, case_name: str, retry_count: int) -> dict[str, Any]:
    raw = run_gh(
        ["api", f"repos/{REPOSITORY}/actions/runs/{run_id}/artifacts"],
        retry_count,
    )
    artifacts = json.loads(raw).get("artifacts", [])
    matches = [
        row
        for row in artifacts
        if not row.get("expired", False)
        and str(row.get("name", "")).startswith("relaxtime-phase-guided-transport-results-")
        and case_name in str(row.get("name", ""))
    ]
    if len(matches) != 1:
        raise ValueError(
            f"run {run_id}: expected one unexpired result artifact, found {len(matches)}"
        )
    return matches[0]


def download_artifact(
    run_id: int,
    artifact_name: str,
    download_root: Path,
    retry_count: int,
) -> None:
    run_dir = download_root / str(run_id)
    existing_scans = list(run_dir.rglob("phase_guided_transport_scan.csv")) if run_dir.exists() else []
    if len(existing_scans) == 1:
        return
    if run_dir.exists():
        raise FileExistsError(
            f"partial or ambiguous download directory already exists: {run_dir}"
        )
    run_dir.parent.mkdir(parents=True, exist_ok=True)
    run_gh(
        [
            "run",
            "download",
            str(run_id),
            "--name",
            artifact_name,
            "--dir",
            str(run_dir),
        ],
        retry_count,
    )


def main() -> int:
    args = parse_args()
    runs = list_runs(args.case_name, args.retry_count)
    selected, duplicates = select_unique_runs(runs)
    if len(selected) != args.expected_count:
        raise SystemExit(
            f"expected {args.expected_count} unique shards, found {len(selected)}"
        )

    invalid = [
        row
        for row in selected
        if row["status"] != "completed"
        or row["conclusion"] != "success"
        or row["headSha"] != args.source_commit
    ]
    if invalid:
        summary = [
            {
                "run_id": row["databaseId"],
                "label": shard_label(str(row["displayTitle"])),
                "status": row["status"],
                "conclusion": row["conclusion"],
                "head_sha": row["headSha"],
            }
            for row in invalid
        ]
        print(json.dumps({"ready": False, "invalid_runs": summary}, indent=2))
        return 2

    manifest: list[dict[str, Any]] = []
    for run in selected:
        run_id = int(run["databaseId"])
        artifact = result_artifact(run_id, args.case_name, args.retry_count)
        if args.download:
            download_artifact(
                run_id,
                str(artifact["name"]),
                args.download_root.resolve(),
                args.retry_count,
            )
        manifest.append(
            {
                "RunId": run_id,
                "Label": shard_label(str(run["displayTitle"])),
                "Status": run["status"],
                "Conclusion": run["conclusion"],
                "ArtifactName": artifact["name"],
                "ArtifactSize": artifact["size_in_bytes"],
                "HeadSha": run["headSha"],
                "Url": run["url"],
                "CreatedAt": run["createdAt"],
            }
        )

    args.manifest_out.parent.mkdir(parents=True, exist_ok=True)
    with args.manifest_out.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
    print(
        json.dumps(
            {
                "ready": True,
                "case_name": args.case_name,
                "source_commit": args.source_commit,
                "run_count": len(manifest),
                "duplicates": duplicates,
                "downloaded": args.download,
                "manifest": str(args.manifest_out),
            },
            indent=2,
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
