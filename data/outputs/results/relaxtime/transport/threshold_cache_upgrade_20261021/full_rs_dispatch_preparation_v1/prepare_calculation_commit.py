"""Refresh this batch's explicit calculation-commit preview without staging it."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import tomllib

ROOT = Path(__file__).resolve().parents[7]
E = "data/outputs/results/relaxtime/transport/threshold_cache_upgrade_20261021"
FIGURES = "data/outputs/figures/relaxtime/transport/threshold_cache_upgrade_20261021"
HERE = Path(__file__).resolve().parent
LEDGER = "config/governance/task_tracks.toml"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def git(*args, env=None, data=None):
    return subprocess.check_output(["git", "-C", str(ROOT), *args], input=data, env=env)


def main():
    original = read(HERE / "calculation_scope_preview.json")
    base = git("rev-parse", "HEAD").decode().strip()
    require(base == original["head_base"], "HEAD changed since the reviewed scope")
    require(not git("diff", "--cached", "--name-only"), "real staging area must be empty")
    real_index = Path(git("rev-parse", "--git-path", "index").decode().strip())
    if not real_index.is_absolute():
        real_index = ROOT / real_index
    index_before = digest(real_index)
    old_frozen = {p: checksum for p, checksum in original["paths"].items()
                  if checksum and (p.startswith((E + "/", FIGURES + "/", "tests/baselines/")))}
    for relative, checksum in old_frozen.items():
        require(digest(ROOT / relative) == checksum, f"retained evidence changed: {relative}")
    paths = set(original["paths"])
    paths.update(("benchmark/relaxtime/transport_runner_timing.jl",
                  "tests/unit/relaxtime/test_transport_runner_timing.jl",
                  "scripts/analysis/relaxtime/review_transport_runner_pilot.py",
                  "tests/unit/python/test_transport_runner_pilot.py",
                  "docs/analysis/relaxtime/README.md"))
    for name in ("receive_reviewed_evidence.py", "bind_reviewed_inputs.py", "prepare_calculation_commit.py",
                 "prepare_runner_pilot.py", "collect_runner_pilot.py"):
        paths.add(f"{E}/full_rs_dispatch_preparation_v1/{name}")
    additions = {
        "full_rs_mode_b_qualification_setup_v1": ["prepare_source_package.py", "run_remote_qualification.ps1",
            "review_remote_qualification.py", "returned/runner_inputs/runner_inputs.json"],
        "full_rs_runner_timing_preflight_setup_v1": ["package_manifest.json", "control_manifest.json", "overlay.zip",
            "run_preflight.ps1", "review_preflight.py", "review_preflight_v2.py", "review_preflight_v2_manifest.json",
            "returned/runtime_timing.toml"],
        "full_rs_runner_review_checks_v1": ["source.zip", "package_manifest.json", "run_checks.ps1",
            "actionlint.json", "actionlint.log"],
    }
    for case, names in additions.items():
        names += ["returned_evidence.zip", "returned_archive_manifest.json", "local_receipt.json",
                  "returned/reviewed_result.json"]
        for name in names:
            relative = f"{E}/{case}/{name}"
            require((ROOT / relative).is_file(), f"required received evidence is missing: {relative}")
            paths.add(relative)
    qualification = read(ROOT / E / "full_rs_mode_b_qualification_setup_v1/returned/reviewed_result.json")
    timing = read(ROOT / E / "full_rs_runner_timing_preflight_setup_v1/returned/reviewed_result.json")
    reviewer = read(ROOT / E / "full_rs_runner_review_checks_v1/returned/reviewed_result.json")
    require(qualification["mode_b_transport_representatives_qualified"] and timing["timing_contract_passed"] and
            reviewer["passed"], "remote qualification or timing checks have not passed")
    ledger_current = (ROOT / LEDGER).read_text(encoding="utf-8-sig").replace("\r\n", "\n")
    ledger_base = git("show", f"{base}:{LEDGER}").decode("utf-8").replace("\r\n", "\n")
    marker = '\n[[tracks]]\nid = "transport-threshold-upgrade"'
    require(ledger_current.count(marker) == 1 and marker not in ledger_base,
            "unexpected transport ledger layout")
    tail = ledger_current[ledger_current.index(marker):]
    added = tomllib.loads(tail)
    require(len(added["tracks"]) == 1 and all(item["track_id"] == "transport-threshold-upgrade"
            for item in added["items"]), "shared ledger tail includes another track")
    ledger = re.sub(r'^updated_at = "[^"]+"', 'updated_at = "2026-10-26"', ledger_base,
                    count=1, flags=re.MULTILINE).rstrip() + "\n" + tail.rstrip() + "\n"
    parsed = tomllib.loads(ledger)
    require(parsed["primary_track"] == tomllib.loads(ledger_base)["primary_track"], "primary track changed")
    version = "v3"
    prefix = HERE / f"calculation-preview-{version}"
    index = prefix.with_suffix(".index")
    require(not index.exists() and not prefix.with_suffix(".json").exists(), "preview already exists")
    environment = os.environ.copy()
    environment["GIT_INDEX_FILE"] = str(index)
    git("read-tree", base, env=environment)
    pathfile = prefix.with_suffix(".paths")
    pathfile.write_bytes(b"\0".join(p.encode("utf-8") for p in sorted(paths)) + b"\0")
    git("add", "-f", "--pathspec-from-file=" + str(pathfile), "--pathspec-file-nul", env=environment)
    ledger_path = HERE / f"ledger_for_calculation_{version}.toml"
    ledger_path.write_bytes(ledger.encode("utf-8"))
    blob = git("hash-object", "-w", "--stdin", data=ledger.encode("utf-8")).decode().strip()
    git("update-index", "--add", "--cacheinfo", f"100644,{blob},{LEDGER}", env=environment)
    tree = git("write-tree", env=environment).decode().strip()
    module_spec = importlib.util.spec_from_file_location("transport_inputs", ROOT / "scripts/analysis/relaxtime/prepare_transport_production_inputs.py")
    helper = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(helper)
    package = read(ROOT / E / "full_rs_mode_b_qualification_setup_v1/package_manifest.json")
    sources = {path: checksum for path, checksum in package["source_files"].items() if helper.is_numerical_source(path)}
    committed_paths = git("ls-tree", "-r", "--name-only", "-z", tree).decode().split("\0")
    require({p for p in committed_paths if helper.is_numerical_source(p)} == set(sources), "numerical file set changed")
    eol_only = []
    for relative, checksum in sources.items():
        current = (ROOT / relative).read_bytes()
        require(hashlib.sha256(current).hexdigest() == checksum, f"qualified source changed: {relative}")
        committed = git("show", f"{tree}:{relative}")
        require(current.replace(b"\r\n", b"\n") == committed.replace(b"\r\n", b"\n"),
                f"preview numerical source differs: {relative}")
        if current != committed:
            eol_only.append(relative)
    frozen = {p for p in paths if p.startswith((E + "/", FIGURES + "/", "tests/baselines/"))}
    for relative in sorted(frozen):
        require(git("show", f"{tree}:{relative}") == (ROOT / relative).read_bytes(),
                f"preview changed frozen bytes: {relative}")
    git("diff", "--cached", "--check", env=environment)
    patch = git("diff", "--cached", "--binary", "--full-index", "--no-textconv", "--no-ext-diff", env=environment)
    prefix.with_suffix(".patch").write_bytes(patch)
    changed = set(git("diff", "--cached", "--name-only", "-z", env=environment).decode().split("\0")) - {""}
    require(changed == paths | {LEDGER}, "preview includes unselected or unchanged paths")
    require(digest(real_index) == index_before and git("rev-parse", "HEAD").decode().strip() == base,
            "real index or HEAD changed")
    report = dict(head_base=base, candidate_tree=tree, changed_paths=len(changed),
                  numerical_files_verified=len(sources), exact_numerical_files=len(sources)-len(eol_only),
                  eol_only_numerical_files=eol_only, frozen_files_verified=len(frozen),
                  prior_frozen_files_preserved=len(old_frozen), patch_sha256=hashlib.sha256(patch).hexdigest(),
                  actual_index_sha256=index_before, actual_index_unchanged=True, actual_head_unchanged=True,
                  commit_created=False, remote_checks_passed=True, primary_track=parsed["primary_track"],
                  paths={p: digest(ROOT / p) if (ROOT / p).is_file() else None for p in sorted(paths)},
                  ledger_sha256=digest(ledger_path))
    prefix.with_suffix(".json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "paths"}, indent=2))


if __name__ == "__main__":
    main()
