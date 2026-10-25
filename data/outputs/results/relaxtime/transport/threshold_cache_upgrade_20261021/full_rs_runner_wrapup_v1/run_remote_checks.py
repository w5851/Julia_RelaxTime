"""Validate the exact documentation and ledger snapshot on quote-remote."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import re
import subprocess
import sys
from urllib.parse import unquote
import zipfile


root = Path(__file__).resolve().parent
assert os.environ.get("COMPUTERNAME", "").upper() == "W"
manifest = json.loads((root / "input_manifest.json").read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    for name, checksum in manifest["files"].items():
        path = (root / name).resolve()
        assert path.is_relative_to(root) and sha(path) == checksum, name


verify()
source = root / "source"
links = []
for name in manifest["documents"]:
    document = source / name
    for target in re.findall(r"\[[^\]]*\]\(([^)]+)\)", document.read_text(encoding="utf-8-sig")):
        if target.startswith(("https://", "http://", "#", "mailto:")):
            continue
        path = (document.parent / unquote(target.strip("<>").split("#", 1)[0])).resolve()
        assert path.is_relative_to(source) and path.exists(), (name, target)
        links.append(dict(document=name, target=path.relative_to(source).as_posix()))
environment = dict(os.environ, JULIA_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
with (root / "ledger.stdout.log").open("xb") as stdout, (root / "ledger.stderr.log").open("xb") as stderr:
    result = subprocess.run([sys.argv[1], "--startup-file=no", "--history-file=no", "--threads=1",
                             "--project=" + str(source), str(root / "run_ledger_snapshot.jl")],
                            cwd=source, env=environment, stdout=stdout, stderr=stderr, check=False)
verify()
report = dict(host="W", ledger_exit_code=result.returncode, links_checked=len(links),
              working_and_staged_ledgers_passed=result.returncode == 0,
              source_files_unchanged=True, verified_files=len(manifest["files"]),
              input_manifest_sha256=sha(root / "input_manifest.json"),
              numerical_tests_rerun=False, numerical_sources_changed=False,
              reason="Only the offline reviewer, tests and current-state records changed; the fixed calculation and existing convergence/regression remain applicable.",
              completed_utc=datetime.now(timezone.utc).isoformat())
with (root / "reviewed_checks.json").open("x", encoding="utf-8") as stream:
    json.dump(report, stream, indent=2)
selected = [root / name for name in ("reviewed_checks.json", "ledger.stdout.log", "ledger.stderr.log", "input_manifest.json")]
with (root / "return_manifest.json").open("x", encoding="utf-8") as stream:
    json.dump({path.name: sha(path) for path in selected}, stream, indent=2)
selected.append(root / "return_manifest.json")
archive = root / "wrapup-checks-return.zip"
with zipfile.ZipFile(archive, "x", compression=zipfile.ZIP_DEFLATED) as bundle:
    for path in selected:
        bundle.write(path, path.name)
print(json.dumps(dict(report=report, archive_sha256=sha(archive), archive_bytes=archive.stat().st_size), indent=2))
sys.exit(result.returncode)
