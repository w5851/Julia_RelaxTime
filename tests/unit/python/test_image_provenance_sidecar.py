from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from scripts.common.provenance_image import (
    verify_image_provenance_sidecar,
    write_image_provenance_sidecar,
)


REQUIRED_SIDECAR_FIELDS = {
    "schema_version",
    "generated_at_utc",
    "image_path",
    "image_sha256",
    "script_path",
    "command",
    "git_commit",
    "julia_version",
    "input_data_hashes",
}


MINIMAL_PNG_BYTES = bytes.fromhex(
    "89504E470D0A1A0A"
    "0000000D49484452000000010000000108060000001F15C489"
    "0000000A49444154789C6360000000020001E527D4A20000000049454E44AE426082"
)
def test_sidecar_includes_required_fields(tmp_path):
    image_path = tmp_path / "sample.png"
    image_path.write_bytes(MINIMAL_PNG_BYTES)
    sidecar_path = tmp_path / "sample.png.provenance.json"

    write_image_provenance_sidecar(
        image_path=image_path,
        sidecar_path=sidecar_path,
        script_path="tests/unit/python/test_image_provenance_sidecar.py",
        command="python -m pytest tests/unit/python/test_image_provenance_sidecar.py -q",
    )

    payload = json.loads(sidecar_path.read_text(encoding="utf-8"))
    assert REQUIRED_SIDECAR_FIELDS.issubset(payload.keys())
    assert isinstance(payload["input_data_hashes"], list)
    assert payload["input_data_hashes"]
    for entry in payload["input_data_hashes"]:
        assert isinstance(entry, dict)
        assert "path" in entry
        assert "sha256" in entry


def test_verify_fails_after_image_tampering(tmp_path):
    image_path = tmp_path / "tampered.png"
    image_path.write_bytes(MINIMAL_PNG_BYTES)
    sidecar_path = tmp_path / "tampered.png.provenance.json"

    write_image_provenance_sidecar(
        image_path=image_path,
        sidecar_path=sidecar_path,
        script_path="tests/unit/python/test_image_provenance_sidecar.py",
        command="python -m pytest tests/unit/python/test_image_provenance_sidecar.py -q",
    )

    image_path.write_bytes(image_path.read_bytes() + b"tamper")

    assert verify_image_provenance_sidecar(
        image_path=image_path,
        sidecar_path=sidecar_path,
    ) is False


def test_verify_fails_on_missing_required_field(tmp_path):
    image_path = tmp_path / "missing.png"
    image_path.write_bytes(MINIMAL_PNG_BYTES)
    sidecar_path = tmp_path / "missing.png.provenance.json"
    sidecar_path.write_text('{"schema_version":"v1"}', encoding="utf-8")

    assert verify_image_provenance_sidecar(
        image_path=image_path,
        sidecar_path=sidecar_path,
    ) is False


def test_cli_verify_pass_and_fail(tmp_path):
    image_path = tmp_path / "cli.png"
    image_path.write_bytes(MINIMAL_PNG_BYTES)

    write_image_provenance_sidecar(
        image_path=image_path,
        script_path="tests/unit/python/test_image_provenance_sidecar.py",
        command="pytest",
    )

    project_root = Path(__file__).resolve().parents[3]
    script = project_root / "scripts" / "plot_scan_csv.py"

    ok_proc = subprocess.run(
        [sys.executable, str(script), "--verify-provenance", str(image_path)],
        capture_output=True,
        text=True,
        cwd=str(project_root),
    )
    assert ok_proc.returncode == 0
    assert "PASS" in ok_proc.stdout

    image_path.write_bytes(image_path.read_bytes() + b"tamper")
    fail_proc = subprocess.run(
        [sys.executable, str(script), "--verify-provenance", str(image_path)],
        capture_output=True,
        text=True,
        cwd=str(project_root),
    )
    assert fail_proc.returncode != 0
    assert "FAIL" in fail_proc.stdout
