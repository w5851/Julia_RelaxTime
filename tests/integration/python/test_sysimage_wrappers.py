"""Exercise launcher policy with disposable fixtures and a fake Julia executable."""

import json
import os
from pathlib import Path
import platform
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[3]
FINGERPRINT = "a" * 64


@pytest.fixture(params=["powershell", "sh"])
def launcher(request, tmp_path):
    repo = tmp_path / "repo"
    scripts = repo / "scripts" / "dev"
    scripts.mkdir(parents=True)
    build = repo / "build"
    build.mkdir()
    family = {"Windows": "windows", "Darwin": "macos"}.get(platform.system(), "linux")
    arch = {"amd64": "x86_64", "arm64": "aarch64"}.get(platform.machine().lower(), platform.machine().lower())
    suffix = {"windows": "dll", "macos": "dylib", "linux": "so"}[family]
    (build / f"JuliaRelaxTime.{suffix}").touch()
    metadata = {"julia_version": "1.12.5", "platform_family": family, "platform_arch": arch,
                "build_inputs_fingerprint": FINGERPRINT, "git_commit": "historical-commit"}
    meta_path = build / "JuliaRelaxTime.sysimage.json"
    fresh_path = tmp_path / "fresh.json"
    fresh_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    log = tmp_path / "args.txt"
    built = tmp_path / "built.txt"
    env = dict(os.environ, TEST_FINGERPRINT=FINGERPRINT, TEST_META_PATH=meta_path.as_posix(),
               TEST_FRESH_META=fresh_path.as_posix(), TEST_ARGS_LOG=log.as_posix(), TEST_BUILT_LOG=built.as_posix())

    if request.param == "powershell":
        shell = shutil.which("pwsh") or shutil.which("powershell")
        if shell is None:
            pytest.skip("PowerShell unavailable")
        wrapper = scripts / "run_with_sysimage.ps1"
        harness = tmp_path / "launch.ps1"
        harness.write_text('''function global:julia {
    $global:LASTEXITCODE = 0
    if ($args[0] -eq '--version') { 'julia version 1.12.5' }
    elseif ($args -contains '--fingerprint') { $env:TEST_FINGERPRINT }
    elseif ($args | Where-Object { $_ -like '*build_sysimage.jl' }) {
        Copy-Item -LiteralPath $env:TEST_FRESH_META -Destination $env:TEST_META_PATH
        Set-Content -LiteralPath $env:TEST_BUILT_LOG -Value 'built'
    } else { Set-Content -LiteralPath $env:TEST_ARGS_LOG -Value ($args -join "`n") }
}
& $env:TEST_WRAPPER -MismatchPolicy $env:TEST_POLICY -JuliaArgs '-e', 'println(42)'
''', encoding="utf-8")
        command = [shell, "-NoProfile", "-File", str(harness)]
        env["TEST_WRAPPER"] = str(wrapper)
    else:
        git_bash = Path("C:/Program Files/Git/bin/bash.exe")
        shell = str(git_bash) if os.name == "nt" and git_bash.is_file() else shutil.which("sh")
        if shell is None:
            pytest.skip("POSIX shell unavailable")
        wrapper = scripts / "run_with_sysimage.sh"
        stub_dir = tmp_path / "bin"
        stub_dir.mkdir()
        stub = stub_dir / "julia"
        stub.write_text('''#!/usr/bin/env sh
case "$*" in
    --version) printf '%s\\n' 'julia version 1.12.5' ;;
    *--fingerprint) printf '%s\\n' "$TEST_FINGERPRINT" ;;
    *build_sysimage.jl*) cp "$TEST_FRESH_META" "$TEST_META_PATH"; printf 'built' > "$TEST_BUILT_LOG" ;;
    *) printf '%s\\n' "$@" > "$TEST_ARGS_LOG" ;;
esac
''', encoding="utf-8", newline="\n")
        stub.chmod(0o755)
        env["PATH"] = str(stub_dir) + os.pathsep + env["PATH"]
        command = [shell, str(wrapper)]
    shutil.copyfile(ROOT / "scripts" / "dev" / wrapper.name, wrapper)

    def run(policy, updates=None, legacy=False):
        current_meta = dict(metadata, **(updates or {}))
        if legacy:
            del current_meta["build_inputs_fingerprint"]
        meta_path.write_text(json.dumps(current_meta, indent=2), encoding="utf-8")
        env["TEST_POLICY"] = policy
        args = [] if request.param == "powershell" else [f"--mismatch-policy={policy}", "-e", "println(42)"]
        result = subprocess.run(command + args, env=env, capture_output=True, text=True, timeout=30)
        return result, log.read_text(encoding="utf-8-sig") if log.exists() else "", built.exists()

    return run


def test_compatible_image_accepts_historical_commit(launcher):
    result, args, rebuilt = launcher("strict")
    assert result.returncode == 0, result.stderr
    assert "--sysimage=" in args
    suffix = {"Windows": "dll", "Darwin": "dylib"}.get(platform.system(), "so")
    assert f"JuliaRelaxTime.{suffix}" in args
    assert not rebuilt


@pytest.mark.parametrize("policy", ["fallback", "strict", "rebuild"])
def test_changed_inputs_follow_policy(launcher, policy):
    result, args, rebuilt = launcher(policy, {"build_inputs_fingerprint": "b" * 64})
    if policy == "strict":
        assert result.returncode != 0
        assert not args
    else:
        assert result.returncode == 0, result.stderr
        assert ("--sysimage=" in args) == (policy == "rebuild")
    assert rebuilt == (policy == "rebuild")


def test_legacy_metadata_requires_migration(launcher):
    result, args, _ = launcher("strict", legacy=True)
    assert result.returncode != 0
    assert "build_inputs_fingerprint" in result.stdout + result.stderr
    assert not args


@pytest.mark.parametrize("updates", [{"julia_version": "1.10.0"}, {"platform_family": "wrong"}, {"platform_arch": "wrong"}])
def test_runtime_compatibility_still_required(launcher, updates):
    result, args, _ = launcher("strict", updates)
    assert result.returncode != 0
    assert not args
