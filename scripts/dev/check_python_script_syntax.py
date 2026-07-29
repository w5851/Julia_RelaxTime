#!/usr/bin/env python3
from __future__ import annotations

import py_compile
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
TARGETS = [
    ROOT / "scripts" / "plot_scan_csv.py",
    ROOT / "scripts" / "analysis" / "collect_pnjl_cep_narrow_pilot.py",
    ROOT / "scripts" / "analysis" / "plot_pnjl_cep_narrow_pilot.py",
    ROOT / "scripts" / "analysis" / "freeze_pnjl_cep_narrow_pilot_v2_windows.py",
    ROOT / "scripts" / "analysis" / "collect_pnjl_cep_narrow_pilot_v2.py",
    ROOT / "scripts" / "analysis" / "plot_pnjl_cep_narrow_pilot_v2.py",
    ROOT / "scripts" / "pnjl",
    ROOT / "scripts" / "relaxtime",
]


def iter_python_files() -> list[Path]:
    files: list[Path] = []
    for target in TARGETS:
        if target.is_dir():
            files.extend(sorted(path for path in target.rglob("*.py") if path.is_file()))
        elif target.is_file():
            files.append(target)
    return files


def main() -> int:
    files = iter_python_files()
    failures: list[tuple[Path, str]] = []

    for path in files:
        try:
            py_compile.compile(str(path), doraise=True)
        except py_compile.PyCompileError as exc:
            failures.append((path, str(exc)))

    if failures:
        print(f"[python-script-syntax] FAILED: {len(failures)} file(s)")
        for path, message in failures:
            try:
                rel_path = path.relative_to(ROOT)
            except ValueError:
                rel_path = path
            print(f" - {rel_path}: {message}")
        return 1

    print(f"[python-script-syntax] OK: checked {len(files)} Python files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
