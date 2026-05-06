#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "data" / "outputs" / "results" / "relaxtime" / "scan"
OUT_ROOT = PROJECT_ROOT / "data" / "outputs" / "results" / "relaxtime" / "meson_density" / "phase_f_stage_v1"

FILES = {
    "stable": SRC_ROOT / "meson_density_scan_208_220_step2.csv",
    "current_bu": SRC_ROOT / "phase_shift_meson_density_scan_208_220_step2.csv",
    "gbu_reference": SRC_ROOT / "phase_shift_meson_density_scan_gbu_reference_208_220_step2.csv",
    "strict_bw_stage2": SRC_ROOT / "strict_bw_meson_density_scan_stage2_208_220_step2_converged.csv",
}


def copy_with_metadata(src: Path, dst: Path) -> int:
    rows = 0
    with src.open("r", encoding="utf-8") as fin, dst.open("w", encoding="utf-8", newline="") as fout:
        for line in fin:
            fout.write(line)
            s = line.strip()
            if s and not s.startswith("#") and "," not in s:
                continue
            if s and not s.startswith("#") and "," in s:
                break
        reader = csv.reader(fin)
        writer = csv.writer(fout)
        for row in reader:
            writer.writerow(row)
            if row:
                rows += 1
    return rows


def main() -> None:
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    manifest_lines = [
        "# meson_density phase-f stage bundle",
        "",
        f"output_dir: `{OUT_ROOT.relative_to(PROJECT_ROOT)}`",
        "",
        "source files:",
    ]
    for key, src in FILES.items():
        dst = OUT_ROOT / src.name
        rows = copy_with_metadata(src, dst)
        manifest_lines.append(f"- `{key}`: `{src.relative_to(PROJECT_ROOT)}` -> `{dst.relative_to(PROJECT_ROOT)}` ({rows} data rows)")

    manifest_lines += [
        "",
        "note:",
        "",
        "- This directory is the first formal meson-density result bundle promoted out of `scan/` for downstream baseline/regression governance.",
        "- It remains a stage-level bundle, not the final paper-production directory.",
    ]
    (OUT_ROOT / "README.md").write_text("\n".join(manifest_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
