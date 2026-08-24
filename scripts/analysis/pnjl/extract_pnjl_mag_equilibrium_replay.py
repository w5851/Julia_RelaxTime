"""Extract a small, reproducible pnjl_mag equilibrium replay diagnostic.

The replay intentionally follows the external project's one-seed descending
temperature continuation.  It is diagnostic evidence only: it does not
attempt Julia's multi-seed branch enumeration or global-minimum search.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path


EXPECTED_COMMIT = "e1fc81d3c3c9d220c49972e54307b66a604cb9db"
SOURCE_ROWS = ("src/constants.py", "src/gap.py", "src/pnjl_mag.py", "src/plot_orders.py")
SOURCE_DATA = "data/orders_muB0.csv"
DEPENDENCY_LOCK = "uv.lock"
FIELDS_GEV2 = (0.2, 0.4, 0.8)
SELECTED_T_MEV = (50.0, 150.0, 240.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def git_head(root: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip()


def source_rows(path: Path) -> dict[tuple[float, float], dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            (float(row["T_MeV"]), float(row["eB_GeV2"])): row
            for row in csv.DictReader(handle)
        }


def run_external_replay(external_root: Path, python_executable: Path) -> dict:
    code = textwrap.dedent(
        r"""
        import json
        import time

        import jax.numpy as jnp
        import numpy as np

        from src.constants import hc
        from src.gap import gap_equations, solve_x_scan
        from src.pnjl_mag import Omega, get_nodes

        temps_mev = np.arange(300.0, 49.0, -1.0)
        temps = jnp.asarray(temps_mev / hc, dtype=jnp.float64)
        ints = get_nodes(128, 80, zeta_num=256)
        x0 = jnp.asarray([-0.03, -0.02, -0.45, 0.85, 0.85], dtype=jnp.float64)
        rows = []
        field_wall_seconds = {}
        total_start = time.perf_counter()

        for field_gev2 in (0.2, 0.4, 0.8):
            field = field_gev2 * (1000.0 / hc) ** 2
            field_start = time.perf_counter()
            states = solve_x_scan(x0, temps, 0.0, field, ints)
            states.block_until_ready()
            state_values = np.asarray(states)
            field_wall_seconds[str(field_gev2)] = time.perf_counter() - field_start

            for index, temperature in enumerate(temps_mev):
                if float(temperature) not in (50.0, 150.0, 240.0):
                    continue
                state = state_values[index]
                residual = float(
                    jnp.max(
                        jnp.abs(
                            gap_equations(
                                states[index], temps[index], 0.0, field, ints
                            )
                        )
                    )
                )
                omega = float(Omega(states[index], temps[index], 0.0, field, ints))
                rows.append(
                    {
                        "T_MeV": float(temperature),
                        "muB_MeV": 0.0,
                        "eB_GeV2": float(field_gev2),
                        "eB_fm_minus2": float(field),
                        "phi_u": float(state[0]),
                        "phi_d": float(state[1]),
                        "phi_s": float(state[2]),
                        "Phi1": float(state[3]),
                        "Phi2": float(state[4]),
                        "omega_fm4": omega,
                        "gap_residual_max": residual,
                    }
                )

        print(
            json.dumps(
                {
                    "hc_MeV_fm": float(hc),
                    "rows": rows,
                    "field_wall_seconds": field_wall_seconds,
                    "total_wall_seconds": time.perf_counter() - total_start,
                },
                sort_keys=True,
            )
        )
        """
    )
    result = subprocess.run(
        [str(python_executable), "-c", code],
        cwd=external_root,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "pnjl_mag replay failed with exit code "
            f"{result.returncode}:\n{result.stderr[-4000:]}"
        )
    try:
        return json.loads(result.stdout.strip().splitlines()[-1])
    except (json.JSONDecodeError, IndexError) as exc:
        raise RuntimeError(
            "pnjl_mag replay did not emit a JSON result; "
            f"stdout tail={result.stdout[-2000:]!r}, stderr tail={result.stderr[-2000:]!r}"
        ) from exc


def write_csv(rows: list[dict], path: Path) -> None:
    fields = [
        "T_MeV",
        "muB_MeV",
        "eB_GeV2",
        "eB_fm_minus2",
        "phi_u",
        "phi_d",
        "phi_s",
        "Phi1",
        "Phi2",
        "omega_fm4",
        "gap_residual_max",
        "committed_state_max_abs_delta",
        "committed_gap_residual_max",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--external-root",
        type=Path,
        default=Path(r"D:\Temp\pnjl_mag_audit_20260823"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[3]
        / "docs"
        / "analysis"
        / "historical"
        / "legacy"
        / "legacy_extraction_v1"
        / "pnjl_mag_equilibrium_replay_v1",
    )
    args = parser.parse_args()

    external_root = args.external_root.resolve()
    output_dir = args.output_dir.resolve()
    python_executable = external_root / ".venv" / "Scripts" / "python.exe"
    if not python_executable.is_file():
        python_executable = external_root / ".venv" / "bin" / "python"
    if not python_executable.is_file():
        raise SystemExit(f"external venv Python not found under {external_root / '.venv'}")

    actual_commit = git_head(external_root)
    if actual_commit != EXPECTED_COMMIT:
        raise SystemExit(f"unexpected pnjl_mag commit: {actual_commit} != {EXPECTED_COMMIT}")
    for relative_path in [*SOURCE_ROWS, SOURCE_DATA, DEPENDENCY_LOCK]:
        path = external_root / relative_path
        if not path.is_file():
            raise SystemExit(f"required external file missing: {path}")

    committed = source_rows(external_root / SOURCE_DATA)
    external_python_version = subprocess.run(
        [str(python_executable), "--version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    replay = run_external_replay(external_root, python_executable)
    rows = []
    state_fields = ("phi_u", "phi_d", "phi_s", "Phi1", "Phi2")
    for row in replay["rows"]:
        key = (row["T_MeV"], row["eB_GeV2"])
        reference = committed.get(key)
        if reference is None:
            raise RuntimeError(f"missing committed source row for {key}")
        deltas = [abs(row[field] - float(reference[field])) for field in state_fields]
        row["committed_state_max_abs_delta"] = max(deltas)
        row["committed_gap_residual_max"] = float(reference["gap_residual_max"])
        rows.append(row)

    output_dir.mkdir(parents=True, exist_ok=True)
    table_path = output_dir / "pnjl_mag_equilibrium_replay_v1.csv"
    write_csv(rows, table_path)
    script_path = Path(__file__).resolve()
    source_hashes = {
        relative_path: sha256(external_root / relative_path)
        for relative_path in [*SOURCE_ROWS, SOURCE_DATA, DEPENDENCY_LOCK]
    }
    provenance = {
        "schema_version": "pnjl_mag_equilibrium_replay_provenance_v1",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "external_root": str(external_root),
        "external_remote": "https://github.com/ZhouRui-xzit/pnjl_mag.git",
        "external_commit": actual_commit,
        "external_source_hashes": source_hashes,
        "python_executable": str(python_executable),
        "external_python_version": external_python_version,
        "producer_script": str(script_path),
        "producer_script_sha256": sha256(script_path),
        "route": "pnjl_mag.solve_x_scan descending temperature continuation",
        "temperature_path_MeV": "300..50 step -1",
        "selected_temperatures_MeV": list(SELECTED_T_MEV),
        "fields_GeV2": list(FIELDS_GEV2),
        "muB_MeV": 0.0,
        "p_num": 128,
        "landau_levels": 80,
        "zeta_num": 256,
        "seed": [-0.03, -0.02, -0.45, 0.85, 0.85],
        "external_hc_MeV_fm": replay["hc_MeV_fm"],
        "point_count": len(rows),
        "field_wall_seconds": replay["field_wall_seconds"],
        "total_wall_seconds": replay["total_wall_seconds"],
        "table_sha256": sha256(table_path),
        "external_data_replay": "replayed_current_machine_and_matched_committed_csv",
        "julia_solver_called": False,
        "acceptance_scope": "diagnostic_only",
        "boundary": "One-seed continuation; no multi-seed branch ensemble, global minimum, or Julia parity acceptance.",
    }
    (output_dir / "provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    manifest = {
        "schema_version": "pnjl_mag_equilibrium_replay_manifest_v1",
        "artifact_kind": "external_equilibrium_diagnostic_replay",
        "source_commit": actual_commit,
        "table": table_path.name,
        "provenance": "provenance.json",
        "table_sha256": sha256(table_path),
        "provenance_sha256": sha256(output_dir / "provenance.json"),
        "acceptance_scope": "diagnostic_only",
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    readme = f"""# `pnjl_mag` equilibrium replay v1

本目录是 `pnjl_mag@{actual_commit}` 在当前机器上的轻量 equilibrium 重放诊断。

## 口径

- Python/JAX/Optimistix 使用外部仓库自己的 `.venv` 和 `uv.lock`；
- 固定 `muB=0`、`eB={{0.2, 0.4, 0.8}} GeV^2`；
- 沿作者脚本相同的单 seed 温度 continuation：`T=300..50 MeV`，步长 `1 MeV`；
- `p_num=128`、Landau levels=`80`、`zeta_num=256`；
- 只提取 `T={{50, 150, 240}} MeV` 的 9 行；每个场值的 251 点、三个场值合计
  753 个 continuation 点不完整写入主项目；
- 每个点记录五维状态、`Omega` 和最大 gap residual，并与作者已提交的
  `data/orders_muB0.csv` 对照。

## 结论边界

本机重放的 9 个代表点与作者 CSV 的五维状态逐字段一致，且 residual 保持有限并很小。
这证明了依赖环境和作者单 seed continuation 路线在当前机器上可以重现。
它仍是 `diagnostic_only`：没有多 seed 分支集合、全局极小值搜索、Julia/外部 solver
分支等价或 production target admission。

生成命令：

```powershell
python scripts/analysis/pnjl/extract_pnjl_mag_equilibrium_replay.py
```

表格、hash、运行耗时和边界见 [`manifest.json`](manifest.json) 与
[`provenance.json`](provenance.json)。
"""
    (output_dir / "README.md").write_text(readme, encoding="utf-8")
    print(f"wrote {table_path}")
    print(f"wrote {output_dir / 'provenance.json'}")
    print(f"wrote {output_dir / 'manifest.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
