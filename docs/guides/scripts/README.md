# 脚本入口清单

本目录记录当前仓库推荐给用户直接运行的稳定脚本入口。

目标不是把所有脚本都列进来，而是明确：

- 哪些脚本是稳定入口
- 哪些脚本只是分析/开发/排障工具

---

## 1. 推荐稳定入口

全量 `run_*.jl` 脚本功能目录见：

- [run_script_catalog.md](run_script_catalog.md)

该目录用于治理与迁移盘点；本 README 仅保留稳定白名单入口。

### 默认运行规则

稳定 CLI 默认优先通过以下 wrapper 启动：

- `scripts/dev/run_with_sysimage.ps1`
- `scripts/dev/run_with_sysimage.sh`

用途：

- 若本机已有可用 sysimage，则自动追加 `--sysimage=...`
- 若本机没有 sysimage，则仍回退到普通 `julia --project=.` 运行
- 默认 mismatch policy 为 `fallback`
- PowerShell wrapper 可配合 `-MismatchPolicy strict|fallback|rebuild`
- POSIX wrapper 可配合 `--mismatch-policy=strict|fallback|rebuild`
- `-BuildIfMissing` / `--build-if-missing` 仍保留，作为 `rebuild` 别名

如需先获取预构建 sysimage：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/bootstrap_sysimage.ps1
```

```bash
sh scripts/dev/bootstrap_sysimage.sh
```

bootstrap 脚本会根据当前平台、架构和 Julia 版本，解析 GitHub Release 资产名并解包到 `build/`。

最小示例（Windows / PowerShell）：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/models/run_unified_scan.jl scan tmu --model_kind=PNJL --T_values=150 --mu_values=0,100 --xi_values=0.0 --output_path=data/outputs/results/tmu_smoke.csv --overwrite=true
```

最小示例（Linux / macOS）：

```bash
sh scripts/dev/run_with_sysimage.sh scripts/models/run_unified_scan.jl scan tmu --model_kind=PNJL --T_values=150 --mu_values=0,100 --xi_values=0.0 --output_path=data/outputs/results/tmu_smoke.csv --overwrite=true
```

### Phase A 默认 wrapper 白名单

以下稳定 CLI 默认应经 wrapper 启动：

| 类别 | 入口 | Windows 默认 | Linux/macOS 默认 |
|---|---|---|---|
| PNJL 相图 | `scripts/pnjl/calculate_phase_structure.jl` | `run_with_sysimage.ps1` | `run_with_sysimage.sh` |
| 统一扫描 | `scripts/models/run_unified_scan.jl` | `run_with_sysimage.ps1` | `run_with_sysimage.sh` |
| 守恒荷易感性 | `scripts/pnjl/run_conserved_charge_susceptibilities.jl` | `run_with_sysimage.ps1` | `run_with_sysimage.sh` |
| Relaxtime 编排 | `scripts/relaxtime/run_relaxtime_orchestrator.jl` | `run_with_sysimage.ps1` | `run_with_sysimage.sh` |
| 输运扫描 | `scripts/relaxtime/run_gap_transport_scan.jl` | `run_with_sysimage.ps1` | `run_with_sysimage.sh` |
| 服务器入口 | `scripts/server/server_full.jl` | `run_with_sysimage.ps1` | `run_with_sysimage.sh` |

### PNJL

- [scripts/pnjl/run_conserved_charge_susceptibilities.jl](../../../scripts/pnjl/run_conserved_charge_susceptibilities.jl)
  - 守恒荷广义磁化率、累积量、`Ssigma`、`kappa_sigma2` 的统一单点/小范围扫描入口
- [scripts/models/run_unified_scan.jl](../../../scripts/models/run_unified_scan.jl)
  - 统一扫描入口（`scan tmu|trho`），覆盖 T-μ / T-ρ 网格扫描与 `Models` 主链扫描治理
- [scripts/pnjl/calculate_phase_structure.jl](../../../scripts/pnjl/calculate_phase_structure.jl)
  - 相图自动化产线入口（扫描 -> 判据 -> 报告），支持模板配置 + CLI 覆盖

#### 相图最小单命令产线（PNJL）

在仓库根目录执行：

```powershell
powershell -ExecutionPolicy Bypass -File scripts/dev/run_with_sysimage.ps1 scripts/pnjl/calculate_phase_structure.jl --model_kind=PNJL --mode=research --T_min=150 --T_max=150 --T_step=10 --rho_min=0.1 --rho_max=0.3 --rho_step=0.1 --solver_backend=legacy --output_dir=<your_output_dir>
```

说明：

- 默认会自动加载 `config/models/pnjl/phase_pipeline_default.toml`。
- 可以使用 `--config=<path/to/phase_pipeline.toml>` 指定自定义模板。
- 可以使用 `--preset=smoke` 快速切换到轻量可复现实验参数（随后仍可用 CLI 显式参数覆盖）。
- CLI 显式参数优先级高于模板（同名键会覆盖）。

最小产物结构（输出目录）：

- `trho_scan.csv`
- `first_order_boundary.csv`
- `spinodal.csv`
- `crossover_line.csv`
- `phase_summary.json`
- `phase_report.md`
- `run_manifest.json`（记录 argv、config_path、config_hash、git_commit、artifact_paths）

`run_manifest.json` 关键字段：

- `preset`：若使用 `--preset=...`，记录最终采用的预设名。
- `effective_config`：记录本次运行的最终有效参数快照（含被 CLI 覆盖后的值）。

### Server

- [scripts/server/server_full.jl](../../../scripts/server/server_full.jl)
  - Web + API 联调用完整入口
- [scripts/server/server.jl](../../../scripts/server/server.jl)
  - 仅 API 入口

### RelaxTime

- [scripts/relaxtime/run_manual_relaxation_scan_workflow.jl](../../../scripts/relaxtime/run_manual_relaxation_scan_workflow.jl)
  - 手动组合产物入口（`cross_section` / `plan_a` / `plan_b`）
  - 支持 `--base-output-dir` 将结果写到隔离目录（默认 `data/outputs`）

`plan_a` / `plan_b` 目录最小溯源产物：

- 扫描 CSV（行级 `run_id`）
- `effective_config.json`（最终有效参数快照）
- `run_manifest.json`（`argv`、`git_commit`、`config_hash`、`artifacts`、`summary`）

其中 `plan_b_merged.csv` 额外包含：

- `source_file`
- `source_T_MeV`

用于将合并行反向定位到温度分片 CSV。

---

## 2. 不应视为稳定用户入口的目录

- `scripts/dev/`
  - 开发期导出、迁移、一次性工具
- `scripts/analysis/`
  - 后处理与分析脚本
- `scripts/debug/`
  - 排障脚本
- `scripts/perf/`
  - 性能探针与分析

这些目录仍然重要，但默认不应在用户指南里作为正式入口推荐。

---

## 3. 命名约定建议

为了降低“入口太多”的感受，建议以后遵守：

- `run_*.jl`
  - 稳定、可文档化、面向用户的入口
- `calculate_*.jl`
  - 计算型脚本，但不一定是稳定 CLI
- `export_*.jl`
  - 基线/导出/开发辅助
- `diagnose_*.jl`, `debug_*.jl`
  - 排障专用

结论：

- 不是脚本数量本身有问题。
- 问题是稳定入口和内部工具没有层级。
- 这份清单的作用就是把层级固定下来。
