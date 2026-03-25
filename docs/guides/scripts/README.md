# 脚本入口清单

本目录记录当前仓库推荐给用户直接运行的稳定脚本入口。

目标不是把所有脚本都列进来，而是明确：

- 哪些脚本是稳定入口
- 哪些脚本只是分析/开发/排障工具

---

## 1. 推荐稳定入口

### PNJL

- [scripts/pnjl/run_conserved_charge_susceptibilities.jl](scripts/pnjl/run_conserved_charge_susceptibilities.jl)
  - 守恒荷广义磁化率、累积量、`Ssigma`、`kappa_sigma2` 的统一单点/小范围扫描入口
- [scripts/pnjl/run_tmu_scan.jl](scripts/pnjl/run_tmu_scan.jl)
  - T-μ 参数空间扫描入口
- [scripts/pnjl/calculate_phase_structure.jl](scripts/pnjl/calculate_phase_structure.jl)
  - 相图自动化产线入口（扫描 -> 判据 -> 报告），支持模板配置 + CLI 覆盖

#### 相图最小单命令产线（PNJL）

在仓库根目录执行：

```powershell
julia --project=. scripts/pnjl/calculate_phase_structure.jl --model_kind=PNJL --mode=research --T_min=150 --T_max=150 --T_step=10 --rho_min=0.1 --rho_max=0.3 --rho_step=0.1 --solver_backend=legacy --output_dir=<your_output_dir>
```

说明：

- 默认会自动加载 `config/models/pnjl/phase_pipeline_default.toml`。
- 可以使用 `--config=<path/to/phase_pipeline.toml>` 指定自定义模板。
- CLI 显式参数优先级高于模板（同名键会覆盖）。

最小产物结构（输出目录）：

- `trho_scan.csv`
- `first_order_boundary.csv`
- `spinodal.csv`
- `crossover_line.csv`
- `phase_summary.json`
- `phase_report.md`
- `run_manifest.json`（记录 argv、config_path、config_hash、git_commit、artifact_paths）

### Server

- [scripts/server/server_full.jl](scripts/server/server_full.jl)
  - Web + API 联调用完整入口
- [scripts/server/server.jl](scripts/server/server.jl)
  - 仅 API 入口

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
