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