# PNJL 兼容层说明

本页为历史兼容层说明，不是当前主线开发入口。

## 当前主线入口

仓库已统一使用 `Models` 与 `src/models/entrypoints.jl` 暴露 PNJL 相关能力：

- 平衡求解：`Models.solve_gap(...)`
- 扫描入口：`Models.run_tmu_scan(...)`、`Models.run_trho_scan(...)`
- workflow 入口：`Models.solve_gap_and_transport(...)`、`Models.solve_gap_and_meson_point(...)`
- 相图产线：`Models.run_phase_pipeline(...)`、`Models.find_cep(...)`

优先阅读：

- `docs/api/models/solver/README.md`
- `docs/api/models/scans/README.md`
- `docs/api/models/workflows/README.md`
- `docs/api/models/phase/README.md`

## 兼容层定位

- `docs/api/pnjl/*` 仅保留历史语义与迁移背景说明。
- 新功能与新调用方不应以 `PNJL` 兼容层作为默认入口。
- 若需回放历史流程，可参考兼容说明页并转译为 `Models` 等价调用。

## 单位约定

- 核心求解内部使用自然单位 `fm^-1`。
- 面向脚本/HTTP 的参数通常使用 MeV 输入，再在入口处换算。

## 迁移建议

若你仍在使用旧调用风格，建议按以下顺序迁移：

1. 将入口 include/using 切换到 `src/models/Models.jl`。
2. 将扫描与 workflow 调用改为 `Models.run_*` / `Models.solve_*`。
3. 对照 `docs/api/models/*` 完成参数与返回结构的口径对齐。
