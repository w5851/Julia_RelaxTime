# PNJL求解器解耦 Wave-E 任务单

> 执行主线说明：本任务单用于 Wave-E 的唯一执行主线（勾选、证据、验收）。
> 设计与实现参考：
> - `docs/superpowers/specs/2026-04-01-pnjl-solver-decoupling-wave-e-design.md`
> - `docs/superpowers/plans/2026-04-01-pnjl-solver-decoupling-wave-e-implementation-plan.md`
> - `docs/superpowers/specs/2026-04-01-pnjl-solver-decoupling-wave-d-design.md`
> - `docs/superpowers/plans/2026-04-01-pnjl-solver-decoupling-wave-d-implementation-plan.md`

## 1. 目标

- [x] 完成 compat 退场，仅保留统一治理层单一路径文件集合。
- [x] 提供单一通用扫描脚本入口并替换现有分散扫描脚本入口角色。
- [x] 将 unified 扫描能力扩展到当前已实现模型族：`PNJL/NJL/RPNJL/PNJLMagnetic/Rotation/GasLiquid`。
- [x] 明确 `pnjl_aniso` 作为参数化模式（`xi/profile`），而非独立 `model_kind`。

## 2. 范围

### 2.1 本期范围

- [x] PR #45 收尾并合并（作为 Wave-E 前置门）。
- [x] 单一通用扫描 CLI（子命令化）落地。
- [x] 模型能力矩阵（model_kind × scan_mode）验证与错误口径统一。
- [x] 分散 legacy 扫描脚本退场（remove/archive，阈值化执行）。
- [x] 迁移映射、任务证据与 docs/api 同步。

### 2.2 非范围

- [x] 不更换底层求解后端。
- [x] 不新增新的模型物理实现。
- [x] 不做与入口治理无关的大规模重构。

## 3. 任务分解

### E0：PR #45 收口前置门

- [x] 检查 PR #45 review/checks 状态并记录。
- [x] 若有 review 缺口，最小修复并避免 scope 扩散。
- [x] 合并 PR #45 并回填证据。

### E1：单一通用扫描脚本入口

- [x] 先补失败 integration 测试（单入口 CLI 契约）。
- [x] 落地统一脚本（建议 `scripts/models/run_unified_scan.jl`），采用子命令：`scan tmu` / `scan trho` / `workflow phase`。
- [x] 保持参数校验和错误提示可追溯、可复现。

### E2：全模型族扫描能力并轨

- [x] 先补失败矩阵测试（model_kind × scan_mode）。
- [x] 支持模型族：`PNJL/NJL/RPNJL/PNJLMagnetic/Rotation/GasLiquid`。
- [x] 对不支持组合提供统一错误口径，不允许 silent fallback。
- [x] 将 `pnjl_aniso` 以参数化方式纳入（`xi/profile`），不新增独立 `model_kind`。
- [x] 统一参数容器契约采用 NamedTuple + 语义分层（`control/background/model/numerics`），见 `docs/architecture/models_param_container_contract.md`。

### E3：compat 退场执行（hard-deprecated -> removed/archived）

- [x] 先补失败测试，约束 legacy 脚本退场行为。
- [x] 在阈值满足后移除或归档：
  - `scripts/pnjl/run_tmu_scan.jl`
  - `scripts/pnjl/run_dense_trho_scan.jl`
  - `scripts/pnjl/run_adaptive_trho_scan.jl`
- [x] 保持 migration 状态可查询，明确 `removed/archived`。

### E4：文档治理与证据同步

- [x] 更新 Wave-E 迁移映射台账：`docs/dev/active/2026-04-01_统一扫描入口兼容退场映射表_Wave-E.md`。
- [x] 回填删除门槛满足证据、风险与回退说明。
- [x] 若稳定公共入口行为变化，同步 `docs/api/`。

### E5：验证矩阵与收口

- [x] 运行定向 Wave-E integration/regression。
- [x] 运行 unit/integration/regression smoke。
- [x] 运行 docs governance checks。
- [x] 运行 `scripts/dev/check_unit_skip_policy.jl`。
- [x] 回填任务单勾选与证据。

## 4. 验收标准

- [x] 统一扫描入口成为唯一稳定脚本路径，分散 legacy 扫描脚本完成退场。
- [x] 模型能力矩阵覆盖当前实现模型族，并具备稳定测试证据。
- [x] `pnjl_aniso` 已以参数化模式纳入统一入口，不再依赖并行脚本范式。
- [x] migration 映射与 docs/api 口径一致，governance 检查通过。

## 5. 验证命令

- [x] `julia --project=. -e 'include("tests/integration/models/test_wavee_unified_scan_cli_smoke.jl")'`
- [x] `julia --project=. -e 'include("tests/regression/pnjl/test_wavee_unified_scan_stability.jl")'`
- [x] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [x] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [x] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [x] `julia --project=. scripts/dev/check_unit_skip_policy.jl`
- [x] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [x] `julia --project=. scripts/dev/check_active_docs_governance.jl`
- [x] PR #45 merge 状态核验（`gh pr view 45 --json state,mergeStateStatus,reviewDecision`）

## 6. DoD

- [x] 任务项与验收项全部勾选。
- [x] 关键验证命令可复现通过。
- [x] Wave-E 迁移台账与 docs/api 同步更新。
- [x] PR #45 收尾状态已记录并与 Wave-E 前置门一致。

## 7. 执行记录

- [x] 2026-04-01：创建 Wave-E 三件套草案（active/spec/plan），目标为“compat 彻底退场 + 单一路径文件集合 + 全模型族统一扫描入口”。
- [x] 2026-04-01：E1 测试先行完成。先新增 `tests/integration/models/test_wavee_unified_scan_cli_smoke.jl` 并验证初次失败（缺失 `scripts/models/run_unified_scan.jl`）；随后最小实现统一 CLI（`scan tmu` / `scan trho` / `workflow phase`）并修复 phase 参数透传；定向验证通过：`julia --project=. -e 'include("tests/integration/models/test_wavee_unified_scan_cli_smoke.jl")'`（7 passed）。
- [x] 2026-04-01：E2 测试先行完成。新增 `tests/integration/models/test_wavee_model_scan_matrix_smoke.jl`，初次失败暴露 `run_tmu_scan` 仅支持 `PNJL/RPNJL`；随后在 `src/models/scans/TmuScan.jl`、`src/models/scans/TrhoScan.jl` 扩展模型族校验为 `PNJL/NJL/RPNJL/PNJLMagnetic/Rotation/GasLiquid`，并统一 `pnjl_aniso` 错误口径（要求使用 `PNJL + profile/xi` 参数化）；定向验证通过：`julia --project=. -e 'include("tests/integration/models/test_wavee_model_scan_matrix_smoke.jl")'`（42 passed），并复跑 `julia --project=. -e 'include("tests/integration/models/test_wavee_unified_scan_cli_smoke.jl")'`（11 passed）。
- [x] 2026-04-01：E3 测试先行完成。新增 `tests/integration/models/test_wavee_compat_retirement_smoke.jl`，初次失败暴露 legacy 脚本仍存在且 migration 状态仍为 `hard_deprecated/compat_adapter`；随后执行退场实现：删除 `scripts/pnjl/run_tmu_scan.jl`、`scripts/pnjl/run_dense_trho_scan.jl`、`scripts/pnjl/run_adaptive_trho_scan.jl`，并在 `src/models/entrypoints.jl` 将对应迁移状态更新为 `removed`、`route=:unified_cli`、`removal_wave=:E`；定向验证通过：`julia --project=. -e 'include("tests/integration/models/test_wavee_compat_retirement_smoke.jl")'`（15 passed），并回归校验 `test_wavec_scan_workflow_parity_smoke.jl` 及 `test_run_tmu_scan_*` 单元兼容测试均通过。
- [x] 2026-04-01：E4 文档治理同步完成。已更新迁移映射台账状态与证据；治理校验通过：`julia --project=. scripts/dev/check_docs_consistency.jl`（OK）、`julia --project=. scripts/dev/check_active_docs_governance.jl`（OK）。`docs/api/` 本步无需增量更新（公共 API 函数签名未变，变更集中在脚本入口治理与迁移状态）。
- [x] 2026-04-01：E5 验证矩阵完成。定向：`test_wavee_unified_scan_cli_smoke.jl`（11 passed）、`test_wavee_unified_scan_stability.jl`（33 passed）；smoke：unit（781 passed）、integration（353 passed）、regression（462 passed, 1 broken，属于既有可选 fixture skip）；governance：`check_unit_skip_policy`（OK）、`check_docs_consistency`（OK）、`check_active_docs_governance`（OK）。
- [x] 2026-04-01：PR #45 状态复核完成：`gh pr view 45 --json state,mergeStateStatus,reviewDecision` 返回 `state=MERGED`。
