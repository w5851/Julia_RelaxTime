# PNJL求解器解耦 Wave-E 任务单

> 执行主线说明：本任务单用于 Wave-E 的唯一执行主线（勾选、证据、验收）。
> 设计与实现参考：
> - `docs/superpowers/specs/2026-04-01-pnjl-solver-decoupling-wave-e-design.md`
> - `docs/superpowers/plans/2026-04-01-pnjl-solver-decoupling-wave-e-implementation-plan.md`
> - `docs/superpowers/specs/2026-04-01-pnjl-solver-decoupling-wave-d-design.md`
> - `docs/superpowers/plans/2026-04-01-pnjl-solver-decoupling-wave-d-implementation-plan.md`

## 1. 目标

- [ ] 完成 compat 退场，仅保留统一治理层单一路径文件集合。
- [ ] 提供单一通用扫描脚本入口并替换现有分散扫描脚本入口角色。
- [ ] 将 unified 扫描能力扩展到当前已实现模型族：`PNJL/NJL/RPNJL/PNJLMagnetic/Rotation/GasLiquid`。
- [ ] 明确 `pnjl_aniso` 作为参数化模式（`xi/profile`），而非独立 `model_kind`。

## 2. 范围

### 2.1 本期范围

- [ ] PR #45 收尾并合并（作为 Wave-E 前置门）。
- [ ] 单一通用扫描 CLI（子命令化）落地。
- [ ] 模型能力矩阵（model_kind × scan_mode）验证与错误口径统一。
- [ ] 分散 legacy 扫描脚本退场（remove/archive，阈值化执行）。
- [ ] 迁移映射、任务证据与 docs/api 同步。

### 2.2 非范围

- [ ] 不更换底层求解后端。
- [ ] 不新增新的模型物理实现。
- [ ] 不做与入口治理无关的大规模重构。

## 3. 任务分解

### E0：PR #45 收口前置门

- [ ] 检查 PR #45 review/checks 状态并记录。
- [ ] 若有 review 缺口，最小修复并避免 scope 扩散。
- [ ] 合并 PR #45 并回填证据。

### E1：单一通用扫描脚本入口

- [ ] 先补失败 integration 测试（单入口 CLI 契约）。
- [ ] 落地统一脚本（建议 `scripts/models/run_unified_scan.jl`），采用子命令：`scan tmu` / `scan trho` / `workflow phase`。
- [ ] 保持参数校验和错误提示可追溯、可复现。

### E2：全模型族扫描能力并轨

- [ ] 先补失败矩阵测试（model_kind × scan_mode）。
- [ ] 支持模型族：`PNJL/NJL/RPNJL/PNJLMagnetic/Rotation/GasLiquid`。
- [ ] 对不支持组合提供统一错误口径，不允许 silent fallback。
- [ ] 将 `pnjl_aniso` 以参数化方式纳入（`xi/profile`），不新增独立 `model_kind`。

### E3：compat 退场执行（hard-deprecated -> removed/archived）

- [ ] 先补失败测试，约束 legacy 脚本退场行为。
- [ ] 在阈值满足后移除或归档：
  - `scripts/pnjl/run_tmu_scan.jl`
  - `scripts/pnjl/run_dense_trho_scan.jl`
  - `scripts/pnjl/run_adaptive_trho_scan.jl`
- [ ] 保持 migration 状态可查询，明确 `removed/archived`。

### E4：文档治理与证据同步

- [ ] 更新 Wave-E 迁移映射台账：`docs/dev/active/2026-04-01_统一扫描入口兼容退场映射表_Wave-E.md`。
- [ ] 回填删除门槛满足证据、风险与回退说明。
- [ ] 若稳定公共入口行为变化，同步 `docs/api/`。

### E5：验证矩阵与收口

- [ ] 运行定向 Wave-E integration/regression。
- [ ] 运行 unit/integration/regression smoke。
- [ ] 运行 docs governance checks。
- [ ] 运行 `scripts/dev/check_unit_skip_policy.jl`。
- [ ] 回填任务单勾选与证据。

## 4. 验收标准

- [ ] 统一扫描入口成为唯一稳定脚本路径，分散 legacy 扫描脚本完成退场。
- [ ] 模型能力矩阵覆盖当前实现模型族，并具备稳定测试证据。
- [ ] `pnjl_aniso` 已以参数化模式纳入统一入口，不再依赖并行脚本范式。
- [ ] migration 映射与 docs/api 口径一致，governance 检查通过。

## 5. 验证命令

- [ ] `julia --project=. -e 'include("tests/integration/models/test_wavee_unified_scan_cli_smoke.jl")'`
- [ ] `julia --project=. -e 'include("tests/regression/pnjl/test_wavee_unified_scan_stability.jl")'`
- [ ] `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] `julia --project=. scripts/dev/check_unit_skip_policy.jl`
- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/check_active_docs_governance.jl`
- [ ] PR #45 merge 状态核验（`gh pr view 45 --json state,mergeStateStatus,reviewDecision`）

## 6. DoD

- [ ] 任务项与验收项全部勾选。
- [ ] 关键验证命令可复现通过。
- [ ] Wave-E 迁移台账与 docs/api 同步更新。
- [ ] PR #45 收尾状态已记录并与 Wave-E 前置门一致。

## 7. 执行记录

- [x] 2026-04-01：创建 Wave-E 三件套草案（active/spec/plan），目标为“compat 彻底退场 + 单一路径文件集合 + 全模型族统一扫描入口”。
