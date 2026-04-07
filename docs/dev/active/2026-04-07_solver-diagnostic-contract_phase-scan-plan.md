# Solver 诊断契约与 Phase 扫描适配任务单

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在不影响现有两类求解能力返回结构的前提下，为 solver 新增“全候选池诊断契约（含可行标注与失败原因）”，并让 phase/scan 优先消费该契约。

**Architecture:** 采用“新增模式，不改旧模式”策略：新增 `diagnostic_level=:none|:summary|:full` 与候选池诊断对象，仅在新求解模式启用；旧路径保持输出兼容。phase 层通过适配器读取诊断对象，不再各处重复拼装分支诊断语义。

**Tech Stack:** Julia 1.10+、Models solver 主链、phase/scans 子系统、现有 Constraint/ProblemSpec 治理模块。

---

## 范围与非目标

- 范围：`src/models/solver/` 诊断契约输出、`src/models/phase/` 与 `src/models/scans/` 适配消费。
- 非目标：
  - 不改动物理公式与数值目标函数。
  - 不删除旧 PM 专用实现（本期仅保留兼容回退）。
  - 不改变现有旧模式返回字段。

## 约束基线（必须满足）

- [ ] 新增能力仅作用于“新求解模式（候选池诊断模式）”，旧模式输出结构不变。
- [ ] `diagnostic_level=:none` 时额外开销可忽略（不构造 full 候选明细）。
- [ ] 诊断对象包含 phase 需要的最小字段：`attempt_origin`、`seed_source`、`hard_constraint_ok`、`failed_constraints`、`endpoint_cause`、`continuity_distance`。
- [ ] phase/scan 读取统一契约后，可在不依赖 PM 专用拼装逻辑的情况下得到 valid/invalid/unknown 判定依据。

## Chunk 1: Solver 诊断契约定义

### Task 1: 定义统一诊断结构与级别开关

**Files:**
- Modify: `src/models/solver/ProblemSpec.jl`
- Modify: `src/models/solver/ConstraintSolver.jl`
- Modify: `src/models/solver/Solver.jl`
- Test: `tests/unit/models/test_solver_diagnostic_contract.jl`（新建）

- [ ] **Step 1.1: 新建契约单测骨架**
  - 覆盖 `:none/:summary/:full` 三档输出行为。

- [ ] **Step 1.2: 写失败断言（旧模式不变）**
  - 调用旧入口，断言返回字段集与历史一致。

- [ ] **Step 1.3: 写失败断言（新模式有诊断对象）**
  - 断言 `full` 包含候选池与失败原因。

- [ ] **Step 1.4: 实现诊断类型与开关**
  - 新增 `diagnostic_level` 参数解析与默认值。

- [ ] **Step 1.5: 将候选池治理信息映射到诊断对象**
  - 包括 `attempt_origin/seed_index/hard_constraint_ok/failed_constraints`。

- [ ] **Step 1.6: 运行单测并确认通过**
  - Run: `julia --project=. -e 'include("tests/unit/models/test_solver_diagnostic_contract.jl")'`

## Chunk 2: Phase/Scan 适配层迁移

### Task 2: Phase/Scan 优先消费 solver 诊断契约

**Files:**
- Modify: `src/models/phase/PMPhaseDiagnostic.jl`
- Modify: `src/models/phase/PMPhaseSeeds.jl`
- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`
- Test: `tests/integration/models/test_phase_solver_diagnostic_adapter_smoke.jl`（新建）

- [ ] **Step 2.1: 新建适配 smoke 测试骨架**

- [ ] **Step 2.2: 写失败断言（FixedMu 诊断字段可读）**

- [ ] **Step 2.3: 写失败断言（Trho 分支诊断可读）**

- [ ] **Step 2.4: 适配 phase 读取新诊断对象**
  - 仅在新模式下启用；旧路径保留回退。

- [ ] **Step 2.5: 适配 scan 读取 summary/full 诊断级别**

- [ ] **Step 2.6: 运行 smoke 并确认通过**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_phase_solver_diagnostic_adapter_smoke.jl")'`

## Chunk 3: 回归验证与收尾

### Task 3: 兼容性与治理检查

**Files:**
- Modify: `docs/dev/active/2026-04-07_solver-diagnostic-contract_phase-scan-plan.md`
- Test: `tests/integration/models/test_solver_auto_backend_semantic_parity.jl`
- Test: `tests/integration/models/test_ad_implicit_contract_smoke.jl`

- [ ] **Step 3.1: 复跑旧模式关键回归**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_solver_auto_backend_semantic_parity.jl")'`

- [ ] **Step 3.2: 复跑 AD 契约 smoke（防回归）**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_ad_implicit_contract_smoke.jl")'`

- [ ] **Step 3.3: 跑文档治理检查**
  - Run: `julia --project=. scripts/dev/check_docs_consistency.jl`
  - Run: `julia --project=. scripts/dev/check_active_docs_governance.jl`

- [ ] **Step 3.4: 回填执行证据与遗留项**

## DoD

- [ ] 新模式可返回候选池诊断对象（含可行标注与失败原因）。
- [ ] phase/scan 已优先消费统一契约，且可回退旧路径。
- [ ] 旧模式输出结构未变化（兼容测试通过）。
- [ ] 关键回归与治理脚本通过。

## 风险与回退

- 风险 1：诊断对象字段膨胀导致调用端耦合上升。
  - 缓解：分级输出，默认 `:none`；在 phase 侧通过适配器隔离。
- 风险 2：phase 迁移过程中诊断语义不一致。
  - 缓解：保留旧路径回退并做并行对照 smoke。
- 回退策略：按任务粒度回退新增适配；保留 solver 主求解路径不动。
