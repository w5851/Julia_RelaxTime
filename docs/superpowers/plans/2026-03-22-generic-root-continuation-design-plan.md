# Generic Root Continuation Engine Design Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在不立即改动生产实现的前提下，产出一份零知识可读的“通用连续性求解引擎”设计与迁移计划，支持后续 PR 分阶段落地。

**Architecture:** 采用 80/20 策略：抽象通用根求解内核（continuation 状态机、quality gate、fallback pipeline、诊断输出），并通过领域适配器承载具体物理问题（残差定义、物理性判据、分支语义）。首批覆盖 Meson 混合态与 PNJL Gap 两条主线。

**Tech Stack:** Julia 1.10+, NLsolve, 现有 `src/models/solver/*` 与 `src/models/workflows/*` 模块，Test stdlib，docs/dev + docs/superpowers 计划文档体系。

---

## 范围与边界（零知识前提）

### 明确要解决

- [ ] 统一“连续性求解”骨架能力：
  - continuation（上一点 seed 追踪）
  - quality gate（收敛 + 残差 + 可选领域约束）
  - fallback（newton -> trust_region -> 多初值）
  - 统一 root diagnostics（quality tag、attempt trace）
- [ ] 支持不同维度变量（2D、5D、8D）与不同扫描自变量（T、mu、rho、混合参数）。
- [ ] 允许分支类型差异：单分支、light-heavy 双分支、多分支字典。

### 明确不在本计划内解决

- [ ] 不承诺“一套 API 无配置适配所有未来物理模型”。
- [ ] 不在本 PR 改写全部历史求解路径。
- [ ] 不在本 PR 完成混合态两阶段求根（该属于后续实现任务）。

---

## 现状盘点（为后续实现提供上下文）

### 已有可复用能力

**Files:**
- `src/models/solver/ImplicitSolver.jl`
- `src/models/solver/SeedStrategies.jl`
- `src/models/workflows/MesonMassWorkflow.jl`

- [ ] `ImplicitSolver` 已具备：
  - `_nlsolve_with_tr_fallback`（方法回退）
  - `residual_norm_max` 门禁
  - 候选解选择逻辑 `_choose_candidate`
- [ ] `SeedStrategies` 已具备：
  - `ContinuitySeed`、`HybridContinuitySeed`、`PhaseAwareContinuitySeed`
- [ ] `MesonMassWorkflow` 已具备（近期新增）：
  - `meson_seed_state`
  - `meson_root_policy`
  - `root_quality`

### 当前重复/割裂点

- [ ] Gap 求解和 Meson 求解各自维护 fallback/quality 逻辑，语义相近但实现分散。
- [ ] 领域判据与求解引擎耦合（例如 Gap 侧依赖 omega/thermo；Meson 侧没有对应量）。
- [ ] 诊断字段命名与质量等级未统一，跨模块难以聚合分析。

---

## 目标抽象（建议接口草案）

### Task 1: 定义最小通用接口（仅设计，不实现）

**Files:**
- Create: `docs/dev/active/2026-03-22_generic_root_engine_design.md`

- [x] 给出核心类型草案：
  - `RootProblemSpec`
  - `RootPolicy`
  - `ContinuationState`
  - `RootAttempt` / `RootDiagnostics`
- [x] 给出最小执行入口草案：
  - `solve_root_with_policy(spec, x0; policy, continuation_state)`
  - `solve_root_continuation(scan_points, spec_factory; policy, tracker)`
- [x] 给出适配器接口草案：
  - `residual!(F, x, ctx)`
  - `domain_quality(x, meta)`（可选）
  - `branch_key(x, ctx)`（可选，用于多分支追踪）

---

## 迁移策略（分阶段）

### Task 2: Phase A 迁移（Meson 优先）

**Files:**
- Modify (future): `src/models/workflows/MesonMassWorkflow.jl`
- Add (future): `src/models/solver/GenericRootEngine.jl`

- [ ] 将 Meson 当前 `_solve_meson_mass_with_policy` 的通用部分迁入 `GenericRootEngine`。
- [ ] 保留 Meson 专属部分在 workflow：
  - 混合态分支语义（light/heavy）
  - threshold/gap 计算与输出格式
- [ ] 验证迁移后行为等价（至少保持现有测试全绿）。

### Task 3: Phase B 迁移（Gap 求解接入）

**Files:**
- Modify (future): `src/models/solver/ImplicitSolver.jl`
- Modify (future): `src/models/solver/SeedStrategies.jl`（若需轻量接口适配）

- [ ] 把 `_nlsolve_with_tr_fallback` 的通用部分改为调用 `GenericRootEngine`。
- [ ] 通过 adapter 注入 Gap 特有判据（omega/physicality/thermo finite）。
- [ ] 确保 FixedMu / FixedRho / FixedAsymmetricRho 回归不退化。

### Task 4: Phase C 统一诊断输出

**Files:**
- Modify (future): `src/models/workflows/MesonMassWorkflow.jl`
- Modify (future): `src/models/solver/ImplicitSolver.jl`
- Add (future): `docs/api/...`（如该输出成为稳定公共字段）

- [ ] 统一质量标签枚举：`good/fallback/degraded/bad`。
- [ ] 统一 attempt trace 结构（方法、初值来源、residual、收敛标志）。
- [ ] 输出兼容层：避免破坏现有调用方。

---

## 风险与防护

- [ ] **风险：** 抽象过度导致 API 复杂、使用方学习成本高。
  - **防护：** 先最小接口，优先覆盖现有 2 条主线，不提前泛化。
- [ ] **风险：** 行为回归（尤其相变附近/混合态近简并点）。
  - **防护：** 建立“异常点回归清单”，作为迁移门禁测试。
- [ ] **风险：** 性能回退（额外层级、诊断开销）。
  - **防护：** policy 中提供 diagnostics 开关，默认轻量。

---

## 验收标准（计划级）

- [ ] 文档明确回答：
  - 什么是通用，什么是领域适配
  - 哪些路径先迁，哪些后迁
  - 如何保证不破坏现有验证基线
- [ ] 文档能被“零上下文工程师”直接用于下一 PR 实施。
- [ ] 明确列出“本 PR 不实现代码”的边界，避免范围漂移。

---

## 下一 PR 执行清单（预告）

- [ ] 新建 `GenericRootEngine.jl`（最小骨架 + 单元测试）
- [ ] 先迁 Meson，再迁 Gap
- [ ] 回归：
  - `tests/validation/relaxtime/test_mixed_meson_root_quality_continuation.jl`
  - 关键 legacy/literature 验证栈
- [ ] 文档同步：
  - `docs/dev/active/*`
  - 如涉及稳定接口，更新 `docs/api/*`
