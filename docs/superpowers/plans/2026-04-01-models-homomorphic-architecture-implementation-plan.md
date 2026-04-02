# Models Homomorphic Architecture Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在 Plan-B 稳定后，完成 7 模型边界/骨架同构化，并在最终阶段与 Plan-B 共同收口到“单一新架构主线、无 compat 残留”。

**Architecture:** 先执行 Gate A 前置检查，确认 dimension-agnostic 主线已稳定；随后进入 Gate B 做目录/API/测试模板同构；最后执行 Gate C 做 program-level 联合回归与治理收口。坚持“接口同构、内核异构”的最小公共抽象策略。

**Tech Stack:** Julia 1.10+, `Models` registry/factory, solver/scan/workflow entrypoints, unit/integration/regression + governance scripts.

---

## Program Gates

- Gate A (Prerequisite): Plan-B 已达成（schema 主线稳定 + 旧架构移除）
- Gate B (This plan core): 7 模型同构化（目录/API/测试模板）
- Gate C (Program final): A+B 联合收口（全模型回归 + docs/governance + 无 compat）

## Preflight Checklist (Before Chunk 0)

- [ ] Gate A 完成证据可追溯（命令、输出、结论在 Program 主任务单有记录）。
- [ ] 7 模型目录现状快照已记录（便于对照“骨架同构”改造前后差异）。
- [ ] 外部契约影响清单已冻结（对外 API/导出/产物字段）。
- [ ] 明确“禁止抽象清单”初版并作为代码评审门禁。
- [ ] 迁移执行在隔离分支/工作树中进行。

## Scope

- 模型范围：`NJL/NJL2/PNJL/PNJLMagnetic/RPNJL/Rotation/GasLiquid`
- 改造范围：接口层、目录骨架、测试模板、文档与验收矩阵
- 非范围：物理方程重写、数值策略统一化

## File Responsibility Map

- Modify: `src/models/Models.jl`
  - 责任：统一模型注册/导出契约可发现性，减少散点入口。
- Modify: `src/models/abstract_model.jl`
  - 责任：固化最小公共抽象接口与文档契约。
- Modify: model adapters under `src/models/**/<ModelName>Model.jl`
  - 责任：对齐接口签名与 schema/result contract。
- Modify/Create: `src/models/**/workflows/*.jl`
  - 责任：统一 workflow adapter 形态。
- Create: `tests/unit/models/test_model_interface_homomorphism.jl`
  - 责任：7 模型接口一致性门禁。
- Create: `tests/integration/models/test_model_skeleton_homomorphism_smoke.jl`
  - 责任：7 模型主链路 smoke 验证。
- Create: `tests/regression/models/test_model_homomorphic_regression.jl`
  - 责任：7 模型最小基线回归。
- Modify: `docs/api/models/*`
  - 责任：记录同构化后的统一边界契约。

## Chunk 0: Gate A Precheck (Blocking)

### Task 0: 验证 Plan-B 前置门

**Files:**
- No code changes (verification only)

- [ ] **Step 1: 运行 Gate A 前置检查**
  - Run:
    - `julia --project=. -e 'include("src/constants/Constants_PNJL.jl"); include("src/models/Models.jl"); println(Models.registered_model_kinds())'`
    - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`

- [ ] **Step 2: 记录 Gate A 判定结果**
  - 若未达标：停止本计划执行，回到 Plan-B 补齐。

- [ ] **Step 3: 运行文档治理基线检查**
  - Run:
    - `julia --project=. scripts/dev/check_docs_consistency.jl`
    - `julia --project=. scripts/dev/check_active_docs_governance.jl`

## Chunk 1: Interface Homomorphism

### Task 1: 统一 7 模型最小公共接口签名

**Files:**
- Modify: `src/models/abstract_model.jl`
- Modify: model adapters in `src/models/**/<ModelName>Model.jl`
- Create: `tests/unit/models/test_model_interface_homomorphism.jl`

- [ ] **Step 1: 写失败 unit 测试（7 模型接口一致性）**
  - 覆盖：`solve_gap/calculate_mass_vec/number_densities/gap_state_dim` 的最小契约可调用性。

- [ ] **Step 2: 运行测试确认失败**
  - Run: `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_model_interface_homomorphism.jl"; include("tests/unit/runtests.jl")'`

- [ ] **Step 3: 对齐接口签名与返回契约**
  - 在适配器层统一签名与错误语义。

- [ ] **Step 4: 复跑测试转绿**
  - Run: same as Step 2

- [ ] **Step 5: 提交（小步）**
  - Message style suggestion: `refactor(models): align minimal interface contracts across all model families`

## Chunk 2: Skeleton Homomorphism

### Task 2: 目录骨架与 workflow adapter 同构化

**Files:**
- Modify/Create: `src/models/**/workflows/*.jl`
- Modify: `src/models/Models.jl`
- Create: `tests/integration/models/test_model_skeleton_homomorphism_smoke.jl`

- [ ] **Step 1: 写失败 integration 测试（骨架/链路）**
  - 覆盖 7 模型至少一条 solver/scan/workflow 主链路可用。

- [ ] **Step 2: 运行测试确认失败**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_model_skeleton_homomorphism_smoke.jl")'`

- [ ] **Step 3: 对齐目录骨架与 workflow 薄适配结构**
  - 仅调整组织与边界，不改动物理核心公式。

- [ ] **Step 4: 运行 integration smoke**
  - Run: `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`

- [ ] **Step 5: 提交（小步）**
  - Message style suggestion: `refactor(models/workflows): standardize model skeleton adapters and routing`

## Chunk 3: Regression + Anti-Over-Abstraction Guard

### Task 3: 回归门禁与“禁止抽象清单”

**Files:**
- Create: `tests/regression/models/test_model_homomorphic_regression.jl`
- Modify: `docs/api/models/README.md` (if exists) or nearest API overview docs

- [ ] **Step 1: 写失败 regression 测试（7 模型最小基线）**
  - 每模型 1-2 点，校验可收敛/可产物/关键量范围。

- [ ] **Step 2: 运行 regression smoke 确认失败**
  - Run: `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`

- [ ] **Step 3: 实现最小修复并补充文档“禁止抽象清单”**
  - 禁止项示例：
    - flag-heavy 单函数覆盖全部模型内核
    - 以统一为名删除模型特有物理语义

- [ ] **Step 4: 回归与治理检查**
  - Run:
    - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
    - `julia --project=. scripts/dev/check_docs_consistency.jl`

- [ ] **Step 5: 提交（小步）**
  - Message style suggestion: `test(models): add homomorphic regression matrix across seven model families`

## Chunk 4: Gate C Program Finalization

### Task 4: 联合收口（A+B）

**Files:**
- Modify: docs/specs/plans index references if needed

- [ ] **Step 1: 全套验证**
  - Run:
    - `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
    - `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
    - `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
    - `julia --project=. scripts/dev/check_docs_consistency.jl`
    - `julia --project=. scripts/dev/check_active_docs_governance.jl`

- [ ] **Step 2: compat 残留扫描**
  - 检查是否仍存在旧架构回退入口与临时 shim。

- [ ] **Step 3: 收口提交**
  - Message style suggestion: `refactor(models): finalize homomorphic architecture and remove residual compat branches`

## Done Definition

- Gate A/B/C 全部满足。
- 7 模型在同一边界契约与测试模板下稳定通过。
- 旧架构与 compat 无残留。
- 文档、计划、验证证据一致。
