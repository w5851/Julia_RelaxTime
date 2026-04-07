# Solver Diagnostic Contract Phase Scan Implementation Plan

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

## 诊断契约规范（实现前先冻结）

### 顶层开关

- `diagnostic_level::Symbol = :none`，允许值：`:none | :summary | :full`。
- `:none`：不构造候选池明细，不增加旧模式字段。
- `:summary`：仅构造每次求解的聚合诊断（不含全量候选逐条明细）。
- `:full`：构造聚合诊断 + 候选池逐条诊断明细。

### 诊断对象最小字段（phase/scan 依赖）

- `attempt_origin::Symbol`：候选来源，如 `:seed`, `:warm_start`, `:fallback`。
- `seed_source::Union{Symbol,Nothing}`：种子来源；无种子时 `nothing`。
- `hard_constraint_ok::Union{Bool,Nothing}`：硬约束判定；未评估时 `nothing`。
- `failed_constraints::Vector{Symbol}`：失败约束名列表；无失败时 `Symbol[]`。
- `endpoint_cause::Union{Symbol,Nothing}`：终止原因，如 `:converged`, `:max_iter`, `:nan_guard`。
- `continuity_distance::Union{Float64,Nothing}`：连续性距离；不可计算时 `nothing`。

### valid/invalid/unknown 判定边界（适配器统一语义）

- `valid`：`hard_constraint_ok == true` 且 `failed_constraints` 为空，且无致命 `endpoint_cause`。
- `invalid`：`hard_constraint_ok == false` 或 `failed_constraints` 非空。
- `unknown`：诊断信息不足（例如 `hard_constraint_ok === nothing`）或仅部分字段可得。

## 约束基线（必须满足）

- [x] 新增能力仅作用于“新求解模式（候选池诊断模式）”，旧模式输出结构不变。
- [x] `diagnostic_level=:none` 时额外开销可忽略（不构造 full 候选明细）。
- [x] 诊断对象至少包含 phase 需要的最小字段。
- [x] phase/scan 读取统一契约后，可在不依赖 PM 专用拼装逻辑的情况下得到 valid/invalid/unknown 判定依据。

## Chunk 1: Solver 诊断契约定义

### Task 1: 定义统一诊断结构与级别开关

**Files:**
- Modify: `src/models/solver/ProblemSpec.jl`
- Modify: `src/models/solver/ConstraintSolver.jl`
- Modify: `src/models/solver/Solver.jl`
- Test: `tests/unit/models/test_solver_diagnostic_contract.jl`（新建）

> 注：本 Task 代码块中的函数名为伪代码占位。落地时请替换为 `src/models/entrypoints.jl` 与 `Models` 暴露的真实入口函数名。

- [x] **Step 1.1: 新建单测文件并写测试集骨架**

```julia
using Test

@testset "solver diagnostic contract" begin
    @testset "legacy mode parity" begin end
    @testset "diagnostic level none" begin end
    @testset "diagnostic level summary" begin end
    @testset "diagnostic level full" begin end
end
```

- [x] **Step 1.2: 先写失败断言（旧模式返回字段不变）**

```julia
result_legacy = solve_problem(problem; mode=:legacy)
@test haskey(result_legacy, :diagnostic) == false
@test Set(keys(result_legacy)) == expected_legacy_keys
```

- [x] **Step 1.3: 运行单测确认失败（红灯）**
  - Run: `julia --project=. -e 'include("tests/unit/models/test_solver_diagnostic_contract.jl")'`
  - Expected: `FAIL`，报错点来自 `:summary/:full` 相关断言或新字段缺失。

- [x] **Step 1.4: 在 ProblemSpec 增加 `diagnostic_level` 参数解析与校验**

```julia
const DIAGNOSTIC_LEVELS = (:none, :summary, :full)

diagnostic_level = get(kwargs, :diagnostic_level, :none)
diagnostic_level in DIAGNOSTIC_LEVELS ||
    throw(ArgumentError("invalid diagnostic_level: $(diagnostic_level)"))
```

- [x] **Step 1.5: 在 Solver/ConstraintSolver 实现诊断对象映射**

```julia
diag = (
    attempt_origin = attempt_origin,
    seed_source = seed_source,
    hard_constraint_ok = hard_constraint_ok,
    failed_constraints = failed_constraints,
    endpoint_cause = endpoint_cause,
    continuity_distance = continuity_distance,
)
```

- [x] **Step 1.6: 写 `summary/full` 行为断言（先写到足够严格）**

```julia
result_summary = solve_problem(problem; mode=:new, diagnostic_level=:summary)
@test haskey(result_summary, :diagnostic)
@test !haskey(result_summary[:diagnostic], :candidates)

result_full = solve_problem(problem; mode=:new, diagnostic_level=:full)
@test haskey(result_full[:diagnostic], :candidates)
@test result_full[:diagnostic][:candidates] isa AbstractVector
```

- [x] **Step 1.7: 运行单测确认通过（绿灯）**
  - Run: `julia --project=. -e 'include("tests/unit/models/test_solver_diagnostic_contract.jl")'`
  - Expected: `PASS`，`solver diagnostic contract` 全部通过。

- [ ] **Step 1.8: 提交 Task 1（小步提交）**（未执行：本轮未收到用户提交指令）
  - Run: `git log -10 --oneline`
  - Expected: 可识别最近提交前缀风格（如 `feat:`/`refactor:`/`docs:`），用于本次提交消息对齐。
  - Run: `git add tests/unit/models/test_solver_diagnostic_contract.jl src/models/solver/ProblemSpec.jl src/models/solver/ConstraintSolver.jl src/models/solver/Solver.jl`
  - Run: `git commit -m "feat: add solver diagnostic contract levels for new mode"`
  - Expected: commit 成功，且仅包含 Task 1 相关文件。

## Chunk 2: Phase/Scan 适配层迁移

### Task 2: Phase/Scan 优先消费 solver 诊断契约

**Files:**
- Modify: `src/models/phase/PMPhaseDiagnostic.jl`
- Modify: `src/models/phase/PMPhaseSeeds.jl`
- Modify: `src/models/scans/TmuScan.jl`
- Modify: `src/models/scans/TrhoScan.jl`
- Test: `tests/integration/models/test_phase_solver_diagnostic_adapter_smoke.jl`（新建）

> 注：本 Task 代码块中的函数名为伪代码占位。落地时请替换为 `src/models/entrypoints.jl` 与 `Models` 暴露的真实入口函数名。

- [x] **Step 2.1: 新建 integration smoke 文件与测试集骨架**

```julia
using Test

@testset "phase solver diagnostic adapter smoke" begin
    @testset "fixed-mu diagnostic readable" begin end
    @testset "trho diagnostic readable" begin end
end
```

- [x] **Step 2.2: 写失败断言（FixedMu 可读新诊断字段）**

```julia
diag = run_fixedmu_case(; mode=:new, diagnostic_level=:summary)
@test haskey(diag, :attempt_origin)
@test haskey(diag, :hard_constraint_ok)
@test haskey(diag, :failed_constraints)
```

- [x] **Step 2.3: 写失败断言（Trho 可读并可判定 unknown）**

```julia
diag = run_trho_case(; mode=:new, diagnostic_level=:summary)
status = infer_phase_status(diag)
@test status in (:valid, :invalid, :unknown)
```

- [ ] **Step 2.4: 运行 integration smoke 确认失败（红灯）**（未单独执行：实现已先落地，后续仅做绿灯验证）
  - Run: `julia --project=. -e 'include("tests/integration/models/test_phase_solver_diagnostic_adapter_smoke.jl")'`
  - Expected: `FAIL`，失败点为 phase/scan 尚未接入统一诊断对象。

- [x] **Step 2.5: 在 phase 侧实现统一适配入口（保留旧路径回退）**

```julia
if has_solver_diagnostic(result)
    return adapt_solver_diagnostic(result[:diagnostic])
else
    return adapt_legacy_pm_diagnostic(result)
end
```

- [x] **Step 2.6: 在 scan 侧按 `summary/full` 读取诊断并透传**
  - `summary`：只读取聚合状态，不假设存在 `candidates`。
  - `full`：允许读取 `candidates`，但调用端不得强依赖其长度语义。

- [x] **Step 2.7: 运行 integration smoke 确认通过（绿灯）**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_phase_solver_diagnostic_adapter_smoke.jl")'`
  - Expected: `PASS`，FixedMu 与 Trho 两组 smoke 通过。

- [ ] **Step 2.8: 提交 Task 2（小步提交）**（未执行：本轮未收到用户提交指令）
  - Run: `git log -10 --oneline`
  - Expected: 可识别最近提交前缀风格，并为本次提交选择同类前缀。
  - Run: `git add tests/integration/models/test_phase_solver_diagnostic_adapter_smoke.jl src/models/phase/PMPhaseDiagnostic.jl src/models/phase/PMPhaseSeeds.jl src/models/scans/TmuScan.jl src/models/scans/TrhoScan.jl`
  - Run: `git commit -m "refactor: route phase scan diagnostics through solver contract"`
  - Expected: commit 成功，且仅包含 Task 2 相关文件。

## Chunk 3: 回归验证与收尾

### Task 3: 兼容性与治理检查

**Files:**
- Modify: `docs/dev/active/2026-04-07_solver-diagnostic-contract_phase-scan-plan.md`
- Test: `tests/integration/models/test_solver_auto_backend_semantic_parity.jl`
- Test: `tests/integration/models/test_ad_implicit_contract_smoke.jl`

- [x] **Step 3.1: 复跑旧模式关键回归（兼容性）**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_solver_auto_backend_semantic_parity.jl")'`
  - Expected: `PASS`，旧模式语义一致。

- [x] **Step 3.2: 复跑 AD 契约 smoke（防回归）**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_ad_implicit_contract_smoke.jl")'`
  - Expected: `PASS`，AD 隐式契约不受影响。

- [x] **Step 3.3: 跑文档治理检查（active 文档合规）**
  - Run: `julia --project=. scripts/dev/check_docs_consistency.jl`
  - Expected: `PASS`，无 docs consistency 错误。

- [x] **Step 3.4: 跑 active docs 治理检查**
  - Run: `julia --project=. scripts/dev/check_active_docs_governance.jl`
  - Expected: `PASS`，当前任务单格式合规。

- [x] **Step 3.5: 回填执行证据与遗留项（只写事实，不写推测）**
  - 在本任务单末尾追加：测试命令、通过时间、失败遗留及后续 owner。

- [ ] **Step 3.6: 提交 Task 3（小步提交）**（未执行：本轮未收到用户提交指令）
  - Run: `git log -10 --oneline`
  - Expected: 可识别最近提交前缀风格，并确保文档提交使用历史存在的前缀样式。
  - Run: `git add docs/dev/active/2026-04-07_solver-diagnostic-contract_phase-scan-plan.md`
  - Run: `git commit -m "docs: tighten execution plan for diagnostic contract rollout"`
  - Expected: commit 成功，且仅包含文档与必要治理修订。

## DoD

- [x] 新模式可返回候选池诊断对象（含可行标注与失败原因）。
- [x] phase/scan 已优先消费统一契约，且可回退旧路径。
- [x] 旧模式输出结构未变化（兼容测试通过）。
- [x] 关键回归与治理脚本通过。
- [x] 执行记录包含每个关键命令及其 Expected/Actual。

## 风险与回退

- 风险 1：诊断对象字段膨胀导致调用端耦合上升。
  - 缓解：分级输出，默认 `:none`；在 phase 侧通过适配器隔离。
- 风险 2：phase 迁移过程中诊断语义不一致。
  - 缓解：保留旧路径回退并做并行对照 smoke。
- 风险 3：`unknown` 判定被误用为 `invalid`。
  - 缓解：适配器集中推导状态，并在测试中显式覆盖 `unknown` 分支。
- 回退策略：按任务粒度回退新增适配；保留 solver 主求解路径不动。

## Execution Log（模板）

> 执行阶段请按 append-only 方式记录；每条命令单独一行，避免事后重写。

| Time (UTC+8) | Task/Step | Command | Expected | Actual | Result |
| --- | --- | --- | --- | --- | --- |
| YYYY-MM-DD HH:MM | 1.3 | `julia --project=. -e 'include("tests/unit/models/test_solver_diagnostic_contract.jl")'` | FAIL（红灯） | 填写实际输出摘要 | PASS/FAIL |
| YYYY-MM-DD HH:MM | 1.7 | `julia --project=. -e 'include("tests/unit/models/test_solver_diagnostic_contract.jl")'` | PASS（绿灯） | 填写实际输出摘要 | PASS/FAIL |
| YYYY-MM-DD HH:MM | 2.4 | `julia --project=. -e 'include("tests/integration/models/test_phase_solver_diagnostic_adapter_smoke.jl")'` | FAIL（红灯） | 填写实际输出摘要 | PASS/FAIL |
| YYYY-MM-DD HH:MM | 2.7 | `julia --project=. -e 'include("tests/integration/models/test_phase_solver_diagnostic_adapter_smoke.jl")'` | PASS（绿灯） | 填写实际输出摘要 | PASS/FAIL |
| YYYY-MM-DD HH:MM | 3.1 | `julia --project=. -e 'include("tests/integration/models/test_solver_auto_backend_semantic_parity.jl")'` | PASS | 填写实际输出摘要 | PASS/FAIL |
| YYYY-MM-DD HH:MM | 3.2 | `julia --project=. -e 'include("tests/integration/models/test_ad_implicit_contract_smoke.jl")'` | PASS | 填写实际输出摘要 | PASS/FAIL |
| YYYY-MM-DD HH:MM | 3.3 | `julia --project=. scripts/dev/check_docs_consistency.jl` | PASS | 填写实际输出摘要 | PASS/FAIL |
| YYYY-MM-DD HH:MM | 3.4 | `julia --project=. scripts/dev/check_active_docs_governance.jl` | PASS | 填写实际输出摘要 | PASS/FAIL |

| 2026-04-07 22:20 | 1.3 | `julia --project=. -e 'include("tests/unit/models/test_solver_diagnostic_contract.jl")'` | FAIL（红灯） | `:summary/:full` 诊断字段断言失败（`haskey(result,:diagnostic)`） | PASS |
| 2026-04-07 22:28 | 1.7 | `julia --project=. -e 'include("tests/unit/models/test_solver_diagnostic_contract.jl")'` | PASS（绿灯） | `solver diagnostic contract | 14/14` | PASS |
| 2026-04-07 22:31 | 2.4 | `julia --project=. -e 'include("tests/integration/models/test_phase_solver_diagnostic_adapter_smoke.jl")'` | FAIL（红灯） | 本轮未单独执行（实现完成后直接绿灯验证） | FAIL |
| 2026-04-07 22:33 | 2.7 | `julia --project=. -e 'include("tests/integration/models/test_phase_solver_diagnostic_adapter_smoke.jl")'` | PASS（绿灯） | `phase solver diagnostic adapter smoke | 4/4` | PASS |
| 2026-04-07 22:35 | 3.1 | `julia --project=. -e 'include("tests/integration/models/test_solver_auto_backend_semantic_parity.jl")'` | PASS | 首次因新增关键字默认值缺失报错；修复后 `solver auto backend semantic parity | 4/4` 与 `problemspec parity guard | 8/8` | PASS |
| 2026-04-07 22:36 | 3.2 | `julia --project=. -e 'include("tests/integration/models/test_ad_implicit_contract_smoke.jl")'` | PASS | `models AD implicit contract smoke | 10/10` | PASS |
| 2026-04-07 22:37 | 3.3 | `julia --project=. scripts/dev/check_docs_consistency.jl` | PASS | `[docs-consistency] OK` | PASS |
| 2026-04-07 22:38 | 3.4 | `julia --project=. scripts/dev/check_active_docs_governance.jl` | PASS | `[active-docs-governance] OK` | PASS |

遗留项：
- Task 1/2/3 的提交步骤（1.8、2.8、3.6）未执行，原因：本轮用户未要求创建 git commit。
