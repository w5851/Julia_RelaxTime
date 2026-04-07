---
title: ModeSchema 与统一残差内核实施任务单
archived: true
original: docs/dev/active/2026-04-07_mode-schema_unified-residual-kernel_followup-plan.md
archived_date: 2026-04-07
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# ModeSchema 与统一残差内核实施任务单

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在不改动物理公式的前提下，将 `models` 中 solver/implicit 残差路径收口为单一 base 内核，并以 schema 驱动命名输入与向量输入转换。

**Architecture:** 采用“单 base 残差 + mode 约束块 + schema 适配层”架构。base 残差统一表达为 `F_base(x_state, mu_vec, params)`，模式差异通过附加约束拼接形成完整残差。用户层 `NamedTuple` 与求解层 `Vector` 的互转集中在 schema adapter，不允许在核心残差内散落字段映射逻辑。

**Tech Stack:** Julia 1.10+、ForwardDiff、StaticArrays、NLsolve、ImplicitDifferentiation、现有 `Models` 求解器模块。

---

## 范围与非目标

- 范围：`src/models/solver/` 中残差构建、schema 适配、implicit 适配与对应测试。
- 非目标：
  - 不改动物理模型公式与参数。
  - 不一次性重写所有 mode 路径。
  - 不做目录重构（本任务不迁移文件目录）。

## 约束基线（必须满足）

- [x] 唯一残差入口：base 残差只保留一个实现源。
- [x] 类型泛化：base 残差保持 Dual 友好，不在核心路径强转 `Float64`。
- [x] 映射集中：`NamedTuple <-> Vector` 仅在 schema adapter 层完成。
- [x] 线性化一致：solver 与 implicit 对同一 base 残差做线性化。
- [x] 容器边界化：`Vector/SVector` 差异仅在适配边界处理。

## Chunk 1: Base 残差单入口化

### Task 1: 抽取 `gap_core_residual!` 并建立最小回归

**Files:**
- Modify: `src/models/solver/Conditions.jl`
- Test: `tests/unit/models/test_gap_core_residual.jl`（新建）

- [x] **Step 1.1: 新建测试文件骨架**
  - 在 `tests/unit/models/test_gap_core_residual.jl` 创建 `@testset "gap_core_residual parity"`。

- [x] **Step 1.2: 准备最小测试输入**
  - 在测试内构造 `model = Models.create_model(:PNJL)`、`x_state`、`mu_vec`、`params`（低开销网格）。

- [x] **Step 1.3: 写首个失败断言**
  - 调用未来接口 `gap_core_residual!` 与现有 `gap_conditions` 比较逐分量 `isapprox`。

- [x] **Step 1.4: 运行单测并确认失败**
  - Run: `julia --project=. -e 'include("tests/unit/models/test_gap_core_residual.jl")'`
  - Expected: FAIL（未定义或断言不通过）。

- [x] **Step 1.5: 在 `Conditions.jl` 声明并导出 `gap_core_residual!`**
  - 先补最小函数签名与 `export`，确保测试可调用。

- [x] **Step 1.6: 用现有 `gap_conditions` 逻辑填充实现**
  - 函数仅做 core block 计算，不拼接 mode 约束。

- [x] **Step 1.7: 处理输出容器写入**
  - 约定 `F` 为外部分配向量，函数内部原地写入。

- [x] **Step 1.8: 再跑单测确认通过**
  - Run: `julia --project=. -e 'include("tests/unit/models/test_gap_core_residual.jl")'`
  - Expected: PASS。

- [ ] **Step 1.9: 记录变更并提交**
  - `git add src/models/solver/Conditions.jl tests/unit/models/test_gap_core_residual.jl`
  - `git commit -m "refactor(models/solver): extract single gap core residual entry"`

### Task 2: 让 `build_residual!` 委托 base 残差

**Files:**
- Modify: `src/models/solver/Conditions.jl`
- Test: `tests/unit/models/test_gap_core_residual.jl`

- [x] **Step 2.1: 在同一测试文件新增第二个测试集**
  - 名称建议：`@testset "build_residual delegates core block"`。

- [x] **Step 2.2: 写失败断言（FixedMu）**
  - 调用 `build_residual!(FixedMu...)` 得到 `F` 后，比较 `F[1:5]` 与 `gap_core_residual!`。

- [x] **Step 2.3: 写失败断言（FixedRho）**
  - 同上，验证前 5 维。

- [ ] **Step 2.4: 运行并确认失败**
  - Run: `julia --project=. -e 'include("tests/unit/models/test_gap_core_residual.jl")'`

- [x] **Step 2.5: 修改 `build_residual!(::FixedMu, ...)` 委托 core**

- [x] **Step 2.6: 修改 `build_residual!(::FixedRho, ...)` 委托 core**

- [x] **Step 2.7: 按同样方式处理 `FixedEntropy/FixedSigma/FixedAsymmetricRho`**

- [x] **Step 2.8: 运行测试确认通过**
  - Run: `julia --project=. -e 'include("tests/unit/models/test_gap_core_residual.jl")'`

- [ ] **Step 2.9: 提交**
  - `git add src/models/solver/Conditions.jl tests/unit/models/test_gap_core_residual.jl`
  - `git commit -m "refactor(models/solver): route mode residual builders through gap core"`

## Chunk 2: implicit/solver 线性化语义对齐

### Task 3: implicit 条件路径对齐同一 base 残差

**Files:**
- Modify: `src/models/solver/ImplicitGapLegacy.jl`
- Modify: `src/models/solver/ImplicitAdapters.jl`
- Test: `tests/integration/models/test_ad_implicit_contract_smoke.jl`
- Test: `tests/integration/models/test_solver_implicit_residual_parity_smoke.jl`（新建）

- [x] **Step 3.1: 新建 parity smoke 测试文件**
  - 文件：`tests/integration/models/test_solver_implicit_residual_parity_smoke.jl`。

- [x] **Step 3.2: 写失败断言（fixed-mu）**
  - 固定 `θ/x`，比较 solver residual 与 implicit conditions residual。

- [x] **Step 3.3: 写失败断言（flavor-mu，可选最小一例）**

- [ ] **Step 3.4: 运行新测试并确认失败**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_solver_implicit_residual_parity_smoke.jl")'`

- [x] **Step 3.5: 修改 `ImplicitGapLegacy.jl` 的 `conditions` 委托共享 core**

- [x] **Step 3.6: 修改 `ImplicitAdapters.jl` 的 NJL 条件路径委托共享 core**

- [x] **Step 3.7: 运行 AD 契约 smoke**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_ad_implicit_contract_smoke.jl")'`

- [x] **Step 3.8: 运行 parity smoke**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_solver_implicit_residual_parity_smoke.jl")'`

- [ ] **Step 3.9: 提交**
  - `git add src/models/solver/ImplicitGapLegacy.jl src/models/solver/ImplicitAdapters.jl tests/integration/models/test_ad_implicit_contract_smoke.jl tests/integration/models/test_solver_implicit_residual_parity_smoke.jl`
  - `git commit -m "refactor(models/solver): align implicit conditions with shared gap core residual"`

## Chunk 3: schema/accessor 收口位置语义

### Task 4: 引入 mode accessor，减少 `x[i]` 扩散

**Files:**
- Modify: `src/models/solver/SchemaAdapter.jl`
- Modify: `src/models/solver/Conditions.jl`
- Test: `tests/unit/models/test_mode_accessor_schema.jl`（新建）

- [x] **Step 4.1: 新建 accessor/schema 单测文件**
  - 文件：`tests/unit/models/test_mode_accessor_schema.jl`。

- [x] **Step 4.2: 写失败断言（FixedMu 索引映射）**

- [x] **Step 4.3: 写失败断言（FixedRho 索引映射）**

- [x] **Step 4.4: 写失败断言（维度错配抛 `ArgumentError`）**

- [ ] **Step 4.5: 运行并确认失败**
  - Run: `julia --project=. -e 'include("tests/unit/models/test_mode_accessor_schema.jl")'`

- [x] **Step 4.6: 在 `SchemaAdapter.jl` 增加最小 accessor API**
  - 例如：`state_view(schema, x)`、`mu_view(schema, x)` 或等价接口。

- [x] **Step 4.7: 在 `Conditions.jl` 接入 `FixedMu` 路径**

- [x] **Step 4.8: 在 `Conditions.jl` 接入 `FixedRho` 路径**

- [x] **Step 4.9: 运行单测确认通过**
  - Run: `julia --project=. -e 'include("tests/unit/models/test_mode_accessor_schema.jl")'`

- [ ] **Step 4.10: 提交**
  - `git add src/models/solver/SchemaAdapter.jl src/models/solver/Conditions.jl tests/unit/models/test_mode_accessor_schema.jl`
  - `git commit -m "refactor(models/solver): centralize mode index semantics via schema accessors"`

## Chunk 4: 闭环验证与文档

### Task 5: 回归、治理检查与任务单收尾

**Files:**
- Modify: `docs/dev/active/2026-04-07_mode-schema_unified-residual-kernel_followup-plan.md`
- Test: `tests/unit/models/test_implicit_gap.jl`
- Test: `tests/integration/models/test_ad_implicit_contract_smoke.jl`
- Test: `tests/integration/models/test_solver_auto_backend_semantic_parity.jl`

- [x] **Step 5.1: 跑 unit 关键用例**
  - `julia --project=. -e 'include("tests/unit/models/test_implicit_gap.jl")'`

- [x] **Step 5.2: 跑 integration AD 契约用例**
  - `julia --project=. -e 'include("tests/integration/models/test_ad_implicit_contract_smoke.jl")'`

- [x] **Step 5.3: 跑 integration 语义等价用例**
  - `julia --project=. -e 'include("tests/integration/models/test_solver_auto_backend_semantic_parity.jl")'`

- [x] **Step 5.4: 跑文档一致性检查**
  - `julia --project=. scripts/dev/check_docs_consistency.jl`

- [x] **Step 5.5: 跑 active 文档治理检查**
  - `julia --project=. scripts/dev/check_active_docs_governance.jl`

- [x] **Step 5.6: 回填任务单完成状态**
  - 将已完成 Step 打勾，并为失败项记录阻塞原因。

- [x] **Step 5.7: 写验收记录小节**
  - 记录测试命令、结果、关键风险与遗留项。

- [ ] **Step 5.8: 提交文档收尾**
  - `git add docs/dev/active/2026-04-07_mode-schema_unified-residual-kernel_followup-plan.md`
  - `git commit -m "docs(dev): finalize modeschema unified residual task sheet with acceptance status"`

## 完成定义（DoD）

- [x] solver 与 implicit 都通过同一 base 残差入口构建核心残差。
- [x] `FixedMu` / `FixedRho` 至少一条主路径完成 schema/accessor 化。
- [x] 新增一致性测试能够稳定复现并防回归。
- [x] 现有关键测试不退化。
- [x] 本任务单状态与证据已更新完整。

## 风险与回滚

- 风险 1：数值路径收口导致微小浮点偏差。
  - 缓解：使用 `rtol/atol` 回归阈值，保留 parity smoke。
- 风险 2：AD 类型在新入口出现意外强转。
  - 缓解：增加 Dual smoke，禁止在 base 残差中 `Float64(...)`。
- 回滚策略：按任务提交粒度回滚单个 commit，不做跨任务混合回退。

## 验收记录（2026-04-07）

### 测试与治理命令

- `julia --project=. -e 'include("tests/unit/models/test_gap_core_residual.jl")'`：PASS
- `julia --project=. -e 'include("tests/unit/models/test_mode_accessor_schema.jl")'`：PASS
- `julia --project=. -e 'include("tests/integration/models/test_solver_implicit_residual_parity_smoke.jl")'`：PASS
- `julia --project=. -e 'include("tests/integration/models/test_ad_implicit_contract_smoke.jl")'`：PASS
- `julia --project=. -e 'include("tests/integration/models/test_models_implicitdiff_njl2_smoke.jl")'`：PASS
- `julia --project=. -e 'include("tests/unit/models/test_implicit_gap.jl")'`：PASS
- `julia --project=. -e 'include("tests/integration/models/test_solver_auto_backend_semantic_parity.jl")'`：PASS
- `julia --project=. scripts/dev/check_docs_consistency.jl`：PASS
- `julia --project=. scripts/dev/check_active_docs_governance.jl`：PASS

### 关键改动

- 新增统一 core 残差入口：`src/models/solver/Conditions.jl` 中 `gap_core_residual!`，并由各 mode `build_residual!` 统一委托 core block。
- implicit 条件路径对齐 core：`src/models/solver/ImplicitGapLegacy.jl`（fixed-mu / flavor-mu adapters）。
- schema accessor 收口：`src/models/solver/SchemaAdapter.jl` 新增 `state_view/mu_view`，并在 `src/models/solver/Conditions.jl` 接入切片。
- 新增回归/契约测试：
  - `tests/unit/models/test_gap_core_residual.jl`
  - `tests/unit/models/test_mode_accessor_schema.jl`
  - `tests/integration/models/test_solver_implicit_residual_parity_smoke.jl`

### 未完成项与原因

- Step 2.4 / 3.4 / 4.5（“先失败再修复”）未完整保留失败证据：相关路径在现有主线已接近目标形态，新增断言初次运行未稳定复现预期失败。
- NJL 2/3 维 fallback 收口已在“追加任务单（NJL fallback 收口，2026-04-07）”完成：implicit 条件路径改为通过 shared core 统一入口处理 2/3/5 维。

---

## 追加任务单（NJL fallback 收口，2026-04-07）

> 目标：按本文档开头 Goal/Architecture 的原口径，补齐 NJL2/NJL3 在 implicit 路径中的 shared core 收口，移除“仅 5 维走 shared core”的实现门槛。

### 范围与原则

- 范围：`src/models/solver/Conditions.jl`、`src/models/solver/ImplicitAdapters.jl`、`src/models/solver/ImplicitGapLegacy.jl` 与对应 unit/integration 测试。
- 原则 1（不改公式）：2/3 维数值表达保持与现有 `gap_residual` 一致，仅调整统一入口与装配路径。
- 原则 2（单内核）：`gap_core_residual!` 升级为 2/3/5 维统一残差入口。
- 原则 3（映射收口）：`NamedTuple <-> Vector` 与 `x/mu` 切片继续仅在 schema/accessor 层。
- 原则 4（Dual 友好）：shared core 路径禁止 `Float64(...)` 强转。

### Chunk A: Conditions 单内核扩展到 2/3/5 维

### Task A1: 扩展 `gap_core_residual!` 维度分发

**Files:**
- Modify: `src/models/solver/Conditions.jl`
- Test: `tests/unit/models/test_gap_core_residual_njl_parity.jl`（新建）

- [x] **Step A1.1: 新建 NJL parity 单测骨架**
  - 新建 `@testset "gap_core_residual njl parity"`。

- [x] **Step A1.2: 写失败断言（NJL2, dim=2）**
  - 比较 `gap_core_residual!` 与 `gap_residual` 分量一致性。

- [x] **Step A1.3: 写失败断言（NJL, dim=3）**
  - 比较 `gap_core_residual!` 与 `gap_residual` 分量一致性。

- [ ] **Step A1.4: 运行并确认失败**
  - Run: `julia --project=. -e 'include("tests/unit/models/test_gap_core_residual_njl_parity.jl")'`

- [x] **Step A1.5: 在 `Conditions.jl` 增加 2/3/5 维统一入口**
  - 保留现有 5 维 fast-path；2/3 维通过 `gap_residual(_get_model(params.model_kind), ...)` 计算并写回 `F`。

- [x] **Step A1.6: 统一长度校验与错误信息**
  - `length(F) == length(x_state)`；错误信息包含期望维度与实际维度。

- [x] **Step A1.7: 重跑 NJL parity 单测并通过**
  - Run: `julia --project=. -e 'include("tests/unit/models/test_gap_core_residual_njl_parity.jl")'`

### Task A2: 扩展 `_gap_conditions_dynamic` 支持 2/3/5

**Files:**
- Modify: `src/models/solver/Conditions.jl`
- Test: `tests/unit/models/test_gap_core_residual_njl_parity.jl`

- [x] **Step A2.1: 将 dynamic 维度守卫从“5/3固定”改为“state in (2,3,5) 且 mu_dim=3”**

- [x] **Step A2.2: `FixedRho/FixedEntropy/FixedSigma/FixedAsymmetricRho` 改为统一调用 `gap_core_residual!` 结果**

- [x] **Step A2.3: 回归已有关联单测**
  - Run: `julia --project=. -e 'include("tests/unit/models/test_gap_core_residual.jl")'`
  - Run: `julia --project=. -e 'include("tests/unit/models/test_mode_accessor_schema.jl")'`

### Chunk B: ImplicitAdapters 去除 NJL fallback 门槛

### Task B1: 去掉 `state_n == 5` 限制并统一 shared core

**Files:**
- Modify: `src/models/solver/ImplicitAdapters.jl`
- Test: `tests/integration/models/test_solver_implicit_residual_parity_njl_smoke.jl`（新建）

- [x] **Step B1.1: 新建 NJL implicit parity smoke 测试骨架**

- [x] **Step B1.2: 写失败断言（NJL2）**
  - `build_njl_problem(...).conditions` 与 solver-side `gap_core_residual!` 一致。

- [x] **Step B1.3: 写失败断言（NJL）**
  - 同上，3 维版本。

- [ ] **Step B1.4: 运行并确认失败**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_solver_implicit_residual_parity_njl_smoke.jl")'`

- [x] **Step B1.5: 修改 `_implicit_conditions_with_shared_core`**
  - 删除 `state_n == 5` 门槛；改为 `length(x) == gap_state_dim(model)` 且 `length(θ) in (2, 4)` 统一走 shared core。

- [x] **Step B1.6: 仅保留真正异常输入 fallback/报错路径**
  - 例如 `θ` 维度非法或 `x` 维度不匹配时给出显式 `ArgumentError`。

- [x] **Step B1.7: 重跑 NJL implicit parity smoke 并通过**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_solver_implicit_residual_parity_njl_smoke.jl")'`

### Chunk C: 契约回归与文档收口

### Task C1: AD/语义回归

**Files:**
- Test: `tests/integration/models/test_models_implicitdiff_njl2_smoke.jl`
- Test: `tests/integration/models/test_ad_implicit_contract_smoke.jl`
- Test: `tests/integration/models/test_solver_implicit_residual_parity_smoke.jl`

- [x] **Step C1.1: 复跑 NJL2 implicit differentiation smoke**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_models_implicitdiff_njl2_smoke.jl")'`

- [x] **Step C1.2: 复跑 AD implicit contract smoke**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_ad_implicit_contract_smoke.jl")'`

- [x] **Step C1.3: 复跑 PNJL parity smoke（防回归）**
  - Run: `julia --project=. -e 'include("tests/integration/models/test_solver_implicit_residual_parity_smoke.jl")'`

### Task C2: 注释与任务单更新

**Files:**
- Modify: `src/models/solver/ImplicitGapLegacy.jl`
- Modify: `docs/dev/active/2026-04-07_mode-schema_unified-residual-kernel_followup-plan.md`

- [x] **Step C2.1: 更新 `ImplicitGapLegacy.jl` 相关注释**
  - 明确 NJL 路径已接入 shared core，不再描述为 fallback-only。

- [x] **Step C2.2: 回填本追加任务单执行状态与证据**
  - 记录新增测试命令与 PASS/FAIL。

- [x] **Step C2.3: 治理检查**
  - Run: `julia --project=. scripts/dev/check_docs_consistency.jl`
  - Run: `julia --project=. scripts/dev/check_active_docs_governance.jl`

### 追加任务执行证据（2026-04-07）

- `julia --project=. -e 'include("tests/unit/models/test_gap_core_residual_njl_parity.jl")'`：PASS
- `julia --project=. -e 'include("tests/integration/models/test_solver_implicit_residual_parity_njl_smoke.jl")'`：PASS
- `julia --project=. -e 'include("tests/integration/models/test_models_implicitdiff_njl2_smoke.jl")'`：PASS
- `julia --project=. -e 'include("tests/integration/models/test_ad_implicit_contract_smoke.jl")'`：PASS
- `julia --project=. -e 'include("tests/integration/models/test_solver_implicit_residual_parity_smoke.jl")'`：PASS
- `julia --project=. -e 'include("tests/unit/models/test_gap_core_residual.jl")'`：PASS
- `julia --project=. -e 'include("tests/unit/models/test_mode_accessor_schema.jl")'`：PASS
- `julia --project=. -e 'include("tests/unit/models/test_implicit_gap.jl")'`：PASS
- `julia --project=. -e 'include("tests/integration/models/test_solver_auto_backend_semantic_parity.jl")'`：PASS
- `julia --project=. scripts/dev/check_docs_consistency.jl`：PASS
- `julia --project=. scripts/dev/check_active_docs_governance.jl`：PASS
