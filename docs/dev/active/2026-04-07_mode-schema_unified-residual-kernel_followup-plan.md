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
- NJL 2/3 维路径当前仍通过 fallback `gap_residual` 兼容：shared core 优先覆盖 5 维 PNJL（θ 维度 2/4）路径，避免错误地把低维状态强行映射到 5 维。
