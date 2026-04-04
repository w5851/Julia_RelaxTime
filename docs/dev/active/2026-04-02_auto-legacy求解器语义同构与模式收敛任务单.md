# auto/legacy 求解器语义同构与模式收敛任务单

## 1. 背景与目标

本任务单聚焦你指出的核心问题：`legacy` 与 `auto`（models 路径）在同一物理语义下出现行为差异。目标不是继续打补丁维持双轨，而是把“约束问题定义”收敛为单一可组合框架。

本期目标：

- [ ] 统一 `FixedRho` 在 `legacy` 与 `models` 下的数值问题定义（语义同构，不要求实现细节完全同构）。
- [ ] 建立“基础驻点约束 + 约束组件拼接”的通用残差组装机制，去除维度写死假设。
- [ ] 给出求解器模式收敛结论：收敛为“两类语义模式”还是“一个求解器下两种大模式”。
- [ ] 明确 continuity 仅作为策略层（seed/branch tracking），不再作为语义模式。
- [ ] 最终目标：完全移除旧求解实现与旧路由特判（不保留长期 legacy 双轨）。

---

## 2. 范围与非范围

### 2.1 本期范围

- [ ] `src/models/solver` 与 `src/models/constraint_solver.jl` 的约束组装与求解主链收敛。
- [ ] `FixedMu/FixedRho/FixedAsymmetricRho/FixedEntropy/FixedSigma` 的约束组件抽象统一。
- [ ] 扫描入口 `TmuScan/TrhoScan/DualBranchScan` 的 backend 选择口径与行为对齐。
- [ ] 覆盖 unit/integration/regression 的最小验证矩阵与 benchmark 冒烟证据。

### 2.2 非范围

- [ ] 不新增新的模型物理项（仅重构求解/约束工程结构）。
- [ ] 不在本期引入新的扫描类型或数据格式。
- [ ] 不扩大到与该问题无关的大规模目录重排。

---

## 3. 现状诊断（冻结）

- [x] `legacy` FixedRho：当前为联立残差求解主链（等价 8 维语义）。
- [x] `models` FixedRho：当前为外层 `μ` 一维根 + 内层 gap 解的分层链路。
- [x] 两条链路虽物理目标一致，但数值问题定义不同，导致收敛域与分支选择差异。
- [x] `rho_norm` 归一化曾存在实现偏差（已修复），属于差异放大因子而非唯一根因。
- [x] `auto` 对 PNJL 的路由已临时回退到 `legacy`，属于稳定性兜底，不是终态方案。

---

## 4. 设计约束（必须满足）

- [ ] 禁止写死“5/3/8”维度切片作为通用逻辑前提。
- [ ] 约束系统必须由组件拼接：
  - [ ] 基础驻点约束（`∂Ω/∂x_state = 0`）
  - [ ] 化学势关系约束（如 `μ_u=μ_d=μ_s`、或不对称关系）
  - [ ] 宏观目标约束（如 `n_B`、`s`、`σ`）
- [ ] 组件数量与维度由 mode/schema 决定，不由文件内常量硬编码。
- [ ] continuity 归类为“求解策略插件”，不得承担语义定义职责。

---

## 5. 任务分解（可勾选）

### A0：约束组件抽象层（先测试后实现）

- [x] 新增失败单测：约束组件可声明输出维度、依赖变量域（state/mu/thermo）。
- [x] 定义约束组件契约（建议：`constraint_dim`, `eval_constraint!`, `constraint_name`）。
- [x] 将现有 `Fixed*` mode 映射为组件列表，而非写死残差函数。
- [x] 验证组件拼接后 residual 总维度一致性（静态检查 + 运行时检查）。

### A1：ProblemSpec 主链收敛

- [x] 新增失败测试：`build_problem_spec(mode)` 生成的 residual 能覆盖 `FixedRho/FixedAsymmetricRho`。
- [x] 将 `ProblemSpec` 从占位结构升级为执行主链（residual/forward_solve/postprocess 同步接线）。
- [x] 让 `solve_constraint` 统一走 `ProblemSpec`，仅保留必要适配层。
- [x] 对旧特化路径加迁移标记，禁止继续扩散新逻辑。
  - [x] 增量：`solve_constraint` 已支持 `FixedEntropy` / `FixedSigma` / `FixedAsymmetricRho` 的 `problem_spec` override（保留兼容回退）。

### A2：FixedRho 同构化改造

- [x] 新增失败回归：同点同 seed 下 `legacy` 与 `models` 的 residual 口径一致（容差内）。
- [x] 将 `models FixedRho` 改为组件拼接后的联立约束主链（而非外1维内5维特化主链）。
- [x] 保留可控 fallback，但 fallback 仅做数值稳健兜底，不改变语义定义。
- [x] 补充 branch 选择与候选治理规则（pressure/physicality/residual 的统一优先级）。

### A3：扫描入口与 backend 语义收敛

- [x] 新增失败 integration：`solver_backend=:auto` 与显式 backend 在语义一致时结果一致。
- [x] 统一 `auto` 解释为“按 model capability 选实现，不改语义定义”。
- [x] 去除 PNJL 上仅因历史差异导致的永久路由特判（以阶段开关受控下线）。

### A4：验证、文档与治理

- [x] unit：覆盖组件拼接、维度检查、FixedRho/AsymRho 等价性。
- [x] integration：覆盖 T-μ/T-ρ 扫描行为与 auto/backend 一致性。
- [x] regression：覆盖关键固定点与首阶附近点的分支稳定性。
- [x] benchmark：`scripts/perf/pnjl/scan_perf.jl` 收敛率基线不退化。
- [x] 同步 `docs/api/models/solver/*` 与相关迁移文档口径。

---

## 6. 架构决策讨论：两类模式 vs 一个求解器两种大模式

### 6.1 备选方案

- [ ] 方案 S1：两个独立求解器（ground-state solver / manifold solver）。
- [x] 方案 S2：一个统一求解器内的两种语义模式（推荐）。

### 6.2 推荐结论（当前）

- [x] 推荐 S2：单一求解器框架 + 两种语义模式。
- [x] 语义模式 M1：`GroundState`（物理约束下 Ω 极小优选）。
- [x] 语义模式 M2：`ConstrainedManifold`（物理约束下返回稳定/亚稳可行分支集）。
- [x] continuity 定位：策略层（seed continuity / branch tracking / warm start），不是 M1/M2 之一。

### 6.3 选择 S2 的原因

- [x] 避免两套求解器重复维护 residual/thermo/postprocess 逻辑。
- [x] 可共享同一约束组件拼接内核，降低 `legacy` 与 `auto` 漂移概率。
- [x] 在首阶附近可由同一内核输出候选集，再按模式进行 selector 分流。
- [x] 更利于将来扩展非三味化学势关系，不会绑死在“8 维 PNJL 特例”。

---

## 7. 验收标准

- [ ] `FixedRho` 在 `legacy/models/auto` 下语义一致，差异仅限数值容差内。
- [ ] 约束残差构造不依赖固定维度常量，组件可组合并通过维度一致性检查。
- [ ] `auto` 不再作为“行为差异入口”，仅作为实现路由入口。
- [ ] M1/M2 模式边界清晰，continuity 被验证为策略层能力。
- [ ] unit/integration/regression/governance 与 benchmark 冒烟通过。
- [ ] 旧求解实现完成下线：无 `ImplicitSolver` 直接业务依赖、无 PNJL 特判回退主路径。

### 7.1 PR-0 准入硬门槛（R1 前置冻结）

- [x] FixedRho 精度护栏（关键点：`T=90/110/130 MeV`, `rho*=0.2/0.6/1.0`）：
  - [x] `residual_norm <= 1e-6`（目标 `<= 1e-8`）
  - [x] `abs(rho_norm - rho_target) <= 1e-6`
- [x] 多初值候选治理固定排序（可重复）：
  - [x] `hard_constraint_ok` 优先
  - [x] `residual_norm` 最小优先
  - [x] `pressure` 最大优先
  - [x] `seed_index` 最小优先（缺省回退到候选枚举顺序）
- [x] `auto/models/legacy` 同语义对照矩阵通过：
  - [x] 收敛一致性
  - [x] 关键热力学量容差内一致
  - [x] 分支标签/选择口径一致（当前链路无显式分支标签，按同候选治理口径核对）

### 7.2 PR-0 停止线

- [x] 任一关键点未满足 FixedRho 精度护栏。（未触发）
- [x] 任一对照点 `auto/models/legacy` 语义不一致。（未触发）
- [x] regression smoke 相比上一基线出现退化。（未触发）

触发任一停止线：禁止进入 R1，先修复并补齐证据。

---

## 8. 验证命令（草案）

- [ ] `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_constraint_components.jl,models/test_problem_spec_modes.jl"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_constraint_solver.jl"; include("tests/unit/runtests.jl")'`
- [ ] `julia --project=. -e 'include("tests/integration/models/test_dimension_agnostic_scan_smoke.jl")'`
- [ ] `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- [ ] `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- [ ] `julia --project=. scripts/perf/pnjl/scan_perf.jl`
- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/check_active_docs_governance.jl`

---

## 9. DoD

- [ ] 任务项与验收项全部勾选。
- [ ] 核心差异问题（legacy vs auto）有可复现实验前后对照证据。
- [ ] 架构模式（M1/M2）在代码与文档中均有单一口径。
- [ ] 无新增长期双轨与隐式 fallback 语义漂移点。
- [ ] 旧实现物理移除完成（含代码、导出、文档与测试口径收敛）。

---

## 10. 可执行实现蓝图（文件级）

### 10.1 目标架构（单求解器 + 双语义模式）

- [ ] 新增 `SolveSemanticMode`：`GroundState`、`ConstrainedManifold`。
- [ ] 约束改为组件拼接：`Stationarity` + `MuRelations` + `MacroTargets`。
- [ ] `ConstraintMode` 只表达业务约束目标；不再直接绑定固定维度残差实现。
- [ ] continuity/multiseed/branch tracking 统一归入 strategy 插件层。

### 10.2 文件改动清单（建议顺序）

#### B0：约束组件层（新增）

- [ ] 新增 `src/models/solver/ConstraintComponents.jl`
  - [ ] `AbstractConstraintComponent`
  - [ ] `constraint_dim(component, schema)`
  - [ ] `eval_constraint!(F, offset, component, state_ctx)`
  - [ ] 内置组件：
    - [ ] `StationarityComponent`
    - [ ] `EqualMuComponent` / `LinearMuRelationComponent`
    - [ ] `FixedBaryonDensityComponent`
    - [ ] `FixedEntropyComponent`
    - [ ] `FixedSigmaComponent`
    - [ ] `AsymmetricDensityComponent`

#### B1：组装与契约层（扩展）

- [ ] 扩展 `src/models/solver/ProblemSpec.jl`
  - [ ] `ConstraintAssembly`
  - [ ] `build_constraint_assembly(mode, schema; kwargs...)`
  - [ ] `build_residual_from_assembly!(assembly, ...)`
  - [ ] `solver_spec(mode; semantic, strategy, selector, backend)`

- [ ] 扩展 `src/models/solver/ConstraintModes.jl`
  - [ ] mode -> 默认 `MuRelations` 映射函数（不写死三味相等）
  - [ ] 增加 mode 级参数校验入口（防止无效组合静默通过）

#### B2：执行主链（收敛）

- [ ] 重构 `src/models/solver/Solver.jl`
  - [ ] `solve_constraint` 统一走 `ProblemSpec/SolverSpec` 主链
  - [ ] 对外保留兼容签名，但内部不再分流到旧特化实现

- [ ] 重构 `src/models/constraint_solver.jl`
  - [ ] 逐步移除 `_solve_constraint_fixedrho` 外1维/内5维主路径
  - [ ] 改为“联立 residual + 统一候选治理/fallback”
  - [ ] 将 fallback 下沉为 backend 细节，不改变语义定义

#### B3：扫描入口与路由（收敛）

- [ ] 修改 `src/models/scans/TmuScan.jl`
  - [ ] `solver_backend=:auto` 仅选 backend，不改方程语义
  - [x] 新增 `semantic_mode` 透传（默认 `:ground_state`）

- [ ] 修改 `src/models/scans/TrhoScan.jl`
  - [ ] 删除 PNJL 的历史特判路由（阶段开关保护后下线）
  - [x] 支持 `semantic_mode=:constrained_manifold` 返回候选集/标签

- [x] 修改 `src/models/scans/DualBranchScan.jl`
  - [x] 用 `semantic_mode + selector` 显式表达双分支扫描行为

#### B4：旧实现下线（最终目标）

- [ ] 清理 `src/models/solver/ImplicitSolver.jl` 的业务入口角色（先冻结、后删除）。
  - [x] 阶段收敛：扫描层与脚本层不再直接依赖 `ImplicitSolver` 私有符号；统一走 `Models` 公共转发。
  - [ ] 终态收口：删除 `ImplicitSolver` 业务入口能力，仅保留必要底层数值内核或彻底移除。
- [ ] 移除 `Models.jl` 中对旧路径的导出/接线依赖。
  - [x] 阶段收敛：已移除对 `ImplicitSolver` 的顶层模块导入。
  - [ ] 终态收口：清理 legacy 相关公共导出与接线残留。
- [ ] 移除临时兼容开关与 legacy-only fallback 分支。
  - [x] 阶段收敛：`use_problem_spec=false` 已显式门禁（`allow_legacy_path=true`）、参数互斥和扫描入口透传拦截。
  - [ ] 终态收口：删除 `use_problem_spec` / `allow_legacy_path` / `warn_on_legacy_path` 兼容参数。
- [ ] 文档统一更新：不再暴露 legacy solver 作为推荐或默认。
  - [x] 阶段收敛：`ConstraintModes.md` 已明确 ProblemSpec 主链默认与 legacy 仅兼容路径。
  - [ ] 终态收口：同步 `docs/api/models/*` 与 `docs/api/data_contracts.md` 全量去 legacy 推荐口径。

##### B4-7 下线窗口与准入标准（收口计划）

- [ ] W1（当前→下一次 PR 合并前）：冻结新增 legacy 调用面
  - [x] `src/` 新增入口已禁止透传 legacy 开关（扫描/solver 主链收紧完成）。
  - [ ] CI 侧新增“legacy 参数泄漏检查”脚本或测试门禁（防回归）。
- [ ] W2（下一里程碑）：移除用户侧兼容开关
  - [ ] 删除 `solve_constraint` 中 `use_problem_spec=false` 分支与相关参数。
  - [ ] 保留迁移说明（breaking change）并给出等价调用方式。
- [ ] W3（R4）：物理删除 legacy 业务路径
  - [ ] 清理 `ImplicitSolver` 业务入口/导出与扫描层 legacy backend 特判。
  - [ ] 回归矩阵通过（unit/integration/regression + 文档治理）后归档任务单。

### 10.3 测试新增与迁移清单

#### Unit（新增）

- [ ] `tests/unit/models/test_constraint_components.jl`
  - [ ] 组件维度声明与拼接顺序稳定
  - [ ] 组件依赖域校验（state/mu/thermo）

- [ ] `tests/unit/models/test_problem_spec_semantic_modes.jl`
  - [ ] `GroundState` 与 `ConstrainedManifold` 返回契约
  - [ ] selector 排序优先级可重复

- [ ] `tests/unit/models/test_mu_relations_component.jl`
  - [ ] 相等关系/线性关系/不对称关系约束正确

#### Integration（新增）

- [ ] `tests/integration/models/test_solver_auto_backend_semantic_parity.jl`
  - [ ] 同语义下 `auto` 与显式 backend 结果一致（容差内）

- [ ] `tests/integration/pnjl/test_trho_scan_semantic_modes_smoke.jl`
  - [ ] 同一扫描点验证 `GroundState` 单解与 `Manifold` 多解关系

#### Regression（新增）

- [ ] `tests/regression/models/test_fixedrho_semantic_equivalence_regression.jl`
  - [ ] 关键固定点 legacy-era 基线对齐（迁移期）

- [ ] `tests/regression/models/test_firstorder_manifold_branch_stability.jl`
  - [ ] 首阶附近分支可重复识别与排序稳定

### 10.4 阶段里程碑（含移除旧实现）

- [ ] R1：组件层与 ProblemSpec 主链可运行，旧路径仍保留只读镜像。
- [ ] R2：`FixedRho` 主路径切换到新链，`auto` 与显式 backend 语义一致。
- [ ] R3：扫描层语义参数化完成，PNJL 历史特判移除。
- [ ] R4：旧实现彻底删除（代码/导出/文档/测试全部收敛）。

### 10.5 回滚与安全阈值

- [ ] 回滚策略：仅允许回滚到“新链+旧链并存但新链默认”的前一里程碑，不回到永久 legacy 默认。
- [ ] 停止线：
  - [ ] `scan_perf.jl` 收敛率显著退化
  - [ ] regression 首阶分支不稳定
  - [ ] docs/governance 检查失败
- [ ] 任一停止线触发时，冻结下一里程碑推进，先修复当前阶段。

---

## 11. 执行记录（append-only）

- [x] 2026-04-02：基于当前 PR #47 讨论，确认差异根因是“同语义的两套数值问题定义”未收敛，而非物理目标主动分叉。
- [x] 2026-04-02：冻结本任务单，明确下一阶段以“组件拼接约束 + 单一求解器双模式（M1/M2）”为主线推进。
- [x] 2026-04-02：完成 PR-0 硬门槛固化与护栏测试接线。
  - 新增 `tests/regression/models/test_fixedrho_precision_guard_regression.jl`，覆盖 `T=90/110/130 MeV` + `rho*=0.2/0.6/1.0` 关键点。
  - 新增 `tests/unit/models/test_multiseed_candidate_governance.jl`，冻结 `hard_constraint > residual > pressure > seed_index` 规则。
  - 新增 `tests/integration/models/test_solver_backend_semantic_parity_guard.jl`，对照 `auto/models/legacy` 同语义一致性。
  - `tests/regression/runtests.jl` 与 `tests/integration/runtests.jl` 已接线新护栏测试。
- [x] 2026-04-02：PR-0 验证矩阵执行通过（R1 准入通过）。
  - unit smoke：`781/781` 通过。
  - integration smoke：`394/394` 通过。
  - regression smoke：`512 pass, 1 broken(optional fixture)`，与既有可选跳过口径一致。
  - 文档治理：`check_docs_consistency.jl` 与 `check_active_docs_governance.jl` 均通过。
- [x] 2026-04-02：创建 R1 开发分支 `feat/models-solver-semantic-convergence-r1`，按 stacked PR 方式基于 PR #48 开启下一轮功能开发（PR #49）。
- [x] 2026-04-02：R1 A0/A1 首批完成（PR #49）。
  - 新增 `src/models/solver/ConstraintComponents.jl`，建立 `AbstractConstraintComponent` 契约及默认 mode->components 映射。
  - `build_problem_spec` 接线 `conditions` 与 `forward_solve`，并通过 `tests/unit/models/test_constraint_components.jl` / `test_problem_spec_contract.jl`。
- [x] 2026-04-02：R1 A2 深化完成（PR #49）。
  - `FixedRho` 的 `ProblemSpec.forward_solve` 升级为显式候选池 + 统一 hard constraints + selector 治理。
  - `solve_constraint(FixedRho)` 支持 `problem_spec` 优先链路；不传时保持兼容旧链路。
  - 新增对照：`tests/integration/models/test_problem_spec_fixedrho_forwardsolve_smoke.jl`、`tests/regression/models/test_problem_spec_fixedrho_parity_regression.jl`。
- [x] 2026-04-02：R1 A3 完成（PR #49）。
  - `TmuScan/TrhoScan` 中 PNJL `solver_backend=:auto` 默认路由改为 `:models`；新增受控开关 `auto_pnjl_backend` 支持回退 `:legacy`。
  - 新增 `tests/integration/models/test_auto_backend_models_routing_smoke.jl`，并复核 `test_solver_backend_semantic_parity_guard.jl`。
- [x] 2026-04-02：完成 A4 验证矩阵（除 API 文档同步项）。
  - unit smoke：`781/781` 通过。
  - integration smoke：`402/402` 通过。
  - regression smoke：`530 pass, 1 broken(optional fixture)`，与既有可选跳过口径一致。
  - benchmark smoke：`scripts/perf/pnjl/scan_perf.jl` 收敛率 100%，无退化证据。
  - 治理检查：`check_docs_consistency.jl` 与 `check_active_docs_governance.jl` 通过。
- [x] 2026-04-02：完成 `docs/api/models/solver/*` 文档口径同步。
  - 更新 `README.md`、`CoreConcepts.md`、`ConstraintModes.md`，补充 `ProblemSpec` 与 `ConstraintComponents` 的 R1 主链语义。
  - 通过脚本重生成 `docs/api/models/solver/generated/Exports.md`（`generate_api_export_index.jl`）。
- [x] 2026-04-02：继续推进 B1 增量收敛（PR #50）。
  - `solve_constraint` 对 `FixedEntropy` / `FixedSigma` / `FixedAsymmetricRho` 已支持 `problem_spec` override（保持旧路径兼容回退）。
  - 新增 integration 冒烟 `tests/integration/models/test_problem_spec_other_modes_forwardsolve_smoke.jl` 并接入 integration smoke。
- [x] 2026-04-02：继续推进 B3 路由收敛（PR #50）。
  - `DualBranchScan` 已引入 `auto_pnjl_backend` 参数，`solver_backend=:auto` 下与 `TmuScan/TrhoScan` 对齐：`thermo_backend=:models` 固定走 `:models`，否则由 `auto_pnjl_backend` 显式决定（默认 `:models`）。
  - 新增 `tests/unit/models/test_dual_branch_scan.jl` 路由规则覆盖，冻结 `auto` 路由行为。
- [x] 2026-04-02：继续推进 B2/B3 收敛（PR #50）。
  - `solve_constraint` 对 `problem_spec` 新增类型校验，`FixedEntropy/FixedSigma/FixedAsymmetricRho` 统一复用 `problem_spec` override 主链适配逻辑。
  - `ProblemSpec` 非 `FixedRho` governed forward_solve 增加 `semantic_mode` selector 路由（`:ground_state` -> 压强优先，`:constrained_manifold` -> 残差优先），并支持自定义 `selector`。
  - `TrhoScan` 新增 `semantic_mode`/`selector` 透传到 models `ProblemSpec` 主链；补充 `tests/unit/models/test_trho_scan.jl` 与 `tests/integration/pnjl/test_trho_scan_solver_backend_models_smoke.jl` 覆盖参数校验与冒烟链路。
- [x] 2026-04-03：继续推进 B3 扫描语义口径对齐（PR #50）。
  - `TmuScan` 新增 `semantic_mode` 与 `selector` 参数并接入显式校验：当前 `FixedMu` 语义仅允许 `semantic_mode=:ground_state` 且 `selector=nothing`，避免 `auto` 路径产生隐式语义漂移。
  - 新增 `tests/unit/models/test_tmu_scan.jl` 参数校验覆盖，并在 `tests/integration/pnjl/test_tmu_scan_smoke.jl` 增加非法 `semantic_mode` 冒烟断言。
  - integration smoke 全量复跑通过：`427/427`。
- [x] 2026-04-03：继续推进 B3 配置入口收敛（PR #50）。
  - `ScanConfig.TmuScanConfig` / `ScanConfig.TrhoScanConfig` 新增 `auto_pnjl_backend`、`semantic_mode`、`selector` 字段并接入 `scan_kwargs` 透传，避免 kwargs 与配置对象语义能力不一致。
  - 在 `tests/integration/pnjl/test_scan_config_equivalence_smoke.jl` 增加语义参数透传冒烟，验证 config 入口可驱动 `TmuScan`/`TrhoScan` 新参数链路。
  - integration smoke 复跑通过：`427/427`。
- [x] 2026-04-03：继续推进 B3→B4 依赖收敛（PR #50）。
  - `TmuScan/TrhoScan/DualBranchScan/ScanCommon/ScanResultFinalize` 不再直接依赖 `ImplicitSolver` 子模块类型或入口，统一改为从 `Main.Models` 公共接口导入 `solve`/`SolverResult`（`TrhoScan` 同步使用 `solve_weighted_block_fallback` 公共转发）。
  - `Solver.jl` 新增 `solve_weighted_block_fallback` 公共转发并在 `Models.jl` 导出，减少扫描层对旧求解子模块的编排依赖。
  - integration smoke 复跑通过：`427/427`。
- [x] 2026-04-03：继续推进 B4 去旧化接口收敛（PR #50）。
  - 新增 `Models.is_physical_solution` 公共转发（由 `Solver.jl` 暴露），将 `scripts/relaxtime/run_gap_transport_scan.jl` 中对 `Main.Models.ImplicitSolver._default_is_physical_solution` 的直接依赖替换为公共 API。
  - 该改动进一步压缩业务脚本对 `ImplicitSolver` 私有符号的耦合，为后续 `ImplicitSolver` 入口冻结/下线清理提供前置条件。
  - integration smoke 复跑通过：`427/427`。
- [x] 2026-04-03：继续推进 B4 模块边界收敛（PR #50）。
  - `Models.jl` 移除对 `ImplicitSolver` 的顶层模块导入（保留 include 与 `Solver.jl` 内部转发），减少外层模块对旧求解子模块命名空间的直接暴露。
  - 补充回归：`tests/unit/pnjl/test_solver_implicit.jl` 全集与 integration smoke 通过，验证边界收敛未破坏既有求解能力。
- [x] 2026-04-03：继续推进 B4 私有符号去耦（PR #50）。
  - `ImplicitSolver` 将物理性判据函数标准化为导出符号 `default_is_physical_solution`，并保留 `_default_is_physical_solution` 兼容别名。
  - `Solver.jl` 的 `is_physical_solution` 公共转发改为依赖导出符号，避免继续绑定 `ImplicitSolver` 私有命名约定。
  - 回归验证通过：`tests/unit/pnjl/test_solver_implicit.jl`、`tests/unit/models/test_scan_result_finalize.jl`、integration smoke `427/427`。
- [x] 2026-04-03：继续推进 B4 内部默认口径收敛（PR #50）。
  - `ImplicitSolver` 内部调用与参数默认中的物理性判据统一切换为 `default_is_physical_solution`，`_default_is_physical_solution` 仅保留兼容别名，不再作为主路径默认符号。
  - 该变更确保“公共导出符号”为代码主语义来源，进一步降低后续私有符号清理的回归风险。
  - 回归验证通过：`tests/unit/pnjl/test_solver_implicit.jl`、`tests/unit/models/test_problem_spec_contract.jl`、integration smoke `427/427`。
- [x] 2026-04-03：继续推进 B4 私有兼容符号清理（PR #50）。
  - 移除 `ImplicitSolver` 中 `_default_is_physical_solution` 兼容别名，确认代码库内不再存在该私有符号引用。
  - 回归验证通过：`tests/unit/pnjl/test_solver_implicit.jl`、`tests/unit/models/test_problem_spec_contract.jl`、integration smoke `427/427`。
- [x] 2026-04-03：继续推进 B4 兼容回退治理（PR #50）。
  - `Solver.solve_constraint` 的 `use_problem_spec=false` 兼容回退新增显式告警（可通过 `warn_on_legacy_path=false` 静默），用于标注该路径为过渡能力并避免新逻辑继续依赖旧链。
  - 在 `tests/unit/models/test_problem_spec_contract.jl` 增加 legacy 回退告警/静默与参数校验覆盖，确保治理行为可测试。
  - integration smoke 复跑通过：`427/427`。
- [x] 2026-04-03：继续推进 B4 兼容开关收紧（PR #50）。
  - `Solver.solve_constraint` 对 `use_problem_spec=false` 增加显式门禁：仅当 `allow_legacy_path=true` 时允许进入 legacy 回退路径；否则抛 `ArgumentError`，避免新调用方继续绑定旧链。
  - 保留并复用 `warn_on_legacy_path` 告警治理能力，仅在显式放行 legacy 路径时触发告警/静默策略。
  - `tests/unit/models/test_problem_spec_contract.jl` 新增门禁行为与 `allow_legacy_path` 参数校验覆盖，确保“默认禁用、显式放行”语义可回归。
- [x] 2026-04-03：继续推进 B4 测试侧兼容最小化（PR #50）。
  - 完成 `use_problem_spec=false` 全仓扫描与分类：`src/` 生产入口无残余 legacy 调用，仅保留 `Solver.jl` 门禁实现与单测/文档治理记录。
  - `tests/unit/models/test_problem_spec_contract.jl` 将 legacy 兼容断言收敛为最小必要集合：门禁拒绝、显式放行并告警、显式放行并静默、参数类型校验；去除重复 legacy 调用路径。
  - 回归验证通过：`julia --project=. -e 'include("tests/unit/models/test_problem_spec_contract.jl")'`（`106/106`）。
- [x] 2026-04-03：继续推进 B4-3 兼容参数暴露面收缩（PR #50）。
  - `solve_constraint` 在 `use_problem_spec=true`（默认）时，显式拒绝 `allow_legacy_path` / `warn_on_legacy_path`，将 legacy 治理参数限定在兼容回退分支，避免主链调用保留无效开关。
  - `tests/unit/models/test_problem_spec_contract.jl` 新增“默认主链禁止 legacy 参数”断言覆盖，确保参数可见性收缩可回归。
  - 同文件回归验证通过：`julia --project=. -e 'include("tests/unit/models/test_problem_spec_contract.jl")'`（`108/108`）。
- [x] 2026-04-03：继续推进 B4-4 兼容参数互斥收敛（PR #50）。
  - `solve_constraint` 新增互斥守卫：`use_problem_spec=false` 时禁止同时传 `problem_spec`，避免 legacy 回退与 `ProblemSpec` 主链契约发生混合语义。
  - `tests/unit/models/test_problem_spec_contract.jl` 新增对应断言，固定“legacy 回退 vs 主链契约”参数互斥规则。
  - 同文件回归验证通过：`julia --project=. -e 'include("tests/unit/models/test_problem_spec_contract.jl")'`（`109/109`）。
- [x] 2026-04-03：继续推进 B4-5 扫描入口 legacy 参数下线准备（PR #50）。
  - `TrhoScan/TmuScan/DualBranchScan` 在 models 求解路径新增 `_reject_legacy_solver_kwargs` 守卫，显式拒绝 `use_problem_spec` / `allow_legacy_path` / `warn_on_legacy_path` / `problem_spec` 透传，避免扫描入口绕过主链治理边界。
  - 三个扫描模块单测新增“models 路径禁止 legacy solver 开关”覆盖，固定外部可见面的收缩行为。
  - 回归验证通过：`tests/unit/models/test_trho_scan.jl`（`21/21`）、`tests/unit/models/test_tmu_scan.jl`（`19/19`）、`tests/unit/models/test_dual_branch_scan.jl`（`34/34`）。
- [x] 2026-04-03：继续推进 B4-6 集成冒烟回归准入（PR #50）。
  - 执行 integration smoke 全量回归：`julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`。
  - 回归通过：`427/427`，确认 B4-3~B4-5 的 legacy 参数边界收紧未引入集成侧回归。
- [x] 2026-04-03：继续推进 B4-7 收口计划建模（PR #50）。
  - 将 B4 主目标拆分为“阶段收敛 vs 终态收口”两层任务，补齐下线窗口 W1/W2/W3 与对应准入标准，避免后续合并阶段出现范围漂移。
  - 明确当前已完成边界（主链默认 + legacy 门禁 + 参数互斥 + 扫描入口拦截）与剩余终态动作（删除兼容参数、清理 legacy 业务入口、文档全量去推荐口径）。
- [x] 2026-04-03：继续推进 W1 门禁落地（PR #50）。
  - 新增 `scripts/dev/check_legacy_solver_switch_leakage.jl`，对 `src/**/*.jl` 与 `tests/**/*.jl` 扫描 `use_problem_spec` / `allow_legacy_path` / `warn_on_legacy_path`，并将可出现位置收敛到“治理边界 allowlist”。
  - 新增 CI workflow `.github/workflows/legacy-switch-governance.yml`，在 PR/push 上执行该门禁，阻断 legacy 参数外溢回归。
  - 根据扫描层“拦截型守卫”实现现状，将 allowlist 扩展为 `Solver.jl` + `TrhoScan/TmuScan/DualBranchScan` 守卫实现及对应 unit tests，确保“仅治理边界可见”而非“全仓任意可见”。
  - 验证通过：`scripts/dev/check_legacy_solver_switch_leakage.jl`（OK）；`tests/unit/models/test_problem_spec_contract.jl`（`109/109`）、`test_trho_scan.jl`（`21/21`）、`test_tmu_scan.jl`（`19/19`）、`test_dual_branch_scan.jl`（`34/34`）。
- [x] 2026-04-03：继续推进 W2 兼容参数移除（PR #50）。
  - `solve_constraint` 删除 `use_problem_spec=false` legacy 回退分支与 `allow_legacy_path` / `warn_on_legacy_path` 兼容参数；`FixedRho/FixedEntropy/FixedSigma/FixedAsymmetricRho` 统一固定为 `ProblemSpec` 主链。
  - `tests/unit/models/test_problem_spec_contract.jl` 删除 legacy 回退行为断言，改为“兼容参数已移除并抛错”断言。
  - `docs/api/models/solver/ConstraintModes.md` 同步更新为“兼容参数已移除、problem_spec 仍可显式覆盖”的最终口径。
  - 验证通过：`scripts/dev/check_legacy_solver_switch_leakage.jl`（OK）；`tests/unit/models/test_problem_spec_contract.jl`（`103/103`）、`test_trho_scan.jl`（`21/21`）、`test_tmu_scan.jl`（`19/19`）、`test_dual_branch_scan.jl`（`34/34`）；integration smoke（`427/427`）。
- [x] 2026-04-03：进入 W3 准备（ImplicitSolver 业务入口物理下线前评估，PR #50）。
  - 影响面盘点：`src/models/scans/TrhoScan.jl`、`src/models/scans/TmuScan.jl`、`src/models/scans/DualBranchScan.jl` 仍在 `solver_backend=:legacy` 分支调用 `solve(...)`；`src/models/phase/PMPhaseSeeds.jl` 与 `src/models/phase/PMPhaseDiagnostic.jl` 仍保留 `solver_backend=:legacy` 路径。
  - 公共 API 依赖盘点：`Models.jl` 仍导出 `create_implicit_solver` / `solve_with_derivatives` / `solve_with_root_diagnostics` / `solve_weighted_block_fallback`，且对应 unit/integration/regression 覆盖仍活跃，说明 W3 需按“先冻结业务入口、后拆分底层能力”推进。
  - W3 建议拆分：
    - W3-A：冻结 legacy backend 可达面（扫描/phase 入口禁用 `:legacy`，保留受控过渡错误提示）；
    - W3-B：将 `ImplicitSolver` 中仅用于业务编排的入口从 `Models` 导出面剥离，保留必要数值内核；
    - W3-C：完成导出/API 文档与回归矩阵收口后，再执行物理删除。
- [x] 2026-04-03：执行 W3-A（冻结 legacy backend 可达面，PR #50）。
  - `TmuScan/TrhoScan/DualBranchScan` 运行路径移除 legacy backend 执行分支：调用层仅接受 models 执行（保留 `:auto` 路由语义，但运行时若解析到 legacy 将抛显式 `ArgumentError`）。
  - `PMPhaseSeeds` 与 `PMPhaseDiagnostic` 移除 legacy backend 执行分支，改为对 `solver_backend=:legacy` 明确报错，冻结 phase 链路对旧业务入口的可达性。
  - 单测与冒烟验证通过：`test_tmu_scan.jl`（`20/20`）、`test_trho_scan.jl`（`22/22`）、`test_dual_branch_scan.jl`（`34/34`）、`test_pm_phase_diagnostic.jl`（全 testset 通过）、integration `test_tmu_scan_smoke.jl`、`test_trho_scan_solver_backend_models_smoke.jl`、`test_scan_solver_boundary_error_smoke.jl`。
- [x] 2026-04-03：执行 W3-B（导出面剥离首批，PR #50）。
  - `Models` 公共导出面移除 `create_implicit_solver`、`solve_with_root_diagnostics`、`solve_weighted_block_fallback`，保留 `solve_with_derivatives` 作为稳定导数入口；对应内部能力仍保留在 solver 子模块供受控内部路径使用。
  - `tests/unit/models/test_solver_dimension_agnostic.jl` 增加“上述 3 个入口不再出现在 Models 顶层导出”断言，冻结导出面收敛结果。
  - API 文档同步：`docs/api/models/solver/README.md`、`CoreConcepts.md`、`ImplicitSolvers.md` 移除 `create_implicit_solver` 作为公开推荐入口的描述。
- [x] 2026-04-03：执行 W3-C 收口（导出索引与治理校验，PR #50）。
  - 使用 `generate_api_export_index.jl` 重生成导出索引：`docs/api/generated/models/ModelsExportIndex.md` 与 `docs/api/models/solver/generated/Exports.md`，确保导出面剥离后索引与真实导出一致。
  - 治理检查通过：`scripts/dev/check_docs_consistency.jl`、`scripts/dev/check_active_docs_governance.jl`。
- [x] 2026-04-03：下一阶段启动（W4 入口一致性收口，PR #50）。
  - 在 W3-B 导出面收敛后修复残余调用：`tests/integration/models/test_solver_config_isolation_smoke.jl` 与 `tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl` 的 `create_implicit_solver` 调用统一切换为 `create_pnjl_implicit_solver`。
  - 同步收口 Phase CLI preset：`scripts/pnjl/calculate_phase_structure.jl` 的 `--preset=smoke` 默认 backend 改为 `:models`，并将 CLI 参数校验收紧为 `models|auto`，避免 W3-A 后 preset 仍回落到 legacy 配置。
  - 回归验证通过：上述 2 个 integration 定向测试通过；`tests/integration/models/test_phase_cli_smoke.jl` 通过；integration smoke `427/427`。
- [x] 2026-04-03：继续推进 W4 深化（phase/integration legacy 参数口径统一，PR #50）。
  - integration/phase 相关测试中 `solver_backend=:legacy` 统一切换为 `:models`（`test_phase_pipeline_smoke.jl`、`test_phase_artifacts_promotion_smoke.jl`、`test_models_unified_entrypoints_smoke.jl`、`test_pm_phase_diagnostic_smoke.jl`）；对应 unit `test_phase_pipeline.jl` 与 `test_pm_phase_diagnostic.jl` 同步收口。
  - 该批次确保 W3-A 后“legacy backend 已冻结”语义在 phase 与统一入口测试层不再出现反向示例，减少后续 W4/W5 回归噪声。
  - 回归验证通过：`test_phase_pipeline.jl`、`test_pm_phase_diagnostic.jl`、`test_phase_pipeline_smoke.jl`、`test_pm_phase_diagnostic_smoke.jl`，integration smoke `427/427`。
- [x] 2026-04-03：继续推进 W4 深化（CLI smoke 口径与 models backend 一致，PR #50）。
  - `test_wavee_unified_scan_cli_smoke.jl` 中 tmu/trho/phase 三条 CLI 烟测命令统一改为 `--solver_backend=models`，对齐 W3-A 冻结后运行事实。
  - `test_phase_cli_direct_smoke.jl` 内 direct smoke、默认模板运行、manifest 写入与 preset 断言统一收口到 models（含 config fixture 中 `solver_backend` 从 `legacy` 改为 `models`）。
  - 回归验证通过：`test_wavee_unified_scan_cli_smoke.jl`、`test_phase_cli_direct_smoke.jl`。
- [x] 2026-04-03：继续推进 W4 深化（models 单测遗留口径清理，PR #50）。
  - `test_pm_phase_diagnostic.jl` 中 Maxwell 参考 wiring 数据夹具的 `solver_backend` 字段及同文件相关调用参数由 `:legacy` 统一改为 `:models`，避免语义冻结后测试数据仍反映旧后端。
  - 验证通过：`tests/unit/models/test_pm_phase_diagnostic.jl`。
- [x] 2026-04-03：继续推进 W4 深化（scripts 层首批清理，PR #50）。
  - `scripts/pnjl/diagnose_pm_phase.jl` 默认 backend 改为 `:models`，并在参数解析阶段显式拒绝非 models 取值，帮助脚本行为与 PM 诊断主链路冻结语义一致。
  - `scripts/dev/export_phase_pipeline_regression_baseline.jl` 将 phase pipeline baseline 导出默认 backend 切换为 `:models`。
  - `scripts/analysis/relaxtime_xi_fluctuation_study.jl` 移除 equilibrium 的 legacy fallback 路径，统一使用 models backend。
  - 验证通过：`tests/integration/models/test_pm_phase_diagnostic_smoke.jl`（覆盖 diagnose_pm_phase CLI 路径），以及 `scripts/dev/export_phase_pipeline_regression_baseline.jl --help`。
- [x] 2026-04-03：继续推进 W4 深化（relaxtime scripts fallback 清理，PR #50）。
  - `scripts/relaxtime/run_gap_transport_scan.jl` 移除 `solve_legacy_equilibrium` 分支与 metadata 中 legacy fallback 描述；当 models 路径无有效解时改为显式抛错并由单点容错逻辑记录。
  - `scripts/relaxtime/scan_relaxation_times_vs_T.jl` 与 `scripts/relaxtime/scan_transport_vs_xi_T150_muB800.jl` 去除 models→legacy 回退块，统一仅走 models equilibrium backend。
  - 同步修复 `scan_transport_vs_xi_T150_muB800.jl` 的 `Models` include 依赖，确保脚本入口可独立加载。
  - 验证通过：`scripts/relaxtime/run_gap_transport_scan.jl --help`、`scripts/relaxtime/scan_relaxation_times_vs_T.jl --help`、`scripts/relaxtime/scan_transport_vs_xi_T150_muB800.jl --help`。
- [x] 2026-04-03：继续推进 W4 深化（dev baseline export 脚本收口，PR #50）。
  - `scripts/dev/export_transport_fixedpoint_baseline.jl` 默认 backend 改为 `:models`，并将 `--backend` 校验收紧为仅接受 models。
  - 移除脚本内 legacy 分支导出路径，统一以 models backend 生成 fixed-point baseline。
  - 增补 `Models.jl` include，确保脚本可稳定解析 `Main.Models.NLsolveGapSolver` 依赖。
  - 验证通过：`scripts/dev/export_transport_fixedpoint_baseline.jl --help`。
- [x] 2026-04-03：继续推进 W4 深化（scripts/dev remaining legacy 口径清理，PR #50）。
  - `scripts/dev/export_pnjl_scan_fixedpoint_baseline.jl` 中 tmu/trho baseline 导出默认后端由 `thermo_backend=:legacy` 改为 `:models`。
  - `scripts/dev/profile_paramtypes_hotspots.jl` 中 `_transport_inputs_from_equilibrium` 的热力学后端参数由 `:legacy` 改为 `:models`，与当前主线口径保持一致。
  - 验证通过：`scripts/dev/export_pnjl_scan_fixedpoint_baseline.jl --help`。
