# Solver 五大痛点治理 PR-D 任务单（Seed SoT + Residual Spine）

## 1. 目标

- [ ] 建立种子单一事实来源（SoT），移除 `GapSolver.jl` 的重复常量与重复候选列表。
- [ ] 建立残差单一主路径（Residual Spine），避免相同约束在多处各写一套。

## 2. 范围

### 2.1 本期范围

- [ ] D1.1 统一 PNJL 默认种子常量定义到 `SeedStrategies`（仅保留一处定义）。
- [ ] D1.2 `GapSolver` 的多种子求解改为消费 `SeedStrategies.get_all_seeds(MultiSeed(), ...)`（或等价统一入口）。
- [ ] D1.3 删除/停用 `GapSolver` 内重复种子常量与 `_PNJL_MULTI_SEEDS`。
- [ ] D2.1 确立残差主入口（建议：`gap_core_residual!` 为核心残差引擎）。
- [ ] D2.2 `ConstraintSolver` 中 mode-specific `residual_fn!` 优先组合核心残差入口，而非重写同构逻辑。
- [x] D2.3 明确 `gap_residual` 在 NJL/PNJL 下的角色边界（是前端 API 还是核心引擎），并文档化。

### 2.2 非范围

- [ ] 不改变求解策略（Newton/TR/多 seed）默认行为优先级。
- [ ] 不在本 PR 内调整 selector 的业务偏好（pressure/residual 优先级）。

## 3. 设计约束

- [ ] 外部 API 行为保持兼容：`solve_gap`/`solve` 调用方式不变。
- [ ] 在相同 seed 输入下，数值漂移仅允许在公差内，且需可解释。
- [ ] 新增 helper 需避免循环依赖（`SeedStrategies` 不反向依赖 `GapSolver`）。

## 3.1 深层治理判据（防止只做表层整理）

- [ ] 必须删除 `GapSolver` 中默认 seed 常量副本，而不是仅转发/别名保留副本。
- [ ] 必须收敛到“一个默认 seed 目录 + 一个消费入口”，禁止 `GapSolver`/`WeightedFallback` 再定义隐式默认池。
- [ ] 必须将残差主链定义为“可被复用的唯一实现”，并为旧路径加退役注记或删除。
- [ ] 必须用测试证明“跨入口同点同序候选一致”，而不只验证单一路径可运行。

## 3.2 唯一实现点声明（本 PR 结束时）

- [ ] Seed SoT：`SeedStrategies`（含默认候选顺序、去重与扩展语义）。
- [ ] Residual Spine：`Conditions.gap_core_residual!`（或最终命名的等价单点）。

## 4. 实施任务（可勾选）

### 4.1 Seed SoT

- [ ] 建立 `seed_catalog(mode/context)`（命名可调整）统一返回默认候选序列。
- [ ] `GapSolver.solve_gap(::AbstractPNJLModel, ...)` 多种子分支改为使用 catalog。
- [ ] `WeightedFallback` 的 `MultiSeed` 候选获取也复用 catalog（若当前已有复用，补契约测试即可）。
- [ ] 增加“同点不同入口 seed 顺序一致性”单测。

### 4.2 Residual Spine

- [ ] 整理 `Conditions.build_residual!` 与 `gap_core_residual!` 的职责说明与调用关系。
- [ ] 将 `ConstraintSolver` 的重复 gap 残差拼装改为调用统一入口。
- [ ] 为 FixedMu/FixedRho/FixedEntropy/FixedSigma/FixedAsymmetricRho 增加残差等价性测试（旧路径 vs 新路径）。

## 5. 验证清单

- [ ] 单测：新增 `tests/unit/models/test_seed_catalog_contract.jl`（或等价文件）。
- [ ] 单测：新增 `tests/unit/models/test_residual_spine_contract.jl`（或等价文件）。
- [ ] 回归：
  - `models/test_solver_attempt_engine_convergence_regression.jl`
  - `pnjl/test_constraint_selection_regression.jl`
- [ ] 针对 5 类 mode 的代表点做 before/after 快照对比（记录 residual_norm/selection_reason/seed_index）。
- [ ] 新增“唯一实现点守护测试”：扫描断言 `GapSolver` 不再定义默认 PNJL seed 常量。
- [ ] 新增“残差主链守护测试”：主要调用链经过统一 residual 入口（可用契约/spy 方式验证）。

## 6. 风险与缓解

- 风险：seed 顺序统一后触发边缘点分支切换。
  - 缓解：固定候选顺序与去重键；保留 explicit seed 覆盖通道。
- 风险：残差路径收敛导致某些模式容差边界变紧或变松。
  - 缓解：保留回归基线点，必要时分模式设定兼容阈值并记录理由。

## 7. PR-D DoD

- [ ] `GapSolver` 内不再维护 PNJL 默认种子常量副本。
- [ ] 残差主路径被明确且在主要调用链唯一化。
- [ ] 关键 unit/regression 通过，且无不可解释漂移。
- [ ] 代码审计可证明“并行 seed/residual 通道数量下降”，并附删除路径清单。

## 8. 进展记录（2026-04-08）

- [x] D1.1 已完成：新增 `SeedStrategies.default_multiseed_candidates_5()` 作为默认 5 维候选单一来源。
- [x] D1.2 已完成：`GapSolver.solve_gap(::AbstractPNJLModel, ...)` 多种子路径改为消费 `default_multiseed_candidates_5()`。
- [x] D1.3 已完成：删除 `GapSolver` 中 `_HADRON_SEED_5` 等重复常量与 `_PNJL_MULTI_SEEDS`。
- [x] 4.1 seed catalog 已落地（当前为 `default_multiseed_candidates_5()` 命名）。
- [x] 4.1 GapSolver 多种子分支已切换到 catalog。
- [x] 4.1 增加守护测试：`tests/unit/pnjl/test_solver_seed_strategies.jl` 新增 `Seed SoT guard`，覆盖
  - catalog 合约
  - `GapSolver` 不再持有默认种子常量副本。
- [ ] 4.1 `WeightedFallback` 显式改用 catalog（当前通过 `MultiSeed()` 间接复用，后续可做显式收敛）。
- [ ] D2 Residual Spine 尚未开始。
- [x] D2.1 已完成：新增 `ConstraintSolver._gap_norm_from_state(...)`，统一以 `gap_core_residual!` 计算 gap 范数。
- [x] D2.2 已完成：`ConstraintSolver` 内 FixedMu/FixedRho/FixedEntropy/FixedSigma/FixedAsymmetricRho 的 gap-norm 计算均改为复用 `_gap_norm_from_state`，移除局部 `gap_residual(...)` 直算路径。
- [x] D2.3 已完成：在 `docs/api/models/solver/ConstraintModes.md` 补充角色边界说明：
  - `gap_residual`：模型级公开残差入口（API/适配层）；
  - `gap_core_residual!`：solver 主链统一残差内核；
  - `ConstraintSolver` 的 gap-norm 统一复用主链内核。

### 本轮验证

- [x] `julia --project=. -e 'include("tests/unit/pnjl/test_solver_seed_strategies.jl")'`
- [x] `julia --project=. -e 'include("tests/unit/models/test_gap_solver.jl")'`
- [x] `julia --project=. -e 'include("tests/unit/models/test_solver.jl")'`
- [x] `julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_solver_attempt_engine_convergence_regression.jl;pnjl/test_constraint_selection_regression.jl"; include("tests/regression/runtests.jl")'`
- [x] `julia --project=. -e 'include("tests/unit/models/test_constraint_solver.jl")'`
- [x] `julia --project=. -e 'include("tests/unit/models/test_gap_core_residual.jl")'`
- [x] `julia --project=. -e 'include("tests/unit/models/test_gap_core_residual_njl_parity.jl")'`

### 删除重复路径清单（本轮）

- [x] 删除 `ConstraintSolver.jl` 中 6 处 `gap_residual(model, ...)` 直接计算 gap-norm 的并行路径。
- [x] 统一替换为 `_gap_norm_from_state -> gap_core_residual!` 主链。
