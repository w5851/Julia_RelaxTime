# Models 求解器主链收敛与 ImplicitSolver 退役任务单

更新日期：2026-04-04

目标：
- 在不改变物理结果与回归基线的前提下，将求解器能力收敛到 models 主链；
- 最终移除 `src/models/solver/ImplicitSolver.jl` 的运行时依赖。

---

## 0. 当前共识（会话记录）

- `solve_constraint` 主链已存在（`Solver.jl -> ProblemSpec/constraint_solver.jl`）。
- 扫描链路默认使用 models 路径，但 `Trho` 的 rescue 仍通过 `_solve_weighted_block_fallback` 间接触达 `ImplicitSolver`。
- 隐式求导主实现已在 `src/models/implicit_gap.jl`（`create_pnjl_implicit_solver` / `solve_pnjl_with_derivatives`）。
- `ImplicitSolver.solve_with_derivatives` 当前是转发壳。

---

## 1. 迁移步骤（最小风险序列）

### Step 1：断开求解主入口对 `ImplicitSolver` 的硬依赖

- [ ] 将 `src/models/solver/Solver.jl` 中 `solve(...)` / `solve_multi(...)` 从 `ImplicitSolver.*` 转发，改为统一走 `solve_constraint(...)` 主链（必要时保留参数适配层）。
- [ ] 将 `is_physical_solution` 收敛为 models 主链可见规则，不再经 `ImplicitSolver.default_is_physical_solution` 转发。

验收（最小集）：
- [ ] `tests/unit/models/test_constraint_solver.jl`
- [ ] `tests/unit/models/test_problem_spec_contract.jl`
- [ ] `tests/integration/models/test_solver_auto_backend_semantic_parity.jl`

### Step 2：将 weighted fallback 平移到 models 主链

- [ ] 将 `solve_weighted_block_fallback` 从 `ImplicitSolver` 迁移到 models 侧（建议新建 `src/models/solver/WeightedFallback.jl` 或并入 `constraint_solver.jl`）。
- [ ] `src/models/scans/TrhoScan.jl` 直接调用 models 侧 fallback 能力，不再经 `Solver.jl` 私有桥 `_solve_weighted_block_fallback`。
- [ ] 删除 `Solver.jl` 中 `_solve_weighted_block_fallback` 桥接。

验收（最小集）：
- [ ] `tests/integration/pnjl/test_trho_scan_semantic_modes_smoke.jl`
- [ ] `tests/integration/pnjl/test_trho_scan_solver_backend_models_smoke.jl`
- [ ] 对应 regression smoke

### Step 3：收敛隐式求导入口，仅保留 `implicit_gap.jl` 主实现

- [ ] 移除 `ImplicitSolver.solve_with_derivatives` 对外导出（或降级 internal）。
- [ ] 统一对外入口为 `Models.solve_pnjl_with_derivatives` 与 `Models.create_pnjl_implicit_solver`。
- [ ] 清理 `Solver.jl` 中与导数相关的重复包装（保留一个稳定入口）。

验收（最小集）：
- [ ] `tests/unit/models/test_implicit_gap.jl`
- [ ] `tests/integration/models/test_models_implicitdiff_flavor_mu_smoke.jl`
- [ ] `tests/unit/pnjl/test_solver_implicit.jl`

### Step 4：移除 `ImplicitSolver` 模块并收口导出/文档

- [ ] 删除 `src/models/solver/ImplicitSolver.jl`。
- [ ] 清理 `Models.jl` 中 include/export 与所有 `ImplicitSolver.` 调用残留。
- [ ] 更新 `docs/api/models/solver/*` 与 generated 导出索引。

验收（全量收口）：
- [ ] unit/integration/regression smoke 全通过
- [ ] `julia --project=. scripts/dev/check_docs_consistency.jl`
- [ ] `julia --project=. scripts/dev/check_legacy_solver_switch_leakage.jl`

---

## 2. DoD（完成判据）

- [ ] 全仓代码调用中无 `ImplicitSolver.`（仅允许历史文档文字）。
- [ ] scan / constraint / derivative 三类能力均通过 models 主链到达。
- [ ] PR 证据包含：旧入口移除、主链稳定、测试闭环。

---

## 2.5 实现计划（基于当前讨论，待执行）

说明：本计划先聚焦“能力收敛与维度契约化”，不改变物理含义与数值语义。

### Task 1：维度契约最小 API 落地（先统一元信息入口）

目标：把“5/3/8”从散落硬编码收敛为统一访问接口。

Files:
- Modify: `src/models/solver/StateSchema.jl`
- Modify: `src/models/solver/ConstraintModes.jl`
- Modify: `src/models/solver/ConstraintComponents.jl`
- Test: `tests/unit/models/test_solver_dimension_agnostic.jl`

- [ ] 增加维度访问函数（例如 `state_dim(schema)`, `mu_dim(mode)`, `total_dim(mode, schema)`）。
- [ ] 在不改行为前提下，为 PNJL/RPNJL 约束模式返回与现状一致的维度。
- [ ] 在单测中加入维度契约断言，覆盖 FixedMu/FixedRho/FixedEntropy/FixedSigma/FixedAsymmetricRho。

验证命令：
- `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver_dimension_agnostic.jl"; include("tests/unit/runtests.jl")'`

### Task 2：`Conditions.jl` 去硬编码（A 档第一优先）

目标：让条件函数只依赖 schema+mode，不直接切片 `x[1:5]` / `x[6:8]`。

Files:
- Modify: `src/models/solver/Conditions.jl`
- Test: `tests/unit/models/test_problem_spec_contract.jl`
- Test: `tests/unit/models/test_constraint_solver.jl`

- [x] 让 `build_conditions(::FixedRho, params)` 成为 schema 版 wrapper（内部统一走 schema 路径）。
- [x] 对 `FixedAsymmetricRho/FixedEntropy/FixedSigma` 同步做 schema 化切片。
- [x] 将 `_gap_conditions_dynamic` 从“固定 5/3”改为“由 schema/mode 约束验证”，并保留清晰错误信息。

验证命令：
- `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_problem_spec_contract.jl,models/test_constraint_solver.jl"; include("tests/unit/runtests.jl")'`

### Task 3：`constraint_solver.jl` pack/unpack 统一化（A 档第二优先）

目标：把 candidate / solution 的固定维度构造替换为维度契约 API。

Files:
- Modify: `src/models/constraint_solver.jl`
- Modify: `src/models/solver/StateSchema.jl`（如需补 helper）
- Test: `tests/unit/models/test_constraint_solver.jl`

- [x] 引入统一 `pack_solution(...)` / `unpack_solution(...)` / `empty_candidate(...)` helper。
- [x] 替换 `_solve_constraint_fixedrho/_fixedentropy/_fixedsigma/_fixedasymrho` 中 `SVector{5}/SVector{3}/Float64[...]` 散点拼装。
- [x] 保证返回字段与当前契约完全一致。

验证命令：
- `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_constraint_solver.jl,models/test_problem_spec_contract.jl"; include("tests/unit/runtests.jl")'`

### Task 4：联合方程组试点（FixedRho）

目标：验证“展平求解”在单模式上的可行性与稳定性。

Files:
- Modify: `src/models/constraint_solver.jl`
- Modify: `src/models/solver/ProblemSpec.jl`
- Add/Modify: `tests/integration/models/test_solver_auto_backend_semantic_parity.jl`
- Add/Modify: `tests/regression/models/test_problem_spec_fixedrho_parity_regression.jl`

- [x] 新增 `FixedRho` 联合残差构造（未知量与方程组一次性求解）。
- [x] 保留旧路径开关用于 A/B 对比（默认不切换行为）。
- [x] 输出对比证据：收敛率、残差、耗时、结果偏差（容差内）。

验证命令：
- `julia --project=. -e 'ENV["INTEGRATION_FILES"]="models/test_solver_auto_backend_semantic_parity.jl"; include("tests/integration/runtests.jl")'`
- `julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_problem_spec_fixedrho_parity_regression.jl"; include("tests/regression/runtests.jl")'`

### Task 5：扫描链路与 fallback 对齐（Step 2 衔接）

目标：把 `Trho` rescue 从 `Solver.jl` 私有桥迁走。

Files:
- Modify: `src/models/scans/TrhoScan.jl`
- Modify/Add: `src/models/solver/WeightedFallback.jl`（或并入 `constraint_solver.jl`）
- Modify: `src/models/solver/Solver.jl`
- Test: `tests/integration/pnjl/test_trho_scan_semantic_modes_smoke.jl`

- [x] 提供 models 侧明确 fallback 能力入口。
- [x] `TrhoScan` 直接调用该入口，不经 `_solve_weighted_block_fallback`。
- [x] 删除 `Solver.jl` 私有桥并更新对应测试。

验证命令：
- `julia --project=. -e 'ENV["INTEGRATION_FILES"]="pnjl/test_trho_scan_semantic_modes_smoke.jl,pnjl/test_trho_scan_solver_backend_models_smoke.jl"; include("tests/integration/runtests.jl")'`

### Task 6：收敛导数入口 + 预备退役 `ImplicitSolver`

目标：导数能力仅保留 `implicit_gap.jl` 主实现对外入口。

Files:
- Modify: `src/models/solver/Solver.jl`
- Modify: `src/models/solver/ImplicitSolver.jl`
- Modify: `src/models/Models.jl`
- Test: `tests/unit/models/test_implicit_gap.jl`

- [x] 去掉 `ImplicitSolver.solve_with_derivatives` 对外可见性（或降 internal）。
- [x] 对外 API 统一为 `Models.solve_pnjl_with_derivatives` 与 `Models.create_pnjl_implicit_solver`。
- [x] 补文档与测试口径一致性。

验证命令：
- `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_implicit_gap.jl,pnjl/test_solver_implicit.jl"; include("tests/unit/runtests.jl")'`

### Task 7：全量收口与退役检查

目标：确保迁移后可安全删除 `ImplicitSolver.jl`（若本轮覆盖完整则执行删除，否则形成下一轮 gate）。

Files:
- Modify: `docs/api/models/solver/*`
- Modify: `docs/api/generated/models/ModelsExportIndex.md`（脚本生成）

- [x] 跑全套 smoke + docs/governance。
- [x] 生成“剩余 Implicit 依赖点”报告（必须为空才可删文件）。
- [x] 若不为空，记录 blockers 与下一轮拆分项。

验证命令：
- `julia --project=. -e 'ENV["UNIT_PROFILE"]="smoke"; include("tests/unit/runtests.jl")'`
- `julia --project=. -e 'ENV["INTEGRATION_PROFILE"]="smoke"; include("tests/integration/runtests.jl")'`
- `julia --project=. -e 'ENV["REGRESSION_PROFILE"]="smoke"; include("tests/regression/runtests.jl")'`
- `julia --project=. scripts/dev/check_docs_consistency.jl`
- `julia --project=. scripts/dev/check_legacy_solver_switch_leakage.jl`

---

## 3. 执行记录（append-only）

- 2026-04-04：根据会话共识建立本任务单，记录 Step 1~4 最小风险迁移序列与 DoD。
- 2026-04-04：补充“联合方程组展平求解”优化路线共识：
  1) 先抽“按模式生成联合残差方程组”的统一接口，不立即替换全量求解器；
  2) 先在 `FixedRho` 做展平联合求解试点，与现有路径对比收敛/耗时/稳健性；
  3) 试点稳定后再推广到 `FixedEntropy/FixedSigma/FixedAsymmetricRho`，最后删除外层-内层分治路径。
- 2026-04-04：完成 Task2（`Conditions.jl` schema+mode 去硬编码）并补充单测：
  1) `build_conditions(::FixedRho, params)` 改为 schema wrapper；
  2) 为 `FixedAsymmetricRho/FixedEntropy/FixedSigma` 新增 schema 版 `build_conditions`；
  3) `_gap_conditions_dynamic` 增加 schema/mode 维度契约校验与清晰报错；
  4) 新增 non-rho schema 路径 parity 与 mismatch 报错单测。
- 2026-04-04：Task3 进行中：
  1) 在 `constraint_solver.jl` 新增 `_pack_solution/_unpack_solution/_empty_candidate`；
  2) 已将多个模式返回结果中的 `Float64[...]` 散点拼装替换为 `_pack_solution(...)`；
  3) 已补充 helper 单测；下一步继续收敛内部 `SVector`/切片构造。
- 2026-04-04：Task3 收口完成：
  1) `_solve_constraint_fixedrho/_fixedentropy/_fixedsigma/_fixedasymrho` 返回值统一走 `_pack_solution(...)`；
  2) 引入 `_to_state_svec/_to_mu_svec/_mass_from_state` 收敛 `SVector` 与 mass 构造散点；
  3) 目标验证通过：
     - `models/test_problem_spec_contract.jl` + `models/test_constraint_solver.jl`（155/155）；
     - `tests/integration/models/test_solver_auto_backend_semantic_parity.jl`（direct include 4/4）。
- 2026-04-04：Task4 启动（阶段化落地，先 gate 后联合求解）：
  1) 在 `ProblemSpec` FixedRho forward_solve 增加 `fixedrho_joint_solve::Bool` 参数校验；
  2) 作为 A/B 开关占位，当前默认和开启时均走现有路径，并回传 `fixedrho_joint_solve_requested/fixedrho_joint_solve_active` 元信息；
  3) 新增/更新测试覆盖 flag 语义与回归口径：
     - `tests/unit/models/test_problem_spec_contract.jl`；
     - `tests/regression/models/test_problem_spec_fixedrho_parity_regression.jl`；
  4) 验证通过：
     - unit 目标集（160/160）；
     - regression 目标集（27/27）。
- 2026-04-04：Task4 第二阶段完成（FixedRho 联合求解试点）：
  1) 在 `ProblemSpec` 中新增 `_fixedrho_joint_problem_spec_forward_solve`，基于 `Conditions.build_residual!(FixedRho, ...)` 做 8 维一次性 NLsolve；
  2) `fixedrho_joint_solve=true` 时优先尝试联合求解；若残差不可接受或失败，自动回退到既有 `_solve_constraint_fixedrho`（保序保稳）；
  3) 返回元信息扩展：`fixedrho_joint_solve_requested / fixedrho_joint_solve_active / fixedrho_joint_fallback`；
  4) 新增与更新测试：
     - `tests/unit/models/test_problem_spec_contract.jl`（joint+fallback 语义覆盖）；
     - `tests/regression/models/test_problem_spec_fixedrho_parity_regression.jl`（joint flag 回归口径）。
  5) 验证结果：
     - unit 目标集：`165/165`；
     - regression 目标集：`33/33`。
- 2026-04-04：Task5 完成（Trho fallback 链路对齐）：
  1) 新增 models 侧 `src/models/solver/WeightedFallback.jl`，导出 `solve_weighted_block_fallback`；
  2) `src/models/scans/TrhoScan.jl` 改为直接调用 `Main.Models.solve_weighted_block_fallback`；
  3) 删除 `src/models/solver/Solver.jl` 私有桥 `_solve_weighted_block_fallback`；
  4) `src/models/Models.jl` 增加 include/export；
  5) 更新单测口径（`test_solver_dimension_agnostic.jl`）；
  6) 验证通过：
     - unit: `models/test_solver_dimension_agnostic.jl,models/test_problem_spec_contract.jl`（163/163）；
     - integration direct file:
       - `tests/integration/pnjl/test_trho_scan_semantic_modes_smoke.jl`（6/6）
       - `tests/integration/pnjl/test_trho_scan_solver_backend_models_smoke.jl`（22/22 分组通过）。
- 2026-04-04：Task6 完成（导数入口收敛）：
  1) `src/models/solver/ImplicitSolver.jl` 移除 `solve_with_derivatives` export 与实现；
  2) `src/models/solver/Solver.jl` 中 `solve_with_derivatives` 改为转发 `solve_pnjl_with_derivatives`（models 主实现）；
  3) 更新口径测试：`tests/unit/models/test_solver_dimension_agnostic.jl` 增加 `:solve_with_derivatives ∉ names(Models.ImplicitSolver)`；
  4) 验证通过：
     - unit 目标集：`models/test_implicit_gap.jl,pnjl/test_solver_implicit.jl,models/test_solver_dimension_agnostic.jl`（123/123）；
     - integration direct file：`tests/integration/models/test_models_native_solver_phase1_smoke.jl`（11/11）。
- 2026-04-04：Task7 收口检查完成：
  1) Smoke 验证：
     - unit smoke：`781/781`；
     - integration smoke：`437/437`（首次运行因工具超时中断，重跑 with longer timeout 通过）；
     - regression smoke：`423/424`，其中 1 个 Broken 为已声明可选跳过（tau xi probe fixture 缺失）；
  2) Governance：
     - `scripts/dev/check_docs_consistency.jl` => OK；
     - `scripts/dev/check_legacy_solver_switch_leakage.jl` => OK；
  3) “剩余 Implicit 依赖点”清点（`src/`）：
     - `src/models/solver/Solver.jl` 仍有 `ImplicitSolver.*` 依赖：
       - `SolverResult = ImplicitSolver.SolverResult`
       - `solve/solve_multi` 多个入口仍转发 `ImplicitSolver.solve`/`solve_multi`
       - `is_physical_solution` 仍转发 `ImplicitSolver.default_is_physical_solution`
       - `src/models/Models.jl` 仍 include `solver/ImplicitSolver.jl`
  4) 结论：当前 **不可删除** `src/models/solver/ImplicitSolver.jl`。
  5) 下一轮 blockers（建议拆分）：
     - B1：将 `solve/solve_multi` 从 `ImplicitSolver` 转发迁移到 models 主链（对应 Step 1 尾项）。
     - B2：收敛 `SolverResult` 与 `is_physical_solution` 到 models 原生实现，解除对 `ImplicitSolver` 类型/函数依赖。
     - B3：完成 B1/B2 后再执行模块删除与 include/export 收口。
- 2026-04-04：B1 继续推进（分段迁移）：
  1) `Solver.jl` 中 `solve(model, FixedMu, ...)` 与 `solve_multi(model, FixedMu, ...)` 已迁移到 models 主链（`solve_constraint + selector`）；
  2) `FixedRho/FixedEntropy/FixedSigma/FixedAsymmetricRho` 在本轮尝试全量迁移时引入质量回归（integration smoke 失败），已回滚其迁移部分并保持旧行为；
  3) 为防止递归路径，`constraint_solver.jl` 中 `FixedRho` 失败回退已改为显式调用 `Main.Models.ImplicitSolver.solve(...)`（避免通过 `Main.Models.solve(...)` 再入新分发）；
  4) 验证通过：
     - `tests/unit/pnjl/test_solver_implicit.jl` 全通过；
     - `tests/unit/models/test_solver.jl` 全通过；
     - `tests/integration/models/test_models_native_solver_phase1_smoke.jl` 全通过；
     - `tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl` 全通过。
  5) 当前 B1 状态：**部分完成（FixedMu 主链化完成，其余模式待下一轮拆分迁移）**。
- 2026-04-04：B1-2 尝试记录（FixedRho 单独迁移实验）：
  1) 尝试将 `solve(model, FixedRho, ...)` / `solve_multi(model, FixedRho, ...)` 迁移到 models 主链；
  2) 回归观察到 `fixedrho semantic equivalence` 口径偏移（legacy 与 models 数值误差超阈，且 legacy convergence 断言不稳定）；
  3) 已回滚该迁移，恢复 FixedRho 走 `ImplicitSolver` 路径，确保当前门禁稳定全绿；
  4) 仍保留已稳定项：
     - `solve(model, FixedMu, ...)` 与 `solve_multi(model, FixedMu, ...)` 主链化；
     - `constraint_solver` 中 FixedRho fallback 明确调用 `Main.Models.ImplicitSolver.solve(...)`，避免递归。
- 2026-04-04：B2 推进完成（类型与物理判据去耦）：
  1) `src/models/solver/Solver.jl` 新增本地 `SolverResult` 结构定义，解除 `const SolverResult = ImplicitSolver.SolverResult` 依赖；
  2) 新增 `_coerce_solver_result(...)`，将仍由 `ImplicitSolver` 返回的结果统一转为 models 本地 `SolverResult`；
  3) `is_physical_solution(...)` 在 `Solver.jl` 本地实现，移除对 `ImplicitSolver.default_is_physical_solution` 的转发；
  4) 依赖检查：`Solver.jl` 内已无 `ImplicitSolver.SolverResult` / `ImplicitSolver.default_is_physical_solution` 引用；
  5) 验证通过：
     - unit：`pnjl/test_solver_implicit.jl,models/test_solver_dimension_agnostic.jl,models/test_constraint_solver.jl`（161/161）；
     - integration：`test_models_native_solver_phase1_smoke.jl`（11/11），`test_solver_constraints_models_backend_smoke.jl`（43/43）。
- 2026-04-04：B3 探索结论（本轮）
  1) 目标：清理 `solve/solve_multi` 对 `ImplicitSolver` 的剩余转发（非 FixedMu 模式）。
  2) 结果：对 `FixedRho/FixedEntropy/FixedSigma/FixedAsymmetricRho` 迁移到 models 主链后，`test_solver_constraints_models_backend_smoke.jl` 出现系统性回归（收敛状态与约束量显著偏移）。
  3) 处理：已回滚该段迁移，恢复上述模式继续走 `ImplicitSolver` 路径，确保门禁稳定；仅保留 B1/B2 已验证稳定改动。
  4) 当前可执行下一步建议：
     - 先做“语义桥接层”拆分（为 non-FixedMu 模式补齐与 `ImplicitSolver.solve` 一致的参数语义），再逐模式迁移；
     - 避免直接把 `solve_constraint` 默认语义替换到 `solve` 旧入口。
- 2026-04-04：语义桥接层第一步已落地（non-FixedMu）
  1) 在 `Solver.jl` 新增 `_resolve_nonfixedmu_bridge(...)` 与 `_resolve_solver_model(...)`：
     - 统一解析 `seed_guess/seed_candidates/semantic_mode/selector/rho0/model_kind`；
     - 当检测到 problem-spec 语义参数时切换到 `solve_constraint` 路径；
     - 默认（无桥接参数）仍走 `ImplicitSolver.solve`，保持旧行为稳定。
  2) 新增单测覆盖桥接参数可用性：`tests/unit/models/test_solver.jl`。
  3) 验证通过：
     - unit：`models/test_solver.jl,pnjl/test_solver_implicit.jl`（75/75）；
     - integration：`tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl`（43/43）。
- 2026-04-04：语义桥接层第二步已落地（non-FixedMu solve_multi）
  1) `solve_multi(model, Union{FixedRho,FixedAsymmetricRho}, ...)` 增加桥接语义：
     - 当携带 `seed_guess/seed_candidates/semantic_mode/selector/problem_spec` 等语义参数时，走 `solve_constraint` 多种子筛选路径；
     - 默认无桥接参数时仍走 `ImplicitSolver.solve_multi`，保持旧行为。
  2) 新增单测覆盖：`tests/unit/models/test_solver.jl` 增加 non-FixedMu `solve_multi` 桥接参数测试。
  3) 验证通过：
     - unit：`models/test_solver.jl,pnjl/test_solver_implicit.jl`（77/77）；
     - integration：
       - `tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl`（43/43）
       - `tests/integration/models/test_models_native_solver_phase1_smoke.jl`（11/11）。
- 2026-04-04：B3 主迁移进展（当前轮）
  1) 已完成：
     - non-FixedMu `solve` / `solve_multi` 的语义桥接路径可用并稳定（显式语义参数时走主链）；
     - `solve_multi(FixedRho/FixedAsymmetricRho)` 默认路径改为“多种子 + 调 `solve` 聚合筛选”，不再直接依赖 `ImplicitSolver.solve_multi` 语义。
  2) 保守回退策略：
     - `solve` 默认路径（无语义桥接参数）对 `FixedRho/FixedEntropy/FixedSigma/FixedAsymmetricRho` 仍保留 `ImplicitSolver.solve`，以避免已确认的约束回归。
  3) 当前剩余 direct forwarding（`Solver.jl`）：
     - `ImplicitSolver.solve(mode, T_fm; ...)`（默认 non-FixedMu path）
     - `ImplicitSolver.solve_multi(mode, T_fm; ...)`（保护性尾路径，当前主分支不走）
  4) 下一步建议：
     - 为 `FixedRho/FixedEntropy/FixedSigma` 建立“默认语义等价桥接”后，再替换默认路径；
     - 最后再处理 `FixedAsymmetricRho` 默认路径与 weighted-fallback 协同。

- 2026-04-05：B3 默认路径收口（solve_multi 尾转发移除）
  1) 已完成：
     - `Solver.jl` 中 non-FixedMu `solve_multi` 保护性尾路径已移除，不再调用 `ImplicitSolver.solve_multi(mode, T_fm; ...)`；
     - 该分支改为显式 `ArgumentError("unsupported mode for solve_multi bridge: ...")`，避免静默回落与语义漂移。
  2) 当前 direct forwarding 余量：
     - `ImplicitSolver.solve(mode, T_fm; ...)`（默认 non-FixedMu `solve` 路径，保守保留）。
  3) 验证通过：
     - `julia --project=. -e 'include("tests/unit/models/test_solver.jl")'`
     - `julia --project=. -e 'include("tests/unit/pnjl/test_solver_implicit.jl")'`
     - `julia --project=. -e 'include("tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl")'`
     - `julia --project=. -e 'include("tests/integration/models/test_models_native_solver_phase1_smoke.jl")'`
  4) 下一步建议（保持低风险顺序）：
     - 为 `FixedEntropy/FixedSigma` 先补“默认语义等价桥接”的最小验证集（点位+多 seed 选择一致性）；
     - 再尝试替换 `solve` 默认 non-FixedMu 对 `ImplicitSolver.solve` 的最后一处转发。
