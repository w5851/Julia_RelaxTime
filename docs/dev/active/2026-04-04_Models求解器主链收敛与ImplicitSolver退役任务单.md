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

- 2026-04-05：B3 默认路径替换（Entropy/Sigma）
  1) 红灯验证（先证伪）：
     - 新增 `tests/unit/models/test_solver.jl` 中两组 parity 用例：
       - `FixedEntropy default/bridge semantic parity (single point)`
       - `FixedSigma default/bridge semantic parity (single point)`
     - 在替换前，测试稳定失败：默认 `solve` 可收敛而桥接 `solve_constraint` 路径不收敛，确认语义差距真实存在。
  2) 绿灯修复（最小改动）：
     - 在 `constraint_solver.jl` 的 `_solve_constraint_fixedentropy/_solve_constraint_fixedsigma` 中增加 PNJL 兜底：
       - 当 models 外层约束求解未达收敛时，显式回退 `Main.Models.ImplicitSolver.solve(...)`；
       - 回填为统一返回结构，保证 `solve_constraint` 调用方语义稳定。
     - 在 `Solver.jl` 中将 `FixedEntropy/FixedSigma` 默认 `solve` 路径切至 ProblemSpec 主链（不再走 `ImplicitSolver.solve` 直转发）。
  3) 验证通过：
     - `julia --project=. -e 'include("tests/unit/models/test_solver.jl")'`（17/17）
     - `julia --project=. -e 'include("tests/unit/pnjl/test_solver_implicit.jl")'`
     - `julia --project=. -e 'include("tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl")'`（43/43）
     - `julia --project=. -e 'include("tests/integration/models/test_models_native_solver_phase1_smoke.jl")'`（11/11）
  4) 当前 direct forwarding 余量（`Solver.jl`）：
     - `ImplicitSolver.solve(mode, T_fm; ...)`：仅保留在默认 `FixedRho/FixedAsymmetricRho` 路径。

- 2026-04-05：B3 收口完成（Rho/Asym 默认路径）
  1) 红灯验证（先证伪）：
     - 在 `tests/unit/models/test_solver.jl` 新增两组 parity 用例：
       - `FixedRho default/bridge semantic parity (single point)`
       - `FixedAsymmetricRho default/bridge semantic parity (single point)`
     - 初始失败，确认默认路径与桥接主链在该两模式上仍存在语义差距。
  2) 绿灯修复（最小收口）：
     - `constraint_solver.jl`：
       - `_solve_constraint_fixedasymrho` 增加 PNJL 兜底（失败时显式回退 `Main.Models.ImplicitSolver.solve(...)`）；
       - 保持 `_solve_constraint_fixedrho` 既有兜底路径；
       - 结合上一轮 Entropy/Sigma 兜底，形成 non-FixedMu 主链可用闭环。
     - `Solver.jl`：
       - 默认 `solve(model, mode::Union{FixedRho,FixedAsymmetricRho,FixedEntropy,FixedSigma}, ...)` 统一改走 ProblemSpec 主链；
       - 删除 `solve` 中 non-FixedMu 对 `ImplicitSolver.solve` 的最后 direct forwarding 分支。
  3) 验证通过：
     - `julia --project=. -e 'include("tests/unit/models/test_solver.jl")'`（25/25）
     - `julia --project=. -e 'include("tests/unit/pnjl/test_solver_implicit.jl")'`
     - `julia --project=. -e 'include("tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl")'`（43/43）
     - `julia --project=. -e 'include("tests/integration/models/test_models_native_solver_phase1_smoke.jl")'`（11/11）
  4) 状态结论：
     - `Solver.jl` 已无 `ImplicitSolver.solve`/`solve_multi` direct forwarding；
     - 剩余 `ImplicitSolver` 依赖已下沉到 `constraint_solver` 的兼容兜底层，用于维持语义与回归稳定。

- 2026-04-05：Step4 强收口尝试与回滚结论（阻塞确认）
  1) 本轮尝试内容：
     - 试图删除 `src/models/solver/ImplicitSolver.jl`，并清理 `Models.jl` include 与测试中对 `Models.ImplicitSolver` 的访问。
  2) 实测结果（失败）：
     - non-FixedMu 默认语义出现系统性回归：`tests/unit/models/test_solver.jl` 中 `FixedEntropy/FixedSigma/FixedRho/FixedAsymmetricRho` parity 大面积失败；
     - 独立复现显示默认 `FixedEntropy` 可退化到不收敛（`converged=false`，残差显著偏大）。
  3) 根因判断：
     - 当前主链对 non-FixedMu 默认语义仍存在历史算法依赖；`constraint_solver` 侧尚未具备与历史路径完全等价的稳定实现；
     - 直接移除 `ImplicitSolver` 会破坏既有默认行为与回归基线。
  4) 处理决定：
     - 已按低风险策略回滚本轮 Step4 强删尝试，恢复到 `cfb7228` 后的稳定语义状态；
     - 保持 B3 已完成成果（`Solver.jl` 无 direct forwarding）不变。
   5) 当前任务单状态：
      - **未闭环**（Step4 / DoD 尚未满足）；
      - 下一轮需单独立项完成 non-FixedMu 默认语义等价重构后，方可真正删除 `ImplicitSolver.jl`。

- 2026-04-05：FixedRho joint 策略重构（按“展平与求解策略解耦 + 非 FixedMu 懒触发”执行）
  1) 根因复核结论：
     - `ImplicitSolver` 在 `FixedRho` 上是 8 维联合方程单层求解（非双层 μ 外层 + gap 内层）；
     - 先前 joint 不稳定主因是“求解策略耦合与默认方法不匹配”（非展平方程本身错误）。
  2) 代码策略调整（`src/models/solver/ProblemSpec.jl`）：
     - 保留 joint 负责展平残差（`build_residual!(FixedRho, ...)`），求解策略统一走 `solve_root_with_policy(...)`；
     - `FixedRho` 引入 staged/lazy attempts：
       - 先 primary 单次尝试；
       - 失败后才触发 method-rescue / seed-rescue；
       - 命中满足 hard-constraints 的候选后提前停止（early-stop）；
     - 方法默认值改为“连续性优先”自适应：
       - `continuity_seed=true` 默认 `:newton`；
       - 否则默认 `:trust_region`；
     - 增补 joint 诊断元信息：
       - `fixedrho_joint_selected_method`
       - `fixedrho_joint_selected_quality`
       - `fixedrho_joint_fallback_used`
       - `fixedrho_joint_attempt_origin`。
  3) 测试口径更新：
     - regression：`tests/regression/models/test_problem_spec_fixedrho_parity_regression.jl`
       - 新增 continuity-like `joint-only` 稳定性回归（滚动 seed）；
       - 保留/强化 joint parity gates（no-fallback 与 with-fallback 分支口径）。
     - unit：`tests/unit/models/test_problem_spec_contract.jl`
       - 同步 FixedRho forward_solve 合同与元信息断言。
     - unit：`tests/unit/models/test_solver.jl`
       - `FixedRho default/bridge semantic parity (single point)` 调整为与默认 seed 语义一致（bridge 使用 `get_seed(...)` 的 mode-extended 初值），避免“手工未扩展 seed”造成非目标分支。
  4) 本轮验证通过：
     - `julia --project=. -e 'include("tests/regression/models/test_problem_spec_fixedrho_parity_regression.jl")'`
     - `julia --project=. -e 'include("tests/unit/models/test_problem_spec_contract.jl")'`
     - `julia --project=. -e 'include("tests/unit/models/test_solver.jl")'`（25/25）
     - `julia --project=. -e 'include("tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl")'`（43/43）。
  5) 当前状态与约束：
     - `FixedMu` 强制多 attempts 选优逻辑保持不变；
     - 非 `FixedMu` 按“懒触发 attempts”推进，优先对齐真实扫描使用场景（continuity seed）；
     - `fixedrho_joint_solve=false` 兼容分支仍保留（便于过渡期对照/回归），尚未进行物理语义层面的全量 legacy 退役。

- 2026-04-05：FixedRho legacy 分支收口（ProblemSpec 侧）
  1) 处理原则：在 `FixedRho` 路径中继续遵循“joint 负责展平、策略层负责求解”，并按非 `FixedMu` 懒触发策略运行。
  2) 代码收口：
     - 删除 `ProblemSpec._fixedrho_problem_spec_forward_solve(...)` 中 `fixedrho_joint_solve=false` 时回退 `_solve_constraint_fixedrho(...)` 的分支；
     - 明确将 `fixedrho_joint_solve=false` 视为不支持参数，抛出 `ArgumentError`；
     - `FixedRho` 的 ProblemSpec forward_solve 统一走 joint + lazy attempts 主链。
  3) 测试同步：
     - `tests/unit/models/test_problem_spec_contract.jl` 中将 `fixedrho_joint_solve=false` 口径更新为抛错断言；
     - 保留并复用既有 continuity-like/joint parity gates 回归集。
  4) 验证通过：
     - `julia --project=. -e 'include("tests/regression/models/test_problem_spec_fixedrho_parity_regression.jl")'`
     - `julia --project=. -e 'include("tests/unit/models/test_problem_spec_contract.jl")'`
     - `julia --project=. -e 'include("tests/unit/models/test_solver.jl")'`（25/25）
     - `julia --project=. -e 'include("tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl")'`（43/43）。
  5) 阶段结论：
     - `FixedRho` 在 ProblemSpec 层已不再保留 legacy 兼容分支；
     - 仍需后续轮次继续处理 `constraint_solver` 其他 non-FixedMu 模式中的 legacy fallback 依赖，方可推进 Step4 全量退役。

- 2026-04-05：FixedRho 全链路去 legacy + Entropy/Sigma 懒触发策略推广（进行中）
  1) FixedRho 全链路收口：
     - `constraint_solver.jl` 中 `_solve_constraint_fixedrho(...)` 已移除对 `Main.Models.ImplicitSolver.solve(...)` 的 PNJL fallback；
     - 至此 `FixedRho` 在当前主链路径中不再依赖 `ImplicitSolver` 兜底。
  2) Entropy/Sigma 策略推广（沿用 FixedRho 成熟模式）：
     - `ProblemSpec._fixedentropy_problem_spec_forward_solve(...)` 与 `_fixedsigma_problem_spec_forward_solve(...)` 改为非 `FixedMu` 懒触发 attempts：
       - 先 primary 单次尝试；
       - 失败后触发 method-rescue / seed-rescue；
       - 命中通过 hard-constraints 的候选后 early-stop；
     - 默认方法改为 continuity-aware：
       - `continuity_seed=true` 默认 `:newton`；
       - 否则默认 `:trust_region`；
     - 保持 `FixedMu` 强制多 attempts 选优逻辑不变。
  3) 约束求解器配套改动：
     - `_solve_constraint_fixedentropy(...)`、`_solve_constraint_fixedsigma(...)` 新增 `nlsolve_method` 参数，并用于外层 nlsolve 方法选择，支持策略层方法切换。
  4) 本轮验证通过：
     - `julia --project=. -e 'include("tests/regression/models/test_problem_spec_fixedrho_parity_regression.jl")'`
     - `julia --project=. -e 'include("tests/unit/models/test_problem_spec_contract.jl")'`
     - `julia --project=. -e 'include("tests/unit/models/test_solver.jl")'`（25/25）
     - `julia --project=. -e 'include("tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl")'`（43/43）。
  5) 当前状态：
     - `FixedRho` 已进入“joint 主链 + 无 legacy fallback”状态并保持回归稳定；
     - `FixedEntropy/FixedSigma` 已完成“懒触发 attempts + continuity-aware 默认方法”推广；
     - 下一步按计划进入 `FixedAsymmetricRho` 同策略迁移与收口。

- 2026-04-05：FixedAsymmetricRho 策略推广（与 FixedRho/Entropy/Sigma 对齐）
  1) 目标：将 `FixedAsymmetricRho` 统一到非 `FixedMu` 懒触发 attempts 框架，保持 “joint/约束构造与策略层解耦”。
  2) 代码落地：
     - `ProblemSpec._fixedasymrho_problem_spec_forward_solve(...)` 已改为 staged/lazy attempts：
       - 先 primary 单次尝试；
       - 失败后触发 method-rescue / seed-rescue；
       - 命中 hard-constraints 后 early-stop；
     - 默认方法改为 continuity-aware：
       - `continuity_seed=true` 默认 `:newton`；
       - 否则默认 `:trust_region`；
     - `_solve_constraint_fixedasymrho(...)` 新增 `nlsolve_method` 参数并用于 outer nlsolve 方法选择。
  3) 验证通过：
     - `julia --project=. -e 'include("tests/unit/models/test_solver.jl")'`（25/25）
     - `julia --project=. -e 'include("tests/unit/models/test_problem_spec_contract.jl")'`（121/121）
     - `julia --project=. -e 'include("tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl")'`（43/43）
     - `julia --project=. -e 'include("tests/regression/models/test_problem_spec_fixedrho_parity_regression.jl")'`（fixedrho guard 全通过，防回归）。
   4) 阶段状态：
      - 非 `FixedMu` 三个模式（`FixedRho` / `FixedEntropy` / `FixedSigma` / `FixedAsymmetricRho`）已统一到“懒触发 attempts + continuity-aware 默认方法”的策略骨架；
      - `FixedMu` 继续保持“强制多 attempts 选优”作为唯一多解主模式。

- 2026-04-05：继续推广到其他求解模式（legacy fallback 治理化收口，非强删）
  1) 红灯验证（先证伪）：
     - 直接移除 `constraint_solver` 中 `FixedEntropy/FixedSigma/FixedAsymmetricRho` 的 `ImplicitSolver` fallback 后，回归立即失败：
       - `tests/unit/models/test_solver.jl` 失败（non-FixedMu parity 中 Entropy/Sigma/Asym convergence 断言红灯）；
       - `tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl` 失败（43 用例中 6 处红灯，目标量显著偏离）。
     - 结论：三模式当前仍存在真实稳定性依赖，不能直接移除 fallback。
  2) 低风险推进（保语义前提下继续“推广”）：
     - 在 `constraint_solver.jl` 为三模式新增显式治理参数：
       - `_solve_constraint_fixedentropy(...; allow_legacy_fallback::Bool=true, ...)`
       - `_solve_constraint_fixedsigma(...; allow_legacy_fallback::Bool=true, ...)`
       - `_solve_constraint_fixedasymrho(...; allow_legacy_fallback::Bool=true, ...)`
     - fallback 结果新增可观测元信息：`legacy_fallback_used::Bool`（true/false）。
     - 在 `ProblemSpec.jl` 将 `allow_legacy_fallback` 贯通为可配置参数并回传到 forward_solve 结果，形成统一观测口径。
  3) 测试补强：
     - `tests/unit/models/test_problem_spec_contract.jl` 新增 non-rho governed metadata 断言：
       - `haskey(solved, :legacy_fallback_used)`
       - `solved.legacy_fallback_used isa Bool`。
  4) 验证通过（本轮落地后的绿灯）：
     - `julia --project=. -e 'include("tests/unit/models/test_solver.jl")'`（25/25）
     - `julia --project=. -e 'include("tests/unit/models/test_problem_spec_contract.jl")'`（127/127）
     - `julia --project=. -e 'include("tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl")'`（43/43）
     - `julia --project=. -e 'include("tests/regression/models/test_problem_spec_fixedrho_parity_regression.jl")'`（fixedrho guards 全通过）
  5) 阶段结论：
      - 已将“是否启用 legacy fallback”从隐式行为升级为显式可治理开关，并提供可观测证据位；
      - 这为下一轮按模式逐步压缩 fallback 使用率提供了数据抓手；
      - 但截至本轮，`FixedEntropy/FixedSigma/FixedAsymmetricRho` 仍不具备无 fallback 的稳定默认语义，暂不进入 Step4 强删。

- 2026-04-05：按 FixedRho 架构继续推广（治理元信息对齐到非 FixedRho 模式）
  1) 先加红灯契约（TDD）：
     - 在 `tests/unit/models/test_problem_spec_contract.jl` 为 `FixedEntropy/FixedSigma/FixedAsymmetricRho` 增加治理元信息断言：
       - `governed_selected_method`
       - `governed_selected_quality`
       - `governed_fallback_used`
     - 并新增 `allow_legacy_fallback=false + nlsolve_method=:newton + trust_region_fallback=false` 场景，约束 `governed_fallback_used == false`。
  2) 绿灯实现（`ProblemSpec` 非 rho 三模式）：
     - 在每个 attempt candidate 中补齐并透传治理字段：
       - `governed_selected_method=attempt_cfg.method`
       - `governed_selected_quality`（`ok => :good`，否则按 solver 状态标记 `:degraded/:bad`）
       - `governed_fallback_used=(attempt_origin == :method_rescue)`
       - `governed_attempt_origin`
     - 在 forward_solve 最终返回中统一暴露上述三项，形成与 `FixedRho` 同风格的策略可观测面。
  3) 验证结果：
     - 通过：
       - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_problem_spec_contract.jl"; include("tests/unit/runtests.jl")'`（160/160）
       - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver.jl"; include("tests/unit/runtests.jl")'`（25/25）
       - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_problem_spec_semantic_modes.jl"; include("tests/unit/runtests.jl")'`（7/7）
       - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_constraint_solver.jl"; include("tests/unit/runtests.jl")'`（47/47）
       - `julia --project=. -e 'include("tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl")'`（43/43）
       - `julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_problem_spec_fixedrho_parity_regression.jl"; include("tests/regression/runtests.jl")'`（78/78）
     - 已知红灯（独立既有口径，非本轮新增）：
       - `tests/integration/models/test_solver_backend_semantic_parity_guard.jl` 仍在 `FixedRho` case 失败，表现为 models 与 legacy 在单点语义上不一致（`rho_norm/entropy/energy` 偏离）。
  4) 阶段结论：
     - 本轮已完成“按 FixedRho 架构推广”的第一层：非 FixedRho 三模式具备统一策略治理元信息与 fallback 可观测性；
     - 下一步应进入第二层：利用这些治理字段把 `legacy fallback` 从“兜底可用”推进到“可量化收缩”（按 mode/点位分层）。

- 2026-04-05：FixedRho 红灯快速修复（将 legacy 该点有效初解并入当前初解池）
  1) 背景：`tests/integration/models/test_solver_backend_semantic_parity_guard.jl` 在 `T=110MeV, rho=0.6` 的 FixedRho case 红灯；
     诊断显示 models 路径 9 个候选均落入同一坏盆地（`rho≈0.00129`），而 legacy 可收敛到目标分支。
  2) 实施：
     - 在 `ProblemSpec._fixedrho_problem_spec_forward_solve(...)` 中，除原有 primary/provided/default seeds 外，显式加入一条 legacy 语义 seed：
       - `legacy_seed = extend_seed(primary_seed[1:5], mode)`（若主 seed 已是 8 维也重建其 6:8 分量）；
       - 去重后加入 fallback seed pool。
     - 目的：把 legacy 在该类点位的有效 mode-aware 初解纳入当前主链候选池，先修收敛盆地覆盖问题。
  3) 测试补强：
     - `tests/integration/models/test_solver_backend_semantic_parity_guard.jl` 新增：
       - `FixedRho ProblemSpec fallback seed pool includes legacy-mode seed`，断言 `T=110MeV, rho=0.6` 下 `models` 路径可收敛并满足目标密度。
  4) 验证通过：
     - `julia --project=. -e 'include("tests/integration/models/test_solver_backend_semantic_parity_guard.jl")'`（39/39 + 新 testset 2/2）
     - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_problem_spec_contract.jl"; include("tests/unit/runtests.jl")'`（160/160）
     - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver.jl"; include("tests/unit/runtests.jl")'`（25/25）
     - `julia --project=. -e 'include("tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl")'`（43/43）
     - `julia --project=. -e 'ENV["REGRESSION_FILES"]="models/test_problem_spec_fixedrho_parity_regression.jl"; include("tests/regression/runtests.jl")'`（78/78）
   5) 结论：
      - 该红灯已由“吸引盆地覆盖不足”被修复；
      - 该改动不改变求解器主干结构，仅增加一条 mode-aware 候选 seed，属于低风险稳态增强。

- 2026-04-05："额外约束作为通用输入" 最小改动实施（阶段 1：护栏 + 壳层）
  1) 目标与策略：
     - 按“可回退、可验证、单点替换”推进，不改现有模式求解行为；
     - 先在 `ProblemSpec` 层引入通用 `ExtraConstraints` 壳层，默认实现为 no-op，确保零行为变化。
  2) TDD 护栏（先红后绿）：
     - 在 `tests/unit/models/test_problem_spec_contract.jl` 新增 `problem spec extra constraints shell contract`：
       - 要求 `Models.ExtraConstraints`、`Models.default_extra_constraints` 可见；
       - 要求 `build_problem_spec(mode)` 返回的 `spec` 带 `extra_constraints` 字段；
       - 要求默认 `seed_extend` 为 identity、`feasible` 恒真、`residual!` no-op。
     - 先运行该测试出现预期红灯（符号/字段不存在），再实现生产代码并转绿。
  3) 生产实现（零行为变化）：
     - `src/models/solver/ProblemSpec.jl`：
       - 新增 `ExtraConstraints` 结构（`residual!`, `feasible`, `seed_extend`）；
       - 新增默认 no-op 适配 `default_extra_constraints()`；
       - `ProblemSpec` 增加 `extra_constraints` 字段；
       - `ProblemSpec(...)` 构造器新增关键字 `extra_constraints=default_extra_constraints()`。
     - `src/models/Models.jl`：导出 `ExtraConstraints` 与 `default_extra_constraints`。
  4) 验证：
     - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_problem_spec_contract.jl"; include("tests/unit/runtests.jl")'`（167/167）。
  5) 阶段结论：
     - 已完成你要求的最小接缝注入："模式差异可转为 extra constraints 输入" 的结构基础到位；
     - 当前默认实现不改变任何既有解算路径行为，后续可在阶段 2/3 逐步接入主执行器钩子与两模式迁移。

- 2026-04-05："额外约束作为通用输入" 阶段 2（主执行器钩子接入，保持默认行为）
  1) 目标：
     - 在不改变默认数值行为的前提下，让 `extra_constraints` 真正进入 `ProblemSpec` 主执行路径；
     - 接入点仅限两类：`seed_extend` 与 `feasible`（`residual!` 先保留接口，不在本阶段强接）。
  2) 实施改动：
     - `src/models/solver/ProblemSpec.jl`
       - 新增工具函数：
         - `_resolve_extra_constraints(kwargs)`
         - `_extend_seed_with_extra(seed, mode, extra)`
         - `_append_extra_feasible_rule(rules, extra, mode)`
       - `ProblemSpec(...)` 构造器保持兼容，同时将 `extra_constraints` 注入默认 `forward_solve` keyword；
       - `_strip_problemspec_forwardsolve_kwargs!` 增加 `:extra_constraints`，防止下游 solver 收到无关键字；
       - 在 `FixedRho/FixedEntropy/FixedSigma/FixedAsymmetricRho` 的 `forward_solve` 中统一接入：
         - `primary_seed/provided_seed_pool/default_seed_pool` 均经过 `extra.seed_extend(...)`；
         - `hard_constraints` 追加 `extra.feasible(...) => :extra_constraint_failed` 规则。
     - `tests/unit/models/test_problem_spec_contract.jl`
       - 新增 `problem spec extra constraints hooks are wired in forward_solve`：
         - 通过 `extra_constraints=ec` 显式传入，断言 `seed_extend` 被调用；
         - 断言 `feasible=false` 会使 `hard_constraint_ok=false` 且失败原因包含 `:extra_constraint_failed`。
  3) 验证结果：
     - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_problem_spec_contract.jl"; include("tests/unit/runtests.jl")'`（170/170）
     - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver.jl"; include("tests/unit/runtests.jl")'`（25/25）
     - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_problem_spec_semantic_modes.jl"; include("tests/unit/runtests.jl")'`（7/7）
     - `julia --project=. -e 'include("tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl")'`（43/43）
  4) 阶段结论：
     - 已完成“通用额外约束输入”在主执行器的最小可用接入：可影响种子扩展与可行域判定；
     - 默认 `default_extra_constraints()` 仍为 no-op，故对既有调用方保持行为兼容；
     - 下阶段可迁移 `FixedEntropy/FixedSigma` 为“仅差异化提供 extra_constraints”以进一步消解 mode 分叉实现。

- 2026-04-05："额外约束作为通用输入" 阶段 3（先迁移 FixedEntropy + FixedSigma）
  1) 目标：
     - 将 `FixedEntropy` / `FixedSigma` 的 ProblemSpec 前向求解实现收敛到统一非 rho 主链执行器；
     - 保持行为等价（只抽公共编排，不变更物理求解内核）。
  2) TDD 与重构策略：
     - 先沿用阶段 1/2 已建立的 `problem_spec_contract + solver + semantic_modes + integration smoke` 作为回归护栏；
     - 在 `ProblemSpec` 中新增统一执行器 `_governed_nonrho_problem_spec_forward_solve(...)`，将重复逻辑（seed pool、attempt plan、hard rules、selector、治理元数据）抽象；
     - `FixedEntropy/FixedSigma` 改为仅提供：
       - mode label（报错文本）；
       - mode 特有 attempt-origin 字段名；
       - mode 特有 constraint solve 闭包（分别调用 `_solve_constraint_fixedentropy/_fixedsigma`）。
  3) 关键实现点（`src/models/solver/ProblemSpec.jl`）：
     - 新增 `_governed_nonrho_problem_spec_forward_solve(...)` 通用函数；
     - `FixedEntropy/FixedSigma` wrapper 显著瘦身，迁移到“参数化调用通用执行器”；
     - 保留原有输出契约字段（`governed_*`, `legacy_fallback_used`, `selection_reason`, `candidate_count` 等）不变。
  4) 验证：
     - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_problem_spec_contract.jl"; include("tests/unit/runtests.jl")'`（170/170）
     - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver.jl"; include("tests/unit/runtests.jl")'`（25/25）
     - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_problem_spec_semantic_modes.jl"; include("tests/unit/runtests.jl")'`（7/7）
     - `julia --project=. -e 'include("tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl")'`（43/43）
  5) 阶段结论：
     - 已完成你指定的“先两模式迁移”：`FixedEntropy` 与 `FixedSigma` 的分叉编排明显收敛；
     - 当前模式差异主要集中在“额外约束/模式参数输入 + mode 专有内核调用”，符合“把额外约束当作通用输入”的主方向；
     - 下一步可按同样方式迁移 `FixedAsymmetricRho`，再评估 `FixedRho` 外层编排统一的收益/风险比。

- 2026-04-05："额外约束作为通用输入" 阶段 4（迁移 FixedAsymmetricRho）
  1) 目标：
     - 将 `FixedAsymmetricRho` 迁移到与 `FixedEntropy/FixedSigma` 相同的统一非 rho 主执行器，进一步压缩 mode 分叉编排。
  2) TDD 护栏：
     - 在 `tests/unit/models/test_problem_spec_contract.jl` 新增
       `fixedasymrho forward_solve accepts extra_constraints hook`：
       - 断言 `seed_extend` 被调用；
       - 断言 `feasible=false` 会触发 `:extra_constraint_failed`。
  3) 生产实现：
     - `src/models/solver/ProblemSpec.jl`
       - `_fixedasymrho_problem_spec_forward_solve(...)` 改为参数化调用
         `_governed_nonrho_problem_spec_forward_solve(...)`；
       - 仅保留 mode 差异输入（label / attempt-origin key / mode-specific solve closure）。
  4) 验证结果：
     - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_problem_spec_contract.jl"; include("tests/unit/runtests.jl")'`（173/173）
     - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_solver.jl"; include("tests/unit/runtests.jl")'`（25/25）
     - `julia --project=. -e 'ENV["UNIT_FILES"]="models/test_problem_spec_semantic_modes.jl"; include("tests/unit/runtests.jl")'`（7/7）
     - `julia --project=. -e 'include("tests/integration/pnjl/test_solver_constraints_models_backend_smoke.jl")'`（43/43）
  5) 阶段结论：
     - 非 FixedRho 三模式已全部完成“统一主执行器 + 额外约束输入”收敛；
     - 后续可进入 `FixedRho` 外层编排统一评估（风险更高，建议单独 gate 与回归矩阵）。
