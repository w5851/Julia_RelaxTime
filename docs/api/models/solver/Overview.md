# Models 核心求解入口总览

本页优先回答“我应该调用哪个入口”，而不是先展开底层实现细节。

## 最小稳定入口

对大多数新调用方，核心入口顺序如下：

1. `Models.create_model`：先选择模型对象
2. `Models.solve_gap`：直接求平衡态平均场状态
3. `Models.solve`：在 PNJL 约束模式下求完整 `SolverResult`
4. `Models.solve_multi`：当你明确需要多初值筛选物理解时使用
5. `Models.coerce_solver_result` / `Models.solver_result_view`：当你需要稳定序列化视图或版本化契约检查时使用

Phase 3 收口说明：`fixedmu_problem_spec_active` 兼容字段与 `fixedmu_use_problem_spec` 开关均已移除，`solve_constraint` 统一走 ProblemSpec 主链。

## 什么时候用 `create_model`

当你明确在 `Models` 多重派发体系下工作时，先构造模型对象：

```julia
model = Models.create_model(:PNJL)
st = Models.solve_gap(model, T_fm, mu_vec)
```

适合场景：

- 你需要直接调用 `omega`、`omega_components`、`model_thermo` 等模型方法
- 你要在 NJL / PNJL / RPNJL 之间切换实现
- 你希望后续工作流沿模型对象接口继续展开

## 什么时候用 `solve_gap`

`Models.solve_gap` 适合“我只想拿到平衡态平均场状态”的场景。

典型输出是 `MeanFieldState` 风格的状态对象，可继续用于：

- `omega` / `omega_components`
- `state_vector`
- 自定义后处理或 TD-first 导数链路

它是最轻量的统一入口，不承担 CSV、扫描或完整 workflow 打包职责。

NJL/PNJL 家族的 `solve_gap(...; iterations=1000)` 将正整数迭代预算传给每次
NLsolve 尝试。FixedMu 的 `solve` / `solve_constraint` 同样透传 `iterations`，
包括 Newton 失败后的 trust-region 尝试；预算不是所有种子的累计上限。
耗尽预算不能绕过残差和物理性判据。未指定时仍为 NLsolve 原默认的 1000；
输运扫描显式使用每次 40（底层 CLI 可用 `--max-iter` 覆盖）。
FixedRho 等其他约束模式继续使用各自既有的联合求解合同。

## 什么时候用 `solve`

`Models.solve` 适合“我已经确定约束模式，希望拿到完整的求解结果载荷”的场景。

典型形式：

```julia
res = Models.solve(Models.FixedMu(), T_fm, mu_fm)
```

适合场景：

- 你要同时拿到 `x_state`、`mu_vec`、`omega`、`pressure`、`rho_norm`、`entropy`、`energy`、`masses`
- 你在固定密度、固定熵等约束模式下工作
- 你需要一个标准化 `SolverResult` 而不是仅状态对象

## 什么时候用 `solve_multi`

`Models.solve_multi` 适合一阶相变附近或多值解区域，需要显式尝试多组初值并按物理判据选解的场景。

模式支持说明：`solve_multi` 当前支持 `FixedMu | FixedRho | FixedAsymmetricRho | FixedMuBConservedCharges`，不覆盖 `FixedEntropy` / `FixedSigma`。

它不是常规入口的默认替代品，而是“我明确知道这里有分支竞争”的增强入口。

## 已实现的约束模式与路径支持

PNJL 约束接口当前有六种模式；下表的未知量数包含五个平均场和需要联立求出的化学势。

| 模式 | 固定条件 | 未知量数 | 通用 PALC 路径 |
| --- | --- | ---: | --- |
| `FixedMu` | T 与共同夸克化学势 μq | 5 | 已接入，固定 T 扫 μq |
| `FixedRho` | T、重子密度，三味化学势相等 | 8 | 已接入，固定 T 扫 ρ |
| `FixedAsymmetricRho` | T、重子密度、u/d 净密度比与 s 味净密度 | 8 | 已接入，固定 T 扫 ρ |
| `FixedMuBConservedCharges` | T、μB、电荷/重子比、净奇异数密度 | 8 | 已接入，固定 mode 中约束沿 T 延拓 |
| `FixedEntropy` | T、熵密度 s，三味化学势相等 | 8 | 已接入，沿 T 延拓；PALC 为显式选项 |
| `FixedSigma` | T、比熵 σ=s/nB，三味化学势相等 | 8 | 已接入，固定比熵沿 T 延拓 |

`FixedEntropy` 固定熵密度；`FixedSigma` 才是保持单位重子熵的等熵轨迹，σ 与 σ 介子无关。
三种约束使用 `ConstraintTPath(mode, T_values)`，温度单位为 fm^-1，继续求五个平均场与
三个化学势的八变量联立方程。标准 PNJL tensor 热积分路径中，`FixedSigma` 和
`FixedMuBConservedCharges` 的 `continuation=:auto` 在至少 21 点时使用 PALC；
`FixedEntropy` 保留 legacy 默认，原因是已测窗口存在等压力的正负化学势分支，应用的
选支语义仍需确定。可用 `continuation=:palc` 显式选择。

`run_isentropic_meson_mass_scan` 的标准长分段已复用 `FixedSigma` PALC 平衡态，直接交给
介子后处理；短路径、RPNJL 和自定义控制保留 legacy 默认。三种模式的独立小规模对照见
[约束路径验证](../../../analysis/pnjl/algorithmic_feasibility/palc_matched_20261108_analysis/README.md)。
`ConstraintTPath` 接受 `:auto|:legacy|:palc`；`single_pass` 是 Tμ/Tρ 扫描的显式选项。

Tμ/Tρ 扫描的 `:legacy` 逐点复用前一点与多初值搜索；`:single_pass` 将每点候选初值展开，
减少重复求解；`:palc` 从锚点双向延拓并在目标点重新求根。
标准 PNJL 多点路径默认使用 `EndpointAnchor()`：从较高参数端取得首个合格根，
温度路径因此从最高温度起步；等温路径从较高 μ 或 ρ 起步。正常路径只建立一个锚点。
首个初猜失败才尝试互补初猜和另一端；延拓仍有缺点时恢复目录搜索。
`SeedPoolAnchor()` 可显式保留中点与两端的完整目录搜索。单点路径、RPNJL 和
`FixedEntropy` 的 PALC 继续默认使用目录策略。
`SingleSeedAnchor(value, seed)` 配合 `PALCContinuation()` 使用同一延拓算法，但只从一个
初值建立锚点并关闭 catalog fallback，适合研究已知解支；它不是第四个 `continuation` 符号。
Tρ 的 `:natural` / `NaturalContinuation()` 使用一次从高密到低密的普通参数延拓，
每点保持完整八变量联立求解。失败时恢复单锚点 PALC；此策略适用于密度路径。
标准 PNJL 的较长 Tρ 扫描默认选择 natural，Tμ 继续默认选择 PALC。

磁场扫描和 NJL/RPNJL 等属于模型或外场变体。`run_magnetic_scan` 的平衡态使用
`FixedMu`，其 Landau 求和残差需要独立验证后才能接入 PALC。

## Models 边界签名与参数合同（Plan-A）

在 `Models` 边界，`solve*` 系列当前稳定公共表面为：

- `solve(model, mode, ...)` 与 `solve(mode, ...)`
- `solve_multi(model, mode, ...)` 与 `solve_multi(mode, ...)`
- `solve_constraint(model, mode, T_fm; ...)`
- `solve_vec(model, mode, theta_vec; ...)`
- `solve_named(model, mode, theta_named; ...)`

关键参数合同：

- `solve_constraint` 统一走 `ProblemSpec` 主链；`use_problem_spec` / `allow_legacy_path` / `warn_on_legacy_path` / `fixedmu_use_problem_spec` / `legacy_fallback_plugin` 已移除。
- `solve_vec`：`FixedMu` 约定 `theta_vec=[T_fm, μ_fm]`（长度 2）；其它模式约定 `theta_vec=[T_fm]`（长度 1）。`FixedMuBConservedCharges` 的 `muB_fm` 固定在 mode 中，不重复放入 `theta_vec`。
- `solve_named`：`FixedMu` 需要 `(:T_fm, :μ_fm)`；其它模式需要 `:T_fm`。

`solve_vec` / `solve_named` 在本页按签名与参数契约层覆盖；逐符号细节以自动生成导出索引为准（见 `docs/api/models/solver/generated/Exports.md`）。

## 用户入口分层建议

### 新使用者

- `create_model`
- `solve_gap`
- `solve`

### 进阶使用者

- `solve_multi`
- `solve_constraint`
- `MeanFieldState` / `state_vector` / `normalize_mu_vec`

### 维护者或算法开发者

- `ConstraintModes`
- `SeedStrategy` 家族
- `build_conditions` / `build_residual!`
- `solve_pnjl_with_derivatives` / `solve_pnjl_with_flavor_mu_derivatives`
- residual problem builders：`build_njl_problem` / `build_pnjl_fixedmu_problem` / `build_pnjl_flavor_mu_problem`

## 与其它主题的关系

- phase 主题、workflow 主题、scans 主题都建立在本主题定义的求解合同之上
- 如果你先理解本主题，再去看 [../workflows/README.md](../workflows/README.md) 或 [../scans/README.md](../scans/README.md)，会更容易看清哪些是“业务流程入口”，哪些是“核心求解公共表面”

## 迁移期契约说明（Plan-B / Gate A）

- `SolverResult` 已从固定 `SVector{5}/SVector{3}` 扩展为向量泛型字段，结果消费侧不应再假设状态维度恒为 5、化学势维度恒为 3。
- scan 结果适配层（`TmuScan` / `TrhoScan`）已按动态向量透传 `x_state` 与 `mu_vec`，为后续 schema-driven 主链路铺路。
- 当需要显式状态布局时，优先使用 `ModelStateSchema` 系列 API（`schema_for_model`、`flatten_state`、`unflatten_state`）。
- `SolverResult` 稳定契约版本字段为 `contract_version=:v1`；上层消费建议统一通过 `solver_result_view` 取公共视图。

## 请求级 solver telemetry

### 通用等温路径

`FixedMuPath(T_fm, mu_values)`、`FixedRhoPath(T_fm, rho_values)` 与
`FixedAsymmetricRhoPath(T_fm, rho_values; ud_ratio_target=0.876, s_target=0)`
接受非空向量；温度为正 fm⁻¹，μ 为 fm⁻¹，ρ 为非负 ρ/ρ₀。路径保留输入顺序和重复点。

```julia
model = Models.create_model(:PNJL)
path = Models.FixedRhoPath(120 / Constants_PNJL.ħc_MeV_fm, collect(0.2:0.1:3.5);
    xi=0.0, p_num=24, t_num=8)
result = Models.solve_path(model, path;
    continuation_strategy=Models.PALCContinuation(), residual_norm_max=1e-6)
points = Models.selected_path_results(model, result)
```

PALC 默认 `anchor_strategy=nothing` 按上述模型和路径选择锚点策略；可显式传入
`EndpointAnchor()`、`SeedPoolAnchor()` 或 `SingleSeedAnchor(value,seed)`。
`EndpointAnchor` 使用双向延拓的全部目标交点；`anchor_values` 仅用于目录策略。
诊断的 `anchor_strategy` 为 `:endpoint_anchor` 或实际触发恢复时的 `:endpoint_catalog_recovery`。
有限起点搜索需要通过对照检查分支覆盖。每个目标的全部合格根交给
`GroundStateBranchPolicy`；相图中的相共存仍由上层 Maxwell 等构造判断。
`selected_path_results` 返回按原路径排列的 `Union{Nothing,SolverResult}` 向量；失败为
`nothing`，成功结果维持 v1 合同。`solve_path` 默认仍为 `SeedContinuation()`。

PALC 支持 PNJL/RPNJL 的 tensor 热积分。路径热力学和残差复算使用相同参数。
通用扫描的 `continuation=:auto` 对 PNJL 的较长 Tρ 路径（至少 21 点）启用 natural，
Tμ 路径（至少 31 点）启用 PALC；RPNJL、短路径、自定义选择器及不支持的求解器选项保留
旧路线。单点 `Models.solve`
和 `solve_multi` 的默认算法独立于该扫描策略。
固定密度复用现有 joint residual 的配置模型；若传入不同于该配置的自定义模型参数，入口
会报 `ArgumentError`，以免残差和后处理使用不同参数。共同 μ 的路径直接使用传入模型。
可传入 `PALCWorkTelemetry()` 读取根求解、延拓残差/Jacobian、起点、校正、回退及失败计数；
对象复用时计数累加，多个路径内保存的累计快照不应再次相加。
其中 `root_time_s` 包含直接求根及候选后处理，`root_candidate_time_s` 是其中的后处理子项；
`continuation_time_s` 记录 BifurcationKit 延拓，`selected_rebuild_time_s` 记录选中结果重构。
后者需将同一个 `path_work` 传给 `selected_path_results`；这些计时不含序列化和独立验证。

`FixedAsymmetricRho` 可显式传入 `fixedasymrho_seed_policy=:single_pass`，并设置
`evaluate_all_attempts=true`。它展开既有外层和内层种子的并集，对每个不同起点进行一次
相同的八变量联立求解，将全部候选交给现有硬约束和 selector；默认 `:nested` 保留原路径。
该选项适用于 `solve`、`solve_multi`、`solve_constraint` 及 SeedContinuation 的 `solve_path`。
`solve_multi` 和默认 SeedContinuation 已评估全部候选；直接 `solve`/`solve_constraint`
需要显式设置 `evaluate_all_attempts=true`。未知策略或提前停止配置抛出 `ArgumentError`。
两策略可能因候选可见性不同而选中不同根，尤其在一阶相变温区；`semantic_mode` 和
自定义 selector 的含义不变。固定密度的压力最大规则与物理相共存判定须分别核对。
新策略全部候选未通过硬约束时返回 `converged=false`；单点 `solve_constraint`
的 `hard_constraint_ok`、`failed_constraints` 和诊断字段保留具体原因。

需要审计 `FixedRho` 或 `FixedAsymmetricRho` equilibrium 工作量时，可显式创建 `Models.SolverWorkTelemetry()`，并将其作为
`work_telemetry=` 传给 `Models.solve`、`Models.solve_multi`、`Models.solve_constraint` 或 `Models.solve_path`。计数器属于单次请求/扫描，
不是全局状态；`SolverResult v1` 的字段和序列化契约不变。`Models.solver_work_snapshot`
返回可写入 evidence 的 `NamedTuple`，包含 residual/Jacobian 调用、Newton/trust-region
迭代、fallback/rescue、异常和 retry 计数；长循环复用对象前可调用
`Models.reset_solver_work!`。未传入 telemetry 时保持原有路径。

`FixedAsymmetricRho` 的 `fixedasymrho_requests` 统计外部 equilibrium 请求，
`governed_attempts` 统计外层候选尝试，`fixedasymrho_inner_solves` 统计每次实际启动的内层
多种子 NLsolve（包括抛出异常的尝试）。`fixedasymrho_residual_evaluations` 统计该内层
残差函数的实际调用，包含 ForwardDiff 求 Jacobian 时的 Dual 调用；
`nlsolve_f_calls/g_calls` 则沿用 NLsolve 自报口径，两者不能直接相加。
后处理残差单独计入 `postprocess_residual_calls`；`nonconverged_attempts` 包含内层
物理/残差检查未通过和异常退出，最终请求是否成功仍由外层治理及返回结果确定。
底层计数累加器 `record_solver_request!`、`record_governed_attempt!`、
`record_nlsolve_work!`、`record_postprocess_residual!`、`record_attempt_outcome!`、
`record_solver_exception!` 和 `record_scan_retry!` 主要供隔离 pilot/测试使用，
不改变 `SolverResult` 字段。
- 诊断契约稳定公共视图与内部调试视图已分层；上层默认应消费 `to_public_namedtuple` / `coerce_solver_diagnostic_public_view`。

## 导数契约骨架（Phase-1 / Issue #79）

solver/diff 当前提供最小稳定骨架，先冻结“上下文 + 目标 + 参数 + Jacobian 形状”约定，再在后续阶段接入真实导数数值实现：

- `build_thermo_diff_context(result; mode, model, theta, spec_override=nothing)`
  - 构建 `ThermoDiffContext`；`theta` 采用 NamedTuple 口径（首批支持 `T_fm` / `mu_fm(兼容 μ_fm)` / `xi`，并归一化到 `mu_fm`）。
- `ParamSpec([:T_fm, :mu_fm, :xi])`
  - 约束外参导数输入维度与顺序；非法参数名会抛 `ArgumentError`。
- `diff_target(:pressure | :entropy | :rho_norm | :energy)`
  - 目标注册入口；未知目标名抛 `ArgumentError`。
- `jacobian(ctx, target, params)` / `jacobian(ctx, targets, params)`
  - 始终返回矩阵；单目标单参数固定 `1x1`，多目标单参数返回 `N x 1`。
  - 当 `ParamSpec` 包含多个参数时，要求 backend 返回 `N x P` 矩阵（`N=目标数`, `P=参数数`）。

Issue #80 已接入默认数值 Jacobian（当前主路径聚焦 `FixedMu`），并保持 `N x P` 形状契约与输入校验语义。

## 试点统一导数服务（Issue #81）

在不改动 `ProblemSpec` 主契约前提下，solver/diff 提供两条试点服务入口，供脚本链与分析链以统一方式接入导数：

- `build_pilot_diff_context(result; mode, model, theta, spec_override=nothing, jacobian_backend=nothing)`
  - 对 `build_thermo_diff_context` 的薄封装；统一 `mu_fm`/`μ_fm` 别名语义。
- `eval_pilot_derivatives(ctx; target_names, param_names)`
  - 统一执行 `DiffTarget + ParamSpec + jacobian`，并返回 `jacobian` 与 `by_name` 命名映射。

该层定位为“试点评估服务”，用于验证可维护性、稳定性与性能边界；是否回收进 `ProblemSpec` 由后续 Phase-2 决策确定。
