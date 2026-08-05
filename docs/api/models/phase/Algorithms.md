# 相图算法核心

本文档面向需要理解相图主流程判据、数据流与算法边界的维护者，不把底层实现细节等同于公开稳定入口。

## 总体数据流

相图主流程可以概括为：

1. 先由 `run_trho_scan` 生成按温度切片组织的 `μ(ρ)` 曲线
2. 对每条曲线做 `detect_s_shape`
3. 对存在 S-shape 的曲线做 `maxwell_construction`
4. 跨温度切片执行 `find_cep`
5. 若启用 crossover 计算，则额外执行 `detect_crossover` 或 `scan_crossover_line`
6. 将结果封装为 `PhasePipelineResult`
7. 通过 `build_phase_artifacts` 写出工件，并可选调用 `promote_phase_artifacts` 晋升

上层编排主要位于 [src/models/phase/PhasePipeline.jl](src/models/phase/PhasePipeline.jl) 与 [src/models/phase/ProductionPhasePipeline.jl](src/models/phase/ProductionPhasePipeline.jl)，结果类型定义位于 [src/models/phase/PhaseTypes.jl](src/models/phase/PhaseTypes.jl)。

## 双入口分工

当前相图主题有两条公开主流程：

- `run_phase_pipeline`
  - 通用 / research 入口
  - 保留插值 CEP、direct CEP、crossover 和策略比较等较灵活的研究向路径
- `run_production_phase_pipeline`
  - production / baseline 入口
  - 使用显式温度区间 `T_start:T_end`、固定 `dT_initial`、unknown budget 和非插值 Maxwell 收口逻辑

两条入口共享底层 `run_trho_scan`、S-shape、Maxwell 与 artifact 治理能力，但在温度扫掠编排和 CEP 收口策略上职责不同。

## 一阶相变判据

### S-shape 检测

`detect_s_shape` 对单条 `μ(ρ)` 曲线检查导数符号是否出现 `正 → 负 → 正` 的变化序列，用于判定一阶相变区域是否存在。

它负责提供：

- 是否存在 S-shape
- 强子侧与夸克侧 spinodal 位置
- 导数符号变化的诊断信息

该步骤是 Maxwell 构造与 CEP 搜索的前置条件。细节见 [PhaseTransition.md](docs/api/models/phase/PhaseTransition.md)。

### Maxwell 等面积构造

`maxwell_construction` 在检测到 S-shape 后，搜索满足等面积条件的共存化学势与两侧共存密度。

它输出的是一阶相变边界的核心几何量：

- `mu_transition`
- `rho_hadron`
- `rho_quark`
- `area_residual`
- `converged`

在 pipeline 中，这一步直接决定 `first_order_boundary` 与 `spinodal` 两类工件的内容质量。

#### Production 容差合同

`run_production_phase_pipeline` 不把 Maxwell 二分的停止条件与跨层几何证书分开硬编码。
它从当前生效的 acceptance tolerances 派生内部 solver tolerance：

`tol_solver = 0.1 × minimum(active acceptance tolerances)`。

`area_tol_good` 始终生效；启用 rho geometry 或 temperature geometry 时，分别加入
`rho_maxwell_area_tol` 或 `temperature_maxwell_area_tol`。派生值写入
`config_snapshot`、`config_hash` 和 diagnostics，便于复现。`maxwell_construction` 的公共
默认值仍保持兼容；production 调用显式传递派生值，并复用同一次 Maxwell 结果填充分类和
共存密度，避免同一曲线的二次求根。外层 `PhaseGeometryTolerances.maxwell_area` 仍然保留，
它负责 coarse/fine 网格之间的离散收敛证书，而不是重复替代内部求根停止条件。

公共核心还要求等面积候选具有唯一的三交点拓扑。endpoint-limit hybrid 只在该候选唯一、
右外支已采样且左外交点落在首个正密度单元时补 `rho=0.003125` anchor 和最多 12 个
局部二分点；失败时保留 `ambiguous_near_critical`，不会把 endpoint 诊断降级成 monotone。

## CEP 搜索

`find_cep` 基于温度切片上的 Maxwell 有效性变化来定位 CEP 证据边界，并输出 [src/models/phase/PhaseTypes.jl](src/models/phase/PhaseTypes.jl#L1) 中的 `CEPResult`。它不再把 ambiguous midpoint 当作 CEP 单点。

当前应重点理解的不是函数签名，而是两类策略边界：

- `:interpolate`
  - 更接近历史插值路径
  - 更依赖相邻温度切片具备可比较或可映射的 `rho` 网格
- `:direct`
  - 直接在温度轴上评估 Maxwell 有效性
  - 不强制要求统一 `rho` 网格
  - 更适合与自适应 `rho` 加密联合使用

`CEPResult` 的关键诊断字段包括：

公开常量 `Models.CEP_RESULT_STATUSES` 列出允许的结果状态：
`:resolved`、`:ambiguous`、`:not_found`。

- `result_status`：`:resolved`、`:ambiguous` 或 `:not_found`
- `T_cep_MeV`
- `mu_cep_MeV`：历史兼容字段，表示夸克化学势 `mu_q`，不是重子化学势 `mu_B`
- `T_bracket_low_MeV` / `T_bracket_high_MeV`
- `bracket_width_T_MeV`
- `T_last_first_order_MeV` / `mu_last_first_order_MeV`
- `T_first_monotone_MeV`
- `ambiguity_width_T_MeV`
- `temperature_resolution_target_MeV`
- `eval_count`
- `unknown_count`
- `uncertainty_T_MeV`
- `reason`
- `method`

相图模块内部的 `mu_MeV` / `mu_transition_MeV` / `mu_cep_MeV` 口径均为 `mu_q`。面向 `mu_B` 的外层扫描或输运脚本需要显式使用 `mu_B = 3mu_q`；新生成的 CEP reference CSV 与 phase artifact 会同时给出 `muq_cep_MeV` 和 `muB_cep_MeV` 以避免歧义。这些字段主要服务于回归治理、失败解释与策略比较，而不只是展示一个 `T_cep_MeV` 数值。

### production 收口路径

`run_production_phase_pipeline` 不把 CEP 搜索抽象成单次“对已有曲线集合调用 `find_cep`”，而是显式执行温度扫掠与 bracket 收缩：

1. 根据 `ProductionPipelineConfig` 生成温度网格
2. 对每个温度切片执行 Maxwell 与 rho-geometry 分类，状态为 `confirmed_first_order`、`confirmed_monotone` 或 `ambiguous_near_critical`
3. 将扫掠结果收集到 `FirstOrderSweepResult`
4. 从最后一个 `confirmed_first_order` 与随后的 `confirmed_monotone` 温度切片构造证据区间
5. 独立收缩低、高两端；ambiguous 温度不被强制归入任一侧

只有独立的同点临界估计产生 `(T, μ)` 时，`result_status=:resolved` 且旧单点字段才有限。
普通 bracket 路径返回 `result_status=:ambiguous`，单点字段为空，
`[T_last_first_order_MeV, T_first_monotone_MeV]` 是实际 ambiguity interval。
`temperature_resolution_target_MeV` 只是端点数值搜索目标，不是物理误差上限。
CLI 与 dense-reference builder 仍接受旧的 `cep_tol`/`--cep-tol`，同时提供
`temperature_resolution_target_MeV`（或 `--temperature-resolution-target`）明确别名；
两者同时给出时以后者为准。

`confirmed_monotone` 至少需要两层 rho 网格都稳定报告 `no_s_shape`；单层无 S 形、Maxwell 失败、面积灰区、geometry 未收敛和 solver/no-curve 失败均保持 ambiguous。`unknown_budget` 只影响停止和诊断，不改变切片物理标签。
其中 `unknown_count` 只统计原始 solver/classification `:unknown`，不把
`weak-S`、Maxwell 灰区或其他 `ambiguous_near_critical` 计入预算。预算超限后
停止继续的温度/端点 refinement，并在 `reason` 与 diagnostics 中记录
`unknown_budget_exceeded`；已取得的证据区间原样保留。

显式设置 `rho_geometry_convergence=false` 时表示调用方选择旧的单层诊断/兼容口径：
通过 Maxwell 的一阶切片仍可作为边界候选输出，以保持既有边界消费者可用；但单层
`no_s_shape` 永远不会被提升为 `confirmed_monotone`，因此该配置不会伪造高侧 CEP
证据。默认 production 配置保持严格的两层 rho-geometry gate。

这条路径的重点不是更“通用”，而是更适合 production/baseline 口径的可解释性与稳定收口。

## Phase 网格几何量收敛

production 路径的 `rho` refinement 不再只在分类为 `unknown` 时触发。启用
`rho_geometry_convergence=true` 后，每个温度切片至少比较一层嵌套粗细 `rho` 网格，并同时约束：

- Maxwell 共存化学势与两侧 spinodal 化学势的位置误差（MeV）；
- 共存密度与 spinodal 密度误差；
- 粗细层的 Maxwell area residual。

达到最大细化层仍未满足几何量门限的 `valid` 点会降为 `unknown`，不能仅凭粗网格上出现
S-shape 就作为正式一阶相变点。温度自适应只处理两端均为 `valid` 的一阶相线区间：求解中点，
比较中点结果与端点线性插值，并按同一组位置量、密度量和面积残差门限继续细化。
`valid/invalid` 转换区由 CEP bracket 负责，不用相线插值代替 CEP 分类。

跨 `xi` 的 dense reference 使用同样的中点插值误差定义；每个待审计区间必须先计算一次中点，
因此即使相线近似线性，启用 adaptive xi 后也至少会把锚点间距减半一层。所有 rho/T/xi 诊断记录写入
`phase_grid_convergence.csv`（dense reference 为带 tag 的同名 sibling），诊断未收敛本身不自动等价于物理异常。

Action 生产面把该递归拆成分层的一 xi job：锚点与第一层中点并行，后续层只计算上一层未收敛区间的中点。
层间评估仍调用同一个 `_phase_result_midpoint_error`，不是在 workflow 中另定义一套误差公式。为使层间重建保留
Maxwell 面积门限，新的 dense boundary CSV 以附加列记录 `area_residual` 与 `converged`；旧 reference 不回填、
也不因该附加合同被覆盖。各单 xi shard 导出的 rho/T convergence record 同时显式写入所属 `xi`，避免跨 shard
确定性 merge 时把不同物理切片误判为同一诊断键。

每个温度的全部粗细层原始扫描保留在 `production_eval/`；聚合 `trho_scan.csv` 只写该温度最终采用的
一层，避免同一 `(T,rho)` 因审计层级不同而出现无标识重复行。

## Crossover 支路

当温度高于 CEP 或未形成一阶相变时，pipeline 可转而计算 crossover line。

相图主题中主要涉及两个能力：

- `detect_crossover`
- `scan_crossover_line`

现有实现支持至少两类判据：

- `:peak`
- `:inflection`

它们的职责是从序参量导数结构中给出 crossover 温度，而不是替代一阶相变的 Maxwell 判据。细节见 [Crossover.md](docs/api/models/phase/Crossover.md)。

crossover 温区上限由 `crossover_T_max_MeV` 显式控制；`NaN` 表示继承主 phase 的 `T_end` 或
`T_grid` 最大值。实现不再隐藏截断到 `220 MeV`。crossover 基态和导数路径使用调用方实际传入的
`p_num` 与 `t_num` 口径。phase 主扫描的 `iterations` 连同实际 `p_num/t_num` 一起进入配置快照和 config hash；
`crossover_mu0_only=true` 时实际采样只包含 `mu_q=0`。本次不新增 crossover 的独立 gap-iteration 参数。

## 自适应 rho 加密

自适应 `rho` 加密不是独立的主流程，而是为了提高 CEP 与 Maxwell 稳定性而存在的数据层辅助。

相关核心接口：

- `AdaptiveRhoConfig`
- `suggest_refinement_points`
- `merge_rho_values`

它们的职责是：

- 检测 `|Δμ/Δρ|` 过小、需要补样的区段
- 生成新的 `rho` 候选点
- 与现有网格做受控合并与去重

更详细的使用背景见 [AdaptiveRhoRefinement.md](docs/api/models/phase/AdaptiveRhoRefinement.md)。

### rho-support cascade（显式 opt-in）

`Models.RhoSupportConfig`、`Models.RhoSupportPrior`、`Models.RhoSupportAssessment` 和
`Models.analyze_rho_support_cascade` 是补点路由层；
local cubic 只用于定位 support window，最终 S-shape 仍必须由真实采样点形成
稳定的 `+ -> - -> +` topology。production cascade 使用两个嵌套 rho 层：
`rho_support_cascade` 默认 `0.05 -> 0.025`，每个温度累计 targeted 点上限为 `12`。
两层均通过 Maxwell 与 C2 geometry gate 才输出 `confirmed_first_order`，两层均
稳定 `no_s_shape` 才输出 `confirmed_monotone`，solver failure、Maxwell/area/geometry
不收敛统一保留 `ambiguous_near_critical`。

补点 session 是请求作用域的 exact `(T,xi,rho)` Float64 cache。已有成功点按同一
`(T,xi)` 下最近 rho 作为首选 seed；等距时遵循扫描方向。失败点也缓存，后续层级
不会静默重复求解。`SolverWorkTelemetry` 只在调用方显式提供时累加，不改变
`SolverResult v1` 或默认 uniform path。

### rho-support hybrid（五级验证链）

`rho_support_hybrid` 在 Stage A 的 cascade 结果之外复用同一请求作用域 cache，执行
完整 Stage B `0.0125 -> 0.00625` dense 层。若 Stage A/B 的状态或 geometry 未闭合，
则从 Stage-B 曲线的两个 μ 极值向外扫描，取首个满足 `mu < mu_low` 与 `mu > mu_high`
的已有采样点作为离散 guard；相等点跳过，不插值、不二分、不使用固定 padding。
Stage C 的分类曲线是完整 Stage-B 全域曲线与 guard 内按 Stage-B 曲率、面积贡献和
Maxwell 特征排序选出的 `0.003125` 点的并集。缺少任一 guard、多 S topology 或
局部点不足时保留 `ambiguous_near_critical`；Stage C 不产生新的 monotone 证书。
`RhoHybridVerificationConfig`、manifest 和 diagnostics 会记录 guard rule、比较 epsilon、
local step、target cap、support μ 极值、采用阶段和 verification 点数。它还固定记录
`candidate_policy=:unique_three_crossing_topology_v1` 与
`endpoint_policy=:bounded_zero_density_v1`；后者只在唯一 endpoint-dependent 候选且右外点
存在时启用历史零密度端点证书。新的显式
`endpoint_policy=:three_crossing_endpoint_local_v2` 不再要求右支越过两个 μ 极值：它要求
右 Maxwell 交点被两个实际 Stage-B 外支采样点 bracket，并始终保留完整 Stage-B 曲线；
Stage-C 只在左 Maxwell 交点的 active bracket 内加入 midpoint。左 bracket 下界保持为零并
收缩到端点预算时，内部诊断证书为 `endpoint_limited_first_order`；下界变为正值且最后两级
geometry 均通过时为 `endpoint_local_geometry_first_order`。这两个名称只描述内部证书，
对外仍映射到既有三态 `confirmed_first_order` / `confirmed_monotone` /
`ambiguous_near_critical`，不新增用户可调容差。

## 工件治理层

相图主题不止是算法判据，还包括工件与基线治理：

- `build_phase_artifacts` 负责把数值结果写成 CSV、JSON、Markdown
- `resolve_phase_output_target` 负责决定 processed/reference/outputs 目标路径
- `promote_phase_artifacts` 负责 gate 校验与 reference 晋升

因此，相图 pipeline 的真正闭环是“算法结果 + 结构化产物 + 晋升治理”，而不是只停留在相变点数值本身。

production 入口同样复用这套工件治理，但会在 `diagnostics` 和 `config_snapshot` 中附加 production 特有字段，例如：

- `mode`
- `dT_initial`
- `unknown_budget`
- `first_point_fallback`
- `forced_invalid_count`（兼容字段；三态路径不再把 unknown 强制改写为 invalid，正式产物中应为 `0`）
- `sweep_temperatures_MeV`
- `sweep_statuses`
- `rho_geometry_convergence` 及其位置量、密度量、面积残差门限
- `adaptive_temperature` 及其中点细化门限
- `crossover_mu0_only`
- `crossover_T_max_MeV`

## 当前边界

- 本页解释的是算法关系与职责分工，不承诺所有底层函数都作为长期稳定公开入口
- 面向外部调用时，仍应优先回到 [Overview.md](docs/api/models/phase/Overview.md) 中列出的 `Models.*` 入口
- 旧路径算法页已降级为兼容跳转说明；当前新目录页已成为算法细节的主版本
