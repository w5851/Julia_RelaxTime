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

## CEP 搜索

`find_cep` 基于温度切片上的 Maxwell 有效性变化来估计 CEP，并输出 [src/models/phase/PhaseTypes.jl](src/models/phase/PhaseTypes.jl#L1) 中的 `CEPResult`。

当前应重点理解的不是函数签名，而是两类策略边界：

- `:interpolate`
  - 更接近历史插值路径
  - 更依赖相邻温度切片具备可比较或可映射的 `rho` 网格
- `:direct`
  - 直接在温度轴上评估 Maxwell 有效性
  - 不强制要求统一 `rho` 网格
  - 更适合与自适应 `rho` 加密联合使用

`CEPResult` 的关键诊断字段包括：

- `T_cep_MeV`
- `mu_cep_MeV`：历史兼容字段，表示夸克化学势 `mu_q`，不是重子化学势 `mu_B`
- `T_bracket_low_MeV` / `T_bracket_high_MeV`
- `bracket_width_T_MeV`
- `eval_count`
- `unknown_count`
- `uncertainty_T_MeV`
- `reason`
- `method`

相图模块内部的 `mu_MeV` / `mu_transition_MeV` / `mu_cep_MeV` 口径均为 `mu_q`。面向 `mu_B` 的外层扫描或输运脚本需要显式使用 `mu_B = 3mu_q`；新生成的 CEP reference CSV 与 phase artifact 会同时给出 `muq_cep_MeV` 和 `muB_cep_MeV` 以避免歧义。这些字段主要服务于回归治理、失败解释与策略比较，而不只是展示一个 `T_cep_MeV` 数值。

### production 收口路径

`run_production_phase_pipeline` 不把 CEP 搜索抽象成单次“对已有曲线集合调用 `find_cep`”，而是显式执行温度扫掠与 bracket 收缩：

1. 根据 `ProductionPipelineConfig` 生成温度网格
2. 对每个温度切片执行 Maxwell 有效性分类，状态为 `valid`、`unknown` 或 `invalid`
3. 将扫掠结果收集到 `FirstOrderSweepResult`
4. 从最后一个 `valid` 与随后的 `invalid` 温度切片构造 production bracket
5. 在 bracket 内继续二分，得到 production CEP

`uncertainty_T_MeV` 表示最终 CEP 温度 bracket 的半宽，完整宽度由
`bracket_width_T_MeV` 单独记录，不能把两者混作同一误差量。

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
本次不新增 crossover 的独立 gap-iteration 参数。

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
- `forced_invalid_count`
- `sweep_temperatures_MeV`
- `sweep_statuses`
- `rho_geometry_convergence` 及其位置量、密度量、面积残差门限
- `adaptive_temperature` 及其中点细化门限
- `crossover_T_max_MeV`

## 当前边界

- 本页解释的是算法关系与职责分工，不承诺所有底层函数都作为长期稳定公开入口
- 面向外部调用时，仍应优先回到 [Overview.md](docs/api/models/phase/Overview.md) 中列出的 `Models.*` 入口
- 旧路径算法页已降级为兼容跳转说明；当前新目录页已成为算法细节的主版本
