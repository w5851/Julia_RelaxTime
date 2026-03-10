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

上层编排主要位于 [src/models/phase/PhasePipeline.jl](src/models/phase/PhasePipeline.jl)，结果类型定义位于 [src/models/phase/PhaseTypes.jl](src/models/phase/PhaseTypes.jl)。

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

- `eval_count`
- `unknown_count`
- `uncertainty_T_MeV`
- `reason`
- `method`

这些字段主要服务于回归治理、失败解释与策略比较，而不只是展示一个 `T_cep_MeV` 数值。

## Crossover 支路

当温度高于 CEP 或未形成一阶相变时，pipeline 可转而计算 crossover line。

相图主题中主要涉及两个能力：

- `detect_crossover`
- `scan_crossover_line`

现有实现支持至少两类判据：

- `:peak`
- `:inflection`

它们的职责是从序参量导数结构中给出 crossover 温度，而不是替代一阶相变的 Maxwell 判据。细节见 [Crossover.md](docs/api/models/phase/Crossover.md)。

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

## 当前边界

- 本页解释的是算法关系与职责分工，不承诺所有底层函数都作为长期稳定公开入口
- 面向外部调用时，仍应优先回到 [Overview.md](docs/api/models/phase/Overview.md) 中列出的 `Models.*` 入口
- 旧路径算法页已降级为兼容跳转说明；当前新目录页已成为算法细节的主版本