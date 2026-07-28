# 相图工作流总览

本文档是相图主题的主入口，优先面向“先跑通相图生成与工件输出”的使用者。

## 何时使用本主题

当你需要以下任一能力时，应优先从 `Models` 聚合入口开始：

- 从 `T-ρ` 扫描构建一阶相变边界与 spinodal
- 按 production/baseline 口径执行高精度一阶相变与 CEP 搜索
- 估计 CEP 并获取诊断信息
- 在 CEP 邻域做 compare-only `P-mu / Omega-mu` 双分支竞争诊断
- 可选生成 crossover line
- 将结果写出为 CSV、JSON、Markdown 报告
- 将已验证的结果晋升到 reference 目录

研究口径下，`run_phase_pipeline` 默认使用 `cep_strategy=:interpolate`，且 crossover 默认方法为 `crossover_method=:peak`。更细的 research 调参应按具体实验目标单独显式传参，而不是把它当作总览页的默认使用方式。

`Models.run_phase_pipeline` 当前同时支持 `mode=:research` 与 `mode=:production` 两条口径。

脚本入口 `scripts/pnjl/calculate_phase_structure.jl` 在未显式传参时默认走 `--mode=production`；当你需要研究口径时，必须显式传 `mode=:research` 或 `--mode=research`。

## 首选公开入口

相图主题当前推荐的统一入口位于 [src/models/Models.jl](src/models/Models.jl#L82) 与 [src/models/entrypoints.jl](src/models/entrypoints.jl#L15)。

优先使用的 API：

- `Models.run_phase_pipeline`
- `Models.run_production_phase_pipeline`
- `Models.find_cep`
- `Models.build_phase_artifacts`
- `Models.resolve_phase_output_target`
- `Models.promote_phase_artifacts`
- `Models.analyze_pm_branch_competition`

结果类型：

- `PhasePipelineResult`
- `CEPResult`
- `FirstOrderSweepResult`
- `ProductionPipelineConfig`
- `PromotionResult`

完整导出基线见 [generated/Exports.md](docs/api/models/phase/generated/Exports.md)。

## 最短工作流

### 通用 / research 入口

```julia
using .Models

tmp = mktempdir()

result = Models.run_phase_pipeline(
    :PNJL;
    mode=:research,
    T_grid=[150.0],
    rho_grid=[0.1, 0.2, 0.3],
    xi=0.0,
    output_dir=tmp,
    profile=:smoke,
    solver_backend=:models,
    p_num=12,
    t_num=4,
    iterations=80,
    rho_geometry_convergence=true,
    adaptive_temperature=true,
    crossover_mu0_only=true,
    crossover_T_max_MeV=240.0,
    promote_reference=false,
)
```

这条链路与 [tests/integration/models/test_phase_pipeline_smoke.jl](tests/integration/models/test_phase_pipeline_smoke.jl) 对齐，适合先验证接口、参数与工件写出是否正常。

### production / baseline 入口

```julia
using .Models

tmp = mktempdir()

result = Models.run_production_phase_pipeline(
    :PNJL;
    T_start=130.0,
    T_end=132.0,
    dT=1.0,
    rho_grid=[0.0, 0.1, 0.2],
    xi=0.0,
    output_dir=tmp,
    profile=:smoke,
    solver_backend=:models,
    p_num=12,
    t_num=4,
    iterations=80,
    promote_reference=false,
)
```

当你希望显式采用 production 的高精度温度扫描、unknown budget 与非插值 CEP 收口逻辑时，优先使用 `Models.run_production_phase_pipeline`。CEP 返回值采用三态合同：`resolved`、`ambiguous`、`not_found`；ambiguous 结果保留最后确认 Maxwell 与首个确认单调温度，不发布温度中点或借用的 Maxwell 化学势。

## Phase 热积分策略

PNJL 标量 phase thermodynamics 入口支持两种显式策略：

- `thermo_quadrature_policy=:tensor_gauss`：保留既有有限动量截断的 `p_num × t_num` Gauss--Legendre 路径，也是兼容默认值；
- `thermo_quadrature_policy=:rs_reduced_adaptive`：仅用于角度依赖完全包含在 RS 分布自变量
  `E_xi = sqrt(m^2 + p^2(1 + xi*cos(theta)^2))` 的标量热项，把角积分约化为测度因子 `1/sqrt(1+xi)`，并在 `[0, Inf)` 上执行一维自适应径向积分。

自适应策略的控制量为 `thermo_quadrature_rtol`、`thermo_quadrature_atol` 与
`thermo_quadrature_maxevals`，最终值会写入 `config_snapshot` 和 config hash。这里的
`E_xi` 只是分布函数自变量，不是新的物理色散关系；chi、Polyakov 势和 vacuum 项不增加角度依赖。该约化不适用于 magnetic、transport 或其他含独立角核的路径。
固定态积分诊断可调用 `Models.PNJLIntegrals.calculate_log_sum_rs_reduced_adaptive_with_error` 同时取得数值和局部求积误差估计。该估计由仓库内 16/32 阶 Gauss--Legendre 双规则的差并经安全因子累加得到，不依赖外部自适应积分库，也不能替代求解级与相线级收敛审计。

当只需要零化学势 crossover 时，设置 `crossover_mu0_only=true`。research 与 production 路径都会只求解
`mu_q=0`，并把这一有效采样策略写入 `config_snapshot` 与 config hash；它不是仅供 manifest 使用的标签。

固定状态热核和直接数密度定义了严格 `T=0` 极限；五变量 PNJL gap/phase solve 因 Polyakov 场在严格零温退化而显式要求 `T>0`。因此正式“全温区”reference 的下限必须是经收敛验证的严格正温，除非后续引入独立零温求解合同。

## 输入与输出约定

- 外部接口的温度、化学势和结果产物统一以 `MeV` 口径表述
- `rho_grid` 是归一化密度网格，数值含义依赖所选模型与扫描配置
- `phase_summary.json` 负责输出 schema、配置快照与诊断统计
- 主要工件包括：
  - `trho_scan.csv`
  - `first_order_boundary.csv`
  - `spinodal.csv`
- `crossover_line.csv`
- `phase_grid_convergence.csv`
- `phase_summary.json`
  - `pm_branch_scan.csv`
  - `pm_phase_summary.json`
  - `pm_vs_maxwell.csv`
  - `phase_report.md`

## 常见使用模式

### 1. 完整跑一次并消费结果

优先调用 `Models.run_phase_pipeline`，然后读取 `PhasePipelineResult` 中的：

- `cep`
- `first_order_boundary`
- `spinodal`
- `crossover_line`
- `artifact_paths`
- `diagnostics`

如果你要显式锁定 production 口径，则改用 `Models.run_production_phase_pipeline`；它仍返回 `PhasePipelineResult`，但 `config_snapshot` 与 `diagnostics` 会包含 production 专有字段，例如 `dT_initial`、`unknown_budget`、`unknown_budget_exhausted`、`first_point_fallback`、兼容保留的 `forced_invalid_count`（三态路径为 `0`）、rho/T 网格误差门限与显式 crossover 温区上限。正式产物应同时消费 `phase_grid_convergence.csv`，不能只检查边界 CSV 是否存在。

### 2. 对已有曲线离线做 CEP 分析

当你已经拥有按温度分组的 `μ(ρ)` 曲线时，直接调用 `Models.find_cep` 即可，不必重新走完整 pipeline。

### 3. 先解析输出路径，再运行与晋升

如果你需要把产物写入受控目录，再决定是否晋升参考数据，先调用 `Models.resolve_phase_output_target`，完成主流程后再调用 `Models.promote_phase_artifacts`。对应的最小验证路径可参考 [tests/integration/models/test_phase_artifacts_promotion_smoke.jl](tests/integration/models/test_phase_artifacts_promotion_smoke.jl)。

### 4. 在 production 模式下查看温度扫掠诊断

当你需要解释 production 路径是如何从 `valid/unknown/invalid` 温度切片收口到 CEP 时，可查看：

- `Models.FirstOrderSweepResult`
- `Models.ProductionPipelineConfig`

这些类型主要服务于 production/baseline 诊断与测试，不是多数用户的第一入口，但它们定义了 production 流程的结果与配置合同。

### 5. 在 CEP 邻域做 `P-mu` 诊断对照

当你不想替换当前 Maxwell 主线，只想在固定温度下比较 hadron-like / quark-like 双分支竞争与 Maxwell 判据是否一致时，优先调用 `Models.analyze_pm_branch_competition`。

该入口当前定位为 compare-only 诊断：

- 输入 `T_values` 与固定 `mu_grid`
- 对每个 `(T, mu)` 同时追踪 hadron / quark 两支
- 输出：
  - `pm_branch_scan.csv`
  - `pm_phase_summary.json`
  - `pm_vs_maxwell.csv`

最小示例见 [PMPhaseDiagnostic.md](docs/api/models/phase/PMPhaseDiagnostic.md)。

## 非首选入口

以下接口虽然重要，但不应作为多数用户的第一站：

- `detect_s_shape`
- `maxwell_construction`
- `detect_crossover`
- `scan_crossover_line`
- `AdaptiveRhoConfig`
- `suggest_refinement_points`
- `merge_rho_values`
- `analyze_pm_branch_competition` 的内部 branch helper 与 artifact helper

这些能力属于算法核心层，关系说明见 [Algorithms.md](docs/api/models/phase/Algorithms.md)。

## 迁移期说明

- 从本阶段起，`docs/api/models/phase/Overview.md` 与 `docs/api/models/phase/README.md` 共同承担相图主题主入口
- 旧 `pnjl` 相图兼容页已退出当前主导航
