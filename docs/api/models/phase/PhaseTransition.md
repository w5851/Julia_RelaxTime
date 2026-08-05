# 相变判据与 Maxwell 构造

本文档描述相图主题中用于识别一阶相变的核心算法，主要实现位于 [src/models/phase/PhaseCore.jl](src/models/phase/PhaseCore.jl#L9) 与 [src/models/phase/PhaseIO.jl](src/models/phase/PhaseIO.jl#L1)。

## 涉及能力

- `SShapeResult`
- `detect_s_shape`
- `MaxwellResult`
- `maxwell_construction`
- `group_curves_by_temperature`

这些接口是 `run_phase_pipeline`、`find_cep` 等上层能力的判据基础，但通常不作为新用户的首选入口。

## S-shape 检测

### `SShapeResult`

`SShapeResult` 定义位于 [src/models/phase/PhaseCore.jl](src/models/phase/PhaseCore.jl#L9)。它用于表达单条 `μ(ρ)` 曲线是否存在一阶相变特征。

关键字段：

- `has_s_shape::Bool`
- `mu_spinodal_hadron`
- `mu_spinodal_quark`
- `rho_spinodal_hadron`
- `rho_spinodal_quark`
- `derivative_sign_changes::Int`

### `detect_s_shape(mu_vals, rho_vals; eps=0.0, min_points=5)`

实现位于 [src/models/phase/PhaseCore.jl](src/models/phase/PhaseCore.jl#L85)。

该函数对离散的 `μ(ρ)` 数据进行排序、求导与符号变化分析，检查是否存在 `正 → 负 → 正` 的导数序列。若成立，则认为该温度切片存在一阶相变区间。

算法步骤：

1. 按 `ρ` 升序整理输入数据
2. 计算离散导数 `dμ/dρ`
3. 忽略绝对值不超过 `eps` 的近零斜率后，检测 `正 → 负 → 正` 的符号变化模式
4. 将第一个 `正 → 负` 极大值记为 hadron-side spinodal，将随后第一个 `负 → 正` 极小值记为 quark-side spinodal
5. 用邻近三点二次插值细化两个 spinodal 极值点位置；若二次项方向或极值位置不可信，则回退到相邻点中点

因此，`mu_spinodal_hadron`/`mu_spinodal_quark` 是固定温度切片上 `μ(ρ)` 曲线的两个离散极值轨迹，不等同于稳定相分支切换点；稳定相切换由 Maxwell 等面积构造的一阶相变线给出。

典型用法：

```julia
mu, rho = get_curve_data(T=70.0)
result = detect_s_shape(mu, rho)

if result.has_s_shape
    println("Spinodal hadron: ρ=$(result.rho_spinodal_hadron), μ=$(result.mu_spinodal_hadron)")
    println("Spinodal quark: ρ=$(result.rho_spinodal_quark), μ=$(result.mu_spinodal_quark)")
end
```

## Maxwell 等面积构造

### `MaxwellResult`

`MaxwellResult` 与 `detect_s_shape` 配套使用，用于表示共存化学势与两侧共存密度的求解结果。

关键字段：

- `converged::Bool`
- `mu_transition`
- `rho_hadron`
- `rho_quark`
- `area_residual`
- `iterations::Int`
- `details`

### `maxwell_construction(mu_vals, rho_vals; kwargs...)`

实现位于 [src/models/phase/PhaseCore.jl](src/models/phase/PhaseCore.jl#L300)。

该函数在检测到 S-shape 后，搜索满足等面积条件的相变化学势与两侧相共存点。

常用参数：

- `min_samples`
- `detect_min_points`
- `detect_eps`
- `candidate_steps`
- `max_iter`
- `tol_area`
- `spinodal_hint`

求解流程：

1. 先做 S-shape 检测或复用传入的 `spinodal_hint`
2. 依据 spinodal 位置估计化学势搜索区间
3. 在每个化学势上枚举全部去重交点；只有恰有三个交点时面积才有效，拓扑间隙会重置前一段面积变号
4. 枚举所有有效面积根；多根或拓扑不唯一返回失败诊断，唯一根才进入结果
5. 用严格二分逼近等面积点；达到 `max_iter` 仍未满足 `tol_area` 时 `converged=false`

`details` 记录 `candidate_count`、`crossing_count`、`candidate_policy`、
`endpoint_dependent`、交点 bracket 和失败原因。因而存在 S 形但三交点/二分证书未闭合时，
上层应保持 numerical ambiguous/unknown，不能据此生成单调证书。

典型用法：

```julia
mu, rho = get_curve_data(T=70.0)
result = maxwell_construction(mu, rho)

if result.converged
    println("相变点: μ=$(result.mu_transition)")
    println("强子相密度: $(result.rho_hadron)")
    println("夸克相密度: $(result.rho_quark)")
end
```

## 温度切片整理

### `group_curves_by_temperature(rows; xi=0.0, tol=1e-6)`

实现位于 [src/models/phase/PhaseIO.jl](src/models/phase/PhaseIO.jl#L1)。

该函数负责从 CSV 行记录中提取指定 `xi` 的样本，并按温度重组为后续相图算法可消费的曲线字典。

返回值为按温度索引的曲线样本集合，通常作为 `find_cep` 或离线相变分析的输入准备步骤。

```julia
using CSV

rows = CSV.File("data/outputs/results/pnjl/scan/trho/trho_scan.csv")
grouped = group_curves_by_temperature(rows; xi=0.0)
```

## 与主流程的关系

- `detect_s_shape` 与 `maxwell_construction` 共同决定 `first_order_boundary` 和 `spinodal` 的质量
- `group_curves_by_temperature` 负责把扫描结果整理成可供 `find_cep` 消费的输入结构
- 如果你只是想跑通完整流程，应优先返回 [Overview.md](docs/api/models/phase/Overview.md)
