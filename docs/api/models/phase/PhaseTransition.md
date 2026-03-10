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
3. 检测符号变化模式
4. 细化两个 spinodal 极值点位置

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

实现位于 [src/models/phase/PhaseCore.jl](src/models/phase/PhaseCore.jl#L274)。

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
3. 搜索面积差符号变化
4. 用二分法逼近等面积点

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