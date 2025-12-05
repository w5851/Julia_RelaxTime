# PhaseTransition 模块（S 形检测）

`src/pnjl/analysis/PhaseTransition.jl` 提供针对 `T-ρ` 扫描数据的判稳工具，用于识别 ρ(μ) 曲线中的 S 形（即出现负压缩率区域）。

## 导出符号

| 名称 | 类型 | 说明 |
|------|------|------|
| `SShapeResult` | 结构体 | 记录是否存在 S 形及自旋odal 点估计 |
| `detect_s_shape` | 函数 | 输入 μ/ρ 序列，检测导数符号变化并返回 `SShapeResult` |
| `group_curves_by_temperature` | 函数 | 解析 CSV 行或 Dict 序列，按温度聚合 `(μ, ρ)` 样本 |

## `SShapeResult`

字段：
- `has_s_shape::Bool`
- `mu_spinodal_low/high::Union{Nothing, Float64}`：两个自旋odal 化学势（μ），若无法估计则为 `nothing`
- `rho_spinodal_low/high::Union{Nothing, Float64}`：对应的密度值
- `derivative_sign_changes::Int`：ρ 对 μ 的导数符号变化次数（用于调试/诊断）

`SShapeResult()` 构造函数默认返回“无 S 形”。

## `detect_s_shape(mu_vals, rho_vals; eps=1e-6, min_points=5)`

- 将数据按 μ 升序排序，计算离散导数 `Δρ/Δμ`；
- 忽略梯度绝对值小于 `eps` 的段；
- 若检测到符号序列 `+ → - → +`，则判定存在 S 形，并用相邻样本的中点近似自旋odal；
- 不满足条件时返回 `SShapeResult()`。

典型用法：
```julia
mu = [-50.0, -20.0, -10.0, 0.0, 5.0, 15.0, 30.0]
rho = [-0.6, -0.25, -0.1, 0.2, 0.35, 0.30, 0.45]
res = detect_s_shape(mu, rho)
@assert res.has_s_shape
```

## `group_curves_by_temperature(rows; xi=0.0, tol=1e-6)`

- `rows` 可为 `Vector{Dict}` 或类似对象，需至少包含字段 `"T_MeV"`、`"xi"`、`"rho"`，以及 `"mu_avg_MeV"` 或 `"mu_MeV"`；
- 根据温度 (MeV) 归类，并过滤出 `ξ` 满足 `|ξ - xi| ≤ tol` 的样本；
- 返回 `Dict{Float64, Vector{Tuple{Float64, Float64}}}`，其中每个样本为 `(μ, ρ)`。

该函数可直接接收由 `TrhoScan` 生成、经 `CSV.read` 或 `CSV.File` 解析后的行。