# PhaseTransition 模块（相变分析）

`src/pnjl/analysis/PhaseTransition.jl` 提供 PNJL 模型的相变分析功能，包括 S 形曲线检测和 Maxwell 等面积构造。

## 导出符号

| 名称 | 类型 | 说明 |
|------|------|------|
| `SShapeResult` | 结构体 | S 形检测结果，包含 spinodal 点信息 |
| `detect_s_shape` | 函数 | 检测 μ(ρ) 曲线的 S 形特征 |
| `MaxwellResult` | 结构体 | Maxwell 构造结果，包含相变点信息 |
| `maxwell_construction` | 函数 | Maxwell 等面积构造，计算相变点 |
| `group_curves_by_temperature` | 函数 | 按温度分组 CSV 数据 |

## S 形检测

### `SShapeResult`

字段：
- `has_s_shape::Bool`：是否检测到 S 形
- `mu_spinodal_hadron::Union{Nothing, Float64}`：μ(ρ) 曲线的局部极大值（强子相侧 spinodal 的 μ）
- `mu_spinodal_quark::Union{Nothing, Float64}`：μ(ρ) 曲线的局部极小值（夸克相侧 spinodal 的 μ）
- `rho_spinodal_hadron::Union{Nothing, Float64}`：强子相侧 spinodal 的 ρ（极大值点）
- `rho_spinodal_quark::Union{Nothing, Float64}`：夸克相侧 spinodal 的 ρ（极小值点）
- `derivative_sign_changes::Int`：导数符号变化次数

### `detect_s_shape(mu_vals, rho_vals; eps=0.0, min_points=5)`

检测 μ(ρ) 曲线是否呈 S 形（导数 dμ/dρ 符号变化：正 → 负 → 正）。

物理意义：
- S 形曲线表示一阶相变区域
- `spinodal_hadron`：μ(ρ) 的局部极大值点，对应强子相的亚稳态边界
- `spinodal_quark`：μ(ρ) 的局部极小值点，对应夸克相的亚稳态边界

算法：
1. 按 ρ 升序排列数据
2. 计算离散导数 dμ/dρ
3. 检测符号变化序列 +1 → -1 → +1
4. 使用二次插值细化极值点位置

```julia
mu, rho = get_curve_data(T=70.0)
result = detect_s_shape(mu, rho)
if result.has_s_shape
    println("Spinodal hadron: ρ=$(result.rho_spinodal_hadron), μ=$(result.mu_spinodal_hadron)")
    println("Spinodal quark: ρ=$(result.rho_spinodal_quark), μ=$(result.mu_spinodal_quark)")
end
```

## Maxwell 构造

### `MaxwellResult`

字段：
- `converged::Bool`：是否收敛
- `mu_transition::Union{Nothing, Float64}`：相变化学势 (MeV)
- `rho_hadron::Union{Nothing, Float64}`：强子相密度 (ρ/ρ₀)
- `rho_quark::Union{Nothing, Float64}`：夸克相密度 (ρ/ρ₀)
- `area_residual::Union{Nothing, Float64}`：面积残差
- `iterations::Int`：迭代次数
- `details::Dict{Symbol, Any}`：详细信息

### `maxwell_construction(mu_vals, rho_vals; kwargs...)`

Maxwell 等面积构造，计算一阶相变的共存化学势和共存密度。

参数：
- `mu_vals`：化学势数组
- `rho_vals`：密度数组
- `min_samples`：最小样本数 (默认 12)
- `detect_min_points`：S 形检测最小点数 (默认 6)
- `detect_eps`：S 形检测斜率阈值 (默认 1e-6)
- `candidate_steps`：初始搜索步数 (默认 64)
- `max_iter`：最大迭代次数 (默认 60)
- `tol_area`：面积容差 (默认 1e-4)
- `spinodal_hint`：预计算的 S 形结果 (可选)

算法：
1. 检测 S 形曲线，获取 spinodal 点
2. 从 spinodal 点估计 μ 搜索区间
3. 在区间内搜索面积差符号变化
4. 二分法精确求解等面积点

```julia
mu, rho = get_curve_data(T=70.0)
result = maxwell_construction(mu, rho)
if result.converged
    println("相变点: μ=$(result.mu_transition) MeV")
    println("强子相密度: $(result.rho_hadron) ρ₀")
    println("夸克相密度: $(result.rho_quark) ρ₀")
end
```

## 工具函数

### `group_curves_by_temperature(rows; xi=0.0, tol=1e-6)`

按温度分组 CSV 数据行，用于下游分析。

参数：
- `rows`：CSV 行的迭代器（每行是 Dict）
- `xi`：要筛选的各向异性参数
- `tol`：xi 匹配容差

返回：
`Dict{Float64, Vector{Tuple{Float64, Float64}}}`：温度 → (μ, ρ) 对列表

```julia
using CSV
rows = CSV.File("data/outputs/results/pnjl/trho_scan.csv")
grouped = group_curves_by_temperature(rows; xi=0.0)
for (T, samples) in grouped
    println("T=$T MeV: $(length(samples)) 个样本")
end
```

## 典型工作流

```julia
using PNJL.PhaseTransition

# 1. 加载数据
rows = CSV.File("trho_scan.csv")
grouped = group_curves_by_temperature(rows; xi=0.0)

# 2. 对每个温度进行分析
for T in sort(collect(keys(grouped)))
    samples = grouped[T]
    mu_vals = [s[1] for s in samples]
    rho_vals = [s[2] for s in samples]
    
    # S 形检测
    s_result = detect_s_shape(mu_vals, rho_vals)
    if !s_result.has_s_shape
        println("T=$T MeV: 无 S 形（crossover 或 CEP 以上）")
        continue
    end
    
    # Maxwell 构造
    m_result = maxwell_construction(mu_vals, rho_vals; spinodal_hint=s_result)
    if m_result.converged
        println("T=$T MeV: μ_c=$(m_result.mu_transition), ρ_h=$(m_result.rho_hadron), ρ_q=$(m_result.rho_quark)")
    end
end
```

## 注意事项

1. **数据质量**：确保 T-ρ 扫描数据足够密集，特别是在 S 形区域
2. **ρ=0 奇异点**：建议使用反向扫描（从 ρ_max 到 0）避免 ρ=0 处的数值问题
3. **CEP 附近**：接近 CEP 时 S 形变弱，可能需要更高的数据精度
