# Crossover 检测 API 文档

## 概述

本模块提供 PNJL 模型中 crossover 相变温度的检测功能。在 T > T_CEP 区域，手征相变和退禁闭相变变为平滑的 crossover，需要通过序参量的导数特征来定义相变温度。

## 物理背景

### 双峰结构

手征凝聚 φ_u 对温度的导数 ∂φ_u/∂T 存在双峰结构：
- **小峰**（~147 MeV）：由退禁闭相变引起，随化学势增加逐渐消失
- **大峰**（~200 MeV）：真正的手征 crossover

因此，检测手征 crossover 需要找到**最大的峰**，而不是第一个峰。

### 多拐点结构

二阶导数 ∂²φ_u/∂T² 存在多个符号变化（拐点）：
- 第一个 + → -（~150 MeV）：小峰的拐点
- 最后一个 + → -（~200 MeV）：大峰的拐点（手征 crossover）

## 检测方法

### 峰值法 (`:peak`)

**定义**：找 |∂φ/∂T| 的最大值位置

**算法**：
1. 从高温向低温扫描，计算 |∂φ/∂T|
2. 找所有局部极大值
3. 选择最大的局部极大值
4. 使用**黄金分割法**在峰值附近细化

**黄金分割法**：
- 用于单峰函数的极值搜索
- 每次迭代将搜索区间缩小为原来的 φ = (√5-1)/2 ≈ 0.618 倍
- 收敛速度：线性收敛，每次迭代需要 2 次函数求值

### 拐点法 (`:inflection`)

**定义**：找 ∂²φ/∂T² = 0 的位置（从正变负）

**算法**：
1. 从高温向低温扫描，计算 ∂²φ/∂T²
2. 找第一个 - → + 的符号变化（对应原方向的 + → -）
3. 使用**二分法**细化拐点位置

**二分法**：
- 用于连续函数的零点搜索
- 每次迭代将搜索区间缩小一半
- 收敛速度：线性收敛，每次迭代需要 1 次函数求值

## API 参考

### `detect_crossover`

```julia
detect_crossover(μ_fm, T_range; method=:peak, variable=:phi_u, kwargs...)
```

检测 crossover 温度。

**参数**：
- `μ_fm::Real`: 化学势 (fm⁻¹)
- `T_range::Tuple{Real, Real}`: 温度搜索范围 (T_min, T_max) (fm⁻¹)
- `method::Symbol`: 检测方法
  - `:peak`: 峰值法（找 |∂φ/∂T| 最大值）
  - `:inflection`: 拐点法（找 ∂²φ/∂T² = 0）
- `variable::Symbol`: 检测变量
  - `:phi_u`: 手征凝聚 φ_u（默认）
  - `:Phi`: Polyakov loop Φ
- `xi::Real=0.0`: 各向异性参数
- `n_scan::Int=20`: 初始扫描点数
- `tol::Real=1e-4`: 收敛容差 (fm⁻¹)
- `max_iter::Int=20`: 最大迭代次数
- `p_num::Int=24`: 动量积分节点数
- `t_num::Int=12`: 角度积分节点数

**返回**：`CrossoverResult` 结构体

### `CrossoverResult`

```julia
struct CrossoverResult
    found::Bool                           # 是否找到 crossover
    T_crossover::Union{Nothing, Float64}  # crossover 温度 (fm⁻¹)
    method::Symbol                        # 使用的方法
    derivative_value::Union{Nothing, Float64}  # 在 crossover 点的导数值
    iterations::Int                       # 迭代次数
    details::Dict{Symbol, Any}            # 详细信息
end
```

**details 字段**：
- `:T_range`: 搜索范围
- `:n_scan`: 扫描点数
- `:variable_index`: 变量索引
- `:n_peaks_found` / `:n_inflections_found`: 找到的峰/拐点数量
- `:all_peaks` / `:all_inflections`: 所有峰/拐点的位置
- `:scan_data`: 扫描数据

### `scan_crossover_line`

```julia
scan_crossover_line(mu_range, T_range; method=:peak, variable=:phi_u, kwargs...)
```

扫描 crossover 线。

**参数**：
- `mu_range::Tuple{Real, Real, Int}`: 化学势范围 (μ_min, μ_max, n_points) (fm⁻¹)
- `T_range::Tuple{Real, Real}`: 温度搜索范围 (T_min, T_max) (fm⁻¹)
- 其他参数同 `detect_crossover`

**返回**：`Vector{NamedTuple}` 包含 `(mu_fm, T_crossover_fm, converged, derivative)`

## 使用示例

```julia
using PNJL

const hbarc = 197.327  # MeV·fm

# 单点检测
μ_fm = 0.0 / hbarc
T_range_fm = (100.0 / hbarc, 250.0 / hbarc)

# 峰值法
result_peak = detect_crossover(μ_fm, T_range_fm; method=:peak, variable=:phi_u)
println("T_crossover (peak) = ", result_peak.T_crossover * hbarc, " MeV")

# 拐点法
result_infl = detect_crossover(μ_fm, T_range_fm; method=:inflection, variable=:phi_u)
println("T_crossover (inflection) = ", result_infl.T_crossover * hbarc, " MeV")

# 扫描 crossover 线
mu_range_fm = (0.0 / hbarc, 200.0 / hbarc, 10)
results = scan_crossover_line(mu_range_fm, T_range_fm; method=:inflection)
for r in results
    println("μ = ", r.mu_fm * hbarc, " MeV, T_c = ", r.T_crossover_fm * hbarc, " MeV")
end
```

## 性能对比

| 方法 | 每次迭代计算量 | 收敛速度 | 典型耗时 |
|------|--------------|---------|---------|
| 峰值法 | 2 次导数计算 | 线性 (φ ≈ 0.618) | ~4-6 s |
| 拐点法 | 1 次导数计算 | 线性 (0.5) | ~1-4 s |

**建议**：拐点法更快，推荐用于批量计算。

## 验证结果

与文献值对比（ξ = 0）：

| 变量 | 本模块结果 | 文献值 |
|------|-----------|--------|
| T_φ (手征) | ~200 MeV | 200 MeV |
| T_Φ (退禁闭) | ~177 MeV | 171 MeV |

## 相关文档

- [PhaseTransition.md](./PhaseTransition.md) - 一阶相变分析
- [PNJL_Phase_Diagram_TODO.md](../../dev/PNJL_Phase_Diagram_TODO.md) - 开发计划
