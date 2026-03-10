# Crossover 检测

本文档描述相图主题中的 crossover 支路。相关实现位于 [src/models/phase/CrossoverLine.jl](src/models/phase/CrossoverLine.jl#L1)。

## 适用场景

当某一温度区间不再形成一阶相变时，相图流程可通过序参量导数结构估计 crossover 温度。当前核心接口为：

- `CrossoverResult`
- `detect_crossover`
- `scan_crossover_line`

## 物理背景

在高温区域，手征相变和退禁闭相变会表现为平滑 crossover。此时主问题不再是求 Maxwell 等面积，而是从序参量导数的峰值或拐点中选出代表性的 crossover 温度。

实现中主要支持两类判据：

- `:peak`
- `:inflection`

它们都服务于 `run_phase_pipeline` 的可选 crossover 分支，而不是替代一阶相变判据。

## `CrossoverResult`

结构体定义位于 [src/models/phase/CrossoverLine.jl](src/models/phase/CrossoverLine.jl#L1)。

关键字段：

- `found::Bool`
- `T_crossover`
- `method::Symbol`
- `derivative_value`
- `iterations::Int`
- `details`

## `detect_crossover`

实现位于 [src/models/phase/CrossoverLine.jl](src/models/phase/CrossoverLine.jl#L209)。

```julia
detect_crossover(μ_fm, T_range; method=:peak, variable=:phi_u, kwargs...)
```

常用参数：

- `μ_fm`
- `T_range`
- `method`
- `variable`
- `xi`
- `n_scan`
- `tol`
- `max_iter`
- `p_num`
- `t_num`

支持的检测变量通常包括：

- `:phi_u`
- `:Phi`

### 峰值法 `:peak`

实现内部会搜索 `|∂φ/∂T|` 的主要峰值，并在候选区间内继续细化，适合直接从导数强度定位 crossover 温度。

### 拐点法 `:inflection`

实现内部会搜索 `∂²φ/∂T² = 0` 的有效拐点，并过滤高温区的虚假拐点。该方法通常更适合批量计算。

## `scan_crossover_line`

实现位于 [src/models/phase/CrossoverLine.jl](src/models/phase/CrossoverLine.jl#L253)。

```julia
scan_crossover_line(mu_range, T_range; method=:peak, variable=:phi_u, kwargs...)
```

该函数沿给定化学势区间批量调用 `detect_crossover`，返回 crossover 线的离散采样结果。

## 示例

```julia
const hbarc = 197.327

μ_fm = 0.0 / hbarc
T_range_fm = (100.0 / hbarc, 250.0 / hbarc)

result_peak = detect_crossover(μ_fm, T_range_fm; method=:peak, variable=:phi_u)
result_infl = detect_crossover(μ_fm, T_range_fm; method=:inflection, variable=:phi_u)

mu_range_fm = (0.0 / hbarc, 200.0 / hbarc, 10)
results = scan_crossover_line(mu_range_fm, T_range_fm; method=:inflection)
```

## 与主流程的关系

- 当 `run_phase_pipeline` 启用 `compute_crossover=true` 时，会进入 crossover 支路
- crossover 结果最终写入 `crossover_line.csv`
- 如果你主要关注一阶相变判据，应优先阅读 [PhaseTransition.md](docs/api/models/phase/PhaseTransition.md)