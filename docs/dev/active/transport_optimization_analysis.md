# 输运系数计算流程优化分析

## 日期：2025-12-28

## 当前计算流程

```
1. PNJL 能隙方程求解 (ThermoDerivatives.jl)
   └─ solve_equilibrium_mu() → 质量、Polyakov 环、热力学量

2. 热力学导数计算 (ThermoDerivatives.jl)
   └─ bulk_viscosity_coefficients() → v_n², ∂μB/∂T|σ, ∂M/∂T, ∂M/∂μB

3. 弛豫时间计算 (RelaxationTime.jl)
   └─ relaxation_times() → τ_u, τ_d, τ_s, τ_ubar, τ_dbar, τ_sbar

4. 输运系数积分 (TransportCoefficients.jl)
   └─ shear_viscosity(), bulk_viscosity_isentropic(), electric_conductivity()
```

## 优化建议

### 1. 代码清晰度优化

#### 问题 A: TransportCoefficients.jl 中存在大量重复代码

`shear_viscosity`, `electric_conductivity`, `bulk_viscosity`, `bulk_viscosity_isentropic` 四个函数有相似的循环结构。

**建议**: 提取通用积分框架

```julia
function integrate_transport_kernel(
    kernel_func::Function,  # (p, flavor, cosθ) -> kernel_value
    quark_params, thermo_params;
    p_nodes, p_max, cos_nodes, ...
)
    # 统一的积分循环
end
```

#### 问题 B: 内部函数定义在循环外部

`flavor_index`, `compute_B`, `one_flavor_pair_contrib` 等函数在每次调用时重新定义。

**建议**: 将这些函数提升为模块级函数或使用闭包缓存

### 2. 性能优化

#### 问题 A: ThermoDerivatives.jl 中重复求解能隙方程 ⭐ 高优先级

`bulk_viscosity_coefficients` 内部多次调用 `IMPLICIT_SOLVER`（计算 s, n, 导数时各调用一次）。

**解决方案**: 重构为一次求解，通过 `ForwardDiff.jacobian` 一并计算所有需要的热力学导数

优点：
- 无额外缓存开销
- 代码更简洁
- 更符合物理逻辑

**性能测试结果** (2025-12-28):
- 优化后：24 ms
- 单独调用总时间：208 ms
- **加速比：8.67x**

**状态**: ✅ 已实现 (2025-12-28)

#### 问题 B: Symbol 比较开销

`species in (:u, :d, :s)` 等 Symbol 比较在热循环中有开销。

**建议**: 使用整数索引代替 Symbol

```julia
const QUARK_U = 1
const QUARK_D = 2
const QUARK_S = 3
const ANTIQUARK_U = 4
```

**状态**: 待实现（低优先级）

#### 问题 C: 动态分派开销

`distribution_for_species` 中的条件分支在每次调用时执行。

**建议**: 使用函数指针数组或预编译的分派表

**状态**: 待实现（低优先级）

### 3. 稳定性优化

#### 问题 A: 缺少输入验证

函数没有检查输入参数的有效性（如 T > 0, μ 范围等）。

**建议**: 添加参数验证

```julia
function shear_viscosity(quark_params, thermo_params; ...)
    @assert thermo_params.T > 0 "Temperature must be positive"
    @assert all(m -> m >= 0, values(quark_params.m)) "Masses must be non-negative"
    ...
end
```

**状态**: 待实现（中优先级）

#### 问题 B: 数值稳定性

在极端参数下（如 T→0 或 μ→∞），可能出现数值问题。

**建议**: 添加数值保护

```julia
E = max(sqrt(p^2 + m^2), eps(Float64))
ff = clamp(f * (1 - f), 0.0, 0.25)
```

**状态**: 待实现（中优先级）

### 4. 接口优化

#### 问题 A: 参数传递繁琐

每个函数都需要传递大量参数（`p_nodes`, `p_max`, `cos_nodes` 等）。

**建议**: 使用配置结构体

```julia
struct IntegrationConfig
    p_nodes::Int
    p_max::Float64
    cos_nodes::Int
    p_grid::Union{Nothing,Vector{Float64}}
    p_w::Union{Nothing,Vector{Float64}}
end

const DEFAULT_CONFIG = IntegrationConfig(64, 10.0, 32, nothing, nothing)
```

**状态**: 待实现（中优先级）

#### 问题 B: 两种 bulk_viscosity 函数容易混淆

`bulk_viscosity` 和 `bulk_viscosity_isentropic` 名称相似但公式不同。

**建议**: 
- 方案1: 更明确的命名 `bulk_viscosity_thermodynamic()` vs `bulk_viscosity_isentropic()`
- 方案2: 统一接口 `bulk_viscosity(...; formula=:isentropic)`

**状态**: 待讨论

### 5. 文档和测试优化

#### 问题 A: 缺少完整的 API 文档

**建议**: 为每个公开函数添加完整的 docstring

**状态**: 待实现（低优先级）

#### 问题 B: 缺少集成测试

**建议**: 添加端到端测试

**状态**: 待实现（低优先级）

## 优先级总结

| 优先级 | 优化项 | 预期收益 | 状态 |
|--------|--------|----------|------|
| 高 | 一并计算热力学导数 | 性能提升 8.67x | ✅ 已实现 |
| 中 | 添加配置结构体 | 代码清晰度 | 待实现 |
| 中 | 提取通用积分框架 | 减少代码重复 | 待实现 |
| 中 | 添加输入验证 | 稳定性 | 待实现 |
| 低 | Symbol → 整数索引 | 微小性能提升 | 待实现 |
| 低 | 完善文档 | 可维护性 | 待实现 |

## 更新记录

- 2025-12-28: 创建文档，完成初步分析
- 2025-12-28: 实现"一并计算热力学导数"优化

## 性能测试结果对比

### bulk_viscosity_coefficients 优化前后对比

| 指标 | 优化前 (17:34) | 优化后 (19:10) | 改进 |
|------|---------------|---------------|------|
| 中位时间 | 162.69 ms | 24.09 ms | **6.75x 加速** |
| 最小时间 | 156.98 ms | 23.40 ms | **6.71x 加速** |
| 内存分配 | 415.36 KB | 70.09 KB | **5.93x 减少** |

### 其他函数性能（优化后）

| 函数 | 中位时间 | 内存分配 |
|------|----------|----------|
| v_n_squared | 61.06 ms | 174.00 KB |
| dmuB_dT_sigma | 61.01 ms | 174.12 KB |
| bulk_viscosity_coefficients | 24.09 ms | 70.09 KB |
| thermo_derivatives | 175.36 ms | 508.20 KB |
| compute_B_bracket | 26.73 ns | 0 bytes |

### 关键发现

1. **bulk_viscosity_coefficients 现在比单独调用 v_n_squared 或 dmuB_dT_sigma 还快**
   - 因为优化后只调用一次 IMPLICIT_SOLVER + 一次 jacobian
   - 而 v_n_squared/dmuB_dT_sigma 仍然使用旧的多次调用方式

2. **内存分配大幅减少**
   - 从 415 KB 减少到 70 KB
   - 减少了约 83% 的内存分配

3. **建议后续优化**
   - 可以考虑用类似方法优化 v_n_squared 和 dmuB_dT_sigma
   - 或者直接使用 bulk_viscosity_coefficients 替代它们
