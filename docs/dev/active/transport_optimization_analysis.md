# 输运系数计算流程优化分析

## 日期：2025-12-28

## 2026-01-24 状态复核（完成度评估）

- ✅ 已实现（与文档描述一致）："ThermoDerivatives.jl 中重复求解能隙方程" 的高优先级优化。
   - 现状：`bulk_viscosity_coefficients` 通过一次 `IMPLICIT_SOLVER` 求解拿到基态 `x_base`，再用 `ForwardDiff.jacobian(solve_state, θ)` 统一得到 $dx/dθ$，后续用链式法则拼出所需导数。
   - 代码位置：[src/pnjl/derivatives/ThermoDerivatives.jl](src/pnjl/derivatives/ThermoDerivatives.jl)

- ✅ 已实现：采用配置结构体简化参数传递（Phase 1：TransportCoefficients）。
   - `TransportIntegrationConfig` 已在 [src/relaxtime/TransportCoefficients.jl](src/relaxtime/TransportCoefficients.jl) 引入，并提供参数校验（`p_grid/p_w`、`cos_grid/cos_w` 成对出现且长度一致）。
   - `transport_coefficients/shear_viscosity/electric_conductivity/bulk_viscosity_isentropic` 已支持 `config=TransportIntegrationConfig(...)`；若同时给出 `p_nodes/p_max/...` 等显式关键字，显式关键字会覆盖 `config`（保持向后兼容、便于调参）。
   - 上层工作流 [src/pnjl/workflows/TransportWorkflow.jl](src/pnjl/workflows/TransportWorkflow.jl) 已支持 `transport_config=...` 作为调用侧入口，并兼容旧的 `transport_kwargs` 传参方式。
   - 单测覆盖：
     - [tests/unit/relaxtime/test_transport_coefficients.jl](tests/unit/relaxtime/test_transport_coefficients.jl)
     - [tests/unit/relaxtime/test_transport_workflow.jl](tests/unit/relaxtime/test_transport_workflow.jl)
   - 注：为避免 include 顺序导致“同名结构体来自不同模块实例”的类型不一致，`TransportWorkflow` 复用了 `Main.TransportCoefficients`（细节见该文件头部）。

- ❌ 仍未实现（因此“整份任务”不算完成，不建议归档到 archived）：
   - `TransportCoefficients.jl` 的通用积分框架抽取与去重复
   - 输入参数验证（如 `T>0` 等）与极端参数数值保护
   - `Symbol` → 整数索引、分派表等低优先级性能微优化
   - 更完整的 API 文档与集成测试补齐
   - 体粘滞函数命名/接口统一方案（仍待讨论）

> 备注：项目层面的“参数结构体化”（例如 `QuarkParams/ThermoParams` 的跨模块复用）已单独抽出进度文档，以避免本文件堆积：
> - [docs/dev/active/parameter_structs_migration.md](docs/dev/active/parameter_structs_migration.md)

## 方案记录：配置结构体简化参数（Phase 1：TransportCoefficients）

目标：在不破坏现有关键字参数调用方式的前提下，引入一个“积分配置结构体”，统一承载动量/角度积分的常用参数，减少函数签名中的重复参数，并让上层 workflow 更容易透传配置。

### 设计原则

1. **完全向后兼容**：原有调用 `...; p_nodes=..., p_max=..., cos_nodes=...` 必须保持可用且行为不变。
2. **新增更高层入口**：新增 `config=...` 作为统一入口；调用方可以只传一个结构体。
3. **优先级规则**：若同时提供 `config` 与显式关键字（如 `p_nodes=...`），则 **显式关键字覆盖 config**（便于逐步迁移与调参）。
4. **最小侵入落地**：Phase 1 只覆盖 `src/relaxtime/TransportCoefficients.jl`；后续再评估是否将 `RelaxationTime/AverageScatteringRate` 的积分参数也统一到同一套配置系统。

### 结构体草案

建议命名：`TransportIntegrationConfig`（避免与其它模块的更通用 `IntegrationConfig` 撞名）。

字段（与现有 API 一一对应）：
- `p_nodes::Int`
- `p_max::Float64`
- `p_grid::Union{Nothing,Vector{Float64}}`
- `p_w::Union{Nothing,Vector{Float64}}`
- `cos_nodes::Int`
- `cos_grid::Union{Nothing,Vector{Float64}}`
- `cos_w::Union{Nothing,Vector{Float64}}`

默认值：沿用模块当前默认（`DEFAULT_*_NODES/WEIGHTS` 与 `p_max=10.0`）。

### 校验规则（先从结构体层做起）

- `p_grid` 与 `p_w` 必须成对出现，且长度一致。
- `cos_grid` 与 `cos_w` 必须成对出现，且长度一致。

### API 改造方式

- `shear_viscosity/electric_conductivity/bulk_viscosity_isentropic/transport_coefficients` 新增 keyword：
   - `config::TransportIntegrationConfig=DEFAULT_TRANSPORT_CONFIG`
- 原有关键字参数仍保留，但其默认值改为从 `config` 读取，从而自动实现“显式关键字覆盖 config”。

### 文档与测试同步

- API 文档 `docs/api/relaxtime/transport/TransportCoefficients.md` 增加 `config` 参数说明与示例。
- 单元测试增加：
   - `config` 用法与旧关键字用法等价性测试。
   - “显式关键字覆盖 config”的优先级测试。
   - 结构体参数校验（例如只给 `p_grid` 不给 `p_w` 应报错）。

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
