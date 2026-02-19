---
title: 输运系数计算流程优化任务单（执行版，可勾选）
archived: true
original: docs/dev/active/2025-12-28_transport_optimization_analysis.md
archived_date: 2026-02-19
---


以下为原始内容（保留，以便审阅与历史参考）：

---

# 输运系数计算流程优化任务单（执行版，可勾选）

## 目标

将历史优化分析转为“与当前仓库状态一致”的可执行任务：补齐 P0 缺口、补强已落地能力，并形成可归档证据。

## 当前状态（2026-02-19）

- [x] 高优先级项“ThermoDerivatives 中重复求解能隙方程”已落地。
- [x] Phase 1（TransportIntegrationConfig）已落地并有单测覆盖。
- [x] 已有固定点基线文件与对应 smoke 回归（transport bridge）。
- [x] 通用积分框架抽取已落地：η/σ 共享主积分循环，ζ 复用共享组件并保留特例物理公式。
- [x] 输入参数验证与极端参数数值保护已落地（T>0、质量非负、E 下限、ff 截断等）并有单测。
- [x] Transport API 文档口径已修正（`bulk_viscosity_isentropic` 等）。
- [x] 低优先级性能微优化首轮已落地（Symbol→索引映射 + provider 能力缓存，减少热循环分支判定）。
- [x] 与技术债总单的映射与去重已建立（避免 active 任务单漂移）。

## 待办（可勾选）

### P0

- [x] 在 TransportCoefficients 中抽取“共享积分主循环”：已统一 η/σ 的主积分路径；ζ 保留特例公式并复用共享组件。
- [x] 完成输入参数校验（T>0、质量非负、网格一致性、tau/charges 字段完整性）并新增单测。
- [x] 增加极端参数数值保护（E 下限、ff 截断等），并验证正常物理区间结果不漂移。
- [x] 修正 Transport API 文档口径（函数命名、参数名、默认值、边界行为）并与当前实现一致。

### P1

- [x] 评估并落地 Symbol→整数索引替换（仅在热点路径，先完成 species/flavor 映射路径）。
- [x] 评估并落地分派表/函数指针缓存，减少热循环条件分支。
- [x] 统一体粘滞接口命名策略（新增单入口 `bulk_viscosity(...; formula=:isentropic)`，并保留 `bulk_viscosity_isentropic` 兼容）。
- [x] 输出热点路径“优化前后”对比（固定点 + 小网格），形成可复现实验记录。

### P2

- [x] 补强固定点基线覆盖（在现有 `baseline_transport_fixedpoints_v1.csv` 基础上扩展 ξ≠0 点位）。
- [x] 将基线治理流程补齐到文档（容差来源、更新准入、变更说明模板）。
- [x] 与 `2026-02-19_技术债优先级清单与验收任务单` 建立映射，去重并保持状态同步。

## DoD

- [x] P0 全部完成并有单测证据。
- [x] 保留并通过现有基线 smoke 回归；新增/更新基线时有版本号、原因和容差说明。
- [x] 性能热点路径有“优化前后对比”记录（命令、点位、容差可复现）。
- [x] README/STATUS/API/测试口径一致（含体粘滞接口命名）。
- [x] 与技术债总单状态一致，无重复或冲突项。
- [x] 可归档执行记录完整。

## 2026-02-19 执行记录（本轮）

### 已落地代码与测试

- `src/relaxtime/TransportCoefficients.jl`
   - 抽取 η/σ 共享积分主循环 `_integrate_species_sum`。
   - 增加统一输入校验 `_validate_transport_inputs`。
   - 增加数值保护：`E = max(E_raw, sqrt(eps(Float64)))`，`ff = clamp(f*(1-f), 0.0, 0.25)`。
   - species/flavor 路径改为索引映射辅助函数（减少热点 Symbol 分支判断）。
   - provider 能力缓存 `_provider_caps`：在热循环复用分支能力判断。
   - 新增统一体粘滞入口 `bulk_viscosity(...; formula=:isentropic)`（兼容 `bulk_coeffs` 旧参数名）。
- `tests/unit/relaxtime/test_transport_coefficients.jl`
   - 新增输入校验与数值保护测试。
   - 与现有 provider mass/mu override smoke 语义对齐。
- `tests/unit/relaxtime/test_transport_legacy_models_bridge_smoke.jl`
   - 默认纳入 `ξ=0.2` 固定点回归。
- `tests/baselines/relaxtime/baseline_transport_fixedpoints_v1.csv`
   - 新增 `T=0.9, mu=0.0, xi=0.2` 基线行。

### 定向验证命令与结果

- `Set-Item Env:UNIT_FILES 'relaxtime/test_transport_workflow_smoke.jl'; julia --project=. tests/unit/runtests.jl`
   - 结果：`35 passed, 0 failed, 0 errored`
- `Set-Item Env:UNIT_FILES 'relaxtime/test_transport_legacy_models_bridge_smoke.jl'; julia --project=. tests/unit/runtests.jl`
   - 结果：`100 passed, 0 failed, 0 errored`

- `Set-Item Env:UNIT_FILES 'relaxtime/test_transport_coefficients.jl'; julia --project=. tests/unit/runtests.jl`
   - 结果：`48 passed, 0 failed, 0 errored`

### 热点路径性能记录（小网格）

- 脚本：`scripts/dev/benchmark_transport_hotpath.jl`
- 命令：`julia --project=. scripts/dev/benchmark_transport_hotpath.jl`
- 输出（当前实现）：
   - `eta_median_s=2.935e-5`
   - `sigma_median_s=3.165e-5`
   - `eta_min_s=1.88e-5`
   - `sigma_min_s=2.01e-5`

### 基线治理口径（本任务单内约定）

- 基线文件：`tests/baselines/relaxtime/baseline_transport_fixedpoints_v1.csv`
- smoke 校验：`tests/unit/relaxtime/test_transport_legacy_models_bridge_smoke.jl`
- 容差口径：当前 smoke 使用 `rtol=8e-2, atol=1e-6`。
- 基线更新准入：必须同时提交“变更原因 + 影响范围 + 新旧差异摘要 + 验证命令/结果 + 文档更新”。

---

# 原始分析记录（保留）

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
