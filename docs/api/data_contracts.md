# 数据契约定义

## 概述

本文档定义QCD模型库中所有数据结构的契约，包括参数格式、转换规则、验证schema和版本兼容性策略。

**目标**：确保新旧API之间、不同模型之间的数据交换清晰、类型安全且可验证。

---

## 0. Models 子系统（阶段 0）契约

本节描述当前仓库 `src/models` 下 **models 子系统** 的“最小稳定入口”，用于支撑 legacy → models 的渐进迁移。

### 0.1 统一入口

- `Models.solve_gap(model, T, mu_vec; kwargs...) -> x_state`
- `Models.omega(model, x_state, T, mu_vec; kwargs...) -> Real`
- `Models.omega_components(model, x_state, T, mu_vec; kwargs...) -> NamedTuple`
- （可选扩展）`Models.number_densities(model, x_state, T, mu_vec; kwargs...) -> NamedTuple`

后端约定（迁移期）：

- `solver_backend` 继续支持 `:legacy | :models`；默认保持兼容（不强制切换）。
- `src/models/pnjl_physics/core/EquilibriumFacade.jl` 支持 `solver_backend=:auto`（按 `thermo_backend` 推导），但默认仍为 `:legacy` 以避免行为突变。
- `PNJLModel.solve_gap` 在 `solver_backend=:models` 下支持失败时回退 legacy 的受控策略（对称化学势场景）。

### 0.2 `x_state`（平均场状态）

推荐的规范表示是：

```julia
Models.MeanFieldState(phi::SVector{3}, Phi::Real, PhiBar::Real)
```

字段约定：

- `phi = (φu, φd, φs)`：三味凝聚参数
- `Phi` / `PhiBar`：Polyakov loop 变量（NJL 类通常可取 1 并被忽略）

兼容输入（会被规范化为 `MeanFieldState`）：

- `MeanFieldState` 本身
- `AbstractVector`：长度为 3（只含 φ）或长度 ≥ 5（φ + Φ + Φbar）
- `NamedTuple`：必须包含 `:φ` 或 `:phi`；`Φ/Φbar` 缺省时默认 1

配套工具函数：

- `Models.meanfield_state(x_state)`：规范化为 `MeanFieldState`
- `Models.state_vector(x_state)`：转为 5 维向量 `(φu, φd, φs, Φ, Φbar)`

### 0.3 `mu_vec`（化学势向量）

内部约定为三味化学势向量 `(μu, μd, μs)`。

兼容输入：

- `Real`：按对称情形扩展为 `(μ, μ, μ)`
- `AbstractVector`：必须长度为 3

配套工具函数：

- `Models.normalize_mu_vec(mu_vec)`：统一转换为 `SVector{3}`

### 0.4 rPNJL 最小契约（阶段 6 MVP）

本节定义“先可运行、后增强”的 rPNJL 最小交付口径，目标是以最小改动打通 models 主链，不在 MVP 阶段承诺完整物理细节。

最小入口：

- `Models.create_model(:RPNJL; profile="default", physics_profile="default")`
- `Models.solve_gap(model::AbstractPNJLModel, T, mu_vec; kwargs...)`
- `Models.omega(model::AbstractQCDModel, x_state, T, mu_vec; kwargs...)`
- `Models.omega_components(model::AbstractQCDModel, x_state, T, mu_vec; kwargs...)`

状态与输入语义（与 PNJL 家族保持一致）：

- `x_state`：兼容 `MeanFieldState` / 长度 3 或 ≥5 向量 / 含 `phi` 的 `NamedTuple`。
- `mu_vec`：兼容 `Real`（扩展为三味同值）或长度为 3 的向量。
- `gap_state_dim`：rPNJL 在 MVP 阶段按 PNJL 家族约束为 5。

参数契约（MVP）：

- 继承 PNJL 参数集（`G, K, Λ, m_u0, m_d0, m_s0, T0, a0..`）。
- 预留 rPNJL 扩展参数位：`g1`, `g2`, `kappa`。
- 在“八夸克项/Vandermonde”未启用前，允许 `g1/g2/kappa` 仅保留为配置字段，不改变现有数值路径。

兼容边界：

- MVP 不改变 NJL/PNJL 既有行为与默认结果。
- 容差策略仅允许存在于测试断言层，不作为被测函数 keyword 暴露。
- 若后续启用 `g1/g2/kappa` 进入主计算路径，需在 active 文档记录“公式来源 + 参数来源 + 新旧对比 smoke”。

### 0.5 rPNJL 公式映射契约（阶段 7）

来源文档：`docs/reference/formula/models/rpnjl/rPNJL_core.md`

阶段 7 约定把公式映射到以下实现入口：

- 式(3.31)（能隙方程，含 `g1/g2`）→ `calculate_mass_vec(::RPNJLModel, φ)`
- 式(3.30)（巨热力学势中的凝聚/八夸克势项）→ `calculate_chiral(::RPNJLModel, φ)`
- 式(3.27)–(3.29)（Polyakov 势 + Vandermonde）→ `polyakov_potential(::RPNJLModel, Φ, Φbar, T)`

参数与单位契约：

- `g1/g2`：配置文件使用 `MeV^-8`，在模型加载时统一换算为 `fm^8`（乘以 `hbarc^8`）。
- `kappa`：无量纲，默认来源于 `config/models/rpnjl/<profile>.toml` 的 `[rpnjl].kappa`。
- `T0,a0,a1,a2,b3,b4`：rPNJL profile 可覆盖 PNJL 同名参数；覆盖行为只作用于 `RPNJLModel`。

退化与兼容约束：

- `use_rpnjl_extensions=false` 时，`RPNJLModel` 应退化到 PNJL 参数面（用于 bridge sanity）。
- `use_rpnjl_extensions=true` 时，允许与 PNJL 基线出现受控偏离，但需通过固定点 smoke 追踪。
- 阶段 7 期间禁止修改测试外部容差 API；所有容差仅在测试断言层设置。

### 0.6 FullServer API 契约（HTTP，phase2-v1）

本节定义前端任务中心可依赖的稳定契约基线（版本：`phase2-v1`）。

#### 0.6.1 路由与调用边界

- **同步短任务**（请求完成即返回结果）
  - `POST /compute`
  - `POST /api/modules/pnjl-gap/run`
- **异步扫描任务**（创建任务后轮询）
  - `POST /api/modules/pnjl-scan/jobs`（创建）
  - `GET /api/modules/pnjl-scan/jobs/{job_id}`（状态）
  - `GET /api/modules/pnjl-scan/jobs/{job_id}/result`（结果）
  - `POST /api/modules/pnjl-scan/jobs/{job_id}/cancel`（取消）
- **模块发现**
  - `GET /api/modules`

边界约束：

- 前端不直接调用 `run_*.jl`，仅通过上述 HTTP 契约访问能力。
- 扫描任务统一走异步 jobs 路由；不提供同步扫描接口。

#### 0.6.2 `POST /api/modules/pnjl-gap/run`（同步）

请求：

- `params.T_mev`（必填，兼容别名 `t_mev`）
- `params.mu_mev`（`FixedMu` 模式必填，兼容别名 `mu`）
- `params.rho_target`（可选；提供后进入 `FixedRho` 模式）
- `params.xi`（可选，默认 `0.0`）
- `params.p_num`（可选，默认 `24`）
- `params.t_num`（可选，默认 `12`）
- `params.allow_seed_fallback`（可选，默认 `true`）

成功响应（`200`）：

- `status="ok"`
- `result`：`converged/omega/pressure/rho_norm/entropy/energy/iterations/residual_norm/xi/seed_fallback_used/x_state/mu_vec/masses`

#### 0.6.3 `POST /api/modules/pnjl-scan/jobs`（异步创建）

请求：

- `kind`（必填，`tmu | trho`）
- `params`（可选字典）
  - `mode`（可选，`scan | point`，默认 `scan`）
  - 公共：`xi | xi_values | xi_grid`（三选一，默认 `xi=0.0`）
  - `tmu + scan`：`T_values`、`mu_values`
  - `trho + scan`：`T_values`、`rho_values`
  - `tmu + point`：`T_mev`、`mu_mev`
  - `trho + point`：`T_mev`、`rho_value`
  - 可选治理参数：`max_retries`（默认 `0`，最大 `3`）、`timeout_seconds`

成功响应（`202`）：

- `status="accepted"`
- `job_id`
- `kind`
- `status_url`
- `result_url`
- `queue.position/max_running/max_pending`
- `idempotency.key/replayed/conflict`
- `diagnostics.job_id/kind/job_status`

#### 0.6.4 `GET /api/modules/pnjl-scan/jobs/{job_id}`（状态）

成功响应（`200`）：

- `status="ok"`
- `job_id/kind/job_status`
- `created_at/started_at/ended_at`
- `progress.total/completed/percent`
- `error`（仅失败态暴露用户可读摘要）
- `queue.position/queued/running/max_running/max_pending`
- `policy`（任务策略快照）
- `events`（结构化生命周期事件）
- `governance`（保留策略与清理计数）
- `metrics`（terminal + duration_buckets + queue）
- `diagnostics.job_id/kind/job_status`

状态机：

- 非终态：`queued`、`running`
- 终态：`succeeded`、`failed`、`cancelled`

#### 0.6.5 `GET /api/modules/pnjl-scan/jobs/{job_id}/result`（结果）

- 任务成功（`200`）：返回 `status="ok"`、`job_id`、`job_status`、`result(output_path,stats)`、`metadata(output_exists,output_mtime)`、`diagnostics`
- 任务未成功（`409`）：返回 `JOB_NOT_SUCCEEDED`

#### 0.6.6 `POST /api/modules/pnjl-scan/jobs/{job_id}/cancel`（取消）

- 可取消状态：`queued | running`
- 已终态或不可取消：返回 `409` + `JOB_NOT_CANCELLABLE`
- 成功取消：`200`，`job_status="cancelled"`

#### 0.6.7 统一错误响应最小契约

适用链路（当前）：

- `/compute`
- `/api/modules/pnjl-gap/run`
- `/api/modules/pnjl-scan/jobs*`

对外错误 payload 最小字段：

- `status`：固定为 `"error"`
- `error_code`：稳定错误码
- `error`：面向调用方的摘要信息（不包含内部栈）
- `message_id`：每次错误唯一标识，用于检索日志与追踪

可选字段：

- `diagnostics`、`job_id`、`job_status`、`policy`、`idempotency`

当前稳定错误码（phase2-v1）：

- `INVALID_INPUT`
- `INVALID_REQUEST`
- `QUEUE_FULL`
- `JOB_NOT_FOUND`
- `JOB_NOT_SUCCEEDED`
- `JOB_NOT_CANCELLABLE`
- `IDEMPOTENCY_KEY_CONFLICT`
- `PNJL_SINGLE_POINT_FAILED`
- `COMPUTATION_ERROR`

安全约束：

- 不向外返回 `backtrace` / `catch_backtrace()` / 内部异常全文。
- 内部诊断可在服务端日志或内存状态中保留，但不直接透传到公共响应。

#### 0.6.8 最小契约示例（可复现）

最小请求（创建扫描任务）：

```json
{
  "kind": "tmu",
  "params": {
    "T_values": [150.0],
    "mu_values": [0.0],
    "xi": 0.0,
    "max_retries": 0
  }
}
```

错误请求示例（非法 kind）：

```json
{
  "status": "error",
  "error_code": "INVALID_REQUEST",
  "error": "Missing/invalid kind; expected tmu or trho",
  "message_id": "<uuid>",
  "diagnostics": {
    "job_status": "rejected"
  }
}
```

成功响应示例（创建成功）：

```json
{
  "status": "accepted",
  "job_id": "<uuid>",
  "kind": "tmu",
  "status_url": "/api/modules/pnjl-scan/jobs/<uuid>",
  "result_url": "/api/modules/pnjl-scan/jobs/<uuid>/result",
  "queue": {
    "position": 1,
    "max_running": 2,
    "max_pending": 32
  },
  "idempotency": {
    "key": null,
    "replayed": false,
    "conflict": false
  },
  "diagnostics": {
    "job_id": "<uuid>",
    "kind": "tmu",
    "job_status": "queued"
  }
}
```

---

## 1. 参数类型系统

### 1.1 核心参数类型

QCD模型库支持三种参数表示格式：

```julia
# 1. 结构体（推荐，类型安全）
struct QuarkParams
    m_u::Float64
    m_d::Float64
    m_s::Float64
    mu_u::Float64
    mu_d::Float64
    mu_s::Float64
end

struct ThermoParams
    T::Float64
    Φ::Float64
    Φbar::Float64
    xi::Float64
end

# 2. NamedTuple（轻量，不可变）
quark_params = (m_u=0.3, m_d=0.3, m_s=0.5, mu_u=0.0, mu_d=0.0, mu_s=0.0)
thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, xi=0.0)

# 3. Dict{Symbol, Any}（灵活，用于配置）
model_params = Dict(
    :G => 10.08,
    :K => -39.0,
    :T0 => 0.19,
    :Lambda => 0.6
)
```

### 1.2 参数格式优先级

当函数接受多种格式时，内部处理优先级：

1. **Struct** → 直接使用（零开销）
2. **NamedTuple** → 转换为Struct（编译时优化）
3. **Dict** → 运行时转换（有验证开销）

**实现模式**：
```julia
function my_function(quark_params::Union{QuarkParams, NamedTuple, Dict})
    # 内部统一转换为NamedTuple（最小公分母）
    qp_nt = _normalize_quark_params(quark_params)
    # ... 使用 qp_nt.m_u, qp_nt.mu_u 等
end
```

---

## 2. 参数字段定义

### 2.1 QuarkParams

| 字段 | 类型 | 单位 | 范围 | 描述 | 必需 |
|------|------|------|------|------|------|
| `m_u` | Float64 | fm⁻¹ | [0, 2.0] | Up quark constituent mass | ✅ |
| `m_d` | Float64 | fm⁻¹ | [0, 2.0] | Down quark constituent mass | ✅ |
| `m_s` | Float64 | fm⁻¹ | [0, 2.0] | Strange quark constituent mass | ✅ |
| `mu_u` | Float64 | fm⁻¹ | [0, 2.0] | Up quark chemical potential | ✅ |
| `mu_d` | Float64 | fm⁻¹ | [0, 2.0] | Down quark chemical potential | ✅ |
| `mu_s` | Float64 | fm⁻¹ | [0, 2.0] | Strange quark chemical potential | ✅ |
| `A` | Float64 | fm⁻² | [0, Inf] | Meson mass parameter (optional) | ❌ |

**注意**：
- `A` 字段是可选的，用于介子质量计算
- 如果缺失，可通过 `ensure_quark_params_has_A` 自动计算

### 2.2 ThermoParams

| 字段 | 类型 | 单位 | 范围 | 描述 | 必需 |
|------|------|------|------|------|------|
| `T` | Float64 | fm⁻¹ | [0, 1.0] | Temperature | ✅ |
| `Φ` | Float64 | - | [0, 1.0] | Polyakov loop | ✅ |
| `Φbar` | Float64 | - | [0, 1.0] | Conjugate Polyakov loop | ✅ |
| `xi` | Float64 | - | [-1.0, 1.0] | Anisotropy parameter | ✅ |

**物理约束**：
- `Φ` 和 `Φbar` 通常满足 `Φ ≈ Φbar`（实际计算中可能略有差异）
- `xi = 0` 表示各向同性
- `xi > 0` 表示纵向压缩，`xi < 0` 表示纵向拉伸

### 2.3 ModelParams (模型配置)

不同模型有不同的配置参数：

**PNJL模型**：
| 字段 | 类型 | 单位 | 默认值 | 描述 |
|------|------|------|--------|------|
| `G` | Float64 | fm² | 10.08 | Four-quark coupling |
| `K` | Float64 | fm⁵ | -39.0 | Six-quark coupling |
| `Lambda` | Float64 | fm⁻¹ | 0.6 | UV cutoff |
| `m_u0` | Float64 | fm⁻¹ | 0.0056 | Current up quark mass |
| `m_d0` | Float64 | fm⁻¹ | 0.0056 | Current down quark mass |
| `m_s0` | Float64 | fm⁻¹ | 0.135 | Current strange quark mass |
| `T0` | Float64 | fm⁻¹ | 0.19 | Polyakov potential scale |
| `a0`, `a1`, `a2`, `b3` | Float64 | - | - | Polyakov potential coefficients |

**rPNJL模型**（继承PNJL，额外参数）：
| 字段 | 类型 | 单位 | 默认值 | 描述 |
|------|------|------|--------|------|
| `g1` | Float64 | fm⁸ | TBD | Eight-quark coupling 1 |
| `g2` | Float64 | fm⁸ | TBD | Eight-quark coupling 2 |
| `kappa` | Float64 | - | TBD | Vandermonde coefficient |

---

## 3. 参数转换规则

### 3.1 Struct ↔ NamedTuple

```julia
# Struct → NamedTuple
function as_namedtuple(qp::QuarkParams)
    return (
        m_u = qp.m_u,
        m_d = qp.m_d,
        m_s = qp.m_s,
        mu_u = qp.mu_u,
        mu_d = qp.mu_d,
        mu_s = qp.mu_s,
        A = hasproperty(qp, :A) ? qp.A : missing
    )
end

# NamedTuple → Struct
function QuarkParams(nt::NamedTuple)
    # 必需字段检查
    required = [:m_u, :m_d, :m_s, :mu_u, :mu_d, :mu_s]
    for field in required
        if !haskey(nt, field)
            throw(ArgumentError("Missing required field: $field"))
        end
    end
    
    # 构造（忽略额外字段）
    return QuarkParams(
        nt.m_u, nt.m_d, nt.m_s,
        nt.mu_u, nt.mu_d, nt.mu_s
    )
end
```

### 3.2 Dict ↔ Struct

```julia
# Dict → Struct（带验证）
function QuarkParams(d::Dict{Symbol, <:Any})
    # 验证必需字段
    required = [:m_u, :m_d, :m_s, :mu_u, :mu_d, :mu_s]
    for field in required
        if !haskey(d, field)
            throw(ModelParameterError(field, missing, "Required field missing"))
        end
    end
    
    # 类型转换和验证
    m_u = Float64(d[:m_u])
    validate_parameter(:m_u, m_u, QUARK_CONSTRAINTS)
    # ... 其他字段类似
    
    return QuarkParams(m_u, m_d, m_s, mu_u, mu_d, mu_s)
end

# Struct → Dict
function as_dict(qp::QuarkParams)
    d = Dict{Symbol, Any}(
        :m_u => qp.m_u,
        :m_d => qp.m_d,
        :m_s => qp.m_s,
        :mu_u => qp.mu_u,
        :mu_d => qp.mu_d,
        :mu_s => qp.mu_s
    )
    if hasproperty(qp, :A)
        d[:A] = qp.A
    end
    return d
end
```

### 3.3 单位转换

**MeV ↔ fm⁻¹**：
```julia
const ħc_MeV_fm = 197.3269804  # MeV·fm

# MeV → fm⁻¹
to_natural_units(value_MeV) = value_MeV / ħc_MeV_fm

# fm⁻¹ → MeV
to_MeV(value_fm_inv) = value_fm_inv * ħc_MeV_fm

# 批量转换
function convert_params_to_natural(params_MeV::NamedTuple)
    return (
        T = to_natural_units(params_MeV.T),
        mu_u = to_natural_units(params_MeV.mu_u),
        # ... 其他需要转换的字段
        Φ = params_MeV.Φ,  # 无量纲，不转换
        Φbar = params_MeV.Φbar,
        xi = params_MeV.xi
    )
end
```

---

## 4. 验证Schema

### 4.1 参数约束定义

```julia
# src/models/base/validation.jl

"""参数约束结构"""
struct ParameterConstraint
    min::Float64
    max::Float64
    unit::String
    description::String
    validator::Union{Function, Nothing}  # 自定义验证函数
end

"""QuarkParams约束"""
const QUARK_CONSTRAINTS = Dict{Symbol, ParameterConstraint}(
    :m_u => ParameterConstraint(0.0, 2.0, "fm⁻¹", "Up quark mass", nothing),
    :m_d => ParameterConstraint(0.0, 2.0, "fm⁻¹", "Down quark mass", nothing),
    :m_s => ParameterConstraint(0.0, 2.0, "fm⁻¹", "Strange quark mass", nothing),
    :mu_u => ParameterConstraint(0.0, 2.0, "fm⁻¹", "Up quark chemical potential", nothing),
    :mu_d => ParameterConstraint(0.0, 2.0, "fm⁻¹", "Down quark chemical potential", nothing),
    :mu_s => ParameterConstraint(0.0, 2.0, "fm⁻¹", "Strange quark chemical potential", nothing),
    :A => ParameterConstraint(0.0, Inf, "fm⁻²", "Meson mass parameter", nothing),
)

"""ThermoParams约束"""
const THERMO_CONSTRAINTS = Dict{Symbol, ParameterConstraint}(
    :T => ParameterConstraint(0.0, 1.0, "fm⁻¹", "Temperature", nothing),
    :Φ => ParameterConstraint(0.0, 1.0, "-", "Polyakov loop", nothing),
    :Φbar => ParameterConstraint(0.0, 1.0, "-", "Conjugate Polyakov loop", nothing),
    :xi => ParameterConstraint(-1.0, 1.0, "-", "Anisotropy parameter", nothing),
)
```

### 4.2 验证函数

```julia
"""验证单个参数"""
function validate_parameter(name::Symbol, value, constraints::Dict)
    c = get(constraints, name, nothing)
    if c === nothing
        @warn "No constraints defined for parameter" name
        return true
    end
    
    # 范围检查
    if !(c.min ≤ value ≤ c.max)
        throw(ModelParameterError(
            name, value,
            "Value must be in range [$(c.min), $(c.max)] $(c.unit)"
        ))
    end
    
    # 自定义验证
    if c.validator !== nothing && !c.validator(value)
        throw(ModelParameterError(
            name, value,
            "Failed custom validation"
        ))
    end
    
    return true
end

"""验证QuarkParams"""
function validate_quark_params(qp::Union{QuarkParams, NamedTuple, Dict})
    qp_nt = _normalize_quark_params(qp)
    for field in [:m_u, :m_d, :m_s, :mu_u, :mu_d, :mu_s]
        validate_parameter(field, getproperty(qp_nt, field), QUARK_CONSTRAINTS)
    end
    if haskey(qp_nt, :A) && qp_nt.A !== missing
        validate_parameter(:A, qp_nt.A, QUARK_CONSTRAINTS)
    end
    return true
end

"""验证ThermoParams"""
function validate_thermo_params(tp::Union{ThermoParams, NamedTuple, Dict})
    tp_nt = _normalize_thermo_params(tp)
    for field in [:T, :Φ, :Φbar, :xi]
        validate_parameter(field, getproperty(tp_nt, field), THERMO_CONSTRAINTS)
    end
    return true
end
```

---

## 5. 版本兼容性策略

### 5.1 参数结构版本化

```julia
"""参数结构版本标记"""
struct VersionedParams{V, T}
    version::Symbol  # :v1, :v2, etc.
    data::T
end

# 当前版本
const CURRENT_QUARK_PARAMS_VERSION = :v1
const CURRENT_THERMO_PARAMS_VERSION = :v1

# 版本转换
function upgrade_quark_params(vp::VersionedParams{:v1, QuarkParams})
    # v1 → v2 (假设v2增加了新字段)
    return VersionedParams(:v2, QuarkParamsV2(
        vp.data.m_u, vp.data.m_d, vp.data.m_s,
        vp.data.mu_u, vp.data.mu_d, vp.data.mu_s,
        0.0  # 新字段的默认值
    ))
end
```

### 5.2 向后兼容性保证

**规则**：
1. **不删除字段**：旧字段必须保留或提供默认值
2. **新增字段可选**：新字段必须有合理的默认值
3. **转换函数**：提供旧版本→新版本的转换函数
4. **弃用警告**：使用`@deprecate`标记过时的接口

**示例**：
```julia
# 旧接口（保留，但标记为弃用）
@deprecate old_function(x, y) new_function(x, y, z=default_z)

# 新接口
function new_function(x, y, z=default_z)
    # ... 实现
end
```

### 5.3 序列化格式

**推荐格式**：JSON（人类可读，跨语言）

```julia
using JSON3

# 序列化
function save_params(filename::String, qp::QuarkParams, tp::ThermoParams)
    data = Dict(
        "version" => "1.0",
        "quark_params" => as_dict(qp),
        "thermo_params" => as_dict(tp),
        "timestamp" => string(now())
    )
    open(filename, "w") do io
        JSON3.pretty(io, data)
    end
end

# 反序列化（带版本检查）
function load_params(filename::String)
    data = JSON3.read(read(filename, String))
    
    version = get(data, "version", "unknown")
    if version != "1.0"
        @warn "Loading parameters from different version" file_version=version current_version="1.0"
    end
    
    qp = QuarkParams(Dict{Symbol, Any}(Symbol(k) => v for (k, v) in data["quark_params"]))
    tp = ThermoParams(Dict{Symbol, Any}(Symbol(k) => v for (k, v) in data["thermo_params"]))
    
    return qp, tp
end
```

---

## 6. 新旧API数据契约

### 6.1 旧API（兼容层）

**旧接口签名**：
```julia
# 旧：使用位置参数
function calculate_mass_vec(φ_u, φ_d, φ_s)
    # ...
end

# 旧：使用全局常量
function calculate_U(T, Φ, Φbar)
    # 依赖全局 T0, a0, a1, a2, b3
end
```

**兼容层实现**：
```julia
# src/compatibility/legacy_api.jl

"""兼容旧API：位置参数 → SVector"""
function calculate_mass_vec(φ_u::Real, φ_d::Real, φ_s::Real)
    φ = SVector{3, Float64}(φ_u, φ_d, φ_s)
    model = get_default_pnjl_model()  # 使用默认模型
    return calculate_mass_vec(model, φ)
end

"""兼容旧API：使用全局常量"""
function calculate_U(T::Real, Φ::Real, Φbar::Real)
    model = get_default_pnjl_model()
    return polyakov_potential(model, Φ, Φbar, T)
end
```

### 6.2 新API（推荐）

**新接口签名**：
```julia
# 新：使用模型实例和结构化参数
function calculate_mass_vec(model::AbstractQCDModel, φ::SVector{3, T}) where {T}
    # ...
end

function polyakov_potential(model::AbstractQCDModel, Φ, Φbar, T)
    # 从model.params读取T0, a0等
end
```

### 6.3 迁移路径

**阶段1**：兼容层（当前）
- 旧代码无需修改即可运行
- 新代码推荐使用新API

**阶段2**：弃用警告（下一版本）
- 旧API标记为`@deprecate`
- 文档说明迁移方法

**阶段3**：移除旧API（未来版本）
- 仅保留新API
- 提供迁移工具脚本

---

## 7. 测试要求

### 7.1 转换测试

```julia
@testset "Parameter conversion" begin
    # Struct ↔ NamedTuple
    qp_struct = QuarkParams(0.3, 0.3, 0.5, 0.1, 0.1, 0.1)
    qp_nt = as_namedtuple(qp_struct)
    qp_back = QuarkParams(qp_nt)
    @test qp_back == qp_struct
    
    # Dict ↔ Struct
    qp_dict = as_dict(qp_struct)
    qp_from_dict = QuarkParams(qp_dict)
    @test qp_from_dict == qp_struct
end
```

### 7.2 验证测试

```julia
@testset "Parameter validation" begin
    # 有效参数
    @test_nowarn validate_quark_params(QuarkParams(0.3, 0.3, 0.5, 0.1, 0.1, 0.1))
    
    # 无效参数
    @test_throws ModelParameterError validate_quark_params(
        QuarkParams(-0.1, 0.3, 0.5, 0.1, 0.1, 0.1)  # 负质量
    )
    @test_throws ModelParameterError validate_thermo_params(
        (T=0.15, Φ=1.5, Φbar=0.5, xi=0.0)  # Φ > 1
    )
end
```

### 7.3 兼容性测试

```julia
@testset "API compatibility" begin
    # 旧API应该仍然工作
    result_old = calculate_mass_vec(0.1, 0.1, 0.2)
    
    # 新API应该产生相同结果
    model = create_model(:PNJL)
    φ = SVector{3}(0.1, 0.1, 0.2)
    result_new = calculate_mass_vec(model, φ)
    
    @test result_old ≈ result_new
end
```

---

## 8. 使用示例

### 8.1 推荐用法（新API）

```julia
using QCDModels

# 创建参数
quark_params = QuarkParams(
    m_u = 0.3,
    m_d = 0.3,
    m_s = 0.5,
    mu_u = 0.1,
    mu_d = 0.1,
    mu_s = 0.1
)

thermo_params = ThermoParams(
    T = 0.15,
    Φ = 0.5,
    Φbar = 0.5,
    xi = 0.0
)

# 创建模型
model = create_model(:PNJL, params=Dict(
    :G => 10.08,
    :K => -39.0,
    :Lambda => 0.6
))

# 计算
masses = calculate_mass_vec(model, SVector(0.1, 0.1, 0.2))
U = polyakov_potential(model, thermo_params.Φ, thermo_params.Φbar, thermo_params.T)
```

### 8.2 灵活用法（混合格式）

```julia
# 使用NamedTuple（轻量）
quark_params = (m_u=0.3, m_d=0.3, m_s=0.5, mu_u=0.1, mu_d=0.1, mu_s=0.1)

# 函数自动处理
result = some_function(quark_params)  # 内部会转换
```

### 8.3 配置文件用法

```julia
# 从JSON加载
quark_params, thermo_params = load_params("config.json")

# 验证
validate_quark_params(quark_params)
validate_thermo_params(thermo_params)

# 使用
model = create_model(:PNJL)
result = calculate_thermodynamics(model, quark_params, thermo_params)
```

---

## 9. 检查清单

实现新功能时，请确认：

- [ ] 参数类型支持Struct、NamedTuple、Dict三种格式
- [ ] 提供了格式之间的转换函数
- [ ] 定义了参数约束并实现了验证
- [ ] 单位转换正确（MeV ↔ fm⁻¹）
- [ ] 新增字段有默认值（向后兼容）
- [ ] 编写了转换和验证的测试
- [ ] 文档说明了参数格式和约束
- [ ] 兼容层正确处理旧API调用

---

## 10. 参考资料

- 参数类型API文档：`docs/api/PARAMETER_TYPES_API.md`
- 错误处理指南：`docs/guides/error_handling.md`
- 迁移指南：`docs/guides/PARAMETER_STRUCT_MIGRATION.md`
