# TransportCoefficients 模块 API 文档

本模块实现弛豫时间近似（RTA）下的夸克物质输运系数计算：剪切粘滞系数 $\eta$、体粘滞系数 $\zeta$、电导率 $\sigma$。

## 单位约定
- 自然单位：$\hbar=c=k_B=1$
- 温度/化学势/质量/动量：fm⁻¹
- 弛豫时间 $\tau$：fm
- 电荷：自然单位（$e = \sqrt{4\pi\alpha}$，无量纲）

## 输出单位

| 输运系数 | 量纲 | 单位 |
|---------|------|------|
| 剪切粘滞系数 η | $[E]^3$ | fm⁻³ |
| 体粘滞系数 ζ | $[E]^3$ | fm⁻³ |
| 电导率 σ | $[E]^1$ | fm⁻¹ |

转换到 MeV 单位：$1 \text{ fm}^{-1} = \hbar c \approx 197.327 \text{ MeV}$

## 物理公式

详见 `docs/reference/formula/relaxtime/transport/TransportCoefficients_FromRelaxationTime.md`。

### 剪切粘滞系数 η

**各向同性 (ξ=0)**：
$$
\eta = \frac{1}{15T} \cdot \frac{1}{2\pi^2} \sum_a \int_0^\infty dp \, \frac{p^6}{E_a^2} \left[ d_a \tau_a f_a^0 (1 - f_a^0) + d_{\bar{a}} \tau_{\bar{a}} f_{\bar{a}}^0 (1 - f_{\bar{a}}^0) \right]
$$

**各向异性 (ξ≠0)**：
$$
\eta = \frac{1}{15T} \cdot \frac{1}{4\pi^2} \sum_a \int_0^\infty dp \int_{-1}^{1} d(\cos\theta) \, \frac{p^6}{E_a^2} \left[ d_a \tau_a f_a (1 - f_a) + ... \right]
$$

### 电导率 σ

**各向同性 (ξ=0)**：
$$
\sigma = \frac{1}{3T} \cdot \frac{1}{2\pi^2} \sum_a \int_0^\infty dp \, \frac{p^4 q_a^2}{E_a^2} \left[ d_a \tau_a f_a^0 (1 - f_a^0) + d_{\bar{a}} \tau_{\bar{a}} f_{\bar{a}}^0 (1 - f_{\bar{a}}^0) \right]
$$

**各向异性 (ξ≠0)**：
$$
\sigma = \frac{1}{3T} \cdot \frac{1}{4\pi^2} \sum_a \int_0^\infty dp \int_{-1}^{1} d(\cos\theta) \, \frac{p^4 q_a^2}{E_a^2} \left[ d_a \tau_a f_a (1 - f_a) + ... \right]
$$

### 体粘滞系数 ζ

按公式文档中的表达式实现，需要热力学导数组合 $(\partial P/\partial\varepsilon)_n$ 与 $(\partial P/\partial n)_\varepsilon$ 以及质量导数 $\partial M/\partial T, \partial M/\partial\mu$。

## 积分核说明

| 输运系数 | 积分核（含相空间测度 p²） |
|---------|-------------------------|
| η | $p^6/E^2$ |
| σ | $p^4 q^2/E^2$ |
| ζ | $p^2 \times (\text{复杂表达式})$ |

## 电荷约定（自然单位）

在自然单位制中，精细结构常数 $\alpha = e^2/(4\pi) \approx 1/137$，因此元电荷：
$$e = \sqrt{4\pi\alpha} \approx 0.303$$

夸克电荷（自然单位）：
| 夸克 | 电荷 | 数值 |
|------|------|------|
| u | $(2/3)e$ | ≈ 0.202 |
| d | $(-1/3)e$ | ≈ -0.101 |
| s | $(-1/3)e$ | ≈ -0.101 |

## 模块位置
- 源码：`src/relaxtime/TransportCoefficients.jl`

## 参数结构体（可复用）

`QuarkParams` 与 `ThermoParams` 属于项目级通用参数结构体，定义在 `src/ParameterTypes.jl` 的 `ParameterTypes` 模块中。

- 推荐：`using ParameterTypes` 后直接使用这些类型。
- 兼容：本模块仍 re-export 了 `QuarkParams/ThermoParams`，旧调用无需修改。

---

## API 参考

### `shear_viscosity(quark_params, thermo_params; tau, ...) -> Float64`

计算剪切粘滞系数 $\eta$。

**必需参数**
- `quark_params::NamedTuple`：至少包含 `m=(u,d,s)` 与 `μ=(u,d,s)`（均为 fm⁻¹）。
- `thermo_params::NamedTuple`：至少包含 `T, Φ, Φbar`，可选 `ξ`（默认 0）。
- `tau::NamedTuple`：至少包含 `u,d,s,ubar,dbar,sbar`（fm）。

**可选参数**
- `config`：积分配置结构体 `TransportIntegrationConfig(...)`，用于统一传递积分相关参数。
  - 默认：`DEFAULT_TRANSPORT_CONFIG`
  - 优先级：调用时提供的同名关键字（例如 `p_nodes=...`）会覆盖 `config` 中对应字段。
  - 允许覆盖的字段：`p_nodes/p_max/p_grid/p_w/cos_nodes/cos_grid/cos_w`
- `degeneracy::Float64`：默认 $2N_c = 6$（物理简并度，不属于积分配置）。

---

### `electric_conductivity(quark_params, thermo_params; tau, charges, ...) -> Float64`

计算电导率 $\sigma$。

**额外参数**
- `charges::NamedTuple`：夸克电荷，默认使用自然单位制电荷：
  ```julia
  (u = (2/3)e, d = (-1/3)e, s = (-1/3)e)
  ```
  其中 $e = \sqrt{4\pi\alpha}$。

---

### `bulk_viscosity(quark_params, thermo_params; tau, bulk_coeffs, ...) -> Float64`

计算体粘滞系数 $\zeta$。

**额外参数**
- `bulk_coeffs::NamedTuple`：热力学导数系数，需包含：
  - `dP_depsilon_n`：$(\partial P/\partial\varepsilon)_n$
  - `dP_dn_epsilon`：$(\partial P/\partial n)_\varepsilon$
  - `dM_dT`：质量对温度的导数（3元素数组）
  - `dM_dmu`：质量对化学势的导数（3元素数组）

建议使用 `PNJL.ThermoDerivatives.bulk_derivative_coeffs(T, mu; ...)` 的返回值。

---

### `transport_coefficients(quark_params, thermo_params; tau, bulk_coeffs=nothing, ...) -> NamedTuple`

一次性计算所有输运系数。

**返回值**
```julia
(eta=..., zeta=..., sigma=...)
```

当 `bulk_coeffs=nothing` 时，`zeta` 返回 `NaN`。

---

### `default_charges() -> NamedTuple`

返回默认夸克电荷（自然单位）。

```julia
(u = (2/3)*sqrt(4π*α), d = (-1/3)*sqrt(4π*α), s = (-1/3)*sqrt(4π*α))
```

---

## 示例

```julia
include("src/relaxtime/TransportCoefficients.jl")
using .TransportCoefficients

# 物理参数
quark_params = (
    m = (u=0.3, d=0.3, s=0.5),  # fm⁻¹
    μ = (u=0.2, d=0.2, s=0.2)   # fm⁻¹
)

thermo_params = (
    T = 0.15,      # fm⁻¹
    Φ = 0.5,
    Φbar = 0.5,
    ξ = 0.0        # 各向同性
)

# 弛豫时间
tau = (
    u = 3.0, d = 3.0, s = 3.0,
    ubar = 1.0, dbar = 1.0, sbar = 2.0
)


# 计算剪切粘滞系数（直接用关键字覆盖默认积分配置）
η = shear_viscosity(quark_params, thermo_params; tau=tau, p_max=15.0, p_nodes=64)

# 计算电导率
σ = electric_conductivity(quark_params, thermo_params; tau=tau, p_max=15.0, p_nodes=64)

# 各向异性情况
thermo_params_aniso = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.1)
η_aniso = shear_viscosity(quark_params, thermo_params_aniso; tau=tau, p_max=15.0, p_nodes=64, cos_nodes=32)

# 使用配置结构体统一传参（推荐：把积分配置集中在一个结构体里）
cfg = TransportIntegrationConfig(p_nodes=64, p_max=15.0, cos_nodes=32)
η_cfg = shear_viscosity(quark_params, thermo_params_aniso; tau=tau, config=cfg)


# 也支持把 config 作为位置参数（更简洁），并可用关键字覆盖其中字段
η_cfg2 = shear_viscosity(quark_params, thermo_params_aniso, cfg; tau=tau)
η_cfg3 = shear_viscosity(quark_params, thermo_params_aniso, cfg; tau=tau, p_nodes=96)

# 推荐入口：使用更高层的请求结构体统一封装所有参数
req = TransportRequest(
  quark_params,
  thermo_params_aniso;
  tau=tau,
  integration=cfg,
  # 可选：charges / degeneracy 也可以在这里统一指定
)

η_req = shear_viscosity(req)                 # 使用 req.integration
η_req2 = shear_viscosity(req; p_nodes=96)    # 覆盖 integration 字段
σ_req = electric_conductivity(req)
all = transport_coefficients(req; bulk_coeffs=nothing)
```

## 注意事项

1. **积分核**：积分核中包含相空间测度 $p^2$，因此 η 的积分核是 $p^6/E^2$，σ 的积分核是 $p^4 q^2/E^2$。

2. **各向异性**：当 ξ≠0 时，使用 Romatschke-Strickland 形式的分布函数，需要完整的角度积分。

3. **电荷单位**：默认使用自然单位制电荷（$e = \sqrt{4\pi\alpha}$）。

4. **与C++/Fortran的一致性**：
   - **剪切粘滞系数 η**：Julia与C++/Fortran一致
   - **电导率 σ**：Julia与C++/Fortran一致（使用相同的公式和电荷约定）
