# TransportCoefficients 模块 API 文档

本模块实现弛豫时间近似（RTA）下的夸克物质输运系数计算：剪切粘滞系数 $\eta$、体粘滞系数 $\zeta$、电导率 $\sigma$。

如果你当前更关心 transport provider 的字段契约、`Models.TransportProvider` 与 `prepare_transport_provider` 的桥接语义，建议先阅读 `CoreConcepts.md`，本页聚焦系数公式、积分配置与直接调用方式。

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

权威公式、能量角色和外部文献映射见 [TransportCoefficients_FromRelaxationTime.md](../../../reference/formula/relaxtime/transport/TransportCoefficients_FromRelaxationTime.md)；守恒荷扩散矩阵与热导率见 [KappaSeries_AndLambda_FromLegacyFortran.md](../../../reference/formula/relaxtime/transport/KappaSeries_AndLambda_FromLegacyFortran.md)。

本 API 页采用以下公式语义：

- $E_{\mathrm{kin}}=\sqrt{p^2+m^2}$ 用于 $\eta$、$\sigma$、$\zeta$、$\kappa_{XY}$ 的运动学核，用于 $\zeta$ 的等熵导数组合，并用于 $\kappa_{XY}$ 的 Landau-Lifshitz 投影；
- $E_{\mathrm{dist}}=\sqrt{p^2+m^2+\xi(p\cos\theta)^2}$ 只作为 RS 分布 $f_\xi$ 的自变量；
- $\lambda=\kappa_{BB}[(\epsilon+P)/(n_BT)]^2$ 不直接选择单粒子能量，但由采用上述能量分工的 $\kappa_{BB}$ 重新计算；
- $\zeta$ 的等熵平方核对应 Albright and Kapusta (2016) Eq. (138)，经典统计版本为 Eq. (111)；历史记号“A26”只是无法溯源到外部文献的 legacy Fortran/内部别名；
- 方括号采用 $p^2+3v_n^2T^2E\,\partial_T[(E\mp\mu)/T]_\sigma$ 的加号形式。内部核对稿件中出现的减号已由作者确认为文稿错误。

当前源码已经按上述 $E_{\mathrm{kin}}/E_{\mathrm{dist}}$ 分工实现 $\eta$、$\sigma$、$\zeta$ 和 $\kappa_{XY}$，$\lambda$ 则由修正后的 $\kappa_{BB}$ 派生。该修改不追溯覆盖旧 production；论文输入资格以 `data/outputs/results/relaxtime/transport/phase_guided/production_registry.json` 为准。

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

采用等熵声速平方核，需要 $v_n^2$、$(\partial\mu_B/\partial T)_\sigma$、准粒子质量，以及 $\partial M/\partial T$、$\partial M/\partial\mu_B$。完整公式与外部方程映射见权威公式文档。

## 积分核说明

| 输运系数 | 积分核（含相空间测度 p²） |
|---------|-------------------------|
| η | $p^6/E^2$ |
| σ | $p^4 q^2/E^2$ |
| ζ | $p^2 \times (\text{复杂表达式})$ |
| $\kappa_{XY}$ | $p^4/E^2$，并含 $q_X-n_XE/(\epsilon+P)$ 投影 |

上述 $\kappa_{XY}$ 分母和投影中的 $E$ 均为 $E_{\mathrm{kin}}$；占据分布使用 $E_{\mathrm{dist}}$。$\lambda$ 是 $\kappa_{BB}$ 的派生量，不新增独立动量积分。

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

### `bulk_viscosity_isentropic(quark_params, thermo_params; tau, bulk_coeffs_isentropic, ...) -> Float64`

计算体粘滞系数 $\zeta$。

**额外参数**
- `bulk_coeffs_isentropic::NamedTuple`：等熵体粘滞导数系数，需包含：
  - `v_n_sq`：$v_n^2$
  - `dμB_dT_sigma`：$(\partial\mu_B/\partial T)_\sigma$
  - `masses`：三个味的准粒子质量
  - `dM_dT`：质量对温度的导数（3元素数组）
  - `dM_dμB`：质量对重子化学势的导数（3元素数组）

建议使用统一公开入口 `Models.bulk_viscosity_coefficients(T_fm, mu_fm; ...)` 的返回值。

### `bulk_viscosity(quark_params, thermo_params; formula=:isentropic, ...) -> Float64`

统一体粘滞入口，当前支持：

- `formula=:isentropic`：转发到 `bulk_viscosity_isentropic(...)`

兼容参数：

- `bulk_coeffs_isentropic=...`（推荐）
- `bulk_coeffs=...`（兼容旧参数名）

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
include("../../../../src/relaxtime/TransportCoefficients.jl")
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

## 输入校验与数值保护

- `thermo_params.T` 必须为有限正数；`m_u,m_d,m_s` 必须为有限非负数。
- `tau` 必须包含 `u,d,s,ubar,dbar,sbar` 且每个值有限且非负。
- `TransportIntegrationConfig` 约束：`p_nodes>0`、`cos_nodes>0`、`p_max>0`，自定义网格/权重必须成对、长度一致且有限。
- `electric_conductivity` 的 `charges.u/d/s` 必须为有限数。
- 数值保护：
  - $E_{\mathrm{kin}}$ 与 $E_{\mathrm{dist}}$ 分别采用下限保护：`E = max(E_raw, sqrt(eps(Float64)))`
  - 费米因子采用物理区间截断：`f(1-f) -> clamp(f(1-f), 0.0, 0.25)`

## provider 摘要

本页只保留 transport 计算直接相关的 provider 摘要；完整字段契约、`Models.TransportProvider`、`prepare_transport_provider` 与 `ctx` 语义见 `CoreConcepts.md`。

- transport 计算至少依赖：`energy_from_p`、`quark_distribution`、`antiquark_distribution`
- `energy_from_p` 提供普通在壳 $E_{\mathrm{kin}}$；`energy_from_p_aniso` 若存在，只生成 RS 分布自变量 $E_{\mathrm{dist}}$
- 当 `xi != 0` 时，`prefer_energy_aniso` 只决定优先把 $E_{\mathrm{dist}}$ 传给普通分布接口，还是调用 provider 自带的 `*_distribution_aniso`；它不改变输运核能量
- 若你需要自定义 species 质量/化学势解析，或理解 `prepare_transport_provider(...)` 如何把平衡态结果注入 provider，应转到 `CoreConcepts.md`

## 注意事项

1. **积分核**：积分核中包含相空间测度 $p^2$，因此 η 的积分核是 $p^6/E^2$，σ 的积分核是 $p^4 q^2/E^2$。

2. **各向异性**：当 ξ≠0 时，使用 Romatschke-Strickland 形式的分布函数，需要完整的角度积分。

3. **provider 路径**：当 `ξ≠0` 时，`energy_from_p_aniso` 与 `prefer_energy_aniso` 只控制 RS 分布的能量直通/各向异性分布接口路由；积分核始终使用 `energy_from_p` 提供的普通在壳能量。如果你需要完整字段表、回退顺序或 `prepare_transport_provider(...)` 的桥接细节，请阅读 `CoreConcepts.md`。

4. **电荷单位**：默认使用自然单位制电荷（$e = \sqrt{4\pi\alpha}$）。

5. **与C++/Fortran的一致性**：
   - **剪切粘滞系数 η**：Julia与C++/Fortran一致
   - **电导率 σ**：Julia与C++/Fortran一致（使用相同的公式和电荷约定）
