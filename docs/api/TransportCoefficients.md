# TransportCoefficients 模块 API 文档

本模块实现弛豫时间近似（RTA）下的夸克物质输运系数计算：剪切粘滞系数 $\eta$、体粘滞系数 $\zeta$、电导率 $\sigma$。

## 单位约定
- 自然单位：$\hbar=c=k_B=1$
- 温度/化学势/质量/动量：fm⁻¹
- 弛豫时间 $\tau$：fm

## 物理公式（与 docs/reference/formula/输运系数by弛豫时间.md 对应）
- 剪切粘滞：

$$
\eta = \frac{1}{15T} \sum_a \int \frac{d^3p}{(2\pi)^3}\, \frac{p^4}{E_a^2}\, d_a\,\tau_a\, f_a^0(1-f_a^0) + (a\to \bar a)
$$

- 电导率：

$$
\sigma = \frac{2}{3T} \sum_a \int \frac{d^3p}{(2\pi)^3}\, p^2\, \frac{q_a^2}{E_a^2}\, d_a\,\tau_a\, f_a^0(1-f_a^0) + (a\to \bar a)
$$

- 体粘滞：按同一份公式文档中的表达式直译实现；需要热力学导数组合 $(\partial P/\partial\varepsilon)_n$ 与 $(\partial P/\partial n)_\varepsilon$ 以及质量导数 $\partial M/\partial T,\partial M/\partial\mu$。

## 模块位置
- 源码：src/relaxtime/TransportCoefficients.jl

## API 参考

### `shear_viscosity(quark_params, thermo_params; tau, ...) -> Float64`
计算剪切粘滞系数 $\eta$。

**必需参数**
- `quark_params::NamedTuple`：至少包含 `m=(u,d,s)` 与 `μ=(u,d,s)`（均为 fm⁻¹）。
- `thermo_params::NamedTuple`：至少包含 `T, Φ, Φbar`，可选 `ξ`（默认 0）。
- `tau::NamedTuple`：至少包含 `u,d,s,ubar,dbar,sbar`（fm）。

**可选参数**
- `degeneracy`：默认 $2N_c$（从 Constants_PNJL.N_color 推出）。
- `p_nodes,p_max,p_grid,p_w`：动量积分设置。
- `cos_nodes,cos_grid,cos_w`：当 `ξ≠0` 时的角度积分设置。

### `electric_conductivity(quark_params, thermo_params; tau, charges, ...) -> Float64`
计算电导率 $\sigma$。

- `charges::NamedTuple` 默认 `(; u=2/3, d=-1/3, s=-1/3)`，公式中使用 $q^2$。

### `bulk_viscosity(quark_params, thermo_params; tau, bulk_coeffs, ...) -> Float64`
计算体粘滞系数 $\zeta$。

- `bulk_coeffs::NamedTuple`：建议使用 `PNJL.ThermoDerivatives.bulk_derivative_coeffs(T_fm, mu_fm; ...)` 的返回值。
- `n` 的物理含义：净重子数密度（仓库中 `rho = sum(rho_vec)/3` 的量）。

### `transport_coefficients(quark_params, thermo_params; tau, bulk_coeffs=nothing, ...) -> NamedTuple`
一次性返回 `(eta, zeta, sigma)`。

- 当 `bulk_coeffs=nothing` 时，`zeta` 返回 `NaN`。

## 示例

```julia
include("../../src/relaxtime/RelaxationTime.jl")
include("../../src/relaxtime/TransportCoefficients.jl")
include("../../src/pnjl/analysis/ThermoDerivatives.jl")

using .RelaxationTime
using .TransportCoefficients
using .ThermoDerivatives

quark_params = (m=(u=0.3,d=0.3,s=0.5), μ=(u=0.2,d=0.2,s=0.2))
thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)

# 这里 densities 由调用方提供；示例用占位值
result_tau = relaxation_times(quark_params, thermo_params, (K_σπ=1.0,K_σK=1.0,K_σ=1.0,K_δπ=1.0,K_δK=1.0);
    densities=(u=1.0,d=1.0,s=1.0,ubar=1.0,dbar=1.0,sbar=1.0),
    existing_rates=(uu_to_uu=1.0,ud_to_ud=1.0,us_to_us=1.0,usbar_to_usbar=1.0,uubar_to_uubar=1.0,uubar_to_ddbar=1.0,uubar_to_ssbar=1.0,udbar_to_udbar=1.0,ss_to_ss=1.0,ssbar_to_ssbar=1.0,ssbar_to_uubar=1.0)
)

bulk = bulk_derivative_coeffs(thermo_params.T, quark_params.μ.u; xi=thermo_params.ξ)

coeffs = transport_coefficients(quark_params, thermo_params;
    tau=result_tau.tau,
    bulk_coeffs=bulk,
    p_max=10.0,
)

println(coeffs)
```
