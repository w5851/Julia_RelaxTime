# RelaxationTime 模块 API

弛豫时间计算，基于平均散射率表达式（见《弛豫时间by平均散射率.md》）。一次计算所需散射过程的平均散射率，可复用于六种夸克（d/dbar 直接复用 u/ubar 的结果）。

## 导出符号
- `REQUIRED_PROCESSES`：需要的散射过程符号列表。
- `compute_average_rates(quark_params, thermo_params, K_coeffs; existing_rates, cs_caches, p_nodes, angle_nodes, phi_nodes, p_grid, p_w, cos_grid, cos_w, phi_grid, phi_w, n_sigma_points)`
- `relaxation_rates(densities, rates)`
- `relaxation_times(quark_params, thermo_params, K_coeffs; densities, existing_rates, cs_caches, p_nodes, angle_nodes, phi_nodes, p_grid, p_w, cos_grid, cos_w, phi_grid, phi_w, n_sigma_points)`

## 数学定义
弛豫时间满足
$$\tau_i^{-1} = \sum_j \rho_j \; \bar{w}_{ij}$$
其中 \(\rho_j\) 为粒子数密度，\(\bar{w}_{ij}\) 为平均散射率。

## 参数说明
- `quark_params`：含质量、化学势等的 NamedTuple（与散射率模块一致）。
- `thermo_params`：热力学参数 `(T, Φ, Φbar, ξ)`。
- `K_coeffs`：有效耦合系数 NamedTuple。
- `densities`：外部提供的数密度，必须包含 `:u, :d, :s, :ubar, :dbar, :sbar`。
- `existing_rates`：可选，预先计算好的平均散射率集合（NamedTuple 或 Dict）。若已包含 `REQUIRED_PROCESSES`，则不会触发新的散射率积分。
- `cs_caches`：可选，`Dict{Symbol,CrossSectionCache}`，用于重用截面拟合。
- `p_nodes/angle_nodes/phi_nodes`：高斯节点数量。
- `p_grid/p_w/cos_grid/cos_w/phi_grid/phi_w`：可选自定义积分节点与权重；为 `nothing` 时内部使用默认高斯节点（来自 `AverageScatteringRate`）。
- `n_sigma_points`：总截面 t 积分节点数。

## 行为细节
- `p_grid` 等网格参数未提供（`nothing`）时，会在 `average_scattering_rate` 内部自动调用默认高斯节点设置；不会传递 `nothing` 继续下去。
- `existing_rates` 已覆盖全部 `REQUIRED_PROCESSES` 时，`compute_average_rates` 不会再计算散射率，直接返回传入值。
- d、dbar 的结果直接复用 u、ubar 的计算（同位旋对称）。

## 典型用法
```julia
include("src/relaxtime/RelaxationTime.jl")
using .RelaxationTime

# 预先给出数密度与散射率（避免重复积分）
densities = (u=1.0, d=1.0, s=2.0, ubar=3.0, dbar=3.0, sbar=4.0)
rates = (
    uu_to_uu=1.0, ud_to_ud=2.0, us_to_us=3.0, usbar_to_usbar=5.0,
    uubar_to_uubar=7.0, uubar_to_ddbar=11.0, uubar_to_ssbar=13.0,
    udbar_to_udbar=17.0, ss_to_ss=19.0, ssbar_to_ssbar=23.0, ssbar_to_uubar=29.0,
)

quark_params = (m=(u=0.1,d=0.1,s=0.2), μ=(u=0.0,d=0.0,s=0.0))
thermo_params = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
K_coeffs = (K_σπ=1.0, K_σK=1.0, K_σ=1.0, K_δπ=1.0, K_δK=1.0)

result = relaxation_times(
    quark_params, thermo_params, K_coeffs;
    densities=densities,
    existing_rates=rates,
)

println(result.tau_inv)  # 各味的 τ^-1
println(result.tau)      # 各味的 τ
```
