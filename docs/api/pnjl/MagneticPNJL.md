# Magnetic PNJL API 说明

本文档说明外磁场 PNJL（Landau 能级离散化）在项目中的已实现 API。

## 模块位置

- 核心积分：`src/pnjl/core/MagneticIntegrals.jl`
- 核心热力学：`src/pnjl/core/MagneticThermodynamics.jl`
- 顶层导出：`src/pnjl/PNJL.jl`
- models 适配：`src/models/pnjl/PNJLMagneticModel.jl`

## 单位约定

- `T_fm`, `mu_fm`, 质量、动量：`fm⁻¹`
- `eB_fm2`：`fm⁻²`
- 配置/脚本中常见 `eB_MeV2`，换算为：
  - `eB_fm2 = eB_MeV2 / ħc_MeV_fm^2`

## 核心数据结构

### `MagneticIMCParams`

用于 IMC 耦合 `G(B)` 参数化：

- `a::Float64`
- `b::Float64`
- `c::Float64`
- `d::Float64`
- `Λ_QCD_MeV::Float64`

默认值来源：`docs/reference/formula/models/pnjl_magnetic/PNJL_magnetic_core.md`。

### `MagneticConfig`

磁场求解配置：

- `eB_fm2::Float64`：外磁场强度
- `n_max::Union{Int,Nothing}`：Landau 截断层数（`nothing` 为自动估计）
- `p_num::Int`：`p_z` 积分节点数
- `pz_max::Float64`：`p_z` 积分上界
- `cutoff_N::Int`：平滑截断指数（默认 10）
- `imc::MagneticIMCParams`：IMC 参数

构造函数：

- `default_imc_params()`
- `default_magnetic_config(; eB_fm2=0.0)`

## 积分与能谱 API

- `alpha_n(n)`：Landau 简并度（`n=0 -> 1`, 其余 `2`）
- `energy_landau(mass, pz, n, q_abs, eB)`：
  - `E_{f,n} = sqrt(2n|q_f|eB + pz^2 + M_f^2)`
- `smooth_cutoff(p; Λ, N)`：平滑截断
- `resolve_nmax_from_cutoff(mass, mu, q_abs, eB; Λ)`：自动估计 `n_max`
- `omega0_flavor_landau(...)`：单味真空项 `Ω_f^0`
- `omegat_flavor_landau(...)`：单味热项 `Ω_f^T`
- `density_flavor_landau(...)`：单味密度近似（低温路径）

## 热力学 API

- `coupling_GB(eB_fm2; G0, imc)`：IMC 耦合 `G(B)`
- `calculate_magnetic_omega_components(x_state, mu_vec, T_fm, magnetic)`：返回
  - `chi`, `poly`, `vac`, `therm`, `masses`, `omega`, `n_max`, `G_B`
- `calculate_magnetic_omega(...)`
- `calculate_magnetic_pressure(...)`
- `calculate_magnetic_rho(...)`（中心差分）
- `calculate_magnetic_number_densities(...)`
- `magnetic_nmax_convergence_report(...; delta_n=6, rtol=3e-2)`

### `x_state` 约定

`SVector{5}` 顺序：`(φ_u, φ_d, φ_s, Φ, Φbar)`。

### `mu_vec` 约定

`SVector{3}` 顺序：`(μ_u, μ_d, μ_s)`。

## eB→0 退化行为

当 `|eB_fm2| <= 1e-14` 时，`calculate_magnetic_omega_components` 自动退化到现有 PNJL `calculate_omega` 路径，保持与零磁场主链一致。

## models 侧 API

`create_model(:PNJLMagnetic; eB_fm2=...)` 可构造 `PNJLMagneticModel`。

当前实现策略：

- gap 解复用 `PNJLModel` 的稳定路径；
- Ω 组装与密度输出走磁场实现；
- 用于 workflow/回归链路的最小可用适配。

## 相关脚本

- 单点示例：`scripts/pnjl/run_magnetic_point.jl`
- 固定点 eB 扫描：`scripts/pnjl/run_magnetic_eb_scan.jl`
- nightly 回归门禁：`scripts/pnjl/run_magnetic_fixedpoint_regression.jl`
- baseline 导出：`scripts/dev/export_pnjl_magnetic_fixedpoint_baseline.jl`

## 相关测试

- `tests/unit/pnjl/test_magnetic_energy_levels.jl`
- `tests/unit/pnjl/test_magnetic_omega_components.jl`
- `tests/unit/pnjl/test_magnetic_density.jl`
- `tests/unit/pnjl/test_magnetic_coupling_GB.jl`
- `tests/unit/pnjl/test_magnetic_nmax_convergence.jl`
- `tests/unit/pnjl/test_magnetic_fixedpoint_baseline_smoke.jl`
