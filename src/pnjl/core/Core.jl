"""
    Core

PNJL 模型核心计算模块，包含积分和热力学量计算。

## 子模块
- `Integrals`: 积分计算（真空项、热项）
- `Thermodynamics`: 热力学量计算（Ω, P, s, ε, ρ, M）
"""
module Core

include("Integrals.jl")
include("Thermodynamics.jl")

using .Integrals
using .Thermodynamics

# 重新导出常用函数
export cached_nodes, vacuum_integral, calculate_energy_sum, calculate_log_sum
export DEFAULT_THETA_COUNT, DEFAULT_MOMENTUM_COUNT
export safe_log, calculate_energy_isotropic, calculate_energy_anisotropic

export calculate_mass_vec, calculate_chiral, calculate_U
export calculate_pressure, calculate_omega
export calculate_rho, calculate_thermo, calculate_number_densities
export ρ0

end # module Core
