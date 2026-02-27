"""Models

新架构 models 子系统入口。

用法（建议）：

```julia
include("src/models/Models.jl")
using .Models

m = Models.create_model(:NJL)
```

说明：当前处于重构早期阶段，因此此入口以“最小可用接口”为目标。
"""
module Models

export AbstractQCDModel, AbstractNJLModel, AbstractPNJLModel
export calculate_mass_vec, calculate_chiral
export vacuum_contribution, thermal_contribution, polyakov_potential
export number_densities
export NJLModel
export NJL2Model
export PNJLModel
export PNJLMagneticModel
export RPNJLModel
export LegacyPNJLModel
export LegacyNJLModel
export create_model
export omega, omega_components, grand_potential
export MeanFieldState
export meanfield_state, state_vector
export normalize_mu_vec
export solve_gap
export gap_state_dim, gap_residual
export create_implicit_gap_solver
export create_pnjl_implicit_solver
export solve_pnjl_with_derivatives
export mass_derivatives, thermo_derivatives, bulk_derivative_coeffs
export bulk_viscosity_coefficients, compute_B_bracket
export dP_dT, dP_dmu
export solve_fixedmu_constraint
export solve_fixedrho_constraint
export solve_fixedentropy_constraint
export solve_fixedsigma_constraint
export solve_fixedasymrho_constraint
export transport_provider
export TransportProvider
export prepare_transport_provider
export run_tmu_scan, run_trho_scan
export build_default_rho_grid
export legacy_pnjl_module
export solve_gap_and_transport, solve_transport_from_equilibrium
export solve_gap_and_meson_point

include(joinpath(@__DIR__, "abstract_model.jl"))

# State representation
include(joinpath(@__DIR__, "state.jl"))

# Concrete models
include(joinpath(@__DIR__, "njl", "NJLModel.jl"))
include(joinpath(@__DIR__, "njl", "NJL2Model.jl"))
include(joinpath(@__DIR__, "pnjl", "PNJLDistributions.jl"))
include(joinpath(@__DIR__, "pnjl", "PNJLModel.jl"))
include(joinpath(@__DIR__, "pnjl", "PNJLMagneticModel.jl"))
include(joinpath(@__DIR__, "rpnjl", "RPNJLModel.jl"))

# Factory
include(joinpath(@__DIR__, "factory.jl"))

# Ω assembly
include(joinpath(@__DIR__, "omega.jl"))
include(joinpath(@__DIR__, "gap_solver.jl"))
include(joinpath(@__DIR__, "implicit_gap.jl"))
include(joinpath(@__DIR__, "constraint_solver.jl"))
include(joinpath(@__DIR__, "derivative_entrypoints.jl"))

# Transport provider (distribution/dispersion) for Stage-4 workflow decoupling
include(joinpath(@__DIR__, "transport_provider.jl"))

# Unified scan/workflow entrypoints (Stage C compatibility layer)
include(joinpath(@__DIR__, "entrypoints.jl"))

end # module Models
