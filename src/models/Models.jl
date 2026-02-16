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
export PNJLModel
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
export transport_provider
export TransportProvider
export prepare_transport_provider

include(joinpath(@__DIR__, "abstract_model.jl"))

# State representation
include(joinpath(@__DIR__, "state.jl"))

# Concrete models
include(joinpath(@__DIR__, "njl", "NJLModel.jl"))
include(joinpath(@__DIR__, "pnjl", "PNJLDistributions.jl"))
include(joinpath(@__DIR__, "pnjl", "PNJLModel.jl"))
include(joinpath(@__DIR__, "rpnjl", "RPNJLModel.jl"))

# Factory
include(joinpath(@__DIR__, "factory.jl"))

# Ω assembly
include(joinpath(@__DIR__, "omega.jl"))
include(joinpath(@__DIR__, "gap_solver.jl"))
include(joinpath(@__DIR__, "implicit_gap.jl"))

# Transport provider (distribution/dispersion) for Stage-4 workflow decoupling
include(joinpath(@__DIR__, "transport_provider.jl"))

end # module Models
