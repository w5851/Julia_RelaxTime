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
export model_pressure, model_rho, model_thermo
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
export solve_constraint
export solve, solve_multi, SolverResult
export SeedStrategy, DefaultSeed, MultiSeed, ContinuitySeed, HybridContinuitySeed, PhaseAwareSeed, PhaseAwareContinuitySeed
export get_seed, update!, reset!, get_all_seeds, set_phase!
export HADRON_SEED_5, QUARK_SEED_5, HADRON_SEED_8, QUARK_SEED_8
export build_conditions, build_residual!, GapParams
export create_implicit_solver, solve_with_derivatives
export QUARK_CHARGE_ABS
export alpha_n, energy_landau, smooth_cutoff, resolve_nmax_from_cutoff
export omega0_flavor_landau, omegat_flavor_landau, density_flavor_landau
export MagneticIMCParams, default_imc_params, coupling_GB
export MagneticConfig, default_magnetic_config
export calculate_magnetic_omega_components, calculate_magnetic_omega
export calculate_magnetic_pressure, calculate_magnetic_rho
export calculate_magnetic_number_densities
export magnetic_nmax_convergence_report
export transport_provider
export TransportProvider
export prepare_transport_provider
export run_tmu_scan, run_trho_scan
export build_default_rho_grid
export load_dual_branch_scan!
export pnjl_module
export solve_gap_and_transport, solve_transport_from_equilibrium
export solve_gap_and_meson_point
export run_phase_pipeline, find_cep, build_phase_artifacts
export resolve_phase_output_target, promote_phase_artifacts
export CEPResult, PromotionResult, PhasePipelineResult

include(joinpath(@__DIR__, "abstract_model.jl"))

# State representation
include(joinpath(@__DIR__, "state.jl"))

# Concrete models
include(joinpath(@__DIR__, "njl", "NJLModel.jl"))
include(joinpath(@__DIR__, "njl", "NJL2Model.jl"))
include(joinpath(@__DIR__, "pnjl", "PNJLDistributions.jl"))
include(joinpath(@__DIR__, "pnjl", "PNJLCore.jl"))
include(joinpath(@__DIR__, "pnjl", "PNJLModel.jl"))
include(joinpath(@__DIR__, "pnjl", "PNJLMagneticModel.jl"))
include(joinpath(@__DIR__, "rpnjl", "RPNJLModel.jl"))

# Backward-compatible access path used by some tests/callers:
# Main.Models.PNJLIntegrals.*
const PNJLIntegrals = PNJLCore.PNJLIntegrals

# Factory
include(joinpath(@__DIR__, "factory.jl"))

# Ω assembly
include(joinpath(@__DIR__, "omega.jl"))
include(joinpath(@__DIR__, "thermo_kernel.jl"))
include(joinpath(@__DIR__, "gap_solver.jl"))
include(joinpath(@__DIR__, "implicit_gap.jl"))
include(joinpath(@__DIR__, "constraint_solver.jl"))
include(joinpath(@__DIR__, "solver", "ConstraintModes.jl"))
include(joinpath(@__DIR__, "solver", "SeedStrategies.jl"))
include(joinpath(@__DIR__, "solver", "Conditions.jl"))
include(joinpath(@__DIR__, "solver", "ImplicitSolver.jl"))
include(joinpath(@__DIR__, "solver", "Solver.jl"))
include(joinpath(@__DIR__, "derivatives", "ThermoDerivatives.jl"))
include(joinpath(@__DIR__, "pnjl", "core", "MagneticIntegrals.jl"))
include(joinpath(@__DIR__, "pnjl", "core", "MagneticThermodynamics.jl"))
include(joinpath(@__DIR__, "scans", "ScanCommon.jl"))
include(joinpath(@__DIR__, "scans", "ScanConfig.jl"))
include(joinpath(@__DIR__, "scans", "ScanResultFinalize.jl"))
include(joinpath(@__DIR__, "scans", "TmuScan.jl"))
include(joinpath(@__DIR__, "scans", "TrhoScan.jl"))
include(joinpath(@__DIR__, "scans", "DualBranchScan.jl"))
include(joinpath(@__DIR__, "scans", "AdaptiveRhoRefinement.jl"))
include(joinpath(@__DIR__, "scans", "DualBranchScanEntry.jl"))

using .SeedStrategies
using .Conditions
import .ImplicitSolver
using .ThermoDerivatives
using .MagneticIntegrals
using .MagneticThermodynamics
using .TmuScan
using .TrhoScan

# Transport provider (distribution/dispersion) for Stage-4 workflow decoupling
include(joinpath(@__DIR__, "transport_provider.jl"))

# Unified scan/workflow entrypoints (Stage C compatibility layer)
include(joinpath(@__DIR__, "phase", "PhaseTypes.jl"))
include(joinpath(@__DIR__, "phase", "AdaptiveRhoRefinement.jl"))
include(joinpath(@__DIR__, "phase", "PhaseCore.jl"))
include(joinpath(@__DIR__, "phase", "PhaseIO.jl"))
include(joinpath(@__DIR__, "phase", "CEPDetector.jl"))
include(joinpath(@__DIR__, "phase", "CrossoverLine.jl"))
include(joinpath(@__DIR__, "phase", "PhaseArtifacts.jl"))
include(joinpath(@__DIR__, "phase", "PhasePipeline.jl"))
include(joinpath(@__DIR__, "entrypoints.jl"))

end # module Models
