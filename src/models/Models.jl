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
export ConstraintModes
export SeedStrategy, DefaultSeed, MultiSeed, ContinuitySeed, HybridContinuitySeed, PhaseAwareSeed, PhaseAwareContinuitySeed
export get_seed, update!, reset!, get_all_seeds, set_phase!
export HADRON_SEED_5, QUARK_SEED_5, HADRON_SEED_8, QUARK_SEED_8
export build_conditions, build_residual!, GapParams
export create_implicit_solver, solve_with_derivatives
export ρ0
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
export Integrals, cached_nodes, vacuum_integral, calculate_energy_sum, calculate_number_densities
export calculate_U, calculate_U_derivative_T
export Constants_PNJL
export calculate_omega, calculate_pressure, calculate_rho, calculate_thermo
export TmuScanConfig, TrhoScanConfig

include(joinpath(@__DIR__, "abstract_model.jl"))

# State representation
include(joinpath(@__DIR__, "state.jl"))

# Concrete models
include(joinpath(@__DIR__, "njl", "NJLModel.jl"))
include(joinpath(@__DIR__, "njl", "NJL2Model.jl"))
include(joinpath(@__DIR__, "pnjl_physics", "PNJLDistributions.jl"))
include(joinpath(@__DIR__, "pnjl_physics", "PNJLCore.jl"))
include(joinpath(@__DIR__, "pnjl_physics", "PNJLModel.jl"))
include(joinpath(@__DIR__, "pnjl_physics", "PNJLMagneticModel.jl"))
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
include(joinpath(@__DIR__, "pnjl_physics", "core", "MagneticIntegrals.jl"))
include(joinpath(@__DIR__, "pnjl_physics", "core", "MagneticThermodynamics.jl"))
include(joinpath(@__DIR__, "scans", "ScanCommon.jl"))
include(joinpath(@__DIR__, "scans", "ScanConfig.jl"))
include(joinpath(@__DIR__, "scans", "ScanResultFinalize.jl"))
include(joinpath(@__DIR__, "scans", "TmuScan.jl"))
include(joinpath(@__DIR__, "scans", "TrhoScan.jl"))
include(joinpath(@__DIR__, "scans", "DualBranchScan.jl"))
include(joinpath(@__DIR__, "scans", "DualBranchScanEntry.jl"))

using .SeedStrategies
using .Conditions
import .ImplicitSolver
const ConstraintModes = @__MODULE__
using .ThermoDerivatives
using .MagneticIntegrals
using .MagneticThermodynamics
using .TmuScan
using .TrhoScan

# Transport provider (distribution/dispersion) for Stage-4 workflow decoupling
include(joinpath(@__DIR__, "transport_provider.jl"))

const ρ0 = Main.Constants_PNJL.ρ0_inv_fm3

# Unified scan/workflow entrypoints (Stage C compatibility layer)
include(joinpath(@__DIR__, "phase", "PhaseTypes.jl"))
include(joinpath(@__DIR__, "phase", "AdaptiveRhoRefinement.jl"))
include(joinpath(@__DIR__, "phase", "PhaseCore.jl"))
include(joinpath(@__DIR__, "phase", "PhaseIO.jl"))
include(joinpath(@__DIR__, "phase", "CEPDetector.jl"))
include(joinpath(@__DIR__, "phase", "CrossoverLine.jl"))
include(joinpath(@__DIR__, "phase", "PhaseArtifacts.jl"))
include(joinpath(@__DIR__, "phase", "PhasePipeline.jl"))
include(joinpath(@__DIR__, "workflows", "WorkflowParamAdapters.jl"))
include(joinpath(@__DIR__, "workflows", "TransportWorkflow.jl"))
include(joinpath(@__DIR__, "workflows", "MesonMassWorkflow.jl"))

if !isdefined(Main, :WorkflowParamAdapters)
	@eval Main const WorkflowParamAdapters = $WorkflowParamAdapters
end
if !isdefined(Main, :TransportWorkflow)
	@eval Main const TransportWorkflow = $TransportWorkflow
end
if !isdefined(Main, :MesonMassWorkflow)
	@eval Main const MesonMassWorkflow = $MesonMassWorkflow
end

include(joinpath(@__DIR__, "entrypoints.jl"))

const Integrals = PNJLCore.PNJLIntegrals
const Constants_PNJL = Main.Constants_PNJL
const TmuScanConfig = ScanConfig.TmuScanConfig
const TrhoScanConfig = ScanConfig.TrhoScanConfig
const update! = SeedStrategies.update!
const reset! = SeedStrategies.reset!

@inline function cached_nodes(
	p_num::Int=PNJLCore.DEFAULT_MOMENTUM_COUNT,
	t_num::Int=PNJLCore.DEFAULT_THETA_COUNT,
)
	return PNJLCore.cached_nodes(p_num, t_num)
end

@inline function vacuum_integral(mass)
	TT = typeof(mass)
	Λ = convert(TT, Main.Constants_PNJL.Λ_inv_fm)
	return PNJLCore.vacuum_integral_with_cutoff(convert(TT, mass), Λ)
end

@inline function calculate_energy_sum(masses)
	total = zero(eltype(masses))
	@inbounds for i in 1:length(masses)
		total += vacuum_integral(masses[i])
	end
	return -2 * Main.Constants_PNJL.N_color * total
end

@inline function calculate_mass_vec(x_state::AbstractVector{<:Real}; kwargs...)
	model = create_model(:PNJL)
	φ = SVector{3, Float64}(Float64(x_state[1]), Float64(x_state[2]), Float64(x_state[3]))
	return calculate_mass_vec(model, φ; kwargs...)
end

@inline function calculate_chiral(φ::AbstractVector{<:Real}; kwargs...)
	model = create_model(:PNJL)
	ϕ = SVector{3, Float64}(Float64(φ[1]), Float64(φ[2]), Float64(φ[3]))
	return calculate_chiral(model, ϕ; kwargs...)
end

@inline function calculate_number_densities(
	x_state,
	mu_vec,
	T_fm,
	thermal_nodes=nothing,
	xi=0.0;
	p_num::Int=PNJLCore.DEFAULT_MOMENTUM_COUNT,
	t_num::Int=PNJLCore.DEFAULT_THETA_COUNT,
	kwargs...
)
	model = create_model(:PNJL)
	return number_densities(
		model,
		x_state,
		T_fm,
		mu_vec;
		thermal_nodes=thermal_nodes,
		p_num=p_num,
		t_num=t_num,
		xi=xi,
		kwargs...,
	)
end

@inline function calculate_U(T_fm, Φ, Φbar)
	model = create_model(:PNJL)
	return polyakov_potential(model, Φ, Φbar, T_fm)
end

@inline function calculate_U_derivative_T(T_fm, Φ, Φbar)
	δT = 1e-6
	return (calculate_U(T_fm + δT, Φ, Φbar) - calculate_U(T_fm - δT, Φ, Φbar)) / (2δT)
end

@inline function calculate_omega(x_state, mu_vec, T_fm, thermal_nodes=nothing, xi=0.0; p_num::Int=PNJLCore.DEFAULT_MOMENTUM_COUNT, t_num::Int=PNJLCore.DEFAULT_THETA_COUNT)
	model = create_model(:PNJL)
	return omega(model, x_state, T_fm, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
end

@inline function calculate_pressure(x_state, mu_vec, T_fm, thermal_nodes=nothing, xi=0.0; p_num::Int=PNJLCore.DEFAULT_MOMENTUM_COUNT, t_num::Int=PNJLCore.DEFAULT_THETA_COUNT)
	model = create_model(:PNJL)
	return model_pressure(model, x_state, mu_vec, T_fm; p_num=p_num, t_num=t_num, xi=xi)
end

@inline function calculate_rho(x_state, mu_vec, T_fm, thermal_nodes=nothing, xi=0.0; p_num::Int=PNJLCore.DEFAULT_MOMENTUM_COUNT, t_num::Int=PNJLCore.DEFAULT_THETA_COUNT)
	model = create_model(:PNJL)
	return model_rho(model, x_state, mu_vec, T_fm; p_num=p_num, t_num=t_num, xi=xi)
end

@inline function calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes=nothing, xi=0.0; p_num::Int=PNJLCore.DEFAULT_MOMENTUM_COUNT, t_num::Int=PNJLCore.DEFAULT_THETA_COUNT)
	model = create_model(:PNJL)
	return model_thermo(model, x_state, mu_vec, T_fm; p_num=p_num, t_num=t_num, xi=xi)
end

if !isdefined(Main, :PNJL)
	@eval Main const PNJL = $(@__MODULE__)
end

end # module Models
