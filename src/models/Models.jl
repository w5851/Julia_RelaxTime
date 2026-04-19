"""Models

新架构 models 子系统统一入口。

推荐用法：

```julia
include("src/models/Models.jl")
using .Models

model = Models.create_model(:NJL)
```

该模块暴露面向模型对象的多重派发接口，以及少量显式模块访问器。
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
export RotationModel
export GasLiquidModel
export create_model
export get_cached_model, clear_model_cache!
export register_model!, unregister_model!, registered_model_kinds
export omega, omega_components, grand_potential
export model_pressure, model_rho, model_thermo
export MeanFieldState
export meanfield_state, state_vector
export normalize_mu_vec
export solve_gap
export gap_state_dim, gap_residual
export ImplicitProblem, ImplicitSolverConfig, build_implicit_solver
export build_pnjl_fixedmu_problem, build_pnjl_flavor_mu_problem, build_njl_problem
export create_implicit_gap_solver
export create_flavor_mu_implicit_gap_solver
export create_pnjl_implicit_solver
export derive_vec, derive_named
export build_pnjl_fixedmu_adapters, build_pnjl_flavor_mu_adapters
export solve_pnjl_with_derivatives
export solve_pnjl_with_flavor_mu_derivatives
export mass_derivatives, thermo_derivatives, bulk_derivative_coeffs
export bulk_viscosity_coefficients, compute_B_bracket
export dP_dT, dP_dmu
export flavor_pressure_derivatives, conserved_charge_susceptibility
export chi_BQS, cumulant_BQS
export chi_B, chi1_B, chi2_B, chi3_B, chi4_B
export chi_Q, chi1_Q, chi2_Q, chi3_Q, chi4_Q
export chi_S, chi1_S, chi2_S, chi3_S, chi4_S
export chi11_BQ, chi11_BS, chi11_QS
export cumulant_B, baryon_Ssigma, baryon_kappa_sigma2
export default_hard_constraint_rules
export evaluate_hard_constraints
export select_pressure_max_candidate
export solve_constraint
export solve, solve_multi, SolverResult
export SOLVER_CONTRACT_VERSION_V1, SOLVER_RESULT_REQUIRED_FIELDS
export solver_contract_version, to_namedtuple
export coerce_solver_result, solver_result_view, solver_result_is_success
export solve_vec, solve_named
export ModelStateSchema, schema_for_model, flatten_state, unflatten_state
export ConstraintModes
export state_var_dim, mu_var_dim, solution_dim
export ProblemSpec, build_problem_spec, ExtraConstraints, default_extra_constraints
export AbstractConstraintComponent
export constraint_name, constraint_dim, build_constraint_components, constraint_total_dim
export HardRule, CandidateSelector, build_candidate_context
export evaluate_candidate_success, normalize_governance_candidate, build_seed_pool
export classify_attempt_error, normalize_error_message
export governance_quality_tag, normalize_selector_candidates, execute_governance_selector
export build_governance_candidate
export SolverDiagnosticSummary, SolverDiagnosticCandidate, SolverDiagnosticFull
export SOLVER_DIAGNOSTIC_VERSION_V1, SOLVER_DIAGNOSTIC_PUBLIC_FIELDS, SOLVER_DIAGNOSTIC_INTERNAL_FIELDS
export solver_diagnostic_version, solver_diagnostic_public_fields, solver_diagnostic_internal_fields
export to_public_namedtuple
export coerce_solver_diagnostic_summary, coerce_solver_diagnostic_public_view
export VarSchema, SchemaRegistry, register_schema!, schema_for, validate_schema
export named_to_vec, vec_to_named
export state_view, mu_view
export PrimaryStrategy
export SeedStrategy, DefaultSeed, MultiSeed, HybridContinuitySeed
export get_seed, update!, get_all_seeds
export HADRON_SEED_5, QUARK_SEED_5, HADRON_SEED_8, QUARK_SEED_8
export build_conditions, build_residual!, gap_core_residual!, GapParams
export explicit_residual, explicit_residual!
export solve_weighted_block_fallback
export solve_with_derivatives
export is_physical_solution
export RootProblemSpec, RootPolicy, ContinuationState, RootAttempt, RootDiagnostics, RootSolveResult
export solve_root_with_policy, solve_root_continuation
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
export default_scan_numeric_options, solve_pnjl_point
export auto_phase_hint
export default_momentum_count, default_theta_count
export default_momentum_nodes, default_momentum_weights
export pnjl_module
export solve_gap_and_transport, solve_transport_from_equilibrium
export solve_gap_and_meson_point
export solve_rotation_point
export solve_gas_liquid_point
export transport_workflow_module, meson_workflow_module
export gas_liquid_workflow_module, rotation_workflow_module
export workflow_param_adapters_module
export workflow_module_for
export magnetic_thermodynamics_module
export run_phase_pipeline, run_production_phase_pipeline, find_cep, build_phase_artifacts
export resolve_phase_output_target, promote_phase_artifacts
export CEPResult, FirstOrderSweepResult, ProductionPipelineConfig, PromotionResult, PhasePipelineResult
export AbstractPipelineIOContract, PipelineIOContract
export PipelineProvenance, PipelineSpec, PipelineStage, PipelineContext, PipelineArtifact, StageResult
export persisted_symbol_to_string, persisted_string_to_symbol
export PipelineStageRecord, PipelineRunResult
export run_pipeline, compute_pipeline_config_hash, compute_pipeline_artifact_hash
export build_manifest_extensions
export normalize_adapter_kwargs, classify_pipeline_error
export PHASE_PIPELINE_STAGE_IDS, run_phase_pipeline_via_runner
export workflow_pipeline_stage_ids, scan_pipeline_stage_ids, relaxtime_orchestrator_stage_ids
export run_workflow_pipeline
export run_scan_pipeline
export run_relaxtime_orchestrator_pipeline
export PM_BRANCH_STATUSES, PM_SEED_SOURCES, PM_ENDPOINT_CAUSES, PM_COMPARISON_STATUSES
export PMSeedPair, normalize_pm_seed_pair, pm_next_seed_source
export derive_pm_seed_pair, analyze_pm_branch_competition
export Integrals, cached_nodes, vacuum_integral, calculate_energy_sum, calculate_number_densities
export Constants_PNJL
export TmuScanConfig, TrhoScanConfig
export run_precompile_workload
export run_precompile_profile, run_precompile_capability
export list_precompile_capabilities, list_precompile_profile
export ThermoDiffContext, DiffTarget, ParamSpec
export build_thermo_diff_context, diff_target, jacobian
export build_pilot_diff_context, eval_pilot_derivatives

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
include(joinpath(@__DIR__, "variants", "gas_liquid", "GasLiquidModel.jl"))
include(joinpath(@__DIR__, "variants", "gas_liquid", "workflows", "GasLiquidWorkflow.jl"))
include(joinpath(@__DIR__, "variants", "rotation", "RotationModel.jl"))
include(joinpath(@__DIR__, "variants", "rotation", "workflows", "RotationWorkflow.jl"))

# Backward-compatible access path used by some tests/callers:
# Main.Models.PNJLIntegrals.*
const PNJLIntegrals = PNJLCore.PNJLIntegrals

# Stable accessor contract for workflow-side thermal grid defaults.
@inline default_momentum_count() = PNJLCore.DEFAULT_MOMENTUM_COUNT
@inline default_theta_count() = PNJLCore.DEFAULT_THETA_COUNT
@inline default_momentum_nodes() = PNJLIntegrals.THERMAL_DEFAULT_NODES
@inline default_momentum_weights() = PNJLIntegrals.THERMAL_DEFAULT_WEIGHTS

@inline function cached_nodes(
	p_num::Int=PNJLCore.DEFAULT_MOMENTUM_COUNT,
	t_num::Int=PNJLCore.DEFAULT_THETA_COUNT,
)
	return PNJLCore.cached_nodes(p_num, t_num)
end

# Factory
include(joinpath(@__DIR__, "factory.jl"))

# Ω assembly
include(joinpath(@__DIR__, "omega.jl"))
include(joinpath(@__DIR__, "thermo_kernel.jl"))
include(joinpath(@__DIR__, "solver", "topology.jl"))
include(joinpath(@__DIR__, "derivatives", "ThermoDerivatives.jl"))
include(joinpath(@__DIR__, "derivatives", "ConservedChargeSusceptibilities.jl"))
include(joinpath(@__DIR__, "pnjl_physics", "core", "MagneticIntegrals.jl"))
include(joinpath(@__DIR__, "pnjl_physics", "core", "MagneticThermodynamics.jl"))
include(joinpath(@__DIR__, "scans", "ScanCommon.jl"))
include(joinpath(@__DIR__, "scans", "ScanConfig.jl"))
include(joinpath(@__DIR__, "scans", "ScanResultFinalize.jl"))
include(joinpath(@__DIR__, "scans", "TmuScan.jl"))
include(joinpath(@__DIR__, "scans", "TrhoScan.jl"))

using .SeedStrategies
using .Conditions
const ConstraintModes = @__MODULE__
using .ThermoDerivatives
using .ConservedChargeSusceptibilities
using .MagneticIntegrals
using .MagneticThermodynamics
using .TmuScan
using .TrhoScan

# Transport provider (distribution/dispersion) for Stage-4 workflow decoupling
include(joinpath(@__DIR__, "transport_provider.jl"))

const ρ0 = Main.Constants_PNJL.ρ0_inv_fm3

# Unified scan/workflow entrypoints (Stage C compatibility layer)
include(joinpath(@__DIR__, "phase", "PhaseTypes.jl"))
include(joinpath(@__DIR__, "phase", "PMPhaseTypes.jl"))
include(joinpath(@__DIR__, "phase", "PMPhaseSeeds.jl"))
include(joinpath(@__DIR__, "phase", "PMPhaseDiagnostic.jl"))
include(joinpath(@__DIR__, "phase", "PMPhaseArtifacts.jl"))
include(joinpath(@__DIR__, "phase", "AdaptiveRhoRefinement.jl"))
include(joinpath(@__DIR__, "phase", "PhaseCore.jl"))
include(joinpath(@__DIR__, "phase", "PhaseIO.jl"))
include(joinpath(@__DIR__, "phase", "CEPDetector.jl"))
include(joinpath(@__DIR__, "phase", "CrossoverLine.jl"))
include(joinpath(@__DIR__, "phase", "PhaseArtifacts.jl"))
include(joinpath(@__DIR__, "phase", "PhasePipeline.jl"))
include(joinpath(@__DIR__, "phase", "ProductionPhasePipeline.jl"))
include(joinpath(@__DIR__, "workflow", "PipelineTypes.jl"))
include(joinpath(@__DIR__, "workflow", "io", "ManifestExtensions.jl"))
include(joinpath(@__DIR__, "workflow", "PipelineRunner.jl"))
include(joinpath(@__DIR__, "workflow", "adapters", "CommonAdapterUtils.jl"))
include(joinpath(@__DIR__, "workflow", "adapters", "WorkflowAdapter.jl"))
include(joinpath(@__DIR__, "workflow", "adapters", "ScanAdapter.jl"))
include(joinpath(@__DIR__, "workflow", "adapters", "RelaxtimeOrchestratorAdapter.jl"))
include(joinpath(@__DIR__, "workflow", "StageCatalog.jl"))
include(joinpath(@__DIR__, "workflow", "catalog", "WorkflowCatalog.jl"))
include(joinpath(@__DIR__, "workflow", "catalog", "ScanCatalog.jl"))
include(joinpath(@__DIR__, "workflow", "catalog", "RelaxtimeOrchestratorCatalog.jl"))
include(joinpath(@__DIR__, "workflows", "WorkflowParamAdapters.jl"))
include(joinpath(@__DIR__, "workflows", "TransportWorkflow.jl"))
include(joinpath(@__DIR__, "workflows", "MesonMassWorkflow.jl"))

include(joinpath(@__DIR__, "entrypoints.jl"))
include(joinpath(@__DIR__, "precompile", "registry.jl"))
include(joinpath(@__DIR__, "precompile_workload.jl"))

using .PrecompileRegistry: run_precompile_profile, run_precompile_capability
using .PrecompileRegistry: list_precompile_capabilities, list_precompile_profile

const Integrals = PNJLCore.PNJLIntegrals
const Constants_PNJL = Main.Constants_PNJL
const TmuScanConfig = ScanConfig.TmuScanConfig
const TrhoScanConfig = ScanConfig.TrhoScanConfig
const update! = SeedStrategies.update!

@inline function vacuum_integral(mass)
	TT = typeof(mass)
	Λ = convert(TT, Main.Constants_PNJL.Λ_inv_fm)
	return PNJLCore.vacuum_integral_with_cutoff(convert(TT, mass), Λ)
end

@inline function calculate_energy_sum(masses)
	total = zero(eltype(masses))
	@inbounds for i in eachindex(masses)
		total += vacuum_integral(masses[i])
	end
	return -2 * Main.Constants_PNJL.N_color * total
end

end # module Models
