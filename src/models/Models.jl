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

include(joinpath(@__DIR__, "exports_public.jl"))

# Keep key entrypoint exports explicitly in this file for governance checks.
export run_tmu_scan, run_trho_scan, run_magnetic_scan, build_default_rho_grid
export run_freezeout_fixedmu_scan
export run_meson_mass_path_scan
export run_freezeout_meson_mass_scan
export run_isentropic_meson_mass_scan
export run_freezeout_meson_density_scan
export run_crossover_meson_density_scan
export run_external_path_meson_density_scan
export default_scan_numeric_options, solve_pnjl_point
export auto_phase_hint
export solve_gap_and_transport, solve_transport_from_equilibrium
export solve_gap_and_meson_point, solve_meson_point_from_equilibrium, solve_meson_density_from_meson_point, solve_gap_and_meson_density_point
export solve_strict_bw_meson_density_from_meson_point, solve_gap_and_strict_bw_meson_density_point
export solve_phase_shift_meson_density_from_meson_point, solve_gap_and_phase_shift_meson_density_point
export solve_phase_shift_point_diagnostic_from_meson_point
export solve_phase_shift_derivative_reference_from_meson_point
export solve_meson_thermo_from_meson_point, solve_gap_and_meson_thermo_point
export solve_strict_bw_meson_thermo_from_meson_point, solve_gap_and_strict_bw_meson_thermo_point
export solve_phase_shift_meson_thermo_from_meson_point, solve_gap_and_phase_shift_meson_thermo_point
export build_meson_thermo_contract_row
export solve_gas_liquid_point, solve_rotation_point
export GasLiquidModel
export solve_gas_liquid_rmf_point, run_gas_liquid_tmu_scan, run_gas_liquid_trho_scan
export build_gas_liquid_result_row, build_gas_liquid_manifest
export transport_workflow_module, meson_workflow_module, meson_density_workflow_module, meson_thermo_workflow_module
export rotation_workflow_module, gas_liquid_workflow_module
export workflow_param_adapters_module, workflow_module_for
export run_workflow_pipeline, run_scan_pipeline, run_relaxtime_orchestrator_pipeline
export pnjl_module, magnetic_thermodynamics_module
export run_phase_pipeline, run_production_phase_pipeline, find_cep, build_phase_artifacts
export resolve_phase_output_target, promote_phase_artifacts
export normalize_pm_seed_pair, pm_next_seed_source, derive_pm_seed_pair, analyze_pm_branch_competition

include(joinpath(@__DIR__, "abstract_model.jl"))

# State representation
include(joinpath(@__DIR__, "state.jl"))

# Concrete models
include(joinpath(@__DIR__, "njl", "NJLModel.jl"))
include(joinpath(@__DIR__, "njl", "NJL2Model.jl"))
include(joinpath(@__DIR__, "pnjl_physics", "PNJLDistributions.jl"))
include(joinpath(@__DIR__, "pnjl_physics", "PNJLCore.jl"))
include(joinpath(@__DIR__, "pnjl_physics", "PNJLModel.jl"))
# Magnetic core must be defined inside Models before its adapter is loaded;
# this avoids injecting a second copy into Main during module initialization.
include(joinpath(@__DIR__, "pnjl_physics", "core", "MagneticIntegrals.jl"))
include(joinpath(@__DIR__, "pnjl_physics", "core", "MagneticThermodynamics.jl"))
include(joinpath(@__DIR__, "pnjl_physics", "PNJLMagneticModel.jl"))
include(joinpath(@__DIR__, "rpnjl", "RPNJLModel.jl"))
include(joinpath(@__DIR__, "variants", "gas_liquid", "GasLiquidModel.jl"))
include(joinpath(@__DIR__, "variants", "gas_liquid", "workflows", "GasLiquidWorkflow.jl"))
include(joinpath(@__DIR__, "variants", "rotation", "RotationModel.jl"))
include(joinpath(@__DIR__, "variants", "rotation", "workflows", "RotationWorkflow.jl"))

# Canonical per-model homomorphic anchors
include(joinpath(@__DIR__, "njl", "capabilities.jl"))
include(joinpath(@__DIR__, "njl", "api.jl"))
include(joinpath(@__DIR__, "njl", "adapters", "entrypoint_adapter.jl"))
include(joinpath(@__DIR__, "njl", "workflows", "noop.jl"))

include(joinpath(@__DIR__, "njl2", "capabilities.jl"))
include(joinpath(@__DIR__, "njl2", "api.jl"))
include(joinpath(@__DIR__, "njl2", "adapters", "entrypoint_adapter.jl"))
include(joinpath(@__DIR__, "njl2", "workflows", "noop.jl"))

include(joinpath(@__DIR__, "pnjl", "capabilities.jl"))
include(joinpath(@__DIR__, "pnjl", "api.jl"))
include(joinpath(@__DIR__, "pnjl", "adapters", "entrypoint_adapter.jl"))
include(joinpath(@__DIR__, "pnjl", "workflows", "noop.jl"))

include(joinpath(@__DIR__, "rpnjl", "capabilities.jl"))
include(joinpath(@__DIR__, "rpnjl", "api.jl"))
include(joinpath(@__DIR__, "rpnjl", "adapters", "entrypoint_adapter.jl"))
include(joinpath(@__DIR__, "rpnjl", "workflows", "noop.jl"))

include(joinpath(@__DIR__, "pnjl_magnetic", "capabilities.jl"))
include(joinpath(@__DIR__, "pnjl_magnetic", "api.jl"))
include(joinpath(@__DIR__, "pnjl_magnetic", "adapters", "entrypoint_adapter.jl"))
include(joinpath(@__DIR__, "pnjl_magnetic", "workflows", "noop.jl"))

include(joinpath(@__DIR__, "rotation", "capabilities.jl"))
include(joinpath(@__DIR__, "rotation", "api.jl"))
include(joinpath(@__DIR__, "rotation", "adapters", "entrypoint_adapter.jl"))
include(joinpath(@__DIR__, "rotation", "workflows", "noop.jl"))

include(joinpath(@__DIR__, "gas_liquid", "capabilities.jl"))
include(joinpath(@__DIR__, "gas_liquid", "api.jl"))
include(joinpath(@__DIR__, "gas_liquid", "adapters", "entrypoint_adapter.jl"))
include(joinpath(@__DIR__, "gas_liquid", "workflows", "noop.jl"))

# Backward-compatible access path used by some tests/callers:
# Models.PNJLIntegrals.*
const PNJLIntegrals = PNJLCore.PNJLIntegrals

# Stable accessor contract for workflow-side thermal grid defaults.
@inline default_momentum_count() = PNJLCore.DEFAULT_MOMENTUM_COUNT
@inline default_theta_count() = PNJLCore.DEFAULT_THETA_COUNT
@inline default_momentum_nodes() = PNJLIntegrals.THERMAL_DEFAULT_NODES
@inline default_momentum_weights() = PNJLIntegrals.THERMAL_DEFAULT_WEIGHTS
@inline thermal_p_max_inv_fm(::AbstractQCDModel) = PNJLIntegrals.DEFAULT_THERMAL_P_MAX_INV_FM
@inline thermal_p_max_inv_fm(model::PNJLModel) = model.params.thermal_p_max_inv_fm
@inline thermal_p_max_inv_fm(model::PNJLMagneticModel) = model.base.params.thermal_p_max_inv_fm
@inline thermal_p_max_inv_fm(model::RPNJLModel) = model.base.params.thermal_p_max_inv_fm

@inline function cached_nodes(
	p_num::Int=PNJLCore.DEFAULT_MOMENTUM_COUNT,
	t_num::Int=PNJLCore.DEFAULT_THETA_COUNT,
	;
	p_max_inv_fm::Float64=PNJLIntegrals.DEFAULT_THERMAL_P_MAX_INV_FM,
)
	return PNJLCore.cached_nodes(p_num, t_num; p_max_inv_fm=p_max_inv_fm)
end

# Factory
include(joinpath(@__DIR__, "factory.jl"))

# Ω assembly
include(joinpath(@__DIR__, "omega.jl"))
include(joinpath(@__DIR__, "thermo_kernel.jl"))
include(joinpath(@__DIR__, "solver", "topology.jl"))
include(joinpath(@__DIR__, "derivatives", "TaylorDiffForwardDiffCompat.jl"))
include(joinpath(@__DIR__, "derivatives", "MixedTaylorJets.jl"))
include(joinpath(@__DIR__, "derivatives", "PNJLChiBTaylorDiff.jl"))
include(joinpath(@__DIR__, "derivatives", "ThermoDerivatives.jl"))
include(joinpath(@__DIR__, "derivatives", "HigherOrderDerivatives.jl"))
include(joinpath(@__DIR__, "derivatives", "AbstractSusceptibilityProvider.jl"))
include(joinpath(@__DIR__, "derivatives", "ConservedChargeSusceptibilities.jl"))
include(joinpath(@__DIR__, "scans", "ScanCommon.jl"))
include(joinpath(@__DIR__, "scans", "ScanConfig.jl"))
include(joinpath(@__DIR__, "scans", "ScanResultFinalize.jl"))
include(joinpath(@__DIR__, "scans", "FreezeoutProfiles.jl"))
include(joinpath(@__DIR__, "scans", "FreezeoutPathProfiles.jl"))
include(joinpath(@__DIR__, "scans", "IsentropicPathProfiles.jl"))
include(joinpath(@__DIR__, "scans", "FlavorChemicalProfiles.jl"))
include(joinpath(@__DIR__, "scans", "MesonChemicalProfiles.jl"))
include(joinpath(@__DIR__, "scans", "TmuScan.jl"))
include(joinpath(@__DIR__, "scans", "FreezeoutPathScan.jl"))
include(joinpath(@__DIR__, "scans", "TrhoScan.jl"))
include(joinpath(@__DIR__, "scans", "MagneticScan.jl"))

using .SeedStrategies
using .Conditions
const ConstraintModes = @__MODULE__
using .ThermoDerivatives
using .HigherOrderDerivatives
using .AbstractSusceptibilityProvider
using .ConservedChargeSusceptibilities
using .MagneticIntegrals
using .MagneticThermodynamics
using .FreezeoutProfiles
using .FreezeoutPathProfiles
using .IsentropicPathProfiles
using .FlavorChemicalProfiles
using .MesonChemicalProfiles
using .TmuScan
using .FreezeoutPathScan
using .TrhoScan
using .MagneticScan

# Transport provider (distribution/dispersion) for Stage-4 workflow decoupling
include(joinpath(@__DIR__, "transport_provider.jl"))

const ρ0 = Main.Constants_PNJL.ρ0_inv_fm3

# Unified scan/workflow entrypoints (Stage C compatibility layer)
include(joinpath(@__DIR__, "phase", "PhaseCore.jl"))
include(joinpath(@__DIR__, "phase", "RhoSupportRefinement.jl"))
include(joinpath(@__DIR__, "phase", "PhaseTypes.jl"))
include(joinpath(@__DIR__, "phase", "PhaseGridConvergence.jl"))
include(joinpath(@__DIR__, "phase", "PMPhaseTypes.jl"))
include(joinpath(@__DIR__, "phase", "PMPhaseSeeds.jl"))
include(joinpath(@__DIR__, "phase", "PMPhaseDiagnostic.jl"))
include(joinpath(@__DIR__, "phase", "PMPhaseArtifacts.jl"))
include(joinpath(@__DIR__, "phase", "AdaptiveRhoRefinement.jl"))
include(joinpath(@__DIR__, "phase", "PhaseIO.jl"))
include(joinpath(@__DIR__, "phase", "CEPDetector.jl"))
include(joinpath(@__DIR__, "phase", "CrossoverLine.jl"))
include(joinpath(@__DIR__, "phase", "PhaseArtifacts.jl"))
include(joinpath(@__DIR__, "phase", "PhasePipeline.jl"))
include(joinpath(@__DIR__, "phase", "ProductionPhasePipeline.jl"))
include(joinpath(@__DIR__, "workflow_engine", "PipelineTypes.jl"))
include(joinpath(@__DIR__, "workflow_engine", "io", "ManifestExtensions.jl"))
include(joinpath(@__DIR__, "workflow_engine", "PipelineRunner.jl"))
include(joinpath(@__DIR__, "workflow_engine", "adapters", "CommonAdapterUtils.jl"))
include(joinpath(@__DIR__, "workflow_engine", "adapters", "WorkflowAdapter.jl"))
include(joinpath(@__DIR__, "workflow_engine", "adapters", "ScanAdapter.jl"))
include(joinpath(@__DIR__, "workflow_engine", "adapters", "RelaxtimeOrchestratorAdapter.jl"))
include(joinpath(@__DIR__, "workflow_engine", "StageCatalog.jl"))
include(joinpath(@__DIR__, "workflow_engine", "catalog", "WorkflowCatalog.jl"))
include(joinpath(@__DIR__, "workflow_engine", "catalog", "ScanCatalog.jl"))
include(joinpath(@__DIR__, "workflow_engine", "catalog", "RelaxtimeOrchestratorCatalog.jl"))
include(joinpath(@__DIR__, "workflow_apps", "WorkflowParamAdapters.jl"))
include(joinpath(@__DIR__, "workflow_apps", "TransportWorkflow.jl"))
include(joinpath(@__DIR__, "workflow_apps", "MesonMassWorkflow.jl"))
include(joinpath(@__DIR__, "workflow_apps", "MesonDensityWorkflow.jl"))
include(joinpath(@__DIR__, "workflow_apps", "MesonThermoWorkflow.jl"))
include(joinpath(@__DIR__, "scans", "FreezeoutMesonDensityScan.jl"))
include(joinpath(@__DIR__, "scans", "MesonMassPathScan.jl"))
include(joinpath(@__DIR__, "scans", "CrossoverMesonDensityScan.jl"))
include(joinpath(@__DIR__, "scans", "ExternalPathMesonDensityScan.jl"))

include(joinpath(@__DIR__, "entrypoints.jl"))
include(joinpath(@__DIR__, "precompile", "registry.jl"))
include(joinpath(@__DIR__, "precompile_workload.jl"))

using .PrecompileRegistry: run_precompile_profile, run_precompile_capability
using .PrecompileRegistry: list_precompile_capabilities, list_precompile_profile

const Integrals = PNJLCore.PNJLIntegrals
const Constants_PNJL = Main.Constants_PNJL
const TmuScanConfig = ScanConfig.TmuScanConfig
const TrhoScanConfig = ScanConfig.TrhoScanConfig
const FreezeoutScanConfig = ScanConfig.FreezeoutScanConfig
const RhoSupportConfig = RhoSupportRefinement.RhoSupportConfig
const RhoSupportPrior = RhoSupportRefinement.RhoSupportPrior
const RhoSupportAssessment = RhoSupportRefinement.RhoSupportAssessment
const analyze_rho_support_cascade = RhoSupportRefinement.analyze_rho_support_cascade
const MesonChemicalProfile = MesonChemicalProfiles.MesonChemicalProfile
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
