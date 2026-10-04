# Dependency graph generated: 2026-10-04T10:29:08.886

Run: julia --project=. scripts/dev/gen_deps.jl

阅读入口：[职责与调用方向](dependencies.manual.md) · [依赖规则](dependency_rules.md)。

本图静态解析 include 的字面量、路径常量及 joinpath/normpath/dirname/@__DIR__；不执行源码。
相对 using/import 只显示模块名，尚未解析完整模块身份；条件分支合并展示。无环结果仅适用于已解析的边，不能证明完整运行时依赖无环。

## 未解析或缺失的 include

- `src/models/workflow_engine/adapters/WorkflowAdapter.jl: Base.include(Main, lib_path)`

## 静态依赖图

```mermaid
%%{init: { 'flowchart': { 'nodeSpacing': 40, 'rankSpacing': 60, 'useMaxWidth': false } }}%%
flowchart LR
  subgraph config
    src_config_ConfigLoader_jl[config/ConfigLoader.jl]
  end
  subgraph constants
    src_constants_Constants_PNJL_jl[constants/Constants_PNJL.jl]
    src_constants_TransportConstants_jl[constants/TransportConstants.jl]
  end
  subgraph integration
    src_integration_CauchyPV_jl[integration/CauchyPV.jl]
    src_integration_GaussLegendre_jl[integration/GaussLegendre.jl]
    src_integration_IntervalQuadratureStrategies_jl[integration/IntervalQuadratureStrategies.jl]
    src_integration_PhaseSpaceSampling_jl[integration/PhaseSpaceSampling.jl]
  end
  subgraph models
    src_models_Models_jl[models/Models.jl]
    src_models_abstract_model_jl[models/abstract_model.jl]
    src_models_derivatives_AbstractSusceptibilityProvider_jl[models/derivatives/AbstractSusceptibilityProvider.jl]
    src_models_derivatives_ConservedChargeSusceptibilities_jl[models/derivatives/ConservedChargeSusceptibilities.jl]
    src_models_derivatives_DiffService_jl[models/derivatives/DiffService.jl]
    src_models_derivatives_HigherOrderDerivatives_jl[models/derivatives/HigherOrderDerivatives.jl]
    src_models_derivatives_MixedTaylorJets_jl[models/derivatives/MixedTaylorJets.jl]
    src_models_derivatives_PNJLChiBTaylorDiff_jl[models/derivatives/PNJLChiBTaylorDiff.jl]
    src_models_derivatives_TaylorDiffForwardDiffCompat_jl[models/derivatives/TaylorDiffForwardDiffCompat.jl]
    src_models_derivatives_ThermoDerivatives_jl[models/derivatives/ThermoDerivatives.jl]
    src_models_entrypoints_jl[models/entrypoints.jl]
    src_models_exports_public_jl[models/exports_public.jl]
    src_models_factory_jl[models/factory.jl]
    src_models_gas_liquid_adapters_entrypoint_adapter_jl[models/gas_liquid/adapters/entrypoint_adapter.jl]
    src_models_gas_liquid_api_jl[models/gas_liquid/api.jl]
    src_models_gas_liquid_capabilities_jl[models/gas_liquid/capabilities.jl]
    src_models_njl_NJL2Core_jl[models/njl/NJL2Core.jl]
    src_models_njl_NJL2Model_jl[models/njl/NJL2Model.jl]
    src_models_njl_NJLCore_jl[models/njl/NJLCore.jl]
    src_models_njl_NJLModel_jl[models/njl/NJLModel.jl]
    src_models_njl_adapters_entrypoint_adapter_jl[models/njl/adapters/entrypoint_adapter.jl]
    src_models_njl_api_jl[models/njl/api.jl]
    src_models_njl_capabilities_jl[models/njl/capabilities.jl]
    src_models_njl2_adapters_entrypoint_adapter_jl[models/njl2/adapters/entrypoint_adapter.jl]
    src_models_njl2_api_jl[models/njl2/api.jl]
    src_models_njl2_capabilities_jl[models/njl2/capabilities.jl]
    src_models_omega_jl[models/omega.jl]
    src_models_phase_AdaptiveRhoRefinement_jl[models/phase/AdaptiveRhoRefinement.jl]
    src_models_phase_CEPDetector_jl[models/phase/CEPDetector.jl]
    src_models_phase_CrossoverLine_jl[models/phase/CrossoverLine.jl]
    src_models_phase_PMPhaseArtifacts_jl[models/phase/PMPhaseArtifacts.jl]
    src_models_phase_PMPhaseDiagnostic_jl[models/phase/PMPhaseDiagnostic.jl]
    src_models_phase_PMPhaseSeeds_jl[models/phase/PMPhaseSeeds.jl]
    src_models_phase_PMPhaseTypes_jl[models/phase/PMPhaseTypes.jl]
    src_models_phase_PhaseArtifacts_jl[models/phase/PhaseArtifacts.jl]
    src_models_phase_PhaseCore_jl[models/phase/PhaseCore.jl]
    src_models_phase_PhaseGridConvergence_jl[models/phase/PhaseGridConvergence.jl]
    src_models_phase_PhaseIO_jl[models/phase/PhaseIO.jl]
    src_models_phase_PhasePipeline_jl[models/phase/PhasePipeline.jl]
    src_models_phase_PhaseTypes_jl[models/phase/PhaseTypes.jl]
    src_models_phase_ProductionPhasePipeline_jl[models/phase/ProductionPhasePipeline.jl]
    src_models_phase_RhoSupportRefinement_jl[models/phase/RhoSupportRefinement.jl]
    src_models_pnjl_adapters_entrypoint_adapter_jl[models/pnjl/adapters/entrypoint_adapter.jl]
    src_models_pnjl_api_jl[models/pnjl/api.jl]
    src_models_pnjl_capabilities_jl[models/pnjl/capabilities.jl]
    src_models_pnjl_magnetic_adapters_entrypoint_adapter_jl[models/pnjl_magnetic/adapters/entrypoint_adapter.jl]
    src_models_pnjl_magnetic_api_jl[models/pnjl_magnetic/api.jl]
    src_models_pnjl_magnetic_capabilities_jl[models/pnjl_magnetic/capabilities.jl]
    src_models_pnjl_physics_MagneticGapSolver_jl[models/pnjl_physics/MagneticGapSolver.jl]
    src_models_pnjl_physics_PNJLCore_jl[models/pnjl_physics/PNJLCore.jl]
    src_models_pnjl_physics_PNJLDistributions_jl[models/pnjl_physics/PNJLDistributions.jl]
    src_models_pnjl_physics_PNJLIntegrals_jl[models/pnjl_physics/PNJLIntegrals.jl]
    src_models_pnjl_physics_PNJLMagneticModel_jl[models/pnjl_physics/PNJLMagneticModel.jl]
    src_models_pnjl_physics_PNJLModel_jl[models/pnjl_physics/PNJLModel.jl]
    src_models_pnjl_physics_QuarkDistribution_jl[models/pnjl_physics/QuarkDistribution.jl]
    src_models_pnjl_physics_core_EquilibriumFacade_jl[models/pnjl_physics/core/EquilibriumFacade.jl]
    src_models_pnjl_physics_core_Integrals_jl[models/pnjl_physics/core/Integrals.jl]
    src_models_pnjl_physics_core_MagneticIntegrals_jl[models/pnjl_physics/core/MagneticIntegrals.jl]
    src_models_pnjl_physics_core_MagneticThermodynamics_jl[models/pnjl_physics/core/MagneticThermodynamics.jl]
    src_models_pnjl_physics_core_ModelThermodynamics_jl[models/pnjl_physics/core/ModelThermodynamics.jl]
    src_models_precompile_registry_jl[models/precompile/registry.jl]
    src_models_precompile_workload_jl[models/precompile_workload.jl]
    src_models_rotation_adapters_entrypoint_adapter_jl[models/rotation/adapters/entrypoint_adapter.jl]
    src_models_rotation_api_jl[models/rotation/api.jl]
    src_models_rotation_capabilities_jl[models/rotation/capabilities.jl]
    src_models_rpnjl_RPNJLModel_jl[models/rpnjl/RPNJLModel.jl]
    src_models_rpnjl_adapters_entrypoint_adapter_jl[models/rpnjl/adapters/entrypoint_adapter.jl]
    src_models_rpnjl_api_jl[models/rpnjl/api.jl]
    src_models_rpnjl_capabilities_jl[models/rpnjl/capabilities.jl]
    src_models_scans_CrossoverMesonDensityScan_jl[models/scans/CrossoverMesonDensityScan.jl]
    src_models_scans_ExternalPathMesonDensityScan_jl[models/scans/ExternalPathMesonDensityScan.jl]
    src_models_scans_FlavorChemicalProfiles_jl[models/scans/FlavorChemicalProfiles.jl]
    src_models_scans_FreezeoutMesonDensityScan_jl[models/scans/FreezeoutMesonDensityScan.jl]
    src_models_scans_FreezeoutPathProfiles_jl[models/scans/FreezeoutPathProfiles.jl]
    src_models_scans_FreezeoutPathScan_jl[models/scans/FreezeoutPathScan.jl]
    src_models_scans_FreezeoutProfiles_jl[models/scans/FreezeoutProfiles.jl]
    src_models_scans_IsentropicPathProfiles_jl[models/scans/IsentropicPathProfiles.jl]
    src_models_scans_MagneticScan_jl[models/scans/MagneticScan.jl]
    src_models_scans_MesonChemicalProfiles_jl[models/scans/MesonChemicalProfiles.jl]
    src_models_scans_MesonMassPathScan_jl[models/scans/MesonMassPathScan.jl]
    src_models_scans_ScanCommon_jl[models/scans/ScanCommon.jl]
    src_models_scans_ScanConfig_jl[models/scans/ScanConfig.jl]
    src_models_scans_ScanResultFinalize_jl[models/scans/ScanResultFinalize.jl]
    src_models_scans_TmuScan_jl[models/scans/TmuScan.jl]
    src_models_scans_TrhoScan_jl[models/scans/TrhoScan.jl]
    src_models_solver_ImplicitProblem_jl[models/solver/ImplicitProblem.jl]
    src_models_solver_api_SolverAPI_jl[models/solver/api/SolverAPI.jl]
    src_models_solver_compat_ImplicitAdapters_jl[models/solver/compat/ImplicitAdapters.jl]
    src_models_solver_compat_ImplicitGapLegacy_jl[models/solver/compat/ImplicitGapLegacy.jl]
    src_models_solver_compat_SchemaAdapter_jl[models/solver/compat/SchemaAdapter.jl]
    src_models_solver_config_SolverRuntimeConfig_jl[models/solver/config/SolverRuntimeConfig.jl]
    src_models_solver_config_StateSchema_jl[models/solver/config/StateSchema.jl]
    src_models_solver_diagnostics_SolverDiagnostics_jl[models/solver/diagnostics/SolverDiagnostics.jl]
    src_models_solver_diagnostics_SolverDiagnosticsTypes_jl[models/solver/diagnostics/SolverDiagnosticsTypes.jl]
    src_models_solver_diagnostics_SolverWorkTelemetry_jl[models/solver/diagnostics/SolverWorkTelemetry.jl]
    src_models_solver_diagnostics_ThermoPostprocess_jl[models/solver/diagnostics/ThermoPostprocess.jl]
    src_models_solver_diff_JacobianEngine_jl[models/solver/diff/JacobianEngine.jl]
    src_models_solver_diff_PilotAdapters_jl[models/solver/diff/PilotAdapters.jl]
    src_models_solver_diff_Targets_jl[models/solver/diff/Targets.jl]
    src_models_solver_diff_ThermoDiffContext_jl[models/solver/diff/ThermoDiffContext.jl]
    src_models_solver_governance_CandidateGovernance_jl[models/solver/governance/CandidateGovernance.jl]
    src_models_solver_governance_WeightedFallback_jl[models/solver/governance/WeightedFallback.jl]
    src_models_solver_orchestrator_PrimaryStrategy_jl[models/solver/orchestrator/PrimaryStrategy.jl]
    src_models_solver_orchestrator_ProblemSpecOrchestrator_jl[models/solver/orchestrator/ProblemSpecOrchestrator.jl]
    src_models_solver_orchestrator_SeedStrategies_jl[models/solver/orchestrator/SeedStrategies.jl]
    src_models_solver_path_PathContinuation_jl[models/solver/path/PathContinuation.jl]
    src_models_solver_runtime_ConstraintSolver_jl[models/solver/runtime/ConstraintSolver.jl]
    src_models_solver_runtime_ConstraintSolverCommon_jl[models/solver/runtime/ConstraintSolverCommon.jl]
    src_models_solver_runtime_ConstraintSolverFixedAsymmetricRho_jl[models/solver/runtime/ConstraintSolverFixedAsymmetricRho.jl]
    src_models_solver_runtime_ConstraintSolverFixedEntropy_jl[models/solver/runtime/ConstraintSolverFixedEntropy.jl]
    src_models_solver_runtime_ConstraintSolverFixedMu_jl[models/solver/runtime/ConstraintSolverFixedMu.jl]
    src_models_solver_runtime_ConstraintSolverFixedMuBConservedCharges_jl[models/solver/runtime/ConstraintSolverFixedMuBConservedCharges.jl]
    src_models_solver_runtime_ConstraintSolverFixedRho_jl[models/solver/runtime/ConstraintSolverFixedRho.jl]
    src_models_solver_runtime_ConstraintSolverFixedSigma_jl[models/solver/runtime/ConstraintSolverFixedSigma.jl]
    src_models_solver_runtime_GapSolver_jl[models/solver/runtime/GapSolver.jl]
    src_models_solver_runtime_GenericRootEngine_jl[models/solver/runtime/GenericRootEngine.jl]
    src_models_solver_spec_Conditions_jl[models/solver/spec/Conditions.jl]
    src_models_solver_spec_ConstraintComponents_jl[models/solver/spec/ConstraintComponents.jl]
    src_models_solver_spec_ConstraintModes_jl[models/solver/spec/ConstraintModes.jl]
    src_models_solver_spec_ProblemSpec_jl[models/solver/spec/ProblemSpec.jl]
    src_models_solver_topology_jl[models/solver/topology.jl]
    src_models_state_jl[models/state.jl]
    src_models_thermo_kernel_jl[models/thermo_kernel.jl]
    src_models_transport_provider_jl[models/transport_provider.jl]
    src_models_variants_gas_liquid_GasLiquidModel_jl[models/variants/gas_liquid/GasLiquidModel.jl]
    src_models_variants_gas_liquid_core_EquationSet_jl[models/variants/gas_liquid/core/EquationSet.jl]
    src_models_variants_gas_liquid_core_Thermodynamics_jl[models/variants/gas_liquid/core/Thermodynamics.jl]
    src_models_variants_gas_liquid_workflows_GasLiquidWorkflow_jl[models/variants/gas_liquid/workflows/GasLiquidWorkflow.jl]
    src_models_variants_rotation_RotationModel_jl[models/variants/rotation/RotationModel.jl]
    src_models_variants_rotation_core_RotationThermo_jl[models/variants/rotation/core/RotationThermo.jl]
    src_models_variants_rotation_workflows_RotationWorkflow_jl[models/variants/rotation/workflows/RotationWorkflow.jl]
    src_models_workflow_apps_ChargedGBUResearchWorkflow_jl[models/workflow_apps/ChargedGBUResearchWorkflow.jl]
    src_models_workflow_apps_MesonDensityWorkflow_jl[models/workflow_apps/MesonDensityWorkflow.jl]
    src_models_workflow_apps_MesonMassWorkflow_jl[models/workflow_apps/MesonMassWorkflow.jl]
    src_models_workflow_apps_MesonThermoWorkflow_jl[models/workflow_apps/MesonThermoWorkflow.jl]
    src_models_workflow_apps_TransportWorkflow_jl[models/workflow_apps/TransportWorkflow.jl]
    src_models_workflow_apps_WorkflowParamAdapters_jl[models/workflow_apps/WorkflowParamAdapters.jl]
    src_models_workflow_engine_PipelineRunner_jl[models/workflow_engine/PipelineRunner.jl]
    src_models_workflow_engine_PipelineTypes_jl[models/workflow_engine/PipelineTypes.jl]
    src_models_workflow_engine_StageCatalog_jl[models/workflow_engine/StageCatalog.jl]
    src_models_workflow_engine_adapters_CommonAdapterUtils_jl[models/workflow_engine/adapters/CommonAdapterUtils.jl]
    src_models_workflow_engine_adapters_RelaxtimeOrchestratorAdapter_jl[models/workflow_engine/adapters/RelaxtimeOrchestratorAdapter.jl]
    src_models_workflow_engine_adapters_ScanAdapter_jl[models/workflow_engine/adapters/ScanAdapter.jl]
    src_models_workflow_engine_adapters_WorkflowAdapter_jl[models/workflow_engine/adapters/WorkflowAdapter.jl]
    src_models_workflow_engine_catalog_RelaxtimeOrchestratorCatalog_jl[models/workflow_engine/catalog/RelaxtimeOrchestratorCatalog.jl]
    src_models_workflow_engine_catalog_ScanCatalog_jl[models/workflow_engine/catalog/ScanCatalog.jl]
    src_models_workflow_engine_catalog_WorkflowCatalog_jl[models/workflow_engine/catalog/WorkflowCatalog.jl]
    src_models_workflow_engine_io_ManifestExtensions_jl[models/workflow_engine/io/ManifestExtensions.jl]
  end
  subgraph relaxtime
    src_relaxtime_AFieldBuilder_jl[relaxtime/AFieldBuilder.jl]
    src_relaxtime_AverageScatteringRate_jl[relaxtime/AverageScatteringRate.jl]
    src_relaxtime_BUPhaseGates_jl[relaxtime/BUPhaseGates.jl]
    src_relaxtime_CausalSpectralBubble_jl[relaxtime/CausalSpectralBubble.jl]
    src_relaxtime_ChargedPhaseBackend_jl[relaxtime/ChargedPhaseBackend.jl]
    src_relaxtime_ChargedRPAKernel_jl[relaxtime/ChargedRPAKernel.jl]
    src_relaxtime_ChargedRPAProvider_jl[relaxtime/ChargedRPAProvider.jl]
    src_relaxtime_DifferentialCrossSection_jl[relaxtime/DifferentialCrossSection.jl]
    src_relaxtime_EffectiveCouplings_jl[relaxtime/EffectiveCouplings.jl]
    src_relaxtime_KinematicChecks_jl[relaxtime/KinematicChecks.jl]
    src_relaxtime_MesonDensity_jl[relaxtime/MesonDensity.jl]
    src_relaxtime_MesonInteractionKernel_jl[relaxtime/MesonInteractionKernel.jl]
    src_relaxtime_MesonMass_jl[relaxtime/MesonMass.jl]
    src_relaxtime_MesonPropagator_jl[relaxtime/MesonPropagator.jl]
    src_relaxtime_MesonRPA_jl[relaxtime/MesonRPA.jl]
    src_relaxtime_MesonRPAAdapter_jl[relaxtime/MesonRPAAdapter.jl]
    src_relaxtime_MesonThermodynamics_jl[relaxtime/MesonThermodynamics.jl]
    src_relaxtime_MottTransition_jl[relaxtime/MottTransition.jl]
    src_relaxtime_OneLoopIntegrals_jl[relaxtime/OneLoopIntegrals.jl]
    src_relaxtime_OneLoopIntegralsAniso_jl[relaxtime/OneLoopIntegralsAniso.jl]
    src_relaxtime_PhaseNormalization_jl[relaxtime/PhaseNormalization.jl]
    src_relaxtime_PolarizationAniso_jl[relaxtime/PolarizationAniso.jl]
    src_relaxtime_PolarizationCache_jl[relaxtime/PolarizationCache.jl]
    src_relaxtime_RelaxTime_jl[relaxtime/RelaxTime.jl]
    src_relaxtime_RelaxationTime_jl[relaxtime/RelaxationTime.jl]
    src_relaxtime_ScatteringAmplitude_jl[relaxtime/ScatteringAmplitude.jl]
    src_relaxtime_TotalCrossSection_jl[relaxtime/TotalCrossSection.jl]
    src_relaxtime_TotalPropagator_jl[relaxtime/TotalPropagator.jl]
    src_relaxtime_TransportCoefficients_jl[relaxtime/TransportCoefficients.jl]
    src_relaxtime_TransportCoefficientsValidation_jl[relaxtime/TransportCoefficientsValidation.jl]
  end
  subgraph root
    AFieldBuilder[AFieldBuilder]
    AbstractSusceptibilityProvider[AbstractSusceptibilityProvider]
    AverageScatteringRate[AverageScatteringRate]
    BUPhaseGates[BUPhaseGates]
    ChargedPhaseBackend[ChargedPhaseBackend]
    ChargedRPAKernel[ChargedRPAKernel]
    ChargedRPAProvider[ChargedRPAProvider]
    Conditions[Conditions]
    ConfigLoader[ConfigLoader]
    ConservedChargeSusceptibilities[ConservedChargeSusceptibilities]
    CrossSectionOrchestratedScan[CrossSectionOrchestratedScan]
    CrossoverMesonDensityScan[CrossoverMesonDensityScan]
    DifferentialCrossSection[DifferentialCrossSection]
    EffectiveCouplings[EffectiveCouplings]
    EllipsoidCalculation[EllipsoidCalculation]
    FlavorChemicalProfiles[FlavorChemicalProfiles]
    FrameTransformations[FrameTransformations]
    FreezeoutPathProfiles[FreezeoutPathProfiles]
    FreezeoutPathScan[FreezeoutPathScan]
    FreezeoutProfiles[FreezeoutProfiles]
    FullServerApp[FullServerApp]
    GasLiquidEquationSet[GasLiquidEquationSet]
    GasLiquidThermodynamics[GasLiquidThermodynamics]
    GaussLegendre[GaussLegendre]
    HigherOrderDerivatives[HigherOrderDerivatives]
    IsentropicPathProfiles[IsentropicPathProfiles]
    KinematicChecks[KinematicChecks]
    MagneticIntegrals[MagneticIntegrals]
    MagneticScan[MagneticScan]
    MagneticThermodynamics[MagneticThermodynamics]
    MesonChemicalProfiles[MesonChemicalProfiles]
    MesonDensity[MesonDensity]
    MesonDensityWorkflow[MesonDensityWorkflow]
    MesonInteractionKernel[MesonInteractionKernel]
    MesonMass[MesonMass]
    MesonMassWorkflow[MesonMassWorkflow]
    MesonPropagator[MesonPropagator]
    MesonRPA[MesonRPA]
    MesonRPAAdapter[MesonRPAAdapter]
    MesonThermodynamics[MesonThermodynamics]
    MixedTaylorJets[MixedTaylorJets]
    Models[Models]
    MomentumMapping[MomentumMapping]
    MottTransition[MottTransition]
    OneLoopIntegrals[OneLoopIntegrals]
    OneLoopIntegralsCorrection[OneLoopIntegralsCorrection]
    PNJLChiBTaylorDiff[PNJLChiBTaylorDiff]
    PNJLCore[PNJLCore]
    ParticleSymbols[ParticleSymbols]
    PhaseNormalization[PhaseNormalization]
    PhaseSpaceSampling[PhaseSpaceSampling]
    PolarizationAniso[PolarizationAniso]
    PolarizationCache[PolarizationCache]
    PrecompileRegistry[PrecompileRegistry]
    RelaxationTime[RelaxationTime]
    RotationThermo[RotationThermo]
    ScanCommon[ScanCommon]
    ScanConfig[ScanConfig]
    ScanResultFinalize[ScanResultFinalize]
    ScatteringAmplitude[ScatteringAmplitude]
    SeedStrategies[SeedStrategies]
    ServerWarmup[ServerWarmup]
    TaylorDiffForwardDiffCompat[TaylorDiffForwardDiffCompat]
    ThermoDerivatives[ThermoDerivatives]
    TmuScan[TmuScan]
    TotalCrossSection[TotalCrossSection]
    TotalPropagator[TotalPropagator]
    TransportCoefficients[TransportCoefficients]
    TransportCoefficientsValidation[TransportCoefficientsValidation]
    TransportConstants[TransportConstants]
    TrhoScan[TrhoScan]
    WorkflowConfig[WorkflowConfig]
    WorkflowConfigAudit[WorkflowConfigAudit]
    WorkflowParamAdapters[WorkflowParamAdapters]
    src_Constants_PNJL_jl[Constants_PNJL.jl]
    src_ParameterTypes_jl[ParameterTypes.jl]
    src_QuarkDistribution_jl[QuarkDistribution.jl]
    src_QuarkDistribution_Aniso_jl[QuarkDistribution_Aniso.jl]
  end
  subgraph scripts
    scripts_analysis_relaxtime_causal_gbu_infinite_qgate_jl[scripts/analysis/relaxtime/causal_gbu_infinite_qgate.jl]
    scripts_relaxtime_config_WorkflowConfig_jl[scripts/relaxtime/config/WorkflowConfig.jl]
    scripts_relaxtime_config_WorkflowConfigAudit_jl[scripts/relaxtime/config/WorkflowConfigAudit.jl]
    scripts_relaxtime_workflow_charged_gbu_plot_jl[scripts/relaxtime/workflow/charged_gbu_plot.jl]
    scripts_relaxtime_workflow_cross_section_orchestrated_jl[scripts/relaxtime/workflow/cross_section_orchestrated.jl]
  end
  subgraph simulation
    src_simulation_EllipsoidCalculation_jl[simulation/EllipsoidCalculation.jl]
    src_simulation_FrameTransformations_jl[simulation/FrameTransformations.jl]
    src_simulation_FullServerApp_jl[simulation/FullServerApp.jl]
    src_simulation_HTTPServer_jl[simulation/HTTPServer.jl]
    src_simulation_MomentumMapping_jl[simulation/MomentumMapping.jl]
    src_simulation_ServerLauncher_jl[simulation/ServerLauncher.jl]
    src_simulation_ServerWarmup_jl[simulation/ServerWarmup.jl]
    src_simulation_fullserver_compute_handlers_jl[simulation/fullserver/compute_handlers.jl]
    src_simulation_fullserver_http_utils_jl[simulation/fullserver/http_utils.jl]
    src_simulation_fullserver_pnjl_handlers_jl[simulation/fullserver/pnjl_handlers.jl]
    src_simulation_fullserver_pnjl_scan_jobs_jl[simulation/fullserver/pnjl_scan_jobs.jl]
    src_simulation_fullserver_routing_jl[simulation/fullserver/routing.jl]
    src_simulation_fullserver_script_task_jobs_jl[simulation/fullserver/script_task_jobs.jl]
    src_simulation_fullserver_shared_jl[simulation/fullserver/shared.jl]
    src_simulation_fullserver_transport_handlers_jl[simulation/fullserver/transport_handlers.jl]
  end
  subgraph types
    src_types_ParameterTypes_jl[types/ParameterTypes.jl]
  end
  subgraph utils
    src_utils_ParameterAdapters_jl[utils/ParameterAdapters.jl]
    src_utils_ParticleSymbols_jl[utils/ParticleSymbols.jl]
    src_utils_ValidationUtils_jl[utils/ValidationUtils.jl]
  end
  src_Constants_PNJL_jl --> ConfigLoader
  src_Constants_PNJL_jl --> TransportConstants
  src_Constants_PNJL_jl --> src_config_ConfigLoader_jl
  src_Constants_PNJL_jl --> src_constants_TransportConstants_jl
  src_QuarkDistribution_Aniso_jl --> src_QuarkDistribution_jl
  src_constants_Constants_PNJL_jl --> src_Constants_PNJL_jl
  src_integration_PhaseSpaceSampling_jl --> GaussLegendre
  src_models_Models_jl --> AbstractSusceptibilityProvider
  src_models_Models_jl --> Conditions
  src_models_Models_jl --> ConservedChargeSusceptibilities
  src_models_Models_jl --> FlavorChemicalProfiles
  src_models_Models_jl --> FreezeoutPathProfiles
  src_models_Models_jl --> FreezeoutPathScan
  src_models_Models_jl --> FreezeoutProfiles
  src_models_Models_jl --> HigherOrderDerivatives
  src_models_Models_jl --> IsentropicPathProfiles
  src_models_Models_jl --> MagneticIntegrals
  src_models_Models_jl --> MagneticScan
  src_models_Models_jl --> MagneticThermodynamics
  src_models_Models_jl --> MesonChemicalProfiles
  src_models_Models_jl --> PrecompileRegistry
  src_models_Models_jl --> SeedStrategies
  src_models_Models_jl --> ThermoDerivatives
  src_models_Models_jl --> TmuScan
  src_models_Models_jl --> TrhoScan
  src_models_Models_jl --> src_models_abstract_model_jl
  src_models_Models_jl --> src_models_derivatives_AbstractSusceptibilityProvider_jl
  src_models_Models_jl --> src_models_derivatives_ConservedChargeSusceptibilities_jl
  src_models_Models_jl --> src_models_derivatives_HigherOrderDerivatives_jl
  src_models_Models_jl --> src_models_derivatives_MixedTaylorJets_jl
  src_models_Models_jl --> src_models_derivatives_PNJLChiBTaylorDiff_jl
  src_models_Models_jl --> src_models_derivatives_TaylorDiffForwardDiffCompat_jl
  src_models_Models_jl --> src_models_derivatives_ThermoDerivatives_jl
  src_models_Models_jl --> src_models_entrypoints_jl
  src_models_Models_jl --> src_models_exports_public_jl
  src_models_Models_jl --> src_models_factory_jl
  src_models_Models_jl --> src_models_gas_liquid_adapters_entrypoint_adapter_jl
  src_models_Models_jl --> src_models_gas_liquid_api_jl
  src_models_Models_jl --> src_models_gas_liquid_capabilities_jl
  src_models_Models_jl --> src_models_njl_NJL2Model_jl
  src_models_Models_jl --> src_models_njl_NJLModel_jl
  src_models_Models_jl --> src_models_njl_adapters_entrypoint_adapter_jl
  src_models_Models_jl --> src_models_njl_api_jl
  src_models_Models_jl --> src_models_njl_capabilities_jl
  src_models_Models_jl --> src_models_njl2_adapters_entrypoint_adapter_jl
  src_models_Models_jl --> src_models_njl2_api_jl
  src_models_Models_jl --> src_models_njl2_capabilities_jl
  src_models_Models_jl --> src_models_omega_jl
  src_models_Models_jl --> src_models_phase_AdaptiveRhoRefinement_jl
  src_models_Models_jl --> src_models_phase_CEPDetector_jl
  src_models_Models_jl --> src_models_phase_CrossoverLine_jl
  src_models_Models_jl --> src_models_phase_PMPhaseArtifacts_jl
  src_models_Models_jl --> src_models_phase_PMPhaseDiagnostic_jl
  src_models_Models_jl --> src_models_phase_PMPhaseSeeds_jl
  src_models_Models_jl --> src_models_phase_PMPhaseTypes_jl
  src_models_Models_jl --> src_models_phase_PhaseArtifacts_jl
  src_models_Models_jl --> src_models_phase_PhaseCore_jl
  src_models_Models_jl --> src_models_phase_PhaseGridConvergence_jl
  src_models_Models_jl --> src_models_phase_PhaseIO_jl
  src_models_Models_jl --> src_models_phase_PhasePipeline_jl
  src_models_Models_jl --> src_models_phase_PhaseTypes_jl
  src_models_Models_jl --> src_models_phase_ProductionPhasePipeline_jl
  src_models_Models_jl --> src_models_phase_RhoSupportRefinement_jl
  src_models_Models_jl --> src_models_pnjl_adapters_entrypoint_adapter_jl
  src_models_Models_jl --> src_models_pnjl_api_jl
  src_models_Models_jl --> src_models_pnjl_capabilities_jl
  src_models_Models_jl --> src_models_pnjl_magnetic_adapters_entrypoint_adapter_jl
  src_models_Models_jl --> src_models_pnjl_magnetic_api_jl
  src_models_Models_jl --> src_models_pnjl_magnetic_capabilities_jl
  src_models_Models_jl --> src_models_pnjl_physics_PNJLCore_jl
  src_models_Models_jl --> src_models_pnjl_physics_PNJLDistributions_jl
  src_models_Models_jl --> src_models_pnjl_physics_PNJLMagneticModel_jl
  src_models_Models_jl --> src_models_pnjl_physics_PNJLModel_jl
  src_models_Models_jl --> src_models_pnjl_physics_core_MagneticIntegrals_jl
  src_models_Models_jl --> src_models_pnjl_physics_core_MagneticThermodynamics_jl
  src_models_Models_jl --> src_models_precompile_registry_jl
  src_models_Models_jl --> src_models_precompile_workload_jl
  src_models_Models_jl --> src_models_rotation_adapters_entrypoint_adapter_jl
  src_models_Models_jl --> src_models_rotation_api_jl
  src_models_Models_jl --> src_models_rotation_capabilities_jl
  src_models_Models_jl --> src_models_rpnjl_RPNJLModel_jl
  src_models_Models_jl --> src_models_rpnjl_adapters_entrypoint_adapter_jl
  src_models_Models_jl --> src_models_rpnjl_api_jl
  src_models_Models_jl --> src_models_rpnjl_capabilities_jl
  src_models_Models_jl --> src_models_scans_CrossoverMesonDensityScan_jl
  src_models_Models_jl --> src_models_scans_ExternalPathMesonDensityScan_jl
  src_models_Models_jl --> src_models_scans_FlavorChemicalProfiles_jl
  src_models_Models_jl --> src_models_scans_FreezeoutMesonDensityScan_jl
  src_models_Models_jl --> src_models_scans_FreezeoutPathProfiles_jl
  src_models_Models_jl --> src_models_scans_FreezeoutPathScan_jl
  src_models_Models_jl --> src_models_scans_FreezeoutProfiles_jl
  src_models_Models_jl --> src_models_scans_IsentropicPathProfiles_jl
  src_models_Models_jl --> src_models_scans_MagneticScan_jl
  src_models_Models_jl --> src_models_scans_MesonChemicalProfiles_jl
  src_models_Models_jl --> src_models_scans_MesonMassPathScan_jl
  src_models_Models_jl --> src_models_scans_ScanCommon_jl
  src_models_Models_jl --> src_models_scans_ScanConfig_jl
  src_models_Models_jl --> src_models_scans_ScanResultFinalize_jl
  src_models_Models_jl --> src_models_scans_TmuScan_jl
  src_models_Models_jl --> src_models_scans_TrhoScan_jl
  src_models_Models_jl --> src_models_solver_topology_jl
  src_models_Models_jl --> src_models_state_jl
  src_models_Models_jl --> src_models_thermo_kernel_jl
  src_models_Models_jl --> src_models_transport_provider_jl
  src_models_Models_jl --> src_models_variants_gas_liquid_GasLiquidModel_jl
  src_models_Models_jl --> src_models_variants_gas_liquid_workflows_GasLiquidWorkflow_jl
  src_models_Models_jl --> src_models_variants_rotation_RotationModel_jl
  src_models_Models_jl --> src_models_variants_rotation_workflows_RotationWorkflow_jl
  src_models_Models_jl --> src_models_workflow_apps_MesonDensityWorkflow_jl
  src_models_Models_jl --> src_models_workflow_apps_MesonMassWorkflow_jl
  src_models_Models_jl --> src_models_workflow_apps_MesonThermoWorkflow_jl
  src_models_Models_jl --> src_models_workflow_apps_TransportWorkflow_jl
  src_models_Models_jl --> src_models_workflow_apps_WorkflowParamAdapters_jl
  src_models_Models_jl --> src_models_workflow_engine_PipelineRunner_jl
  src_models_Models_jl --> src_models_workflow_engine_PipelineTypes_jl
  src_models_Models_jl --> src_models_workflow_engine_StageCatalog_jl
  src_models_Models_jl --> src_models_workflow_engine_adapters_CommonAdapterUtils_jl
  src_models_Models_jl --> src_models_workflow_engine_adapters_RelaxtimeOrchestratorAdapter_jl
  src_models_Models_jl --> src_models_workflow_engine_adapters_ScanAdapter_jl
  src_models_Models_jl --> src_models_workflow_engine_adapters_WorkflowAdapter_jl
  src_models_Models_jl --> src_models_workflow_engine_catalog_RelaxtimeOrchestratorCatalog_jl
  src_models_Models_jl --> src_models_workflow_engine_catalog_ScanCatalog_jl
  src_models_Models_jl --> src_models_workflow_engine_catalog_WorkflowCatalog_jl
  src_models_Models_jl --> src_models_workflow_engine_io_ManifestExtensions_jl
  src_models_derivatives_ConservedChargeSusceptibilities_jl --> Models
  src_models_derivatives_ConservedChargeSusceptibilities_jl --> PNJLChiBTaylorDiff
  src_models_derivatives_ConservedChargeSusceptibilities_jl --> PNJLCore
  src_models_derivatives_MixedTaylorJets_jl --> PNJLCore
  src_models_derivatives_PNJLChiBTaylorDiff_jl --> Conditions
  src_models_derivatives_PNJLChiBTaylorDiff_jl --> MixedTaylorJets
  src_models_derivatives_PNJLChiBTaylorDiff_jl --> Models
  src_models_derivatives_PNJLChiBTaylorDiff_jl --> TaylorDiffForwardDiffCompat
  src_models_derivatives_ThermoDerivatives_jl --> Models
  src_models_derivatives_ThermoDerivatives_jl --> PNJLChiBTaylorDiff
  src_models_derivatives_ThermoDerivatives_jl --> PNJLCore
  src_models_derivatives_ThermoDerivatives_jl --> src_constants_Constants_PNJL_jl
  src_models_derivatives_ThermoDerivatives_jl --> src_models_pnjl_physics_core_EquilibriumFacade_jl
  src_models_entrypoints_jl --> src_models_workflow_apps_ChargedGBUResearchWorkflow_jl
  src_models_njl_NJL2Core_jl --> ConfigLoader
  src_models_njl_NJL2Core_jl --> src_config_ConfigLoader_jl
  src_models_njl_NJL2Model_jl --> GaussLegendre
  src_models_njl_NJL2Model_jl --> src_integration_GaussLegendre_jl
  src_models_njl_NJL2Model_jl --> src_models_njl_NJL2Core_jl
  src_models_njl_NJLCore_jl --> ConfigLoader
  src_models_njl_NJLCore_jl --> src_config_ConfigLoader_jl
  src_models_njl_NJLModel_jl --> GaussLegendre
  src_models_njl_NJLModel_jl --> src_integration_GaussLegendre_jl
  src_models_njl_NJLModel_jl --> src_models_njl_NJLCore_jl
  src_models_phase_RhoSupportRefinement_jl --> Models
  src_models_pnjl_physics_PNJLCore_jl --> src_constants_Constants_PNJL_jl
  src_models_pnjl_physics_PNJLCore_jl --> src_models_pnjl_physics_PNJLIntegrals_jl
  src_models_pnjl_physics_PNJLIntegrals_jl --> src_integration_GaussLegendre_jl
  src_models_pnjl_physics_PNJLMagneticModel_jl --> Models
  src_models_pnjl_physics_PNJLMagneticModel_jl --> src_models_pnjl_physics_MagneticGapSolver_jl
  src_models_pnjl_physics_QuarkDistribution_jl --> src_QuarkDistribution_jl
  src_models_pnjl_physics_core_Integrals_jl --> src_constants_Constants_PNJL_jl
  src_models_pnjl_physics_core_Integrals_jl --> src_integration_GaussLegendre_jl
  src_models_pnjl_physics_core_MagneticIntegrals_jl --> src_constants_Constants_PNJL_jl
  src_models_pnjl_physics_core_MagneticIntegrals_jl --> src_integration_GaussLegendre_jl
  src_models_pnjl_physics_core_MagneticThermodynamics_jl --> MagneticIntegrals
  src_models_pnjl_physics_core_MagneticThermodynamics_jl --> src_constants_Constants_PNJL_jl
  src_models_pnjl_physics_core_MagneticThermodynamics_jl --> src_models_pnjl_physics_PNJLCore_jl
  src_models_pnjl_physics_core_MagneticThermodynamics_jl --> src_models_pnjl_physics_core_MagneticIntegrals_jl
  src_models_pnjl_physics_core_ModelThermodynamics_jl --> Models
  src_models_pnjl_physics_core_ModelThermodynamics_jl --> src_constants_Constants_PNJL_jl
  src_models_pnjl_physics_core_ModelThermodynamics_jl --> src_models_Models_jl
  src_models_precompile_registry_jl --> Models
  src_models_rpnjl_RPNJLModel_jl --> src_constants_Constants_PNJL_jl
  src_models_scans_CrossoverMesonDensityScan_jl --> FlavorChemicalProfiles
  src_models_scans_CrossoverMesonDensityScan_jl --> MesonChemicalProfiles
  src_models_scans_CrossoverMesonDensityScan_jl --> MesonDensityWorkflow
  src_models_scans_CrossoverMesonDensityScan_jl --> Models
  src_models_scans_CrossoverMesonDensityScan_jl --> ScanCommon
  src_models_scans_ExternalPathMesonDensityScan_jl --> CrossoverMesonDensityScan
  src_models_scans_ExternalPathMesonDensityScan_jl --> FlavorChemicalProfiles
  src_models_scans_ExternalPathMesonDensityScan_jl --> MesonChemicalProfiles
  src_models_scans_ExternalPathMesonDensityScan_jl --> ScanCommon
  src_models_scans_FlavorChemicalProfiles_jl --> ConfigLoader
  src_models_scans_FlavorChemicalProfiles_jl --> src_config_ConfigLoader_jl
  src_models_scans_FreezeoutMesonDensityScan_jl --> FlavorChemicalProfiles
  src_models_scans_FreezeoutMesonDensityScan_jl --> FreezeoutPathProfiles
  src_models_scans_FreezeoutMesonDensityScan_jl --> FreezeoutProfiles
  src_models_scans_FreezeoutMesonDensityScan_jl --> MesonChemicalProfiles
  src_models_scans_FreezeoutMesonDensityScan_jl --> MesonDensityWorkflow
  src_models_scans_FreezeoutMesonDensityScan_jl --> ScanCommon
  src_models_scans_FreezeoutPathProfiles_jl --> ConfigLoader
  src_models_scans_FreezeoutPathProfiles_jl --> FreezeoutProfiles
  src_models_scans_FreezeoutPathProfiles_jl --> src_config_ConfigLoader_jl
  src_models_scans_FreezeoutPathScan_jl --> FreezeoutPathProfiles
  src_models_scans_FreezeoutPathScan_jl --> FreezeoutProfiles
  src_models_scans_FreezeoutPathScan_jl --> Models
  src_models_scans_FreezeoutPathScan_jl --> ScanCommon
  src_models_scans_FreezeoutPathScan_jl --> ScanConfig
  src_models_scans_FreezeoutPathScan_jl --> SeedStrategies
  src_models_scans_FreezeoutPathScan_jl --> TmuScan
  src_models_scans_FreezeoutProfiles_jl --> ConfigLoader
  src_models_scans_FreezeoutProfiles_jl --> src_config_ConfigLoader_jl
  src_models_scans_MagneticScan_jl --> Models
  src_models_scans_MesonChemicalProfiles_jl --> ConfigLoader
  src_models_scans_MesonChemicalProfiles_jl --> src_config_ConfigLoader_jl
  src_models_scans_MesonMassPathScan_jl --> FreezeoutPathProfiles
  src_models_scans_MesonMassPathScan_jl --> FreezeoutProfiles
  src_models_scans_MesonMassPathScan_jl --> IsentropicPathProfiles
  src_models_scans_MesonMassPathScan_jl --> MesonMassWorkflow
  src_models_scans_MesonMassPathScan_jl --> Models
  src_models_scans_MesonMassPathScan_jl --> ScanCommon
  src_models_scans_ScanCommon_jl --> Models
  src_models_scans_ScanCommon_jl --> SeedStrategies
  src_models_scans_ScanResultFinalize_jl --> Models
  src_models_scans_TmuScan_jl --> Models
  src_models_scans_TmuScan_jl --> ScanCommon
  src_models_scans_TmuScan_jl --> ScanConfig
  src_models_scans_TmuScan_jl --> ScanResultFinalize
  src_models_scans_TmuScan_jl --> SeedStrategies
  src_models_scans_TrhoScan_jl --> Models
  src_models_scans_TrhoScan_jl --> ScanCommon
  src_models_scans_TrhoScan_jl --> ScanConfig
  src_models_scans_TrhoScan_jl --> ScanResultFinalize
  src_models_scans_TrhoScan_jl --> SeedStrategies
  src_models_solver_orchestrator_SeedStrategies_jl --> Models
  src_models_solver_runtime_ConstraintSolver_jl --> src_models_solver_diagnostics_ThermoPostprocess_jl
  src_models_solver_runtime_ConstraintSolver_jl --> src_models_solver_runtime_ConstraintSolverCommon_jl
  src_models_solver_runtime_ConstraintSolver_jl --> src_models_solver_runtime_ConstraintSolverFixedAsymmetricRho_jl
  src_models_solver_runtime_ConstraintSolver_jl --> src_models_solver_runtime_ConstraintSolverFixedEntropy_jl
  src_models_solver_runtime_ConstraintSolver_jl --> src_models_solver_runtime_ConstraintSolverFixedMu_jl
  src_models_solver_runtime_ConstraintSolver_jl --> src_models_solver_runtime_ConstraintSolverFixedMuBConservedCharges_jl
  src_models_solver_runtime_ConstraintSolver_jl --> src_models_solver_runtime_ConstraintSolverFixedRho_jl
  src_models_solver_runtime_ConstraintSolver_jl --> src_models_solver_runtime_ConstraintSolverFixedSigma_jl
  src_models_solver_spec_Conditions_jl --> Models
  src_models_solver_topology_jl --> src_models_derivatives_DiffService_jl
  src_models_solver_topology_jl --> src_models_solver_ImplicitProblem_jl
  src_models_solver_topology_jl --> src_models_solver_api_SolverAPI_jl
  src_models_solver_topology_jl --> src_models_solver_compat_ImplicitAdapters_jl
  src_models_solver_topology_jl --> src_models_solver_compat_ImplicitGapLegacy_jl
  src_models_solver_topology_jl --> src_models_solver_compat_SchemaAdapter_jl
  src_models_solver_topology_jl --> src_models_solver_config_SolverRuntimeConfig_jl
  src_models_solver_topology_jl --> src_models_solver_config_StateSchema_jl
  src_models_solver_topology_jl --> src_models_solver_diagnostics_SolverDiagnostics_jl
  src_models_solver_topology_jl --> src_models_solver_diagnostics_SolverDiagnosticsTypes_jl
  src_models_solver_topology_jl --> src_models_solver_diagnostics_SolverWorkTelemetry_jl
  src_models_solver_topology_jl --> src_models_solver_diff_JacobianEngine_jl
  src_models_solver_topology_jl --> src_models_solver_diff_PilotAdapters_jl
  src_models_solver_topology_jl --> src_models_solver_diff_Targets_jl
  src_models_solver_topology_jl --> src_models_solver_diff_ThermoDiffContext_jl
  src_models_solver_topology_jl --> src_models_solver_governance_CandidateGovernance_jl
  src_models_solver_topology_jl --> src_models_solver_governance_WeightedFallback_jl
  src_models_solver_topology_jl --> src_models_solver_orchestrator_PrimaryStrategy_jl
  src_models_solver_topology_jl --> src_models_solver_orchestrator_ProblemSpecOrchestrator_jl
  src_models_solver_topology_jl --> src_models_solver_orchestrator_SeedStrategies_jl
  src_models_solver_topology_jl --> src_models_solver_path_PathContinuation_jl
  src_models_solver_topology_jl --> src_models_solver_runtime_ConstraintSolver_jl
  src_models_solver_topology_jl --> src_models_solver_runtime_GapSolver_jl
  src_models_solver_topology_jl --> src_models_solver_runtime_GenericRootEngine_jl
  src_models_solver_topology_jl --> src_models_solver_spec_Conditions_jl
  src_models_solver_topology_jl --> src_models_solver_spec_ConstraintComponents_jl
  src_models_solver_topology_jl --> src_models_solver_spec_ConstraintModes_jl
  src_models_solver_topology_jl --> src_models_solver_spec_ProblemSpec_jl
  src_models_variants_gas_liquid_GasLiquidModel_jl --> GasLiquidEquationSet
  src_models_variants_gas_liquid_GasLiquidModel_jl --> GasLiquidThermodynamics
  src_models_variants_gas_liquid_GasLiquidModel_jl --> src_models_variants_gas_liquid_core_EquationSet_jl
  src_models_variants_gas_liquid_GasLiquidModel_jl --> src_models_variants_gas_liquid_core_Thermodynamics_jl
  src_models_variants_gas_liquid_core_EquationSet_jl --> ConfigLoader
  src_models_variants_gas_liquid_core_EquationSet_jl --> GaussLegendre
  src_models_variants_gas_liquid_core_EquationSet_jl --> src_config_ConfigLoader_jl
  src_models_variants_gas_liquid_core_EquationSet_jl --> src_integration_GaussLegendre_jl
  src_models_variants_gas_liquid_core_Thermodynamics_jl --> GasLiquidEquationSet
  src_models_variants_gas_liquid_workflows_GasLiquidWorkflow_jl --> GasLiquidEquationSet
  src_models_variants_gas_liquid_workflows_GasLiquidWorkflow_jl --> GasLiquidThermodynamics
  src_models_variants_gas_liquid_workflows_GasLiquidWorkflow_jl --> Models
  src_models_variants_rotation_RotationModel_jl --> RotationThermo
  src_models_variants_rotation_RotationModel_jl --> src_models_variants_rotation_core_RotationThermo_jl
  src_models_variants_rotation_core_RotationThermo_jl --> ConfigLoader
  src_models_variants_rotation_core_RotationThermo_jl --> src_config_ConfigLoader_jl
  src_models_variants_rotation_workflows_RotationWorkflow_jl --> Models
  src_models_variants_rotation_workflows_RotationWorkflow_jl --> RotationThermo
  src_models_workflow_apps_ChargedGBUResearchWorkflow_jl --> scripts_analysis_relaxtime_causal_gbu_infinite_qgate_jl
  src_models_workflow_apps_ChargedGBUResearchWorkflow_jl --> scripts_relaxtime_workflow_charged_gbu_plot_jl
  src_models_workflow_apps_MesonDensityWorkflow_jl --> MesonMassWorkflow
  src_models_workflow_apps_MesonDensityWorkflow_jl --> WorkflowParamAdapters
  src_models_workflow_apps_MesonDensityWorkflow_jl --> src_relaxtime_RelaxTime_jl
  src_models_workflow_apps_MesonMassWorkflow_jl --> Models
  src_models_workflow_apps_MesonMassWorkflow_jl --> WorkflowParamAdapters
  src_models_workflow_apps_MesonMassWorkflow_jl --> src_models_pnjl_physics_core_EquilibriumFacade_jl
  src_models_workflow_apps_MesonMassWorkflow_jl --> src_models_workflow_apps_WorkflowParamAdapters_jl
  src_models_workflow_apps_MesonMassWorkflow_jl --> src_relaxtime_RelaxTime_jl
  src_models_workflow_apps_MesonMassWorkflow_jl --> src_types_ParameterTypes_jl
  src_models_workflow_apps_MesonThermoWorkflow_jl --> MesonMassWorkflow
  src_models_workflow_apps_MesonThermoWorkflow_jl --> Models
  src_models_workflow_apps_MesonThermoWorkflow_jl --> WorkflowParamAdapters
  src_models_workflow_apps_MesonThermoWorkflow_jl --> src_relaxtime_RelaxTime_jl
  src_models_workflow_apps_TransportWorkflow_jl --> Models
  src_models_workflow_apps_TransportWorkflow_jl --> TransportCoefficients
  src_models_workflow_apps_TransportWorkflow_jl --> WorkflowParamAdapters
  src_models_workflow_apps_TransportWorkflow_jl --> src_Constants_PNJL_jl
  src_models_workflow_apps_TransportWorkflow_jl --> src_models_pnjl_physics_core_EquilibriumFacade_jl
  src_models_workflow_apps_TransportWorkflow_jl --> src_models_workflow_apps_WorkflowParamAdapters_jl
  src_models_workflow_apps_TransportWorkflow_jl --> src_relaxtime_RelaxTime_jl
  src_models_workflow_apps_TransportWorkflow_jl --> src_types_ParameterTypes_jl
  src_models_workflow_apps_WorkflowParamAdapters_jl --> src_types_ParameterTypes_jl
  src_models_workflow_engine_adapters_RelaxtimeOrchestratorAdapter_jl --> CrossSectionOrchestratedScan
  src_models_workflow_engine_adapters_RelaxtimeOrchestratorAdapter_jl --> WorkflowConfig
  src_models_workflow_engine_adapters_RelaxtimeOrchestratorAdapter_jl --> WorkflowConfigAudit
  src_models_workflow_engine_adapters_RelaxtimeOrchestratorAdapter_jl --> scripts_relaxtime_config_WorkflowConfig_jl
  src_models_workflow_engine_adapters_RelaxtimeOrchestratorAdapter_jl --> scripts_relaxtime_config_WorkflowConfigAudit_jl
  src_models_workflow_engine_adapters_RelaxtimeOrchestratorAdapter_jl --> scripts_relaxtime_workflow_cross_section_orchestrated_jl
  src_relaxtime_AFieldBuilder_jl --> GaussLegendre
  src_relaxtime_AFieldBuilder_jl --> OneLoopIntegrals
  src_relaxtime_AFieldBuilder_jl --> OneLoopIntegralsCorrection
  src_relaxtime_AverageScatteringRate_jl --> AFieldBuilder
  src_relaxtime_AverageScatteringRate_jl --> GaussLegendre
  src_relaxtime_AverageScatteringRate_jl --> ParticleSymbols
  src_relaxtime_AverageScatteringRate_jl --> TotalCrossSection
  src_relaxtime_AverageScatteringRate_jl --> src_utils_ValidationUtils_jl
  src_relaxtime_CausalSpectralBubble_jl --> GaussLegendre
  src_relaxtime_CausalSpectralBubble_jl --> OneLoopIntegrals
  src_relaxtime_ChargedPhaseBackend_jl --> BUPhaseGates
  src_relaxtime_ChargedPhaseBackend_jl --> ChargedRPAKernel
  src_relaxtime_ChargedPhaseBackend_jl --> GaussLegendre
  src_relaxtime_ChargedPhaseBackend_jl --> PhaseNormalization
  src_relaxtime_ChargedRPAKernel_jl --> MesonInteractionKernel
  src_relaxtime_ChargedRPAProvider_jl --> ChargedRPAKernel
  src_relaxtime_ChargedRPAProvider_jl --> OneLoopIntegrals
  src_relaxtime_ChargedRPAProvider_jl --> PolarizationAniso
  src_relaxtime_DifferentialCrossSection_jl --> KinematicChecks
  src_relaxtime_DifferentialCrossSection_jl --> src_relaxtime_KinematicChecks_jl
  src_relaxtime_EffectiveCouplings_jl --> OneLoopIntegrals
  src_relaxtime_EffectiveCouplings_jl --> OneLoopIntegralsCorrection
  src_relaxtime_MesonDensity_jl --> AFieldBuilder
  src_relaxtime_MesonDensity_jl --> BUPhaseGates
  src_relaxtime_MesonDensity_jl --> EffectiveCouplings
  src_relaxtime_MesonDensity_jl --> GaussLegendre
  src_relaxtime_MesonDensity_jl --> MesonInteractionKernel
  src_relaxtime_MesonDensity_jl --> MesonMass
  src_relaxtime_MesonDensity_jl --> MesonPropagator
  src_relaxtime_MesonDensity_jl --> PolarizationAniso
  src_relaxtime_MesonMass_jl --> AFieldBuilder
  src_relaxtime_MesonMass_jl --> EffectiveCouplings
  src_relaxtime_MesonMass_jl --> GaussLegendre
  src_relaxtime_MesonMass_jl --> PolarizationAniso
  src_relaxtime_MesonPropagator_jl --> EffectiveCouplings
  src_relaxtime_MesonPropagator_jl --> ParticleSymbols
  src_relaxtime_MesonRPA_jl --> MesonInteractionKernel
  src_relaxtime_MesonRPAAdapter_jl --> AFieldBuilder
  src_relaxtime_MesonRPAAdapter_jl --> GaussLegendre
  src_relaxtime_MesonRPAAdapter_jl --> MesonInteractionKernel
  src_relaxtime_MesonRPAAdapter_jl --> MesonRPA
  src_relaxtime_MesonRPAAdapter_jl --> PolarizationAniso
  src_relaxtime_MesonThermodynamics_jl --> AFieldBuilder
  src_relaxtime_MesonThermodynamics_jl --> GaussLegendre
  src_relaxtime_MesonThermodynamics_jl --> MesonDensity
  src_relaxtime_OneLoopIntegrals_jl --> GaussLegendre
  src_relaxtime_OneLoopIntegrals_jl --> src_integration_IntervalQuadratureStrategies_jl
  src_relaxtime_OneLoopIntegralsAniso_jl --> GaussLegendre
  src_relaxtime_OneLoopIntegralsAniso_jl --> OneLoopIntegrals
  src_relaxtime_OneLoopIntegralsAniso_jl --> src_integration_IntervalQuadratureStrategies_jl
  src_relaxtime_PolarizationAniso_jl --> OneLoopIntegrals
  src_relaxtime_PolarizationAniso_jl --> OneLoopIntegralsCorrection
  src_relaxtime_PolarizationCache_jl --> PolarizationAniso
  src_relaxtime_RelaxTime_jl --> AFieldBuilder
  src_relaxtime_RelaxTime_jl --> AverageScatteringRate
  src_relaxtime_RelaxTime_jl --> BUPhaseGates
  src_relaxtime_RelaxTime_jl --> ChargedPhaseBackend
  src_relaxtime_RelaxTime_jl --> ChargedRPAKernel
  src_relaxtime_RelaxTime_jl --> ChargedRPAProvider
  src_relaxtime_RelaxTime_jl --> DifferentialCrossSection
  src_relaxtime_RelaxTime_jl --> EffectiveCouplings
  src_relaxtime_RelaxTime_jl --> KinematicChecks
  src_relaxtime_RelaxTime_jl --> MesonDensity
  src_relaxtime_RelaxTime_jl --> MesonInteractionKernel
  src_relaxtime_RelaxTime_jl --> MesonMass
  src_relaxtime_RelaxTime_jl --> MesonPropagator
  src_relaxtime_RelaxTime_jl --> MesonRPA
  src_relaxtime_RelaxTime_jl --> MesonRPAAdapter
  src_relaxtime_RelaxTime_jl --> MesonThermodynamics
  src_relaxtime_RelaxTime_jl --> MottTransition
  src_relaxtime_RelaxTime_jl --> OneLoopIntegrals
  src_relaxtime_RelaxTime_jl --> OneLoopIntegralsCorrection
  src_relaxtime_RelaxTime_jl --> PhaseNormalization
  src_relaxtime_RelaxTime_jl --> PolarizationAniso
  src_relaxtime_RelaxTime_jl --> PolarizationCache
  src_relaxtime_RelaxTime_jl --> RelaxationTime
  src_relaxtime_RelaxTime_jl --> ScatteringAmplitude
  src_relaxtime_RelaxTime_jl --> TotalCrossSection
  src_relaxtime_RelaxTime_jl --> TotalPropagator
  src_relaxtime_RelaxTime_jl --> TransportCoefficients
  src_relaxtime_RelaxTime_jl --> src_Constants_PNJL_jl
  src_relaxtime_RelaxTime_jl --> src_ParameterTypes_jl
  src_relaxtime_RelaxTime_jl --> src_QuarkDistribution_jl
  src_relaxtime_RelaxTime_jl --> src_QuarkDistribution_Aniso_jl
  src_relaxtime_RelaxTime_jl --> src_integration_GaussLegendre_jl
  src_relaxtime_RelaxTime_jl --> src_integration_PhaseSpaceSampling_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_AFieldBuilder_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_AverageScatteringRate_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_BUPhaseGates_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_CausalSpectralBubble_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_ChargedPhaseBackend_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_ChargedRPAKernel_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_ChargedRPAProvider_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_DifferentialCrossSection_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_EffectiveCouplings_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_KinematicChecks_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_MesonDensity_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_MesonInteractionKernel_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_MesonMass_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_MesonPropagator_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_MesonRPA_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_MesonRPAAdapter_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_MesonThermodynamics_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_MottTransition_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_OneLoopIntegrals_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_OneLoopIntegralsAniso_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_PhaseNormalization_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_PolarizationAniso_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_PolarizationCache_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_RelaxationTime_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_ScatteringAmplitude_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_TotalCrossSection_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_TotalPropagator_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_TransportCoefficients_jl
  src_relaxtime_RelaxTime_jl --> src_relaxtime_TransportCoefficientsValidation_jl
  src_relaxtime_RelaxTime_jl --> src_utils_ParameterAdapters_jl
  src_relaxtime_RelaxTime_jl --> src_utils_ParticleSymbols_jl
  src_relaxtime_RelaxTime_jl --> src_utils_ValidationUtils_jl
  src_relaxtime_RelaxationTime_jl --> AFieldBuilder
  src_relaxtime_RelaxationTime_jl --> AverageScatteringRate
  src_relaxtime_RelaxationTime_jl --> TotalCrossSection
  src_relaxtime_RelaxationTime_jl --> src_utils_ValidationUtils_jl
  src_relaxtime_ScatteringAmplitude_jl --> ParticleSymbols
  src_relaxtime_ScatteringAmplitude_jl --> TotalPropagator
  src_relaxtime_TotalCrossSection_jl --> DifferentialCrossSection
  src_relaxtime_TotalCrossSection_jl --> GaussLegendre
  src_relaxtime_TotalCrossSection_jl --> KinematicChecks
  src_relaxtime_TotalCrossSection_jl --> OneLoopIntegrals
  src_relaxtime_TotalCrossSection_jl --> ParticleSymbols
  src_relaxtime_TotalCrossSection_jl --> ScatteringAmplitude
  src_relaxtime_TotalPropagator_jl --> KinematicChecks
  src_relaxtime_TotalPropagator_jl --> MesonPropagator
  src_relaxtime_TotalPropagator_jl --> ParticleSymbols
  src_relaxtime_TotalPropagator_jl --> PolarizationCache
  src_relaxtime_TransportCoefficients_jl --> GaussLegendre
  src_relaxtime_TransportCoefficients_jl --> PhaseSpaceSampling
  src_relaxtime_TransportCoefficients_jl --> TransportCoefficientsValidation
  src_simulation_FullServerApp_jl --> MomentumMapping
  src_simulation_FullServerApp_jl --> src_constants_Constants_PNJL_jl
  src_simulation_FullServerApp_jl --> src_models_Models_jl
  src_simulation_FullServerApp_jl --> src_simulation_MomentumMapping_jl
  src_simulation_FullServerApp_jl --> src_simulation_fullserver_compute_handlers_jl
  src_simulation_FullServerApp_jl --> src_simulation_fullserver_http_utils_jl
  src_simulation_FullServerApp_jl --> src_simulation_fullserver_pnjl_handlers_jl
  src_simulation_FullServerApp_jl --> src_simulation_fullserver_pnjl_scan_jobs_jl
  src_simulation_FullServerApp_jl --> src_simulation_fullserver_routing_jl
  src_simulation_FullServerApp_jl --> src_simulation_fullserver_script_task_jobs_jl
  src_simulation_FullServerApp_jl --> src_simulation_fullserver_shared_jl
  src_simulation_FullServerApp_jl --> src_simulation_fullserver_transport_handlers_jl
  src_simulation_HTTPServer_jl --> MomentumMapping
  src_simulation_HTTPServer_jl --> src_simulation_MomentumMapping_jl
  src_simulation_MomentumMapping_jl --> EllipsoidCalculation
  src_simulation_MomentumMapping_jl --> FrameTransformations
  src_simulation_MomentumMapping_jl --> src_simulation_EllipsoidCalculation_jl
  src_simulation_MomentumMapping_jl --> src_simulation_FrameTransformations_jl
  src_simulation_ServerLauncher_jl --> FullServerApp
  src_simulation_ServerLauncher_jl --> ServerWarmup
  src_simulation_ServerLauncher_jl --> src_simulation_FullServerApp_jl
  src_simulation_ServerLauncher_jl --> src_simulation_ServerWarmup_jl
  src_simulation_ServerWarmup_jl --> src_constants_Constants_PNJL_jl
  src_simulation_ServerWarmup_jl --> src_models_Models_jl
  src_types_ParameterTypes_jl --> src_ParameterTypes_jl
  src_utils_ParameterAdapters_jl --> src_types_ParameterTypes_jl
  src_utils_ParticleSymbols_jl --> src_constants_Constants_PNJL_jl
  src_utils_ParticleSymbols_jl --> src_types_ParameterTypes_jl
  src_utils_ParticleSymbols_jl --> src_utils_ParameterAdapters_jl
```

