# 依赖审计报告

生成时间：2026-10-04T10:55:47.535

来源：当前 src/ 源码；直接调用 gen_deps 的静态解析器，不读取旧 dependencies.mmd。

只审计可定位的 include 文件边；相对导入名、宏展开、运行时参数和被加载脚本内部的依赖不在完整覆盖范围。
解析 include 边：256；跨组边：55；已知脚本桥接：5；违规：0；未解析或缺失：1。

## 跨组依赖

- src/models/derivatives/ThermoDerivatives.jl → src/constants/Constants_PNJL.jl
- src/models/njl/NJL2Core.jl → src/config/ConfigLoader.jl
- src/models/njl/NJL2Model.jl → src/integration/GaussLegendre.jl
- src/models/njl/NJLCore.jl → src/config/ConfigLoader.jl
- src/models/njl/NJLModel.jl → src/integration/GaussLegendre.jl
- src/models/pnjl_physics/PNJLCore.jl → src/constants/Constants_PNJL.jl
- src/models/pnjl_physics/PNJLIntegrals.jl → src/integration/GaussLegendre.jl
- src/models/pnjl_physics/QuarkDistribution.jl → src/QuarkDistribution.jl
- src/models/pnjl_physics/core/Integrals.jl → src/constants/Constants_PNJL.jl
- src/models/pnjl_physics/core/Integrals.jl → src/integration/GaussLegendre.jl
- src/models/pnjl_physics/core/MagneticIntegrals.jl → src/constants/Constants_PNJL.jl
- src/models/pnjl_physics/core/MagneticIntegrals.jl → src/integration/GaussLegendre.jl
- src/models/pnjl_physics/core/MagneticThermodynamics.jl → src/constants/Constants_PNJL.jl
- src/models/pnjl_physics/core/ModelThermodynamics.jl → src/constants/Constants_PNJL.jl
- src/models/rpnjl/RPNJLModel.jl → src/constants/Constants_PNJL.jl
- src/models/scans/FlavorChemicalProfiles.jl → src/config/ConfigLoader.jl
- src/models/scans/FreezeoutPathProfiles.jl → src/config/ConfigLoader.jl
- src/models/scans/FreezeoutProfiles.jl → src/config/ConfigLoader.jl
- src/models/scans/MesonChemicalProfiles.jl → src/config/ConfigLoader.jl
- src/models/variants/gas_liquid/core/EquationSet.jl → src/config/ConfigLoader.jl
- src/models/variants/gas_liquid/core/EquationSet.jl → src/integration/GaussLegendre.jl
- src/models/variants/rotation/core/RotationThermo.jl → src/config/ConfigLoader.jl
- src/models/workflow_apps/ChargedGBUResearchWorkflow.jl → scripts/analysis/relaxtime/causal_gbu_infinite_qgate.jl
- src/models/workflow_apps/ChargedGBUResearchWorkflow.jl → scripts/relaxtime/workflow/charged_gbu_plot.jl
- src/models/workflow_apps/MesonDensityWorkflow.jl → src/relaxtime/RelaxTime.jl
- src/models/workflow_apps/MesonMassWorkflow.jl → src/relaxtime/RelaxTime.jl
- src/models/workflow_apps/MesonMassWorkflow.jl → src/types/ParameterTypes.jl
- src/models/workflow_apps/MesonThermoWorkflow.jl → src/relaxtime/RelaxTime.jl
- src/models/workflow_apps/TransportWorkflow.jl → src/Constants_PNJL.jl
- src/models/workflow_apps/TransportWorkflow.jl → src/relaxtime/RelaxTime.jl
- src/models/workflow_apps/TransportWorkflow.jl → src/types/ParameterTypes.jl
- src/models/workflow_apps/WorkflowParamAdapters.jl → src/types/ParameterTypes.jl
- src/models/workflow_engine/adapters/RelaxtimeOrchestratorAdapter.jl → scripts/relaxtime/config/WorkflowConfig.jl
- src/models/workflow_engine/adapters/RelaxtimeOrchestratorAdapter.jl → scripts/relaxtime/config/WorkflowConfigAudit.jl
- src/models/workflow_engine/adapters/RelaxtimeOrchestratorAdapter.jl → scripts/relaxtime/workflow/cross_section_orchestrated.jl
- src/relaxtime/AverageScatteringRate.jl → src/utils/ValidationUtils.jl
- src/relaxtime/OneLoopIntegrals.jl → src/integration/IntervalQuadratureStrategies.jl
- src/relaxtime/OneLoopIntegralsAniso.jl → src/integration/IntervalQuadratureStrategies.jl
- src/relaxtime/RelaxTime.jl → src/Constants_PNJL.jl
- src/relaxtime/RelaxTime.jl → src/ParameterTypes.jl
- src/relaxtime/RelaxTime.jl → src/QuarkDistribution.jl
- src/relaxtime/RelaxTime.jl → src/QuarkDistribution_Aniso.jl
- src/relaxtime/RelaxTime.jl → src/integration/GaussLegendre.jl
- src/relaxtime/RelaxTime.jl → src/integration/PhaseSpaceSampling.jl
- src/relaxtime/RelaxTime.jl → src/utils/ParameterAdapters.jl
- src/relaxtime/RelaxTime.jl → src/utils/ParticleSymbols.jl
- src/relaxtime/RelaxTime.jl → src/utils/ValidationUtils.jl
- src/relaxtime/RelaxationTime.jl → src/utils/ValidationUtils.jl
- src/simulation/FullServerApp.jl → src/constants/Constants_PNJL.jl
- src/simulation/FullServerApp.jl → src/models/Models.jl
- src/simulation/ServerWarmup.jl → src/constants/Constants_PNJL.jl
- src/simulation/ServerWarmup.jl → src/models/Models.jl
- src/utils/ParameterAdapters.jl → src/types/ParameterTypes.jl
- src/utils/ParticleSymbols.jl → src/constants/Constants_PNJL.jl
- src/utils/ParticleSymbols.jl → src/types/ParameterTypes.jl

## 已知脚本桥接（后续独立迁移）

- src/models/workflow_apps/ChargedGBUResearchWorkflow.jl → scripts/analysis/relaxtime/causal_gbu_infinite_qgate.jl：charged GBU 计算内核迁移待独立验证
- src/models/workflow_apps/ChargedGBUResearchWorkflow.jl → scripts/relaxtime/workflow/charged_gbu_plot.jl：现有 charged GBU 绘图桥接
- src/models/workflow_engine/adapters/RelaxtimeOrchestratorAdapter.jl → scripts/relaxtime/config/WorkflowConfig.jl：现有 orchestrator 配置适配
- src/models/workflow_engine/adapters/RelaxtimeOrchestratorAdapter.jl → scripts/relaxtime/config/WorkflowConfigAudit.jl：现有 orchestrator 配置审计适配
- src/models/workflow_engine/adapters/RelaxtimeOrchestratorAdapter.jl → scripts/relaxtime/workflow/cross_section_orchestrated.jl：现有 orchestrator 工作流适配

## 未解析或缺失的 include

- src/models/workflow_engine/adapters/WorkflowAdapter.jl: Base.include(Main, lib_path)

## 已解析 include 的矩阵违规

- 未发现新增违规；这不是完整运行时依赖无环的证明。

规则与边界见 [依赖规则](dependency_rules.md)。
