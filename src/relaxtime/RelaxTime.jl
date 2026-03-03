"""RelaxTime

`src/relaxtime` 的原生模块化入口（迁移脚手架）。

当前阶段策略：
- 保持现有 Main + IncludeOnce 路径不变，避免打断既有 tests/scripts；
- 先提供统一入口模块，后续按批次将子模块接入并替换调用方；
- 迁移完成后，再移除 legacy include-once 机制。
"""
module RelaxTime

include(joinpath(@__DIR__, "OneLoopIntegrals.jl"))
include(joinpath(@__DIR__, "OneLoopIntegralsAniso.jl"))
include(joinpath(@__DIR__, "PolarizationAniso.jl"))
include(joinpath(@__DIR__, "PolarizationCache.jl"))
include(joinpath(@__DIR__, "AFieldBuilder.jl"))
include(joinpath(@__DIR__, "EffectiveCouplings.jl"))
include(joinpath(@__DIR__, "MesonPropagator.jl"))
include(joinpath(@__DIR__, "MesonMass.jl"))
include(joinpath(@__DIR__, "MottTransition.jl"))
include(joinpath(@__DIR__, "DifferentialCrossSection.jl"))
include(joinpath(@__DIR__, "ScatteringAmplitude.jl"))
include(joinpath(@__DIR__, "TotalCrossSection.jl"))
include(joinpath(@__DIR__, "AverageScatteringRate.jl"))
if isdefined(Main, :RelaxationTime)
	const RelaxationTime = Main.RelaxationTime
else
	include(joinpath(@__DIR__, "RelaxationTime.jl"))
end
include(joinpath(@__DIR__, "TransportCoefficients.jl"))

using .OneLoopIntegrals
using .OneLoopIntegralsCorrection
using .PolarizationAniso
using .PolarizationCache
using .AFieldBuilder
using .EffectiveCouplings
using .MesonPropagator
using .MesonMass
using .MottTransition
using .DifferentialCrossSection
using .ScatteringAmplitude
using .TotalCrossSection
using .AverageScatteringRate
using .TransportCoefficients

end # module RelaxTime
