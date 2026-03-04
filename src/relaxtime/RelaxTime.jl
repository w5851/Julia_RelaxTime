"""RelaxTime

`src/relaxtime` 的原生模块入口。

外部依赖（Constants_PNJL / ParameterTypes / PNJLQuarkDistributions / PNJLQuarkDistributions_Aniso）
在 `module RelaxTime` 之前被加载到 Main 作用域，保持与已有 tests/scripts 中
`using Main.Constants_PNJL` 等引用的兼容。
"""

# ── 加载共享外部依赖到 Main 作用域 ──
Base.include(Main, normpath(joinpath(@__DIR__, "..", "ParameterTypes.jl")))
Base.include(Main, normpath(joinpath(@__DIR__, "..", "Constants_PNJL.jl")))
Base.include(Main, normpath(joinpath(@__DIR__, "..", "QuarkDistribution.jl")))
Base.include(Main, normpath(joinpath(@__DIR__, "..", "QuarkDistribution_Aniso.jl")))

module RelaxTime

# ── Integration & utility submodules ──
include(joinpath(@__DIR__, "..", "integration", "GaussLegendre.jl"))
include(joinpath(@__DIR__, "..", "utils", "ParticleSymbols.jl"))
include(joinpath(@__DIR__, "..", "integration", "PhaseSpaceSampling.jl"))

# ── Relaxtime submodules (dependency order) ──
include("OneLoopIntegrals.jl")
include("OneLoopIntegralsAniso.jl")
include("PolarizationAniso.jl")
include("PolarizationCache.jl")
include("AFieldBuilder.jl")
include("EffectiveCouplings.jl")
include("MesonPropagator.jl")
include("TotalPropagator.jl")
include("MesonMass.jl")
include("MottTransition.jl")
include("DifferentialCrossSection.jl")
include("ScatteringAmplitude.jl")
include("TotalCrossSection.jl")
include("AverageScatteringRate.jl")
include("RelaxationTime.jl")
include("TransportCoefficients.jl")

# ── Re-exports ──
using .OneLoopIntegrals
using .OneLoopIntegralsCorrection
using .PolarizationAniso
using .PolarizationCache
using .AFieldBuilder
using .EffectiveCouplings
using .MesonPropagator
using .TotalPropagator
using .MesonMass
using .MottTransition
using .DifferentialCrossSection
using .ScatteringAmplitude
using .TotalCrossSection
using .AverageScatteringRate
using .RelaxationTime
using .TransportCoefficients

end # module RelaxTime

# ── 向后兼容：将子模块提升到 Main 作用域 ──
# 使现有代码中 `using Main.OneLoopIntegrals` 等引用继续工作。
for _name in (:OneLoopIntegrals, :OneLoopIntegralsCorrection, :PolarizationAniso,
              :PolarizationCache, :AFieldBuilder, :EffectiveCouplings,
              :MesonPropagator, :TotalPropagator, :MesonMass, :MottTransition,
              :DifferentialCrossSection, :ScatteringAmplitude, :TotalCrossSection,
              :AverageScatteringRate, :RelaxationTime, :TransportCoefficients,
              :GaussLegendre, :ParticleSymbols)
    if !isdefined(@__MODULE__, _name)
        Core.eval(@__MODULE__, :(const $_name = RelaxTime.$_name))
    end
end
