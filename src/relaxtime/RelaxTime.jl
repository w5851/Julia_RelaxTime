"""RelaxTime

`src/relaxtime` 的原生模块入口。

外部依赖（Constants_PNJL / ParameterTypes / PNJLQuarkDistributions / PNJLQuarkDistributions_Aniso）
在 `module RelaxTime` 之前被加载到 Main 作用域，保持与已有 tests/scripts 中
`using Main.Constants_PNJL` 等引用的兼容。
"""

# ── 加载共享外部依赖到 Main 作用域 ──
if !isdefined(Main, :ParameterTypes)
    Base.include(Main, normpath(joinpath(@__DIR__, "..", "ParameterTypes.jl")))
end
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, normpath(joinpath(@__DIR__, "..", "Constants_PNJL.jl")))
end
if !isdefined(Main, :PNJLQuarkDistributions)
    Base.include(Main, normpath(joinpath(@__DIR__, "..", "QuarkDistribution.jl")))
end
if !isdefined(Main, :PNJLQuarkDistributions_Aniso)
    Base.include(Main, normpath(joinpath(@__DIR__, "..", "QuarkDistribution_Aniso.jl")))
end
if !isdefined(Main, :ParameterAdapters)
    Base.include(Main, normpath(joinpath(@__DIR__, "..", "utils", "ParameterAdapters.jl")))
end
if !isdefined(Main, :ValidationUtils)
    Base.include(Main, normpath(joinpath(@__DIR__, "..", "utils", "ValidationUtils.jl")))
end

module RelaxTime

# Reuse legacy Main modules and publish new modules immediately, so later
# includes resolve the same types, functions and caches through either path.
macro shared_module(name, source)
    return esc(quote
        if isdefined(Main, $(QuoteNode(name)))
            const $name = getfield(Main, $(QuoteNode(name)))
        else
            $source
            Core.eval(Main, Expr(:const, Expr(:(=), $(QuoteNode(name)), $name)))
        end
    end)
end

# ── Integration & utility submodules ──
@shared_module ParameterAdapters include(joinpath(@__DIR__, "..", "utils", "ParameterAdapters.jl"))
@shared_module ValidationUtils include(joinpath(@__DIR__, "..", "utils", "ValidationUtils.jl"))
@shared_module GaussLegendre include(joinpath(@__DIR__, "..", "integration", "GaussLegendre.jl"))
@shared_module ParticleSymbols include(joinpath(@__DIR__, "..", "utils", "ParticleSymbols.jl"))
@shared_module PhaseSpaceSampling include(joinpath(@__DIR__, "..", "integration", "PhaseSpaceSampling.jl"))
@shared_module KinematicChecks include("KinematicChecks.jl")

# ── Relaxtime submodules (dependency order) ──
@shared_module OneLoopIntegrals include("OneLoopIntegrals.jl")
@shared_module CausalSpectralBubble include("CausalSpectralBubble.jl")
@shared_module OneLoopIntegralsCorrection include("OneLoopIntegralsAniso.jl")
@shared_module PolarizationAniso include("PolarizationAniso.jl")
@shared_module PolarizationCache include("PolarizationCache.jl")
@shared_module AFieldBuilder include("AFieldBuilder.jl")
@shared_module EffectiveCouplings include("EffectiveCouplings.jl")
@shared_module MesonInteractionKernel include("MesonInteractionKernel.jl")
@shared_module ChargedRPAKernel include("ChargedRPAKernel.jl")
@shared_module ChargedRPAProvider include("ChargedRPAProvider.jl")
@shared_module MesonRPA include("MesonRPA.jl")
@shared_module MesonRPAAdapter include("MesonRPAAdapter.jl")
@shared_module MesonPropagator include("MesonPropagator.jl")
@shared_module TotalPropagator include("TotalPropagator.jl")
@shared_module MesonMass include("MesonMass.jl")
@shared_module MottTransition include("MottTransition.jl")
@shared_module BUPhaseGates include("BUPhaseGates.jl")
@shared_module PhaseNormalization include("PhaseNormalization.jl")
@shared_module ChargedPhaseBackend include("ChargedPhaseBackend.jl")
@shared_module MesonDensity include("MesonDensity.jl")
@shared_module MesonThermodynamics include("MesonThermodynamics.jl")
@shared_module DifferentialCrossSection include("DifferentialCrossSection.jl")
@shared_module ScatteringAmplitude include("ScatteringAmplitude.jl")
@shared_module TotalCrossSection include("TotalCrossSection.jl")
@shared_module PairLoopQuadrature include("PairLoopQuadrature.jl")
@shared_module AverageScatteringRate include("AverageScatteringRate.jl")
@shared_module RelaxationTime include("RelaxationTime.jl")
@shared_module TransportCoefficientsValidation include("TransportCoefficientsValidation.jl")
@shared_module TransportCoefficients include("TransportCoefficients.jl")

# ── Re-exports ──
using .OneLoopIntegrals
using .OneLoopIntegralsCorrection
using .PolarizationAniso
using .PolarizationCache
using .AFieldBuilder
using .EffectiveCouplings
using .MesonInteractionKernel
using .ChargedRPAKernel
using .ChargedRPAProvider
using .MesonRPA
using .MesonRPAAdapter
using .MesonPropagator
using .TotalPropagator
using .MesonMass
using .MottTransition
using .BUPhaseGates
using .PhaseNormalization
using .ChargedPhaseBackend
using .MesonDensity
using .MesonThermodynamics
using .DifferentialCrossSection
using .ScatteringAmplitude
using .TotalCrossSection
using .AverageScatteringRate
using .RelaxationTime
using .TransportCoefficients
using .KinematicChecks

end # module RelaxTime
