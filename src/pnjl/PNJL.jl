if !isdefined(Main, :Constants_PNJL)
	include(joinpath(@__DIR__, "..", "Constants_PNJL.jl"))
end

module PNJL

using ..Constants_PNJL

include("SeedCache.jl")
include(joinpath("solvers", "AnisoGapSolver.jl"))
include(joinpath("seeds", "TrhoSeedChain.jl"))
include("SinglePointSolver.jl")
include(joinpath("scans", "TmuScan.jl"))
include(joinpath("scans", "TrhoScan.jl"))
include(joinpath("scans", "AdaptiveRhoRefinement.jl"))
include(joinpath("analysis", "PhaseTransition.jl"))
include(joinpath("analysis", "CEPFinder.jl"))
include(joinpath("analysis", "MaxwellRhoMu.jl"))
include(joinpath("analysis", "ThermoDerivatives.jl"))

end # module
