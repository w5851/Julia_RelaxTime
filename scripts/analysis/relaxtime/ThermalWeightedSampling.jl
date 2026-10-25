# Analysis compatibility entrypoint; the qualified algorithm has one implementation.
if !isdefined(Main,:ThermalAverageCandidates)
    include(joinpath(@__DIR__,"ThermalAverageCandidates.jl"))
end
if !isdefined(Main,:ThresholdCacheCandidates)
    include(joinpath(@__DIR__,"ThresholdCacheCandidates.jl"))
end
if !isdefined(Main,:ThermalWeightedSampling)
    const ThermalWeightedSampling=Main.RelaxTime.AverageScatteringRate.ThermalWeightedSampling
end
