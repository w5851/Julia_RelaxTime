if !isdefined(Main, :PNJLQuarkDistributions)
    Base.include(Main, normpath(joinpath(@__DIR__, "..", "..", "QuarkDistribution.jl")))
end
