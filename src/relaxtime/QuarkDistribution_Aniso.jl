if !isdefined(Main, :PNJLQuarkDistributions_Aniso)
    Base.include(Main, normpath(joinpath(@__DIR__, "..", "QuarkDistribution_Aniso.jl")))
end
