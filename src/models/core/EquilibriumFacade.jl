if !isdefined(Main, :EquilibriumFacade)
    include(normpath(joinpath(@__DIR__, "..", "pnjl", "core", "EquilibriumFacade.jl")))
end
