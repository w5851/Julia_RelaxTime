if !isdefined(Main, :EquilibriumFacade)
    include(normpath(joinpath(@__DIR__, "..", "pnjl_physics", "core", "EquilibriumFacade.jl")))
end
