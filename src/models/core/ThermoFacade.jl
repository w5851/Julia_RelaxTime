if !isdefined(Main, :ThermoFacade)
    include(normpath(joinpath(@__DIR__, "..", "pnjl_physics", "core", "ThermoFacade.jl")))
end
