if !isdefined(Main, :ThermoFacade)
    include(normpath(joinpath(@__DIR__, "..", "pnjl", "core", "ThermoFacade.jl")))
end
