if !isdefined(Main, :ParameterTypes)
    Base.include(Main, normpath(joinpath(@__DIR__, "..", "ParameterTypes.jl")))
end
