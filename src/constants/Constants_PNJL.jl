if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, normpath(joinpath(@__DIR__, "..", "Constants_PNJL.jl")))
end
