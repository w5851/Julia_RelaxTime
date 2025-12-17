#!/usr/bin/env julia

include(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "AverageScatteringRate.jl"))

using .AverageScatteringRate

c = AverageScatteringRate.CrossSectionCache(:udbar_to_udbar)
println("ok: rtol=", c.rtol, " max_refine=", c.max_refine)
