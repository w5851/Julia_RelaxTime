#!/usr/bin/env julia
"""
测试粒子解析
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "utils", "ParticleSymbols.jl"))

using .ParticleSymbols: parse_scattering_process, parse_scattering_process_flavors

println("="^80)
println("粒子解析测试")
println("="^80)
println()

processes = [
    :ssbar_to_uubar,
    :uubar_to_ssbar,
    :ssbar_to_ssbar,
    :uubar_to_uubar,
]

for proc in processes
    particles = parse_scattering_process(proc)
    flavors = parse_scattering_process_flavors(proc)
    println("$proc:")
    println("  粒子: $particles")
    println("  味: $flavors")
    println()
end

println("="^80)
