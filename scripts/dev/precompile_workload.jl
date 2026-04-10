#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
using .Models

# Representative compile workload for AD-heavy solver/derivative paths.
Models.run_precompile_profile(:full)

# Light phase CLI run path used in integration smoke/core.
cli_script = joinpath(@__DIR__, "..", "..", "scripts", "pnjl", "calculate_phase_structure.jl")
output_dir = mktempdir()
run(`$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(cli_script) --preset=smoke --iterations=4 --p_num=6 --t_num=3 --output_dir=$(output_dir)`)
