#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(@__DIR__, "..", "..", "src", "models", "Models.jl"))
using .Models

# Representative compile workload for AD-heavy solver/derivative paths.
Models.run_precompile_profile(:full)

# Real scan CLI entrypoints that dominate current cold-start cost.
scan_output_dir = mktempdir()
scan_cli = joinpath(@__DIR__, "..", "..", "scripts", "models", "run_unified_scan.jl")
run(`$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(scan_cli) scan tmu --model_kind=PNJL --T_values=150 --mu_values=0,100 --xi_values=0.0 --output_path=$(joinpath(scan_output_dir, "tmu.csv")) --overwrite=true`)
run(`$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(scan_cli) scan trho --model_kind=PNJL --T_values=150 --rho_values=0.1,0.2 --xi_values=0.0 --output_path=$(joinpath(scan_output_dir, "trho.csv")) --overwrite=true`)

# Light phase CLI run path used in integration smoke/core.
cli_script = joinpath(@__DIR__, "..", "..", "scripts", "pnjl", "calculate_phase_structure.jl")
output_dir = mktempdir()
run(`$(Base.julia_cmd()) --project=$(PROJECT_ROOT) $(cli_script) --preset=smoke --iterations=4 --p_num=6 --t_num=3 --output_dir=$(output_dir)`)
