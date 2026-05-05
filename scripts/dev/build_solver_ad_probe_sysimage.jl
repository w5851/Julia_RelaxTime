#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using PackageCompiler
using Libdl
using Dates
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const BUILD_DIR = joinpath(PROJECT_ROOT, "build", "solver_ad_probe")
const SYSIMAGE_PATH = joinpath(BUILD_DIR, "JuliaRelaxTimeSolverAdProbe." * Libdl.dlext)
const META_PATH = joinpath(BUILD_DIR, "JuliaRelaxTimeSolverAdProbe.sysimage.json")
const PRECOMPILE_SCRIPT = joinpath(BUILD_DIR, "solver_ad_probe_precompile.jl")

mkpath(BUILD_DIR)

workload = """
using Pkg
Pkg.activate(raw\"$(replace(PROJECT_ROOT, "\\" => "\\\\"))\")
include(joinpath(raw\"$(replace(PROJECT_ROOT, "\\" => "\\\\"))\", "src", "models", "Models.jl"))
using .Models
using Main.Constants_PNJL: ħc_MeV_fm

model = Models.create_model(:PNJL)
T_fm = 150.0 / ħc_MeV_fm
mu_fm = 250.0 / ħc_MeV_fm

Models.solve(
    model,
    Models.FixedMu(),
    T_fm,
    mu_fm;
    seed_guess=Models.HADRON_SEED_5,
    xi=0.0,
    p_num=6,
    t_num=3,
    residual_norm_max=1e-4,
    iterations=24,
)
"""

open(PRECOMPILE_SCRIPT, "w") do io
    write(io, workload)
end

create_sysimage(
    ["CSV", "ForwardDiff", "JSON3", "NLsolve", "StaticArrays", "Test", "TestItemRunner"],
    project=PROJECT_ROOT,
    sysimage_path=SYSIMAGE_PATH,
    precompile_execution_file=PRECOMPILE_SCRIPT,
    cpu_target="generic",
)

git_commit = try
    readchomp(`git -C $(PROJECT_ROOT) rev-parse HEAD`)
catch
    "unknown"
end

open(META_PATH, "w") do io
    JSON3.pretty(io, (
        generated_at = Dates.format(now(), Dates.ISODateTimeFormat),
        julia_version = string(VERSION),
        sysimage_path = SYSIMAGE_PATH,
        git_commit = git_commit,
        precompile_script = PRECOMPILE_SCRIPT,
        workload = "models_solve_fixedmu_solver_ad",
    ))
end

println("Built solver-ad probe sysimage: " * SYSIMAGE_PATH)
println("Wrote probe metadata: " * META_PATH)
