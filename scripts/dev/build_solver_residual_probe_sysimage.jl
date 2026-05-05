#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using PackageCompiler
using Libdl
using Dates
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const BUILD_DIR = joinpath(PROJECT_ROOT, "build", "solver_residual_probe")
const SYSIMAGE_PATH = joinpath(BUILD_DIR, "JuliaRelaxTimeSolverResidualProbe." * Libdl.dlext)
const META_PATH = joinpath(BUILD_DIR, "JuliaRelaxTimeSolverResidualProbe.sysimage.json")
const PRECOMPILE_SCRIPT = joinpath(BUILD_DIR, "solver_residual_probe_precompile.jl")

mkpath(BUILD_DIR)

workload = """
using Pkg
Pkg.activate(raw\"$(replace(PROJECT_ROOT, "\\" => "\\\\"))\")
include(joinpath(raw\"$(replace(PROJECT_ROOT, "\\" => "\\\\"))\", "src", "models", "Models.jl"))
using .Models
using ForwardDiff
using Main.Constants_PNJL: ħc_MeV_fm

model = Models.create_model(:PNJL)
T_fm = 150.0 / ħc_MeV_fm
mu_fm = 250.0 / ħc_MeV_fm
residual! = Models.Conditions.build_residual!(
    model,
    Models.FixedMu(),
    T_fm,
    mu_fm;
    xi=0.0,
    p_num=6,
    t_num=3,
)
x0 = Float64.(Models.HADRON_SEED_5)

vec_residual = x -> begin
    F = similar(x)
    residual!(F, x)
    return F
end

ForwardDiff.jacobian(vec_residual, x0)
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
        workload = "conditions_build_residual_forwarddiff_jacobian",
    ))
end

println("Built solver-residual probe sysimage: " * SYSIMAGE_PATH)
println("Wrote probe metadata: " * META_PATH)
