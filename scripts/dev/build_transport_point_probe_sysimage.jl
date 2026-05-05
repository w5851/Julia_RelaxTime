#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using PackageCompiler
using Libdl
using Dates
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const BUILD_DIR = joinpath(PROJECT_ROOT, "build", "transport_point_probe")
const SYSIMAGE_PATH = joinpath(BUILD_DIR, "JuliaRelaxTimeTransportPointProbe." * Libdl.dlext)
const META_PATH = joinpath(BUILD_DIR, "JuliaRelaxTimeTransportPointProbe.sysimage.json")
const PRECOMPILE_SCRIPT = joinpath(BUILD_DIR, "transport_point_probe_precompile.jl")

mkpath(BUILD_DIR)

workload = """
using Pkg
Pkg.activate(raw\"$(replace(PROJECT_ROOT, "\\" => "\\\\"))\")
include(joinpath(raw\"$(replace(PROJECT_ROOT, "\\" => "\\\\"))\", "src", "models", "Models.jl"))
using .Models
Models.run_precompile_capability(:gap_solver_ad; strict=true)
Models.run_precompile_capability(:ad_shape_stabilization; strict=true)
Models.run_precompile_capability(:transport_point_api; strict=true)
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
        workload = "gap_solver_ad + ad_shape_stabilization + transport_point_api",
    ))
end

println("Built transport-point probe sysimage: " * SYSIMAGE_PATH)
println("Wrote probe metadata: " * META_PATH)
