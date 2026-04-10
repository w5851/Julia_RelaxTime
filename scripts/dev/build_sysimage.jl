#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using PackageCompiler
using Libdl

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PRECOMPILE_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "dev", "precompile_workload.jl")
const SYSIMAGE_PATH = joinpath(PROJECT_ROOT, "build", "JuliaRelaxTime." * Libdl.dlext)

mkpath(dirname(SYSIMAGE_PATH))

create_sysimage(
    ["CSV", "ForwardDiff", "JSON3", "NLsolve", "StaticArrays", "Test", "TestItemRunner"],
    project=PROJECT_ROOT,
    sysimage_path=SYSIMAGE_PATH,
    precompile_execution_file=PRECOMPILE_SCRIPT,
    cpu_target="generic",
)

println("Built sysimage: " * SYSIMAGE_PATH)
