#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using PackageCompiler
using Libdl
using Dates
using JSON3

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PRECOMPILE_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "dev", "precompile_workload.jl")
const SYSIMAGE_PATH = joinpath(PROJECT_ROOT, "build", "JuliaRelaxTime." * Libdl.dlext)
const SYSIMAGE_META_PATH = joinpath(PROJECT_ROOT, "build", "JuliaRelaxTime.sysimage.json")

mkpath(dirname(SYSIMAGE_PATH))

create_sysimage(
    ["CSV", "ForwardDiff", "JSON3", "NLsolve", "StaticArrays", "Test", "TestItemRunner"],
    project=PROJECT_ROOT,
    sysimage_path=SYSIMAGE_PATH,
    precompile_execution_file=PRECOMPILE_SCRIPT,
    cpu_target="generic",
)

println("Built sysimage: " * SYSIMAGE_PATH)

git_commit = try
    readchomp(`git -C $(PROJECT_ROOT) rev-parse HEAD`)
catch
    "unknown"
end

open(SYSIMAGE_META_PATH, "w") do io
    JSON3.pretty(io, (
        generated_at = Dates.format(now(), Dates.ISODateTimeFormat),
        julia_version = string(VERSION),
        sysimage_path = SYSIMAGE_PATH,
        git_commit = git_commit,
        precompile_script = PRECOMPILE_SCRIPT,
    ))
end

println("Wrote sysimage metadata: " * SYSIMAGE_META_PATH)
