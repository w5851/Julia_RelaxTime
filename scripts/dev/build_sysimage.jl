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

@inline function platform_family()
    if Sys.iswindows()
        return "windows"
    elseif Sys.isapple()
        return "macos"
    elseif Sys.islinux()
        return "linux"
    end
    return lowercase(String(Sys.KERNEL))
end

@inline function release_archive_format(family::AbstractString)
    return family == "windows" ? "zip" : "tar.gz"
end

@inline function release_asset_name(family::AbstractString, arch::AbstractString, julia_version::AbstractString, archive_format::AbstractString)
    return "JuliaRelaxTime-sysimage-$(family)-$(arch)-julia$(julia_version).$(archive_format)"
end

mkpath(dirname(SYSIMAGE_PATH))

create_sysimage(
    ["CSV", "ForwardDiff", "JSON3", "NLsolve", "StaticArrays", "TaylorDiff", "Test", "TestItemRunner"],
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
    family = platform_family()
    arch = lowercase(String(Sys.ARCH))
    archive_format = release_archive_format(family)
    artifact_basename = basename(SYSIMAGE_PATH)
    JSON3.pretty(io, (
        generated_at = Dates.format(now(), Dates.ISODateTimeFormat),
        julia_version = string(VERSION),
        sysimage_path = SYSIMAGE_PATH,
        git_commit = git_commit,
        precompile_script = PRECOMPILE_SCRIPT,
        platform_family = family,
        platform_os = lowercase(String(Sys.KERNEL)),
        platform_arch = arch,
        artifact_basename = artifact_basename,
        release_asset_name = release_asset_name(family, arch, string(VERSION), archive_format),
        release_archive_format = archive_format,
    ))
end

println("Wrote sysimage metadata: " * SYSIMAGE_META_PATH)
