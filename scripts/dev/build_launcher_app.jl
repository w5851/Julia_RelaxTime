#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using Libdl
using PackageCompiler

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const LAUNCHER_ROOT = joinpath(PROJECT_ROOT, "app", "launcher")
const OUTPUT_ROOT = joinpath(PROJECT_ROOT, "build", "launcher-app")
const ROOT_SYSIMAGE = joinpath(PROJECT_ROOT, "build", "JuliaRelaxTime." * Libdl.dlext)
const LIGHT_SYSIMAGE_BUILD_ARGS = `-O0 --compile=min --strip-metadata --strip-ir`

function stage_log(msg)
    println("[build-launcher-app] " * msg)
    flush(stdout)
end

Pkg.activate(LAUNCHER_ROOT)
Pkg.resolve()
Pkg.instantiate()
Pkg.activate(PROJECT_ROOT)

ctx = PackageCompiler.create_pkg_context(LAUNCHER_ROOT)
ctx.env.pkg === nothing && error("launcher package missing name/uuid")

stage_log("reset output dir")
PackageCompiler.try_rm_dir(OUTPUT_ROOT; force=true)

stage_log("bundle julia libraries")
stdlibs = PackageCompiler.gather_stdlibs_project(ctx)
stdlibs = unique(vcat(stdlibs, map(pkg -> pkg.name, PackageCompiler.stdlibs_in_default_sysimage())))
PackageCompiler.bundle_julia_libraries(OUTPUT_ROOT, stdlibs)

stage_log("bundle julia libexec")
PackageCompiler.bundle_julia_libexec(ctx, OUTPUT_ROOT)

stage_log("bundle julia executable")
PackageCompiler.bundle_julia_executable(OUTPUT_ROOT)

stage_log("bundle artifacts/project/preferences/cert")
PackageCompiler.bundle_artifacts(ctx, OUTPUT_ROOT; include_lazy_artifacts=false)
PackageCompiler.bundle_project(ctx, OUTPUT_ROOT)
PackageCompiler.bundle_preferences(ctx, OUTPUT_ROOT)
PackageCompiler.bundle_cert(OUTPUT_ROOT)

sysimage_path = joinpath(OUTPUT_ROOT, "lib", "julia", "sys." * Libdl.dlext)
mkpath(dirname(sysimage_path))

base_sysimage = isfile(ROOT_SYSIMAGE) ? ROOT_SYSIMAGE : nothing
if base_sysimage === nothing
    stage_log("root sysimage missing; fallback to default Julia sysimage base")
else
    stage_log("use repo sysimage as base: " * base_sysimage)
end

stage_log("create thin launcher sysimage")
PackageCompiler.create_sysimage(
    [ctx.env.pkg.name];
    sysimage_path=sysimage_path,
    project=dirname(ctx.env.project_file),
    incremental=true,
    base_sysimage=base_sysimage,
    filter_stdlibs=false,
    include_transitive_dependencies=false,
    cpu_target=PackageCompiler.default_app_cpu_target(),
    sysimage_build_args=LIGHT_SYSIMAGE_BUILD_ARGS,
)

stage_log("create launcher executable")
PackageCompiler.create_executable_from_sysimg(
    joinpath(OUTPUT_ROOT, "bin", "jrt-launcher"),
    String(PackageCompiler.DEFAULT_EMBEDDING_WRAPPER),
    string(ctx.env.pkg.name, ".", "julia_main"),
)

println("Built launcher app: " * OUTPUT_ROOT)
