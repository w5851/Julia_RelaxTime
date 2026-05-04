#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))

using JSON3
using SHA

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_META_PATH = joinpath(PROJECT_ROOT, "build", "JuliaRelaxTime.sysimage.json")
const DEFAULT_OUTPUT_DIR = joinpath(PROJECT_ROOT, "build", "releases")

Base.@kwdef struct CliConfig
    meta_path::String = DEFAULT_META_PATH
    output_dir::String = DEFAULT_OUTPUT_DIR
    overwrite::Bool = false
end

function parse_bool(value::AbstractString)
    normalized = lowercase(strip(value))
    normalized in ("1", "true", "yes", "on") && return true
    normalized in ("0", "false", "no", "off") && return false
    throw(ArgumentError("unsupported boolean value: $(value)"))
end

function parse_args(argv)::CliConfig
    meta_path = DEFAULT_META_PATH
    output_dir = DEFAULT_OUTPUT_DIR
    overwrite = false
    for arg in argv
        if startswith(arg, "--meta-path=")
            meta_path = normpath(split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--output-dir=")
            output_dir = normpath(split(arg, "="; limit=2)[2])
        elseif startswith(arg, "--overwrite=")
            overwrite = parse_bool(split(arg, "="; limit=2)[2])
        elseif arg == "--overwrite"
            overwrite = true
        elseif arg in ("-h", "--help")
            println("""
            Usage:
              julia --project=. scripts/dev/package_sysimage_release.jl [--meta-path=<path>] [--output-dir=<dir>] [--overwrite]

            Behavior:
              - read build/JuliaRelaxTime.sysimage.json
              - stage JuliaRelaxTime.<ext> and JuliaRelaxTime.sysimage.json
              - create release archive named by metadata.release_asset_name
              - write adjacent .sha256 checksum file
            """)
            exit(0)
        else
            throw(ArgumentError("unsupported argument: $(arg)"))
        end
    end
    return CliConfig(; meta_path, output_dir, overwrite)
end

function load_metadata(meta_path::String)
    isfile(meta_path) || error("sysimage metadata not found at $(meta_path); run scripts/dev/build_sysimage.jl first")
    return JSON3.read(read(meta_path, String))
end

function require_string(meta, key::Symbol)
    hasproperty(meta, key) || error("sysimage metadata missing field: $(key)")
    value = getproperty(meta, key)
    value isa AbstractString || error("sysimage metadata field $(key) must be a string")
    isempty(value) && error("sysimage metadata field $(key) must not be empty")
    return String(value)
end

function stage_payload!(stage_dir::String, sysimage_path::String, meta_path::String)
    mkpath(stage_dir)
    cp(sysimage_path, joinpath(stage_dir, basename(sysimage_path)); force=true)
    cp(meta_path, joinpath(stage_dir, basename(meta_path)); force=true)
end

function create_archive!(archive_path::String, stage_dir::String, archive_format::String)
    if archive_format == "zip"
        stage_pattern = replace(joinpath(stage_dir, "*"), "\\" => "\\\\")
        escaped_archive = replace(archive_path, "\\" => "\\\\")
        run(`powershell -NoProfile -Command "Compress-Archive -Path '$stage_pattern' -DestinationPath '$escaped_archive' -Force"`)
        return
    end
    if archive_format == "tar.gz"
        run(`tar -czf $archive_path -C $stage_dir .`)
        return
    end
    error("unsupported release archive format: $(archive_format)")
end

function write_checksum!(archive_path::String)
    checksum_path = archive_path * ".sha256"
    digest_hex = bytes2hex(open(sha256, archive_path))
    open(checksum_path, "w") do io
        println(io, "$(digest_hex)  $(basename(archive_path))")
    end
    return checksum_path
end

function main(argv)
    cfg = parse_args(argv)
    meta = load_metadata(cfg.meta_path)
    sysimage_path = normpath(require_string(meta, :sysimage_path))
    archive_name = require_string(meta, :release_asset_name)
    archive_format = require_string(meta, :release_archive_format)

    isfile(sysimage_path) || error("sysimage payload not found at $(sysimage_path); run scripts/dev/build_sysimage.jl first")

    mkpath(cfg.output_dir)
    archive_path = joinpath(cfg.output_dir, archive_name)
    checksum_path = archive_path * ".sha256"
    if !cfg.overwrite && (isfile(archive_path) || isfile(checksum_path))
        error("release archive or checksum already exists at $(cfg.output_dir); re-run with --overwrite")
    end

    if cfg.overwrite
        isfile(archive_path) && rm(archive_path; force=true)
        isfile(checksum_path) && rm(checksum_path; force=true)
    end

    temp_root = mktempdir()
    try
        stage_dir = joinpath(temp_root, splitext(archive_name)[1])
        stage_payload!(stage_dir, sysimage_path, cfg.meta_path)
        create_archive!(archive_path, stage_dir, archive_format)
        checksum_path = write_checksum!(archive_path)
        println("Packaged sysimage release archive: " * archive_path)
        println("Packaged sysimage checksum: " * checksum_path)
    finally
        rm(temp_root; recursive=true, force=true)
    end
end

main(ARGS)
