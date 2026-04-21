using Test
using JSON3

if !isdefined(Main, :Models)
    include("../../../src/models/Models.jl")
end

function _first_data_run_id(path::String)
    payload = read(path, String)
    lines = split(payload, '\n')
    header_seen = false
    for raw in lines
        s = strip(raw)
        isempty(s) && continue
        startswith(s, "#") && continue
        if !header_seen
            header_seen = true
            continue
        end
        cols = split(s, ',')
        isempty(cols) && continue
        return strip(cols[1])
    end
    throw(ArgumentError("cross section csv has no data rows: $(path)"))
end

function _write_minimal_cross_section_config(path::String)
    open(path, "w") do io
        write(io, """
schema_version = "v1"

[scan.cross_section]
muB_MeV = 0.0
T_list_MeV = [150.0]
xi_list = [0.0]
processes = ["ud_to_ud"]

[scan.cross_section.energy]
mode = "list"
sqrt_s_list_MeV = [500.0]
""")
    end
    return path
end

function _write_strict_cross_section_config(path::String)
    open(path, "w") do io
        write(io, """
schema_version = "v1"
strict = true

[scan.cross_section]
muB_MeV = 0.0
T_list_MeV = [150.0]
xi_list = [0.0]
processes = ["ud_to_ud"]

[scan.cross_section.energy]
mode = "list"
sqrt_s_list_MeV = [500.0]
""")
    end
    return path
end

@testset "relaxtime orchestrator pipeline runner smoke" begin
    outdir = mktempdir()
    cfg_dir = mktempdir()
    cfg = _write_minimal_cross_section_config(joinpath(cfg_dir, "smoke_cfg.toml"))
    result = Models.run_relaxtime_orchestrator_pipeline(
        :cross_section;
        config_path=cfg,
        outdir=outdir,
        processes=["ud_to_ud"],
        resume=true,
        overwrite=true,
        fail_on_fallback=true,
    )

    @test hasproperty(result, :manifest_path)
    manifest_path = String(result.manifest_path)
    @test isfile(manifest_path)
    @test isfile(String(result.effective_config_path))
    @test isfile(String(result.consumption_report_path))
    @test isfile(String(result.cross_section_path))
    @test hasproperty(result, :resume)
    @test hasproperty(result, :overwrite)
    @test hasproperty(result, :fail_on_fallback)
    @test hasproperty(result, :fallback_used)
    @test result.resume === true
    @test result.overwrite === true
    @test result.fail_on_fallback === true
    @test result.fallback_used === false
    @test dirname(String(result.cross_section_path)) == String(result.outdir)
    @test _first_data_run_id(String(result.cross_section_path)) == String(result.run_id)

    manifest = JSON3.read(read(manifest_path, String))
    @test String(manifest.pipeline.pipeline_family) == "relaxtime_orchestrator"
    @test String(manifest.pipeline.run_id) == String(result.run_id)
end

@testset "relaxtime orchestrator pipeline default output dir" begin
    cfg_dir = mktempdir()
    cfg = _write_minimal_cross_section_config(joinpath(cfg_dir, "default_cfg.toml"))
    result = Models.run_relaxtime_orchestrator_pipeline(
        :cross_section;
        config_path=cfg,
        processes=["ud_to_ud"],
    )

    manifest_path = String(result.manifest_path)
    @test isfile(manifest_path)
    @test dirname(manifest_path) == String(result.outdir)
    @test isfile(String(result.effective_config_path))
    @test isfile(String(result.consumption_report_path))
    @test isfile(String(result.cross_section_path))
    @test hasproperty(result, :resume)
    @test hasproperty(result, :overwrite)
    @test hasproperty(result, :fail_on_fallback)
    @test hasproperty(result, :fallback_used)
    @test result.resume === nothing
    @test result.overwrite === nothing
    @test result.fail_on_fallback === false
    @test result.fallback_used === false
    @test dirname(String(result.cross_section_path)) == String(result.outdir)
    @test _first_data_run_id(String(result.cross_section_path)) == String(result.run_id)

    manifest = JSON3.read(read(manifest_path, String))
    @test String(manifest.pipeline.run_id) == String(result.run_id)
end

@testset "relaxtime orchestrator pipeline strict cross-section keeps transport overrides isolated" begin
    outdir = mktempdir()
    cfg_dir = mktempdir()
    cfg = _write_strict_cross_section_config(joinpath(cfg_dir, "strict_cfg.toml"))

    result = Models.run_relaxtime_orchestrator_pipeline(
        :cross_section;
        config_path=cfg,
        outdir=outdir,
        processes=["ud_to_ud"],
        resume=true,
        overwrite=true,
    )

    report = JSON3.read(read(String(result.consumption_report_path), String))
    overridden = Set(String.(collect(report.overridden_keys)))
    unused = Set(String.(collect(report.unused_keys)))

    @test "scan.cross_section.processes" in overridden
    @test !("scan.transport.resume" in unused)
    @test !("scan.transport.overwrite" in unused)
end
