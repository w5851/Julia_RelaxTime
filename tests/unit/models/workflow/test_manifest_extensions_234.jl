using Dates: DateTime
using JSON3
using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _mk_contract()
    return Models.PipelineIOContract(
        :v1,
        Symbol[],
        Symbol[:x],
        :artifact_v1,
        :manifest_v1,
    )
end

function _mk_spec(stage_ids::Vector{Symbol})
    return Models.PipelineSpec(
        "manifest_extensions_pipeline",
        "v1",
        :PNJL,
        stage_ids,
        (;),
        _mk_contract(),
    )
end

function _mk_ctx()
    return Models.PipelineContext(
        Dict{Symbol, Any}(),
        (
            git_commit="abc123",
            config_hash="cfg001",
            run_id="run-1",
            timestamp=DateTime(2026, 4, 11, 12, 0, 0),
        ),
    )
end

function _ok_stage(id::Symbol; provides::Vector{Symbol}=Symbol[])
    return Models.PipelineStage(
        id,
        Symbol[],
        provides,
        (ctx) -> Models.StageResult(
            Dict{Symbol, Any}(sym => "ok" for sym in provides),
            Models.PipelineArtifact[],
            Dict{Symbol, Float64}(),
        ),
    )
end

@testset "manifest extensions 2/3/4" begin
    @testset "build_manifest_extensions defaults and explicit values" begin
        defaults = Models.build_manifest_extensions((;))
        @test defaults["pipeline_family"] == "generic"
        @test defaults["baseline_suite"] == "none"
        @test defaults["physics_profile"] == "default"
        @test defaults["adapter_version"] == "v1"

        explicit = Models.build_manifest_extensions((
            pipeline_family="phase_scan",
            baseline_suite="regression_core",
            physics_profile="pnjl_2p1",
            adapter_version="v2.3",
        ))
        @test explicit["pipeline_family"] == "phase_scan"
        @test explicit["baseline_suite"] == "regression_core"
        @test explicit["physics_profile"] == "pnjl_2p1"
        @test explicit["adapter_version"] == "v2.3"
    end

    @testset "run_pipeline manifest merges extension fields" begin
        spec = _mk_spec([:a])
        stages = [_ok_stage(:a; provides=[:x])]
        ctx = _mk_ctx()
        ctx.state[:manifest_extensions] = Models.build_manifest_extensions((
            pipeline_family="scan",
            baseline_suite="nightly",
            physics_profile="core",
            adapter_version="v3",
        ))

        manifest_path = joinpath(mktempdir(), "manifest_extensions_234.json")
        result = Models.run_pipeline(spec, stages, ctx; manifest_path=manifest_path)

        @test result.success == true
        payload = JSON3.read(read(manifest_path, String))
        @test String(payload.pipeline.pipeline_family) == "scan"
        @test String(payload.pipeline.baseline_suite) == "nightly"
        @test String(payload.pipeline.physics_profile) == "core"
        @test String(payload.pipeline.adapter_version) == "v3"
    end

    @testset "run_pipeline normalizes repo-internal extension paths without parent segments" begin
        spec = _mk_spec([:a])
        stages = [_ok_stage(:a; provides=[:x])]
        ctx = _mk_ctx()
        abs_diag_path = joinpath(PROJECT_ROOT, "data", "outputs", "results", "diag", "index.json")
        ctx.state[:manifest_extensions] = Models.build_manifest_extensions((
            pipeline_family="scan",
            baseline_suite="nightly",
            physics_profile="core",
            adapter_version="v3",
            diagnostics_index_path=abs_diag_path,
        ))

        manifest_path = joinpath(mktempdir(), "manifest_extensions_path_norm_234.json")
        result = Models.run_pipeline(spec, stages, ctx; manifest_path=manifest_path)

        @test result.success == true
        payload = JSON3.read(read(manifest_path, String))
        diag_path = String(payload.pipeline.diagnostics_index_path)
        @test diag_path == "data/outputs/results/diag/index.json"
        @test !occursin("..", diag_path)
        @test !occursin('\\', diag_path)
    end

    @testset "build_manifest_extensions rejects non-whitelisted keys" begin
        @test_throws ArgumentError Models.build_manifest_extensions((
            pipeline_family="scan",
            unexpected_key="bad",
        ))
    end

    @testset "run_pipeline rejects reserved key injection via manifest_extensions" begin
        spec = _mk_spec([:a])
        stages = [_ok_stage(:a; provides=[:x])]
        ctx = _mk_ctx()
        ctx.state[:manifest_extensions] = Dict{String, Any}(
            "success" => false,
        )

        manifest_path = joinpath(mktempdir(), "manifest_extensions_reserved_injection_234.json")
        @test_throws ArgumentError Models.run_pipeline(spec, stages, ctx; manifest_path=manifest_path)
    end
end
