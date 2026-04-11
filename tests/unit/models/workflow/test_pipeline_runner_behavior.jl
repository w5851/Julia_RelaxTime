using Dates: DateTime
using JSON3
using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

function _mk_contract(; required_outputs::Vector{Symbol}=Symbol[])
    return Models.PipelineIOContract(
        :v1,
        Symbol[],
        required_outputs,
        :artifact_v1,
        :manifest_v1,
    )
end

function _mk_spec(stage_ids::Vector{Symbol}; required_outputs::Vector{Symbol}=Symbol[])
    return Models.PipelineSpec(
        "runner_test_pipeline",
        "v1",
        :PNJL,
        stage_ids,
        (;),
        _mk_contract(; required_outputs=required_outputs),
    )
end

function _mk_ctx()
    return Models.PipelineContext(
        Dict{Symbol, Any}(),
        (
            git_commit="abc123",
            config_hash="cfg001",
            run_id="run-1",
            timestamp=DateTime(2026, 4, 10, 12, 0, 0),
        ),
    )
end

function _ok_stage(id::Symbol; requires::Vector{Symbol}=Symbol[], provides::Vector{Symbol}=Symbol[])
    return Models.PipelineStage(
        id,
        requires,
        provides,
        (ctx) -> Models.StageResult(
            Dict{Symbol, Any}(sym => "ok" for sym in provides),
            Models.PipelineArtifact[],
            Dict{Symbol, Float64}(),
        ),
    )
end

@testset "pipeline runner behavior" begin
    @testset "missing dependency throws before run" begin
        spec = _mk_spec([:a])
        stages = [
            _ok_stage(:a; requires=[:missing_input], provides=[:a_out]),
        ]
        manifest_path = joinpath(mktempdir(), "manifest_missing_dep.json")

        @test_throws ArgumentError Models.run_pipeline(spec, stages, _mk_ctx(); manifest_path=manifest_path)
    end

    @testset "stage provided but not listed in spec throws before run" begin
        spec = _mk_spec([:a])
        stages = [
            _ok_stage(:a; provides=[:x]),
            _ok_stage(:extra; provides=[:y]),
        ]
        manifest_path = joinpath(mktempdir(), "manifest_extra_stage.json")

        @test_throws ArgumentError Models.run_pipeline(spec, stages, _mk_ctx(); manifest_path=manifest_path)
    end

    @testset "duplicate stage id throws before run" begin
        spec = _mk_spec([:a, :b])
        stages = [
            _ok_stage(:dup; provides=[:x]),
            _ok_stage(:dup; provides=[:y]),
        ]
        manifest_path = joinpath(mktempdir(), "manifest_dup_stage.json")

        @test_throws ArgumentError Models.run_pipeline(spec, stages, _mk_ctx(); manifest_path=manifest_path)
    end

    @testset "duplicate provides throws before run" begin
        spec = _mk_spec([:a, :b])
        stages = [
            _ok_stage(:a; provides=[:shared]),
            _ok_stage(:b; provides=[:shared]),
        ]
        manifest_path = joinpath(mktempdir(), "manifest_dup_provides.json")

        @test_throws ArgumentError Models.run_pipeline(spec, stages, _mk_ctx(); manifest_path=manifest_path)
    end

    @testset "cycle detection throws before run" begin
        spec = _mk_spec([:a, :b])
        stages = [
            _ok_stage(:a; requires=[:y], provides=[:x]),
            _ok_stage(:b; requires=[:x], provides=[:y]),
        ]
        manifest_path = joinpath(mktempdir(), "manifest_cycle.json")

        @test_throws ArgumentError Models.run_pipeline(spec, stages, _mk_ctx(); manifest_path=manifest_path)
    end

    @testset "fail-fast marks downstream skipped" begin
        spec = _mk_spec([:a, :b, :c])
        stages = [
            _ok_stage(:a; provides=[:x]),
            Models.PipelineStage(:b, [:x], [:y], (ctx) -> error("stage_b_failed")),
            _ok_stage(:c; requires=[:y], provides=[:z]),
        ]

        manifest_path = joinpath(mktempdir(), "manifest_fail_fast.json")
        result = Models.run_pipeline(spec, stages, _mk_ctx(); manifest_path=manifest_path)

        @test result.success == false
        @test result.failed_stage == :b
        @test result.error_kind == :ErrorException
        @test occursin("stage_b_failed", result.error_msg)
        @test result.completed_stages == [:a]

        by_id = Dict(rec.id => rec for rec in result.stage_records)
        @test by_id[:a].status == :completed
        @test by_id[:b].status == :failed
        @test by_id[:c].status == :skipped
        @test by_id[:c].started_at === nothing
        @test by_id[:c].ended_at === nothing
    end

    @testset "runtime dependency check before stage execution" begin
        spec = _mk_spec([:a, :b, :c])
        stages = [
            _ok_stage(:a; provides=[:x]),
            Models.PipelineStage(:b, [:late_dep], [:y], (ctx) -> Models.StageResult(Dict{Symbol, Any}(:y => 1), Models.PipelineArtifact[], Dict{Symbol, Float64}())),
            _ok_stage(:c; provides=[:late_dep]),
        ]
        manifest_path = joinpath(mktempdir(), "manifest_runtime_dep_guard.json")
        result = Models.run_pipeline(spec, stages, _mk_ctx(); manifest_path=manifest_path)

        @test result.success == false
        @test result.failed_stage == :b
        @test result.error_kind == :ArgumentError
        @test occursin("missing runtime dependencies", result.error_msg)
    end

    @testset "required_outputs only checked on success path" begin
        spec_success = _mk_spec([:a]; required_outputs=[:must_exist])
        stages_success = [_ok_stage(:a; provides=[:x])]
        success_manifest = joinpath(mktempdir(), "manifest_required_outputs_success.json")

        success_result = Models.run_pipeline(spec_success, stages_success, _mk_ctx(); manifest_path=success_manifest)
        @test success_result.success == false
        @test success_result.failed_stage === nothing
        @test success_result.error_kind == :missing_required_outputs

        spec_failure = _mk_spec([:a, :b]; required_outputs=[:must_exist])
        stages_failure = [
            _ok_stage(:a; provides=[:x]),
            Models.PipelineStage(:b, [:x], [:y], (ctx) -> error("boom_before_output_check")),
        ]
        failure_manifest = joinpath(mktempdir(), "manifest_required_outputs_failure.json")
        failure_result = Models.run_pipeline(spec_failure, stages_failure, _mk_ctx(); manifest_path=failure_manifest)

        @test failure_result.success == false
        @test failure_result.failed_stage == :b
        @test failure_result.error_kind != :missing_required_outputs
        @test occursin("boom_before_output_check", failure_result.error_msg)
    end

    @testset "manifest written on success and failure with required fields" begin
        spec_success = _mk_spec([:a], required_outputs=[:x])
        stages_success = [_ok_stage(:a; provides=[:x])]
        success_manifest_path = joinpath(mktempdir(), "manifest_success.json")
        success_result = Models.run_pipeline(spec_success, stages_success, _mk_ctx(); manifest_path=success_manifest_path)

        @test success_result.success == true
        @test isfile(success_manifest_path)
        success_manifest = JSON3.read(read(success_manifest_path, String))
        @test haskey(success_manifest, :pipeline)
        @test haskey(success_manifest, :stage_records)
        @test haskey(success_manifest.pipeline, :config_hash)
        @test haskey(success_manifest.pipeline, :artifact_hash)
        @test !isempty(String(success_manifest.pipeline.config_hash))
        @test !isempty(String(success_manifest.pipeline.artifact_hash))

        success_records = success_manifest.stage_records
        @test length(success_records) == 1
        record = success_records[1]
        @test haskey(record, :id)
        @test haskey(record, :status)
        @test haskey(record, :started_at)
        @test haskey(record, :ended_at)
        @test haskey(record, :error_kind)
        @test haskey(record, :error_msg)
        @test occursin(r"Z$", String(record.started_at))
        @test occursin(r"Z$", String(record.ended_at))

        spec_failure = _mk_spec([:a, :b, :c])
        stages_failure = [
            _ok_stage(:a; provides=[:x]),
            Models.PipelineStage(:b, [:x], [:y], (ctx) -> error("manifest_failure")),
            _ok_stage(:c; requires=[:y], provides=[:z]),
        ]
        failure_manifest_path = joinpath(mktempdir(), "manifest_failure.json")
        failure_result = Models.run_pipeline(spec_failure, stages_failure, _mk_ctx(); manifest_path=failure_manifest_path)

        @test failure_result.success == false
        @test isfile(failure_manifest_path)
        failure_manifest = JSON3.read(read(failure_manifest_path, String))
        failure_records = Dict(Symbol(rec.id) => rec for rec in failure_manifest.stage_records)
        @test String(failure_records[:c].status) == "skipped"
        @test failure_records[:c].started_at === nothing
        @test failure_records[:c].ended_at === nothing

        success_manifest_path_2 = joinpath(mktempdir(), "manifest_success_2.json")
        _ = Models.run_pipeline(spec_success, stages_success, _mk_ctx(); manifest_path=success_manifest_path_2)
        success_manifest_2 = JSON3.read(read(success_manifest_path_2, String))
        @test String(success_manifest.pipeline.config_hash) == String(success_manifest_2.pipeline.config_hash)
        @test String(success_manifest.pipeline.artifact_hash) == String(success_manifest_2.pipeline.artifact_hash)

        success_manifest_path_3 = joinpath(mktempdir(), "manifest_success_3.json")
        changed_ctx = _mk_ctx()
        changed_ctx.state[:seed] = 42
        _ = Models.run_pipeline(spec_success, stages_success, changed_ctx; manifest_path=success_manifest_path_3)
        success_manifest_3 = JSON3.read(read(success_manifest_path_3, String))
        @test String(success_manifest.pipeline.artifact_hash) != String(success_manifest_3.pipeline.artifact_hash)
    end

    @testset "artifact hash canonical stability for complex state" begin
        records = [
            Models.PipelineStageRecord(:a, :completed, "2026-04-10T12:00:00.000Z", "2026-04-10T12:00:01.000Z", nothing, nothing),
        ]

        state_dict_a = Dict{Symbol, Any}(
            :cfg => Dict("b" => 2, "a" => 1),
            :tags => Set(["x", "y", "z"]),
            :nested => (alpha=1, beta=[3, 2, 1]),
        )
        state_dict_b = Dict{Symbol, Any}(
            :nested => (beta=[3, 2, 1], alpha=1),
            :tags => Set(["z", "x", "y"]),
            :cfg => Dict("a" => 1, "b" => 2),
        )

        hash_a = Models.compute_pipeline_artifact_hash(state_dict_a, records)
        hash_b = Models.compute_pipeline_artifact_hash(state_dict_b, records)
        @test hash_a == hash_b

        state_symbol = Dict{Symbol, Any}(:key => :a)
        state_string = Dict{Symbol, Any}(:key => "a")
        hash_symbol = Models.compute_pipeline_artifact_hash(state_symbol, records)
        hash_string = Models.compute_pipeline_artifact_hash(state_string, records)
        @test hash_symbol != hash_string
    end
end
