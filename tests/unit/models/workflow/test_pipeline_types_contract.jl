using Dates: DateTime
using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "pipeline types contract" begin
    @testset "type existence and field contracts" begin
        @test isdefined(Models, :PipelineSpec)
        @test isdefined(Models, :PipelineStage)
        @test isdefined(Models, :PipelineContext)
        @test isdefined(Models, :PipelineArtifact)
        @test isdefined(Models, :StageResult)

        @test fieldnames(Models.PipelineSpec) == (
            :name,
            :version,
            :model_kind,
            :stages,
            :params,
            :io_contract,
        )
        @test fieldtype(Models.PipelineSpec, :name) == String
        @test fieldtype(Models.PipelineSpec, :version) == String
        @test fieldtype(Models.PipelineSpec, :model_kind) == Symbol
        @test fieldtype(Models.PipelineSpec, :stages) == Vector{Symbol}
        @test fieldtype(Models.PipelineSpec, :params) == NamedTuple

        @test fieldnames(Models.PipelineStage) == (:id, :requires, :provides, :run!)
        @test fieldtype(Models.PipelineStage, :id) == Symbol
        @test fieldtype(Models.PipelineStage, :requires) == Vector{Symbol}
        @test fieldtype(Models.PipelineStage, :provides) == Vector{Symbol}
        @test fieldtype(Models.PipelineStage, :run!) == Function

        ctx_fields = fieldnames(Models.PipelineContext)
        @test ctx_fields == (:state, :provenance)
        @test fieldtype(Models.PipelineContext, :state) == Dict{Symbol, Any}

        @test fieldnames(Models.PipelineArtifact) == (:path, :hash, :schema_version)
        @test fieldtype(Models.PipelineArtifact, :path) == String
        @test fieldtype(Models.PipelineArtifact, :hash) == String
        @test fieldtype(Models.PipelineArtifact, :schema_version) == String

        @test fieldnames(Models.StageResult) == (:produced, :artifacts, :metrics)
        @test fieldtype(Models.StageResult, :produced) == Dict{Symbol, Any}
        @test fieldtype(Models.StageResult, :artifacts) == Vector{Models.PipelineArtifact}
        @test fieldtype(Models.StageResult, :metrics) == Dict{Symbol, Float64}

        provenance = (
            git_commit="abc123",
            config_hash="cfg001",
            run_id="run-1",
            timestamp=DateTime(2026, 4, 10, 12, 0, 0),
        )
        ctx = Models.PipelineContext(Dict{Symbol, Any}(:x => 1), provenance)
        @test hasproperty(ctx.provenance, :git_commit)
        @test hasproperty(ctx.provenance, :config_hash)
        @test hasproperty(ctx.provenance, :run_id)
        @test hasproperty(ctx.provenance, :timestamp)
    end

    @testset "symbol string normalization" begin
        @test isdefined(Models, :persisted_symbol_to_string)
        @test isdefined(Models, :persisted_string_to_symbol)

        @test Models.persisted_symbol_to_string(:Stage_Init) == "stage_init"
        @test Models.persisted_symbol_to_string(Symbol("phase-scan")) == "phase_scan"

        @test Models.persisted_string_to_symbol(" Stage Init ") == :stage_init
        @test Models.persisted_string_to_symbol("phase-scan") == :phase_scan
        @test Models.persisted_string_to_symbol("phase__scan") == :phase_scan

        key = :transport_stage_1
        @test Models.persisted_string_to_symbol(Models.persisted_symbol_to_string(key)) == key

        @test_throws ArgumentError Models.persisted_string_to_symbol("")
        @test_throws ArgumentError Models.persisted_string_to_symbol("   ")
        @test_throws ArgumentError Models.persisted_string_to_symbol("bad/value")
    end

    @testset "basic constructor validation" begin
        @test_throws ArgumentError Models.PipelineSpec(
            "core",
            "v1",
            :PNJL,
            Symbol[],
            (alpha=1.0,),
            (;),
        )

        @test_throws ArgumentError Models.PipelineStage(
            "load",
            Symbol[:input],
            Symbol[:output],
            (ctx) -> nothing,
        )

        @test_throws ArgumentError Models.PipelineContext(
            Dict{Any, Any}("x" => 1),
            (
                git_commit="abc123",
                config_hash="cfg001",
                run_id="run-1",
                timestamp=DateTime(2026, 4, 10, 12, 0, 0),
            ),
        )
    end
end
