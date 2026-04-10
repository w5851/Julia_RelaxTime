using Dates: DateTime
using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "pipeline types contract" begin
    @testset "type existence and field contracts" begin
        @test isdefined(Models, :AbstractPipelineIOContract)
        @test isdefined(Models, :PipelineIOContract)
        @test isdefined(Models, :PipelineProvenance)
        @test isdefined(Models, :PipelineSpec)
        @test isdefined(Models, :PipelineStage)
        @test isdefined(Models, :PipelineContext)
        @test isdefined(Models, :PipelineArtifact)
        @test isdefined(Models, :StageResult)

        @test fieldnames(Models.PipelineIOContract) == (
            :contract_version,
            :required_inputs,
            :required_outputs,
            :artifact_schema_version,
            :manifest_schema_version,
        )
        @test fieldtype(Models.PipelineIOContract, :contract_version) == Symbol
        @test fieldtype(Models.PipelineIOContract, :required_inputs) == Vector{Symbol}
        @test fieldtype(Models.PipelineIOContract, :required_outputs) == Vector{Symbol}
        @test fieldtype(Models.PipelineIOContract, :artifact_schema_version) == Symbol
        @test fieldtype(Models.PipelineIOContract, :manifest_schema_version) == Symbol

        @test fieldnames(Models.PipelineSpec) == (
            :name,
            :version,
            :model_kind,
            :stages,
            :params,
            :io_contract,
        )
        @test fieldtype(Models.PipelineSpec{NamedTuple{(:alpha,), Tuple{Float64}}, Models.PipelineIOContract}, :name) == String
        @test fieldtype(Models.PipelineSpec{NamedTuple{(:alpha,), Tuple{Float64}}, Models.PipelineIOContract}, :version) == String
        @test fieldtype(Models.PipelineSpec{NamedTuple{(:alpha,), Tuple{Float64}}, Models.PipelineIOContract}, :model_kind) == Symbol
        @test fieldtype(Models.PipelineSpec{NamedTuple{(:alpha,), Tuple{Float64}}, Models.PipelineIOContract}, :stages) == Vector{Symbol}
        @test fieldtype(Models.PipelineSpec{NamedTuple{(:alpha,), Tuple{Float64}}, Models.PipelineIOContract}, :params) == NamedTuple{(:alpha,), Tuple{Float64}}
        @test fieldtype(Models.PipelineSpec{NamedTuple{(:alpha,), Tuple{Float64}}, Models.PipelineIOContract}, :io_contract) == Models.PipelineIOContract

        @test fieldnames(Models.PipelineStage{typeof(identity)}) == (:id, :requires, :provides, :run!)
        @test fieldtype(Models.PipelineStage{typeof(identity)}, :id) == Symbol
        @test fieldtype(Models.PipelineStage{typeof(identity)}, :requires) == Vector{Symbol}
        @test fieldtype(Models.PipelineStage{typeof(identity)}, :provides) == Vector{Symbol}
        @test fieldtype(Models.PipelineStage{typeof(identity)}, :run!) == typeof(identity)

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

        io_contract = Models.PipelineIOContract(
            :v1,
            [:model_kind, :grid],
            [:artifact_paths],
            :artifact_v1,
            :manifest_v1,
        )
        spec = Models.PipelineSpec(
            "phase_pipeline",
            "v1",
            :PNJL,
            [:build_model, :solve_points],
            (alpha=1.0,),
            io_contract,
        )
        @test spec.io_contract isa Models.PipelineIOContract
        @test spec.params isa NamedTuple{(:alpha,), Tuple{Float64}}

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
        io_contract = Models.PipelineIOContract(
            :v1,
            [:model_kind, :grid],
            [:artifact_paths],
            :artifact_v1,
            :manifest_v1,
        )

        @test_throws ArgumentError Models.PipelineSpec(
            "core",
            "v1",
            :PNJL,
            Symbol[],
            (alpha=1.0,),
            io_contract,
        )

        @test_throws ArgumentError Models.PipelineSpec(
            "phase_pipeline",
            "v1",
            :PNJL,
            ["build_model"],
            (alpha=1.0,),
            io_contract,
        )

        @test_throws ArgumentError Models.PipelineSpec(
            "phase_pipeline",
            "v1",
            :PNJL,
            [:build_model],
            (alpha=1.0,),
            Dict(:contract => :v1),
        )

        @test_throws ArgumentError Models.PipelineStage(
            "load",
            Symbol[:input],
            Symbol[:output],
            (ctx) -> nothing,
        )

        @test_throws ArgumentError Models.PipelineStage(
            :load,
            ["input"],
            Symbol[:output],
            (ctx) -> nothing,
        )

        @test_throws ArgumentError Models.PipelineStage(
            :load,
            Symbol[:input],
            ["output"],
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

        @test_throws ArgumentError Models.PipelineContext(
            Dict{Symbol, Any}(:x => 1),
            (
                git_commit="abc123",
                config_hash="cfg001",
                run_id="run-1",
            ),
        )

        @test_throws ArgumentError Models.PipelineContext(
            Dict{Symbol, Any}(:x => 1),
            (
                git_commit="abc123",
                config_hash="cfg001",
                run_id="run-1",
                timestamp="2026-04-10T12:00:00",
            ),
        )

        @test_throws ArgumentError Models.PipelineProvenance("", "cfg001", "run-1", DateTime(2026, 4, 10, 12, 0, 0))
        @test_throws ArgumentError Models.PipelineProvenance("abc123", "", "run-1", DateTime(2026, 4, 10, 12, 0, 0))
        @test_throws ArgumentError Models.PipelineProvenance("abc123", "cfg001", "", DateTime(2026, 4, 10, 12, 0, 0))

        @test_throws ArgumentError Models.PipelineIOContract(
            "v1",
            [:model_kind],
            [:artifact_paths],
            :artifact_v1,
            :manifest_v1,
        )
        @test_throws ArgumentError Models.PipelineIOContract(
            :v1,
            [:model_kind],
            [:artifact_paths],
            "artifact_v1",
            :manifest_v1,
        )
        @test_throws ArgumentError Models.PipelineIOContract(
            :v1,
            [:model_kind],
            [:artifact_paths],
            :artifact_v1,
            "manifest_v1",
        )

        @test_throws ArgumentError Models.StageResult(
            Dict{Symbol, Any}(:x => 1),
            Any[(path="out.csv", hash="abc", schema_version="v1")],
            Dict{Symbol, Any}(:elapsed_ms => 1.0),
        )

        @test_throws ArgumentError Models.StageResult(
            Dict{Symbol, Any}(:x => 1),
            Models.PipelineArtifact[],
            Dict{Symbol, Any}(:elapsed_ms => "bad"),
        )
    end
end
