using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "adapter normalization and error classification" begin
    @testset "alias mapping normalizes to canonical keys" begin
        kwargs = (temperature=120.0, mu=10.0)
        alias_map = Dict("temperature" => :T_MeV)

        normalized = Models.normalize_adapter_kwargs(kwargs, alias_map)

        @test hasproperty(normalized, :T_MeV)
        @test !hasproperty(normalized, :temperature)
        @test normalized.T_MeV == 120.0
        @test normalized.mu == 10.0
    end

    @testset "canonical key takes precedence over alias" begin
        kwargs = (T_MeV=150.0, temperature=120.0)
        alias_map = Dict("temperature" => :T_MeV)

        normalized = Models.normalize_adapter_kwargs(kwargs, alias_map)

        @test normalized.T_MeV == 150.0
        @test !hasproperty(normalized, :temperature)
    end

    @testset "pipeline error classification" begin
        @test Models.classify_pipeline_error(ArgumentError("bad input")) == :input_validation_error
        @test Models.classify_pipeline_error(Base.IOError("disk unavailable", 5)) == :artifact_io_error
        @test Models.classify_pipeline_error(ErrorException("Solver did not CONVERGE in max iterations")) == :numerical_convergence_error
        @test Models.classify_pipeline_error(ErrorException("unknown failure")) == :unexpected_error
    end

    @testset "shared adapter helper contracts" begin
        hash_a = Models.adapter_config_hash(:scan, (a=1, b=2))
        hash_b = Models.adapter_config_hash(:scan, (b=2, a=1))
        hash_c = Models.adapter_config_hash(:workflow, (a=1, b=2))

        @test hash_a == hash_b
        @test hash_a != hash_c
        @test length(hash_a) == 64

        outdir_a = Models.resolve_adapter_output_dir((output_dir="abc",), (:output_dir, :outdir))
        outdir_b = Models.resolve_adapter_output_dir((outdir="xyz",), (:output_dir, :outdir))
        @test outdir_a == "abc"
        @test outdir_b == "xyz"

        outdir_from_path = Models.resolve_adapter_output_dir((output_path=joinpath("tmp", "x.csv"),), (:output_dir, :output_path))
        @test endswith(outdir_from_path, "tmp")
    end
end
