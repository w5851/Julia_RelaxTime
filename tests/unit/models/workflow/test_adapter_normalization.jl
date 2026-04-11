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
end
