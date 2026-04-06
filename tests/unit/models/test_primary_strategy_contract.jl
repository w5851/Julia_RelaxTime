using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const P = Models.pnjl_module()

@testset "PrimaryStrategy contract" begin
    @testset "defaults" begin
        s = Models.PrimaryStrategy()
        @test s.primary_method == :newton
        @test s.use_fallback === true
        @test s.fallback_method == :trust_region
        @test s.use_multiseed === false
        @test s.seed_strategy === nothing
    end

    @testset "kwargs bridge mapping" begin
        s = Models.PrimaryStrategy(
            primary_method=:trust_region,
            use_fallback=false,
            fallback_method=:newton,
            use_multiseed=true,
            seed_strategy=Models.MultiSeed(),
        )

        bridged = P._resolve_primary_strategy_kwargs((; primary_strategy=s))
        @test bridged.nlsolve_method == :trust_region
        @test bridged.trust_region_fallback === false
        @test bridged.fallback_method == :newton
        @test bridged.auto_multiseed_fallback === true
        @test bridged.seed_strategy isa Models.MultiSeed
    end

    @testset "explicit kwargs override strategy fields" begin
        s = Models.PrimaryStrategy(primary_method=:trust_region, use_fallback=false)
        bridged = P._resolve_primary_strategy_kwargs((;
            primary_strategy=s,
            nlsolve_method=:newton,
            trust_region_fallback=true,
        ))

        @test bridged.nlsolve_method == :newton
        @test bridged.trust_region_fallback === true
    end

    @testset "invalid strategy type rejected" begin
        @test_throws ArgumentError P._resolve_primary_strategy_kwargs((; primary_strategy=(a=1,)))
    end
end
