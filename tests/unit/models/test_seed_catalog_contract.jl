using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "seed catalog contract" begin
    @testset "FixedMu catalog order matches default source" begin
        θ = [0.5, 0.0]
        from_catalog = Models.seed_catalog(Models.FixedMu(), θ)
        from_default = Models.default_multiseed_candidates_5()

        @test length(from_catalog) == length(from_default)
        @test all(length.(from_catalog) .== 5)
        @test all(length.(from_default) .== 5)
        for i in eachindex(from_default)
            @test from_catalog[i] == from_default[i]
        end
    end

    @testset "non-FixedMu catalog keeps order and extends" begin
        mode = Models.FixedAsymmetricRho(0.2, 1.0, 0.0)
        θ = [0.5]
        from_catalog = Models.seed_catalog(mode, θ)
        base = Models.default_multiseed_candidates_5()

        @test length(from_catalog) == length(base)
        @test all(length.(from_catalog) .== 8)
        for i in eachindex(base)
            @test from_catalog[i][1:5] == base[i]
        end
    end

    @testset "single-source consumption guards" begin
        gap_solver_path = joinpath(PROJECT_ROOT, "src", "models", "solver", "runtime", "GapSolver.jl")
        weighted_path = joinpath(PROJECT_ROOT, "src", "models", "solver", "governance", "WeightedFallback.jl")

        gap_source = read(gap_solver_path, String)
        weighted_source = read(weighted_path, String)

        @test occursin("seed_catalog(", gap_source)
        @test occursin("seed_catalog(", weighted_source)
        @test !occursin("get_all_seeds(MultiSeed(), [T_fm], mode)", weighted_source)
    end
end
