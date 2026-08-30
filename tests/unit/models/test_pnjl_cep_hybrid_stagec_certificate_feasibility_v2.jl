using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
const V2_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "analysis", "pnjl_cep_hybrid_stagec_certificate_feasibility_v2.jl")
if !isdefined(Main, :StageCCertificateFeasibilityV2)
    include(V2_SCRIPT)
end
const V2 = Main.StageCCertificateFeasibilityV2

@testset "Stage-C certificate feasibility v2 helpers" begin
    make_curve(mu) = (rho=Float64.(collect(0:(length(mu) - 1))), mu=Float64.(mu),
        points=[V2.Point(0.0, 1.0, Float64(i - 1), Float64(value), 0.0) for (i, value) in enumerate(mu)])

    @testset "weak S uses sign topology, not fixed slope margin" begin
        curve = make_curve([0.0, 0.01, 0.009, 0.008, 0.018, 0.028])
        candidates = V2._candidates(curve; level=:fine)
        @test length(candidates) == 1
        @test candidates[1].drop_mu > 0
        @test candidates[1].width >= 2.0
    end

    @testset "monotone curve has no candidate" begin
        curve = make_curve([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
        @test isempty(V2._candidates(curve; level=:fine))
    end

    @testset "candidate pairing distinguishes multiple intervals" begin
        coarse = make_curve([0.0, 1.0, 0.0, -1.0, -2.0, -1.0, 0.0, -1.0, -2.0, -1.0, 0.0])
        fine = make_curve([0.0, 1.0, 0.0, -1.0, -2.0, -1.0, 0.0, -1.0, -2.0, -1.0, 0.0])
        stable, coarse_candidates, fine_candidates = V2._stable_candidate_count(coarse, fine)
        @test length(coarse_candidates) == 2
        @test length(fine_candidates) == 2
        @test length(stable) == 2
    end

    @testset "targeted selection is deterministic and capped" begin
        pool = [V2.Point(0.0, 1.0, i / 128, i, 0.0) for i in 1:512]
        first = V2._select_points(pool, [1.0, 2.0], 12)
        second = V2._select_points(pool, [1.0, 2.0], 12)
        @test first == second
        @test length(first) == 12
        @test length(unique(point.rho for point in first)) == 12
    end

    @test V2.AUTHOR_FIRST_ORDER == Set([(-0.5, 147.0947265625), (0.5, 106.9599609375)])
    @test V2.CONSENSUS_MONOTONE == Set([(-0.5, 147.2197265625), (0.5, 107.0849609375)])
    @test V2.STAGE_C_FINE == 0.003125
    @test V2.CEP_AREA_TOL_GOOD == 1e-4
    @test V2.CEP_AREA_TOL_BAD == 5e-4
    @test V2.MAXWELL_SOLVER_TOL == 5e-6
end
