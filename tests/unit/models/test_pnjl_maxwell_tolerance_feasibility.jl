using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const FEASIBILITY_SCRIPT = joinpath(
    PROJECT_ROOT, "scripts", "analysis", "pnjl_maxwell_tolerance_feasibility.jl",
)
if !isdefined(Main, :MaxwellToleranceFeasibility)
    include(FEASIBILITY_SCRIPT)
end
const FEAS = Main.MaxwellToleranceFeasibility

@testset "Maxwell tolerance feasibility helpers" begin
    @test FEAS.AREA_TOLS == (1e-4, 5e-5, 1e-5, 5e-6)
    @test FEAS.STRICT_SOLVER_TOL == 5e-6
    @test FEAS.OUTER_AREA_GATE == 5e-5

    @testset "CSV rounding deduplicates only formatting duplicates" begin
        points = [
            FEAS.Point(0.0, 20.0, 0.003125, 300.0, 1e-8, 0, "grid"),
            FEAS.Point(0.0, 20.0, 0.0031254, 301.0, 1e-6, 1, "targeted"),
            FEAS.Point(0.0, 20.0, 0.00625, 302.0, 2e-8, 0, "grid"),
        ]
        deduped = FEAS._deduplicate_points(points)
        @test length(deduped) == 2
        @test deduped[1].mu == 300.0
        @test deduped[2].rho == 0.00625
    end

    @testset "layer reconstruction preserves coarse and fine union" begin
        points = [
            FEAS.Point(0.0, 20.0, Float64(i), 300.0 + i, 0.0, i <= 11 ? 0 : 1, "grid")
            for i in 0:14
        ]
        layers = FEAS._layer_curves(points)
        @test layers.coarse !== nothing
        @test layers.fine !== nothing
        @test length(layers.coarse.points) == 12
        @test length(layers.fine.points) == 15
    end

    @testset "candidate identity is independent of residual magnitude" begin
        sres = Models.SShapeResult(true, 2.0, 1.0, 0.5, 1.5, 2)
        valid_a = (status=:valid, has_s_shape=true, sres=sres)
        valid_b = (status=:valid, has_s_shape=true, sres=sres)
        monotone_a = (status=:invalid, has_s_shape=false,
            sres=Models.SShapeResult(false, nothing, nothing, nothing, nothing, 0))
        monotone_b = (status=:invalid, has_s_shape=false,
            sres=Models.SShapeResult(false, nothing, nothing, nothing, nothing, 0))
        @test FEAS._candidate_stable(valid_a, valid_b)
        @test FEAS._candidate_stable(monotone_a, monotone_b)
        @test FEAS._candidate_stable_across(valid_a, valid_b)
        @test !FEAS._candidate_stable(valid_a, monotone_a)
    end
end
