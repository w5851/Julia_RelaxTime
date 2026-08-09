using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
const V2_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "analysis",
    "pnjl_stagec_density_certificate_feasibility_v2.jl")
if !isdefined(Main, :StageCDensityCertificateFeasibilityV2)
    include(V2_SCRIPT)
end
const V2 = Main.StageCDensityCertificateFeasibilityV2

@testset "Stage-C density certificate feasibility v2" begin
    make_point(rho, level) = V2.Point(0.0, 51.0, rho, 300.0 + rho, 0.0,
        "memoized_dense", level)

    @testset "source levels keep Stage-B and Stage-C disjoint" begin
        data = V2.SourceData(
            [make_point(0.0, 0), make_point(0.00625, 1), make_point(0.003125, 1)],
            Any[], Any[], Any[], Any[], 0, "completed")
        fine = V2._points(data, "memoized_dense", (0.0, 51.0), _ -> true; levels=(0, 1))
        pool = V2._points(data, "memoized_dense", (0.0, 51.0),
            rho -> V2._on_grid(rho, V2.STAGE_C_FINE) &&
                !V2._on_grid(rho, V2.STAGE_B_FINE); levels=(1,))
        @test [point.level for point in fine] == [0, 1, 1]
        @test [point.rho for point in pool] == [0.003125]
    end

    @testset "missing candidate interval is not stable" begin
        details = Dict{Symbol, Any}(:candidate_count => 1, :crossing_count => 3)
        result = merge(V2._empty_eval("bisection_failed"),
            (raw_maxwell=Models.MaxwellResult(false, nothing, nothing, nothing, nothing,
                0, details),))
        @test !V2._candidate_stable(result, result)
    end

    @testset "route candidate sequences are deterministic and distinct" begin
        curve = (rho=[0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            mu=[300.0, 301.0, 300.0, 299.0, 300.0, 301.0],
            points=[make_point(rho, 0) for rho in 0.0:1.0:5.0])
        stage_b = V2._evaluate(curve)
        pool = [V2.Point(0.0, 51.0, rho, 300.0 + rho / 10, 0.0,
            "independent_oracle", 1) for rho in 0.1:0.1:4.9]
        balanced = V2._balanced_targets(stage_b, pool)
        @test balanced == V2._balanced_targets(stage_b, pool)
        @test length(unique(point.rho for point in pool)) == length(pool)
        @test V2.ROUTES == (:stage_b_features_v1, :balanced_density_features_v2,
            :geometry_feedback_v2)
    end

    @test V2.STAGE_A_COARSE == 0.05
    @test V2.STAGE_B_FINE == 0.00625
    @test V2.STAGE_C_FINE == 0.003125
    @test V2.AREA_TOL == 5e-5
end
