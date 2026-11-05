using Test

if !isdefined(Main, :_compute_relaxtime_literature_transport_point)
    Base.include(Main, joinpath(@__DIR__, "..", "..", "validation", "relaxtime", "literature_validation_helpers.jl"))
end

module KappaLambdaThermalWeightedRegression
using Test, CSV, JSON3, SHA

const BASELINE = joinpath(@__DIR__, "..", "..", "baselines", "relaxtime", "baseline_kappa_lambda_thermal_weighted_v2.csv")

@testset "Thermal-weighted kappa/lambda fixed-point regression" begin
    rows = collect(CSV.File(BASELINE))
    metadata = JSON3.read(read(splitext(BASELINE)[1] * ".json", String))
    @test bytes2hex(sha256(read(BASELINE))) == metadata.baseline_sha256
    @test length(rows) == 14
    @test length(unique(row.record_id for row in rows)) == 14
    @test Set((row.T_MeV, row.muB_MeV) for row in rows) ==
        Set(((210.0, 0.0), (250.0, 0.0), (210.0, 100.0), (250.0, 100.0)))

    for row in rows
        @testset "$(row.record_id)" begin
            @test row.field in ("kappa_BB", "kappa_QQ", "kappa_SS", "lambda")
            result = Main._compute_relaxtime_literature_transport_point(row.T_MeV, row.muB_MeV)
            @test result.equilibrium.converged
            actual = getproperty(result.transport, Symbol(row.field))
            @test isfinite(actual) && actual > 0
            @test isapprox(actual, row.expected; rtol=1e-6, atol=1e-8)
            # A separate recorded, tighter calculation bounds interpolation and
            # outer/t quadrature sensitivity at this fixed background.
            @test isapprox(actual, row.refined; rtol=1e-3, atol=1e-8)
        end
    end
end
end
