using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const REPLAY_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "analysis",
    "pnjl_cep_hybrid_stagec_tolerance_replay.jl")
if !isdefined(Main, :StageCCertificateFeasibilityToleranceReplay)
    include(REPLAY_SCRIPT)
end

const R = Main

@testset "Stage-C tolerance replay helpers" begin
    row(; rho, muq, xi=-0.5, T_MeV=20.0, method="independent_oracle") = (
        xi=xi, method=method, T_MeV=T_MeV, rho=rho, muq_MeV=muq,
        pressure_fm4=21.0, residual_norm=0.0, iterations=3,
        converged=true, finite=true,
    )

    base = [row(rho=0.0, muq=96.0), row(rho=1.0, muq=100.0)]
    tolerance = [row(rho=0.0, muq=101.0), row(rho=1.0, muq=100.0)]

    identity = R._numeric_identity(base, tolerance)
    @test identity.numeric_identity
    @test identity.rho_zero_common_keys == 1
    @test identity.max_abs_muq_MeV_rho_zero == 5.0
    @test identity.max_abs_muq_MeV_positive_rho == 0.0

    merged = R._merge_curves(base, tolerance)
    @test length(merged) == 2
    @test any(getproperty(item, :rho) == 0.0 &&
        getproperty(item, :muq_MeV) == 101.0 for item in merged)
end
