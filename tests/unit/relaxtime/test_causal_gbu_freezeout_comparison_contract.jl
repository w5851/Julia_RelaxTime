using Test
const _FGC_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(_FGC_ROOT, "scripts", "analysis", "relaxtime", "audit_causal_gbu_freezeout_comparison.jl"))
const testcausalgbufreezeoutcomparisoncontract_FGC = CausalGBUFreezeoutComparison

@testset "GBU freezeout comparison algebra" begin
    rows = NamedTuple[]
    for route in ("direct_finite_q", "q0_lambda_reference"), channel in ("pi_plus", "K_plus")
        for (i, q) in enumerate(testcausalgbufreezeoutcomparisoncontract_FGC.R.gauleg(0.0, 2.0, 8)[1])
            push!(rows, (sqrt_s_NN_GeV=200.0, channel=channel, route=route,
                q_inv_fm=q, bound=1.0 + (route == "q0_lambda_reference"),
                unitary=-0.1, landau=0.02, shell_inv_fm2=2.0 + i/10,
                passed=true))
        end
    end
    out = testcausalgbufreezeoutcomparisoncontract_FGC.aggregate_components(rows; qmax=2.0, nq=8)
    @test length(out) == 4
    @test all(r.passed for r in out)
    @test all(r.total_density_inv_fm3 ≈ r.bound_density_inv_fm3 +
        r.unitary_density_inv_fm3 + r.landau_density_inv_fm3 for r in out)
    @test_throws ErrorException testcausalgbufreezeoutcomparisoncontract_FGC.aggregate_components(rows[1:end-1]; qmax=2.0, nq=8)

    direct = [(sqrt_s_NN_GeV=200.0, channel="pi_plus", density_inv_fm3=2.0),
              (sqrt_s_NN_GeV=200.0, channel="K_plus", density_inv_fm3=1.0),
              (sqrt_s_NN_GeV=200.0, channel="pi_minus", density_inv_fm3=2.5),
              (sqrt_s_NN_GeV=200.0, channel="K_minus", density_inv_fm3=0.5),
              (sqrt_s_NN_GeV=7.7, channel="pi_plus", density_inv_fm3=2.0),
              (sqrt_s_NN_GeV=7.7, channel="K_plus", density_inv_fm3=1.0),
              (sqrt_s_NN_GeV=7.7, channel="pi_minus", density_inv_fm3=2.5),
              (sqrt_s_NN_GeV=7.7, channel="K_minus", density_inv_fm3=0.5),
              (sqrt_s_NN_GeV=3.0, channel="pi_plus", density_inv_fm3=2.0),
              (sqrt_s_NN_GeV=3.0, channel="K_plus", density_inv_fm3=1.0),
              (sqrt_s_NN_GeV=3.0, channel="pi_minus", density_inv_fm3=2.5),
              (sqrt_s_NN_GeV=3.0, channel="K_minus", density_inv_fm3=0.5)]
    reference = [merge(r, (density_inv_fm3=2r.density_inv_fm3,)) for r in direct]
    rr = testcausalgbufreezeoutcomparisoncontract_FGC.ratio_comparison(direct, reference)
    @test length(rr) == 3
    @test all(r.reference_Kplus_over_pi_plus ≈ r.direct_Kplus_over_pi_plus for r in rr)
    @test all(r.Kplus_density_factor_reference_over_direct ≈ 2.0 for r in rr)
end
