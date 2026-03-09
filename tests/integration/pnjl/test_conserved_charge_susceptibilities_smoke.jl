using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const PNJL = Models.pnjl_module()

@testset "PNJL conserved-charge susceptibilities smoke" begin
    cases = (
        (T_fm=0.52, muB_fm=0.00, muQ_fm=0.00, muS_fm=0.00),
        (T_fm=0.56, muB_fm=0.21, muQ_fm=0.04, muS_fm=0.01),
        (T_fm=0.62, muB_fm=0.30, muQ_fm=0.06, muS_fm=0.02),
    )

    for case in cases
        vals = (
            PNJL.chi2_B(case.T_fm, case.muB_fm; xi=0.0, p_num=16, t_num=6),
            PNJL.chi3_B(case.T_fm, case.muB_fm; xi=0.0, p_num=16, t_num=6),
            PNJL.chi4_B(case.T_fm, case.muB_fm; xi=0.0, p_num=16, t_num=6),
            PNJL.chi3_Q(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; xi=0.0, p_num=16, t_num=6),
            PNJL.chi4_S(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; xi=0.0, p_num=16, t_num=6),
            PNJL.chi2_Q(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; xi=0.0, p_num=16, t_num=6),
            PNJL.chi2_S(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; xi=0.0, p_num=16, t_num=6),
            PNJL.chi11_BQ(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; xi=0.0, p_num=16, t_num=6),
            PNJL.chi11_BS(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; xi=0.0, p_num=16, t_num=6),
            PNJL.chi11_QS(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; xi=0.0, p_num=16, t_num=6),
            PNJL.baryon_Ssigma(case.T_fm, case.muB_fm; xi=0.0, p_num=16, t_num=6),
            PNJL.baryon_kappa_sigma2(case.T_fm, case.muB_fm; xi=0.0, p_num=16, t_num=6),
        )
        @test all(isfinite, vals)
    end

    T_fm = 0.55
    muB_fm = 0.24
    muQ_fm = 0.05
    muS_fm = 0.02
    V = 40.0

    chi2 = PNJL.chi2_B(T_fm, muB_fm; xi=0.0, p_num=16, t_num=6)
    chi3 = PNJL.chi3_B(T_fm, muB_fm; xi=0.0, p_num=16, t_num=6)
    chi4 = PNJL.chi4_B(T_fm, muB_fm; xi=0.0, p_num=16, t_num=6)
    c110 = PNJL.cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, V; orders=(1, 1, 0), xi=0.0, p_num=16, t_num=6)
    chi110 = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), xi=0.0, p_num=16, t_num=6)

    @test isapprox(PNJL.baryon_Ssigma(T_fm, muB_fm; xi=0.0, p_num=16, t_num=6), chi3 / chi2; rtol=1e-10, atol=1e-12)
    @test isapprox(PNJL.baryon_kappa_sigma2(T_fm, muB_fm; xi=0.0, p_num=16, t_num=6), chi4 / chi2; rtol=1e-10, atol=1e-12)
    @test isapprox(c110, V * T_fm^3 * chi110; rtol=1e-12, atol=1e-12)
    @test abs(PNJL.chi1_B(0.5, 0.0; xi=0.0, p_num=16, t_num=6)) <= 1e-8
    @test abs(PNJL.chi3_B(0.5, 0.0; xi=0.0, p_num=16, t_num=6)) <= 1e-8
end