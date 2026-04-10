using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const PNJL = Models.pnjl_module()

@testset "PNJL conserved-charge susceptibilities core" begin
    cases = (
        (T_fm=0.52, muB_fm=0.00, muQ_fm=0.00, muS_fm=0.00),
        (T_fm=0.56, muB_fm=0.21, muQ_fm=0.04, muS_fm=0.01),
    )
    kwargs = (; xi=0.0, p_num=8, t_num=4)

    for case in cases
        vals = (
            PNJL.chi2_B(case.T_fm, case.muB_fm; kwargs...),
            PNJL.chi3_B(case.T_fm, case.muB_fm; kwargs...),
            PNJL.chi4_B(case.T_fm, case.muB_fm; kwargs...),
            PNJL.chi3_Q(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; kwargs...),
            PNJL.chi4_S(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; kwargs...),
            PNJL.chi2_Q(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; kwargs...),
            PNJL.chi2_S(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; kwargs...),
            PNJL.chi11_BQ(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; kwargs...),
            PNJL.chi11_BS(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; kwargs...),
            PNJL.chi11_QS(case.T_fm, case.muB_fm, case.muQ_fm, case.muS_fm; kwargs...),
            PNJL.baryon_Ssigma(case.T_fm, case.muB_fm; kwargs...),
            PNJL.baryon_kappa_sigma2(case.T_fm, case.muB_fm; kwargs...),
        )
        @test all(isfinite, vals)
    end

    T_fm = 0.55
    muB_fm = 0.24
    muQ_fm = 0.05
    muS_fm = 0.02
    V = 40.0
    c110 = PNJL.cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, V; orders=(1, 1, 0), kwargs...)
    chi110 = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), kwargs...)
    @test isapprox(c110, V * T_fm^3 * chi110; rtol=1e-12, atol=1e-12)
end
