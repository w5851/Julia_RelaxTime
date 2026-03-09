using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

const _CONSTANTS_PNJL_PATH = normpath(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
const _GAUSS_LEGENDRE_PATH = normpath(joinpath(PROJECT_ROOT, "src", "integration", "GaussLegendre.jl"))
const _MODELS_PATH = normpath(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PNJL_PATH)
end
if !isdefined(Main, :GaussLegendre)
    Base.include(Main, _GAUSS_LEGENDRE_PATH)
end
if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_PATH)
end

using .Models
using Main.Constants_PNJL: ħc_MeV_fm
const PNJL = Models.pnjl_module()

@testset "baryon susceptibilities basic finite" begin
    T_mev = 140.0
    muB_mev = 360.0
    T_fm = T_mev / ħc_MeV_fm
    muB_fm = muB_mev / ħc_MeV_fm

    chi1 = PNJL.chi1_B(T_fm, muB_fm; xi=0.0, p_num=48, t_num=12)
    chi2 = PNJL.chi2_B(T_fm, muB_fm; xi=0.0, p_num=48, t_num=12)
    chi3 = PNJL.chi3_B(T_fm, muB_fm; xi=0.0, p_num=48, t_num=12)
    chi4 = PNJL.chi4_B(T_fm, muB_fm; xi=0.0, p_num=48, t_num=12)

    @test all(isfinite, (chi1, chi2, chi3, chi4))
end

@testset "chi1_B matches n_B over T^3" begin
    T_mev = 130.0
    muB_mev = 300.0
    T_fm = T_mev / ħc_MeV_fm
    muB_fm = muB_mev / ħc_MeV_fm
    muq_fm = muB_fm / 3

    td = PNJL.thermo_derivatives(T_fm, muq_fm; xi=0.0, p_num=48, t_num=12)
    chi1 = PNJL.chi1_B(T_fm, muB_fm; xi=0.0, p_num=48, t_num=12)

    @test isapprox(chi1, td.rho / T_fm^3; atol=1e-8, rtol=1e-7)
end

@testset "baryon cumulant and ratio consistency" begin
    T_fm = 0.55
    muB_fm = 0.9
    V = 125.0

    chi2 = PNJL.chi2_B(T_fm, muB_fm; xi=0.0, p_num=40, t_num=10)
    chi3 = PNJL.chi3_B(T_fm, muB_fm; xi=0.0, p_num=40, t_num=10)
    chi4 = PNJL.chi4_B(T_fm, muB_fm; xi=0.0, p_num=40, t_num=10)

    c4 = PNJL.cumulant_B(T_fm, muB_fm, V; order=4, xi=0.0, p_num=40, t_num=10)
    s_sigma = PNJL.baryon_Ssigma(T_fm, muB_fm; xi=0.0, p_num=40, t_num=10)
    kappa_sigma2 = PNJL.baryon_kappa_sigma2(T_fm, muB_fm; xi=0.0, p_num=40, t_num=10)

    @test isapprox(c4, V * T_fm^3 * chi4; atol=1e-10, rtol=1e-10)
    @test isapprox(s_sigma, chi3 / chi2; atol=1e-10, rtol=1e-10)
    @test isapprox(kappa_sigma2, chi4 / chi2; atol=1e-10, rtol=1e-10)
end

@testset "odd-order baryon susceptibilities are near zero on muB=0 line" begin
    T_fm = 0.5
    muB_fm = 0.0

    chi1 = PNJL.chi1_B(T_fm, muB_fm; xi=0.0, p_num=32, t_num=10)
    chi3 = PNJL.chi3_B(T_fm, muB_fm; xi=0.0, p_num=32, t_num=10)
    @test abs(chi1) <= 1e-8
    @test abs(chi3) <= 1e-8
end

@testset "flavor pressure derivatives finite" begin
    T_fm = 0.55
    mu_vec = SVector(0.10, 0.06, 0.02)

    derivs = PNJL.flavor_pressure_derivatives(T_fm, mu_vec; order=2, xi=0.0, p_num=32, t_num=10)

    @test isfinite(derivs.pressure)
    @test all(isfinite.(derivs.grad_mu))
    @test all(isfinite.(derivs.hessian_mu))
    @test size(derivs.hessian_mu) == (3, 3)
end

@testset "second-order conserved susceptibilities finite" begin
    T_fm = 0.55
    muB_fm = 0.30
    muQ_fm = 0.06
    muS_fm = 0.02

    vals = (
        PNJL.chi2_B(T_fm, muB_fm; xi=0.0, p_num=32, t_num=10),
        PNJL.chi2_Q(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10),
        PNJL.chi2_S(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10),
        PNJL.chi11_BQ(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10),
        PNJL.chi11_BS(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10),
        PNJL.chi11_QS(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10),
    )

    @test all(isfinite, vals)
end

@testset "chi2_B agrees with second-order BQS mapping" begin
    T_fm = 0.60
    muB_fm = 0.24

    chi2_from_baryon = PNJL.chi2_B(T_fm, muB_fm; xi=0.0, p_num=32, t_num=10)
    chi2_from_bqs = PNJL.conserved_charge_susceptibility(T_fm, muB_fm, 0.0, 0.0; orders=(2, 0, 0), xi=0.0, p_num=32, t_num=10)

    @test isapprox(chi2_from_baryon, chi2_from_bqs; rtol=1e-6, atol=1e-8)
end

@testset "systematic chi_BQS API matches convenience wrappers" begin
    T_fm = 0.58
    muB_fm = 0.21
    muQ_fm = 0.04
    muS_fm = 0.01

    @test isapprox(
        PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 2, 0), xi=0.0, p_num=32, t_num=10),
        PNJL.chi2_Q(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10);
        rtol=1e-10,
        atol=1e-12,
    )
    @test isapprox(
        PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), xi=0.0, p_num=32, t_num=10),
        PNJL.chi11_BQ(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=32, t_num=10);
        rtol=1e-10,
        atol=1e-12,
    )
    @test isapprox(
        PNJL.chi_BQS(T_fm, muB_fm, 0.0, 0.0; orders=(2, 0, 0), xi=0.0, p_num=32, t_num=10),
        PNJL.chi2_B(T_fm, muB_fm; xi=0.0, p_num=32, t_num=10);
        rtol=1e-10,
        atol=1e-12,
    )
end

@testset "pure Q and S higher-order susceptibilities are available" begin
    T_fm = 0.57
    muB_fm = 0.18
    muQ_fm = 0.05
    muS_fm = 0.02

    vals = (
        PNJL.chi1_Q(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8),
        PNJL.chi3_Q(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8),
        PNJL.chi4_Q(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8),
        PNJL.chi1_S(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8),
        PNJL.chi3_S(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8),
        PNJL.chi4_S(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8),
    )

    @test all(isfinite, vals)
    @test isapprox(
        PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 3, 0), xi=0.0, p_num=24, t_num=8),
        PNJL.chi3_Q(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8);
        rtol=1e-10,
        atol=1e-12,
    )
    @test isapprox(
        PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 0, 4), xi=0.0, p_num=24, t_num=8),
        PNJL.chi4_S(T_fm, muB_fm, muQ_fm, muS_fm; xi=0.0, p_num=24, t_num=8);
        rtol=1e-10,
        atol=1e-12,
    )
end

@testset "cumulant_BQS wrapper matches VT^3 chi" begin
    T_fm = 0.57
    muB_fm = 0.18
    muQ_fm = 0.05
    muS_fm = 0.02
    V = 64.0

    chi_020 = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(0, 2, 0), xi=0.0, p_num=32, t_num=10)
    c_020 = PNJL.cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, V; orders=(0, 2, 0), xi=0.0, p_num=32, t_num=10)

    chi_110 = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), xi=0.0, p_num=32, t_num=10)
    c_110 = PNJL.cumulant_BQS(T_fm, muB_fm, muQ_fm, muS_fm, V; orders=(1, 1, 0), xi=0.0, p_num=32, t_num=10)

    @test isapprox(c_020, V * T_fm^3 * chi_020; rtol=1e-12, atol=1e-12)
    @test isapprox(c_110, V * T_fm^3 * chi_110; rtol=1e-12, atol=1e-12)
end

@testset "cumulant_B matches cumulant_BQS baryon path" begin
    T_fm = 0.54
    muB_fm = 0.33
    V = 50.0

    c4_b = PNJL.cumulant_B(T_fm, muB_fm, V; order=4, xi=0.0, p_num=32, t_num=10)
    c4_bqs = PNJL.cumulant_BQS(T_fm, muB_fm, 0.0, 0.0, V; orders=(4, 0, 0), xi=0.0, p_num=32, t_num=10)

    @test isapprox(c4_b, c4_bqs; rtol=1e-12, atol=1e-12)
end