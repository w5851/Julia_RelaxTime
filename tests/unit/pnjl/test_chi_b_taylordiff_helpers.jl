using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const _MODELS_PATH = normpath(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

if !isdefined(Main, :Models)
    Base.include(Main, _MODELS_PATH)
end

using .Models

const PNJL = Models.pnjl_module()
const TD = Models.PNJLChiBTaylorDiff
const JT = Models.MixedTaylorJets

@testset "TaylorDiff scalar helpers extract derivatives" begin
    x0 = 0.7
    x = TD.taylor_variable(x0, 4)
    y = x^3 + 2 * x

    @test isapprox(TD.nth_derivative_from_series(y, 0), x0^3 + 2 * x0; rtol=1e-14, atol=1e-14)
    @test isapprox(TD.nth_derivative_from_series(y, 1), 3 * x0^2 + 2; rtol=1e-14, atol=1e-14)
    @test isapprox(TD.nth_derivative_from_series(y, 2), 6 * x0; rtol=1e-14, atol=1e-14)
    @test isapprox(TD.nth_derivative_from_series(y, 3), 6.0; rtol=1e-14, atol=1e-14)
    @test isapprox(TD.nth_derivative_from_series(y, 4), 0.0; rtol=0.0, atol=1e-14)

    @test_throws ArgumentError TD.taylor_variable(x0, 0)
    @test_throws ArgumentError TD.nth_derivative_from_series(y, -1)
end

@testset "MixedTaylorJet polynomial helpers extract mixed derivatives" begin
    x0 = 0.2
    y0 = -0.1
    x = JT.jet_variable(x0, 1, 2, 3)
    y = JT.jet_variable(y0, 2, 2, 3)
    f = x^2 * y + exp(x * y)

    expected = 2.0 + (2y0 + x0 * y0^2) * exp(x0 * y0)
    @test isapprox(TD.nth_derivative_from_series(f, (2, 1)), expected; rtol=1e-14, atol=1e-14)
    @test isapprox(TD.nth_derivative_from_series(f, (0, 0)), x0^2 * y0 + exp(x0 * y0); rtol=1e-14, atol=1e-14)
end

@testset "PNJL gap Taylor series residual is controlled" begin
    model = Models.create_model(:PNJL)
    result = TD._solve_gap_series(
        model,
        0.55,
        0.30;
        order=3,
        xi=0.0,
        p_num=4,
        t_num=2,
        series_residual_tol=1e-7,
    )

    @test result.linear_solve === :refactor_each_order
    @test result.residual_norm <= 1e-7
    @test all(isfinite, result.x_state)

    pressure = TD.pressure_series_B(0.55, 0.30; order=3, xi=0.0, p_num=4, t_num=2)
    @test isfinite(TD.nth_derivative_from_series(pressure, 3))
end

@testset "chi_B backend routing and validation" begin
    T_fm = 0.55
    muB_fm = 0.30
    muQ_fm = 0.06
    muS_fm = 0.02

    kwargs = (; xi=0.0, p_num=4, t_num=2)

    td3 = PNJL.chi_B(T_fm, muB_fm; order=3, derivative_backend=:taylordiff, kwargs...)
    auto3 = PNJL.chi_B(T_fm, muB_fm; order=3, derivative_backend=:auto, kwargs...)
    @test isfinite(td3)
    @test isapprox(auto3, td3; rtol=1e-12, atol=1e-12)
    @test_throws ArgumentError PNJL.chi_B(T_fm, muB_fm; order=3, derivative_backend=:forwarddiff, kwargs...)

    for call in (
        (backend -> PNJL.chi_Q(T_fm, muB_fm, muQ_fm, muS_fm; order=2, derivative_backend=backend, kwargs...)),
        (backend -> PNJL.chi_S(T_fm, muB_fm, muQ_fm, muS_fm; order=2, derivative_backend=backend, kwargs...)),
    )
        td = call(:taylordiff)
        auto = call(:auto)
        @test isfinite(td)
        @test isapprox(auto, td; rtol=1e-12, atol=1e-12)
        @test_throws ArgumentError call(:forwarddiff)
    end

    td_q4 = PNJL.chi_Q(T_fm, muB_fm, muQ_fm, muS_fm; order=4, derivative_backend=:taylordiff, kwargs...)
    auto_q4 = PNJL.chi_Q(T_fm, muB_fm, muQ_fm, muS_fm; order=4, derivative_backend=:auto, kwargs...)
    @test isfinite(td_q4)
    @test isapprox(auto_q4, td_q4; rtol=1e-12, atol=1e-12)

    td_mixed = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), derivative_backend=:taylordiff, kwargs...)
    auto_mixed = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), derivative_backend=:auto, kwargs...)
    jet_mixed = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), derivative_backend=:mixedjet, kwargs...)
    @test isfinite(td_mixed)
    @test isapprox(auto_mixed, td_mixed; rtol=1e-12, atol=1e-12)
    @test isapprox(jet_mixed, td_mixed; rtol=1e-12, atol=1e-12)
    @test_throws ArgumentError PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), derivative_backend=:forwarddiff, kwargs...)
end

@testset "mixed BQS Taylor jet supports higher-order mixed derivatives" begin
    T_fm = 0.57
    muB_fm = 0.18
    muQ_fm = 0.05
    muS_fm = 0.02
    kwargs = (; xi=0.0, p_num=4, t_num=2)

    chi11_auto = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), derivative_backend=:auto, kwargs...)
    chi11_jet = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), derivative_backend=:mixedjet, kwargs...)
    @test isapprox(chi11_auto, chi11_jet; rtol=1e-12, atol=1e-12)
    @test_throws ArgumentError PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(1, 1, 0), derivative_backend=:forwarddiff, kwargs...)

    chi211 = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(2, 1, 1), derivative_backend=:auto, kwargs...)
    chi211_jet = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(2, 1, 1), derivative_backend=:mixedjet, kwargs...)
    @test isfinite(chi211)
    @test isapprox(chi211, chi211_jet; rtol=1e-12, atol=1e-12)
    @test_throws ArgumentError PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(2, 1, 1), derivative_backend=:forwarddiff, kwargs...)

    base = SVector(
        muB_fm / 3 + 2 * muQ_fm / 3,
        muB_fm / 3 - muQ_fm / 3,
        muB_fm / 3 - muQ_fm / 3 - muS_fm,
    )
    b_dir = SVector(1 / 3, 1 / 3, 1 / 3)
    q_dir = SVector(2 / 3, -1 / 3, -1 / 3)
    series_kwargs = (; order=3, kwargs...)
    H_B = TD.nth_derivative_from_series(TD.pressure_series_direction(T_fm, base, b_dir; series_kwargs...), 3)
    H_Q = TD.nth_derivative_from_series(TD.pressure_series_direction(T_fm, base, q_dir; series_kwargs...), 3)
    H_BQ = TD.nth_derivative_from_series(TD.pressure_series_direction(T_fm, base, b_dir + q_dir; series_kwargs...), 3)
    H_BmQ = TD.nth_derivative_from_series(TD.pressure_series_direction(T_fm, base, b_dir - q_dir; series_kwargs...), 3)
    d2B_dQ_from_univariate = (H_BQ - H_BmQ - 2H_Q) / 6
    chi210_jet = PNJL.chi_BQS(T_fm, muB_fm, muQ_fm, muS_fm; orders=(2, 1, 0), derivative_backend=:mixedjet, kwargs...)
    @test isapprox(chi210_jet / T_fm^(3 - 4), d2B_dQ_from_univariate; rtol=1e-10, atol=1e-12)
end
