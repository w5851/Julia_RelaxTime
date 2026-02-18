using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AFieldBuilder.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: G_fm2, K_fm5
using .AFieldBuilder: ensure_quark_params_has_A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .TotalCrossSection: total_cross_section, scan_s_dependence, calculate_all_total_cross_sections

function _build_fixture()
    quark_basic = (
        m=(u=1.52, d=1.52, s=2.50),
        μ=(u=0.30, d=0.30, s=0.00),
    )
    thermo = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)

    quark_with_A = ensure_quark_params_has_A(
        quark_basic,
        thermo;
        p_nodes=12,
        p_max=10.0,
        cos_nodes=4,
        use_aniso=false,
        warn_on_auto=false,
    )

    G_u = calculate_G_from_A(quark_with_A.A.u, quark_with_A.m.u)
    G_s = calculate_G_from_A(quark_with_A.A.s, quark_with_A.m.s)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    return (q=quark_with_A, t=thermo, K=K_coeffs)
end

@testset "TotalCrossSection model-ready fixture smoke" begin
    fx = _build_fixture()
    process = :uu_to_uu
    s0 = 31.0

    σ = total_cross_section(process, s0, fx.q, fx.t, fx.K; n_points=2)
    @test isfinite(σ)
    @test σ ≥ 0.0

    scan = scan_s_dependence([25.0, 31.0, 38.0], process, fx.q, fx.t, fx.K; n_points=2)
    @test length(scan) == 3
    @test all(isfinite, scan)

    allσ = calculate_all_total_cross_sections(s0, fx.q, fx.t, fx.K; n_points=2)
    @test length(keys(allσ)) > 0
    @test all(v -> isfinite(v) || isnan(v), values(allσ))
end
