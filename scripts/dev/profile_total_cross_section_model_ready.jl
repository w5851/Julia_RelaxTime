"""
TotalCrossSection model-ready profiling (Phase G)

目标：
- 使用可复现 fixture 自动构造 quark_params.A 与完整 K_coeffs
- 输出 TCS 关键入口的基线耗时

运行：
  julia --project=. scripts/dev/profile_total_cross_section_model_ready.jl
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AFieldBuilder.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: G_fm2, K_fm5
using .AFieldBuilder: ensure_quark_params_has_A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .TotalCrossSection: total_cross_section, scan_s_dependence, calculate_all_total_cross_sections

@inline function avg_elapsed_ms(f, n::Int)
    t = @elapsed for _ in 1:n
        f()
    end
    return t / n * 1e3
end

function build_model_ready_fixture()
    quark_basic = (
        m=(u=1.52, d=1.52, s=2.50),
        μ=(u=0.30, d=0.30, s=0.00),
    )
    thermo = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)

    quark_with_A = ensure_quark_params_has_A(
        quark_basic,
        thermo;
        p_nodes=16,
        p_max=12.0,
        cos_nodes=4,
        use_aniso=false,
        warn_on_auto=false,
    )

    G_u = calculate_G_from_A(quark_with_A.A.u, quark_with_A.m.u)
    G_s = calculate_G_from_A(quark_with_A.A.s, quark_with_A.m.s)
    K_coeffs = calculate_effective_couplings(G_fm2, K_fm5, G_u, G_s)

    return (quark_params=quark_with_A, thermo_params=thermo, K_coeffs=K_coeffs)
end

function main()
    fixture = build_model_ready_fixture()
    q = fixture.quark_params
    t = fixture.thermo_params
    K = fixture.K_coeffs

    process = :uu_to_uu
    s0 = 31.0
    s_values = [25.0, 31.0, 38.0, 45.0]

    # warm-up
    total_cross_section(process, s0, q, t, K; n_points=2)
    scan_s_dependence(s_values, process, q, t, K; n_points=2)
    calculate_all_total_cross_sections(s0, q, t, K; n_points=2)

    n = 50

    σ_single_acc = Ref(0.0)
    t_single_ms = avg_elapsed_ms(() -> begin
        σ = total_cross_section(process, s0, q, t, K; n_points=2)
        σ_single_acc[] += σ
        nothing
    end, n)

    scan_acc = Ref(0.0)
    t_scan_ms = avg_elapsed_ms(() -> begin
        arr = scan_s_dependence(s_values, process, q, t, K; n_points=2)
        scan_acc[] += sum(arr)
        nothing
    end, n)

    all_acc = Ref(0.0)
    t_all_ms = avg_elapsed_ms(() -> begin
        nt = calculate_all_total_cross_sections(s0, q, t, K; n_points=2)
        all_acc[] += sum(values(nt))
        nothing
    end, n)

    println("=== TotalCrossSection Model-Ready Baseline (ms/call) ===")
    println("total_cross_section(:uu_to_uu, s0): ", round(t_single_ms; digits=4))
    println("scan_s_dependence(4 points):       ", round(t_scan_ms; digits=4))
    println("calculate_all_total_cross_sections:", round(t_all_ms; digits=4))

    println("\n=== Fixture Summary ===")
    println("A(u,d,s) = ", q.A)
    println("K keys include: K123_plus=$(K.K123_plus), K123_minus=$(K.K123_minus), det_K_plus=$(K.det_K_plus)")
end

main()
