using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const TOTAL_CROSS_SECTION_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_total_cross_section_fixedpoints_v1.csv")

if !isdefined(Main, :RelaxTime)
    include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxTime.jl"))
end

using .Constants_PNJL: G_fm2, K_fm5
using .RelaxTime.AFieldBuilder: ensure_quark_params_has_A
using .RelaxTime.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .RelaxTime.TotalCrossSection: total_cross_section

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

function _load_baseline(path::String)
    isfile(path) || error("baseline file not found: $path")
    lines = readlines(path)
    length(lines) >= 2 || error("baseline is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue

        cols = split(s, ',')
        length(cols) == 4 || error("invalid row in baseline: $line")
        process = Symbol(strip(cols[1]))
        sval = parse(Float64, strip(cols[2]))
        n_points = parse(Int, strip(cols[3]))
        sigma = parse(Float64, strip(cols[4]))
        push!(rows, (process=process, s=sval, n_points=n_points, sigma=sigma))
    end

    isempty(rows) && error("baseline has no data rows")
    return rows
end

@testset "Total cross section fixedpoint regression" begin
    points = _load_baseline(TOTAL_CROSS_SECTION_BASELINE_PATH)
    fixture = _build_fixture()
    rtol = 1e-9
    atol = 1e-12

    for pt in points
        @testset "process=$(pt.process) s=$(pt.s) n=$(pt.n_points)" begin
            sigma_fast = total_cross_section(pt.process, pt.s, fixture.q, fixture.t, fixture.K; n_points=pt.n_points, fast_path=true)
            sigma_legacy = total_cross_section(pt.process, pt.s, fixture.q, fixture.t, fixture.K; n_points=pt.n_points, fast_path=false)

            @test isapprox(sigma_fast, pt.sigma; rtol=rtol, atol=atol)
            @test isapprox(sigma_legacy, pt.sigma; rtol=rtol, atol=atol)
            @test isapprox(sigma_fast, sigma_legacy; rtol=rtol, atol=atol)
        end
    end
end