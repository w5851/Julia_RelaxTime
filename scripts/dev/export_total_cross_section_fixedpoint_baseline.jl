#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# Load RelaxTime module (provides GaussLegendre, AFieldBuilder, etc. in correct nesting)
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "RelaxTime.jl"))

using .Constants_PNJL: G_fm2, K_fm5
using .RelaxTime.AFieldBuilder: ensure_quark_params_has_A
using .RelaxTime.EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .RelaxTime.TotalCrossSection: total_cross_section

const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_total_cross_section_fixedpoints_v1.csv")

function parse_args(args::Vector{String})
    output = DEFAULT_OUTPUT
    n_points = 6
    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--output"
            i += 1
            output = args[i]
        elseif arg == "--n-points"
            i += 1
            n_points = parse(Int, args[i])
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/dev/export_total_cross_section_fixedpoint_baseline.jl [--output <csv>] [--n-points <int>]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return (output=output, n_points=n_points)
end

function build_fixture()
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

function main(args::Vector{String})
    opts = parse_args(args)
    fx = build_fixture()

    points = (
        (process=:uu_to_uu, s=25.0),
        (process=:uu_to_uu, s=31.0),
        (process=:ud_to_ud, s=31.0),
        (process=:us_to_us, s=31.0),
        (process=:uubar_to_uubar, s=31.0),
    )

    mkpath(dirname(opts.output))
    open(opts.output, "w") do io
        println(io, "process,s,n_points,sigma")
        for pt in points
            σ = total_cross_section(pt.process, pt.s, fx.q, fx.t, fx.K; n_points=opts.n_points, fast_path=true)
            @printf(io, "%s,%.6f,%d,%.16e\n", String(pt.process), pt.s, opts.n_points, σ)
        end
    end

    println("baseline exported to: " * opts.output)
end

main(ARGS)
