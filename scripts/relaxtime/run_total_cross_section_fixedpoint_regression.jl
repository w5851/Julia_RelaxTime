#!/usr/bin/env julia

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AFieldBuilder.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: G_fm2, K_fm5
using .AFieldBuilder: ensure_quark_params_has_A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .TotalCrossSection: total_cross_section

const DEFAULT_BASELINE = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_total_cross_section_fixedpoints_v1.csv")
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "regression", "total_cross_section_fixedpoint_regression_latest.csv")

function parse_args(args::Vector{String})
    baseline = DEFAULT_BASELINE
    output = DEFAULT_OUTPUT
    rtol = 1e-9
    atol = 1e-12

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--baseline"
            i += 1
            baseline = args[i]
        elseif arg == "--output"
            i += 1
            output = args[i]
        elseif arg == "--rtol"
            i += 1
            rtol = parse(Float64, args[i])
        elseif arg == "--atol"
            i += 1
            atol = parse(Float64, args[i])
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/relaxtime/run_total_cross_section_fixedpoint_regression.jl [--baseline <csv>] [--output <csv>] [--rtol <float>] [--atol <float>]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return (baseline=baseline, output=output, rtol=rtol, atol=atol)
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

function load_baseline(path::String)
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

@inline function relerr(actual::Float64, expected::Float64, atol::Float64)
    return abs(actual - expected) / max(abs(expected), atol)
end

function main(args::Vector{String})
    opts = parse_args(args)
    points = load_baseline(opts.baseline)
    fx = build_fixture()

    mkpath(dirname(opts.output))

    breaches = String[]

    open(opts.output, "w") do io
        println(io, "# baseline=" * opts.baseline)
        println(io, @sprintf("# rtol=%.6e, atol=%.6e", opts.rtol, opts.atol))
        println(io, "process,s,n_points,sigma_fast,sigma_legacy,baseline_sigma,fast_relerr,legacy_relerr,fast_vs_legacy_relerr,pass")

        for pt in points
            σ_fast = total_cross_section(pt.process, pt.s, fx.q, fx.t, fx.K; n_points=pt.n_points, fast_path=true)
            σ_legacy = total_cross_section(pt.process, pt.s, fx.q, fx.t, fx.K; n_points=pt.n_points, fast_path=false)

            err_fast = relerr(σ_fast, pt.sigma, opts.atol)
            err_legacy = relerr(σ_legacy, pt.sigma, opts.atol)
            err_fast_vs_legacy = relerr(σ_fast, σ_legacy, opts.atol)

            pass = isapprox(σ_fast, pt.sigma; rtol=opts.rtol, atol=opts.atol) &&
                   isapprox(σ_legacy, pt.sigma; rtol=opts.rtol, atol=opts.atol) &&
                   isapprox(σ_fast, σ_legacy; rtol=opts.rtol, atol=opts.atol)

            @printf(io, "%s,%.6f,%d,%.16e,%.16e,%.16e,%.6e,%.6e,%.6e,%s\n",
                String(pt.process), pt.s, pt.n_points,
                σ_fast, σ_legacy, pt.sigma,
                err_fast, err_legacy, err_fast_vs_legacy,
                pass ? "true" : "false")

            if !pass
                push!(breaches,
                    @sprintf("process=%s s=%.6f n_points=%d fast=%.16e legacy=%.16e baseline=%.16e",
                        String(pt.process), pt.s, pt.n_points, σ_fast, σ_legacy, pt.sigma)
                )
            end
        end
    end

    println("[total-cross-section-fixedpoint-regression] output: " * opts.output)

    if !isempty(breaches)
        println("[FAIL] fixedpoint regression breaches:")
        for b in breaches
            println("  - " * b)
        end
        error("total cross section fixedpoint regression failed")
    end

    println("[PASS] all fixedpoints are within tolerance")
end

main(ARGS)
