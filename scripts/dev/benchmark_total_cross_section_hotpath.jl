#!/usr/bin/env julia

using Statistics
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "AFieldBuilder.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "EffectiveCouplings.jl"))
include(joinpath(PROJECT_ROOT, "src", "relaxtime", "TotalCrossSection.jl"))

using .Constants_PNJL: G_fm2, K_fm5
using .AFieldBuilder: ensure_quark_params_has_A
using .EffectiveCouplings: calculate_G_from_A, calculate_effective_couplings
using .TotalCrossSection: total_cross_section

const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "data", "outputs", "perf", "relaxtime", "total_cross_section_hotpath_benchmark_latest.csv")

function parse_args(args::Vector{String})
    output = DEFAULT_OUTPUT
    repeats = 6
    n_points = 6
    s = 31.0
    process = :uu_to_uu

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--output"
            i += 1
            output = args[i]
        elseif arg == "--repeats"
            i += 1
            repeats = parse(Int, args[i])
        elseif arg == "--n-points"
            i += 1
            n_points = parse(Int, args[i])
        elseif arg == "--s"
            i += 1
            s = parse(Float64, args[i])
        elseif arg == "--process"
            i += 1
            process = Symbol(args[i])
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/dev/benchmark_total_cross_section_hotpath.jl [--output <csv>] [--repeats <int>] [--n-points <int>] [--s <float>] [--process <symbol>]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    repeats > 0 || error("repeats must be > 0")
    n_points > 0 || error("n_points must be > 0")
    return (output=output, repeats=repeats, n_points=n_points, s=s, process=process)
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

function timed_runs(f, repeats::Int)
    ts = Float64[]
    for _ in 1:repeats
        push!(ts, @elapsed f())
    end
    return ts
end

function write_csv(path::String, legacy::Vector{Float64}, fast::Vector{Float64}, process::Symbol, s::Float64, n_points::Int)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "mode,process,s,n_points,mean_s,median_s,min_s,max_s,repeats")
        @printf(io, "legacy,%s,%.6f,%d,%.9e,%.9e,%.9e,%.9e,%d\n",
            String(process), s, n_points, mean(legacy), median(legacy), minimum(legacy), maximum(legacy), length(legacy))
        @printf(io, "fast,%s,%.6f,%d,%.9e,%.9e,%.9e,%.9e,%d\n",
            String(process), s, n_points, mean(fast), median(fast), minimum(fast), maximum(fast), length(fast))
    end
end

function main(args::Vector{String})
    opts = parse_args(args)
    fx = build_fixture()

    σ_legacy = total_cross_section(opts.process, opts.s, fx.q, fx.t, fx.K; n_points=opts.n_points, fast_path=false)
    σ_fast = total_cross_section(opts.process, opts.s, fx.q, fx.t, fx.K; n_points=opts.n_points, fast_path=true)
    isapprox(σ_fast, σ_legacy; rtol=1e-12, atol=0.0) || error("fast_path changed numeric result")

    # warmup
    total_cross_section(opts.process, opts.s, fx.q, fx.t, fx.K; n_points=opts.n_points, fast_path=false)
    total_cross_section(opts.process, opts.s, fx.q, fx.t, fx.K; n_points=opts.n_points, fast_path=true)

    legacy_times = timed_runs(() -> total_cross_section(opts.process, opts.s, fx.q, fx.t, fx.K; n_points=opts.n_points, fast_path=false), opts.repeats)
    fast_times = timed_runs(() -> total_cross_section(opts.process, opts.s, fx.q, fx.t, fx.K; n_points=opts.n_points, fast_path=true), opts.repeats)

    write_csv(opts.output, legacy_times, fast_times, opts.process, opts.s, opts.n_points)

    speedup = median(legacy_times) / max(median(fast_times), eps(Float64))
    println("[benchmark-total-cross-section-hotpath] done")
    println("  output: " * opts.output)
    println(@sprintf("  legacy_median_s=%.6e", median(legacy_times)))
    println(@sprintf("  fast_median_s=%.6e", median(fast_times)))
    println(@sprintf("  speedup=%.3fx", speedup))
end

main(ARGS)
