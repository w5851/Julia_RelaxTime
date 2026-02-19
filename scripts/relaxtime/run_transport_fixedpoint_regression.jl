#!/usr/bin/env julia

"""
固定点回归：读取 transport baseline，计算 legacy/models 两套实现并逐点比对。

输出：
- 默认写入 `data/outputs/results/relaxtime/regression/transport_fixedpoint_regression_latest.csv`
- 若任一点超出阈值，脚本以非零退出并打印可定位信息

示例：
  julia --project=. scripts/relaxtime/run_transport_fixedpoint_regression.jl
  julia --project=. scripts/relaxtime/run_transport_fixedpoint_regression.jl --rtol 8e-2 --atol 1e-6
"""

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_BASELINE = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_transport_fixedpoints_v1.csv")
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "data", "outputs", "results", "relaxtime", "regression", "transport_fixedpoint_regression_latest.csv")

include(joinpath(PROJECT_ROOT, "src", "pnjl", "workflows", "TransportWorkflow.jl"))
using .TransportWorkflow

struct RegressionOptions
    baseline::String
    output::String
    rtol::Float64
    atol::Float64
    overwrite::Bool
end

function parse_args(args::Vector{String})
    baseline = DEFAULT_BASELINE
    output = DEFAULT_OUTPUT
    rtol = 8e-2
    atol = 1e-6
    overwrite = true

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg == "--baseline"
            i == length(args) && error("missing value for --baseline")
            i += 1
            baseline = args[i]
        elseif arg == "--output"
            i == length(args) && error("missing value for --output")
            i += 1
            output = args[i]
        elseif arg == "--rtol"
            i == length(args) && error("missing value for --rtol")
            i += 1
            rtol = parse(Float64, args[i])
        elseif arg == "--atol"
            i == length(args) && error("missing value for --atol")
            i += 1
            atol = parse(Float64, args[i])
        elseif arg == "--no-overwrite"
            overwrite = false
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. scripts/relaxtime/run_transport_fixedpoint_regression.jl [--baseline <csv>] [--output <csv>] [--rtol <x>] [--atol <x>] [--no-overwrite]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    rtol >= 0 || error("rtol must be >= 0")
    atol >= 0 || error("atol must be >= 0")

    return RegressionOptions(baseline, output, rtol, atol, overwrite)
end

function _point_key(T::Float64, mu::Float64, xi::Float64)
    return @sprintf("%.6f|%.6f|%.6f", T, mu, xi)
end

function load_baseline(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("baseline CSV is empty: $path")

    points = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue

        cols = split(s, ',')
        length(cols) == 6 || error("invalid baseline row: $line")

        T = parse(Float64, strip(cols[1]))
        mu = parse(Float64, strip(cols[2]))
        xi = parse(Float64, strip(cols[3]))
        eta = parse(Float64, strip(cols[4]))
        sigma = parse(Float64, strip(cols[5]))
        zeta = parse(Float64, strip(cols[6]))

        push!(points, (T=T, mu=mu, xi=xi, eta=eta, sigma=sigma, zeta=zeta))
    end

    isempty(points) && error("baseline has no data rows: $path")
    return points
end

function rel_error(actual::Float64, expected::Float64, atol::Float64)
    return abs(actual - expected) / max(abs(expected), atol)
end

function format_point(pt)
    return @sprintf("T=%.6f, mu=%.6f, xi=%.6f", pt.T, pt.mu, pt.xi)
end

function run_one(T::Float64, mu::Float64, xi::Float64, tau, models_solver)
    common_kwargs = (
        xi=xi,
        thermo_backend=:models,
        tau=tau,
        compute_tau=false,
        compute_bulk=true,
        p_num=8,
        t_num=4,
        transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
    )

    res_legacy = TransportWorkflow.solve_gap_and_transport(
        T,
        mu;
        common_kwargs...,
        solver_backend=:legacy,
        solver_kwargs=(iterations=30,),
    )

    res_models = TransportWorkflow.solve_gap_and_transport(
        T,
        mu;
        common_kwargs...,
        solver_backend=:models,
        models_solver=models_solver,
        models_residual_norm_max=1e-4,
        seed_state=collect(res_legacy.equilibrium.x_state),
    )

    return res_legacy, res_models
end

function compare_metric!(breaches::Vector{String}, pt, backend::Symbol, metric::Symbol, actual::Float64, expected::Float64, rtol::Float64, atol::Float64)
    ok = isapprox(actual, expected; rtol=rtol, atol=atol)
    if !ok
        push!(breaches,
            @sprintf(
                "%s | backend=%s | metric=%s | actual=%.16e expected=%.16e abs_diff=%.3e rel_err=%.3e",
                format_point(pt), String(backend), String(metric), actual, expected, abs(actual - expected), rel_error(actual, expected, atol)
            )
        )
    end
    return ok
end

function ensure_output_path(path::String, overwrite::Bool)
    mkpath(dirname(path))
    if isfile(path) && !overwrite
        error("output already exists: $path (use default overwrite behavior or provide a different --output)")
    end
end

function write_header(io)
    println(io, "T,mu,xi,backend,eta,sigma,zeta,baseline_eta,baseline_sigma,baseline_zeta,eta_relerr,sigma_relerr,zeta_relerr,pass")
end

function write_row(io, pt, backend::Symbol, eta::Float64, sigma::Float64, zeta::Float64, atol::Float64, rtol::Float64)
    eta_rel = rel_error(eta, pt.eta, atol)
    sigma_rel = rel_error(sigma, pt.sigma, atol)
    zeta_rel = rel_error(zeta, pt.zeta, atol)
    ok = isapprox(eta, pt.eta; rtol=rtol, atol=atol) &&
         isapprox(sigma, pt.sigma; rtol=rtol, atol=atol) &&
         isapprox(zeta, pt.zeta; rtol=rtol, atol=atol)

    @printf(io,
        "%.6f,%.6f,%.6f,%s,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.6e,%.6e,%.6e,%s\n",
        pt.T, pt.mu, pt.xi, String(backend),
        eta, sigma, zeta,
        pt.eta, pt.sigma, pt.zeta,
        eta_rel, sigma_rel, zeta_rel,
        ok ? "true" : "false"
    )
end

function main(args::Vector{String})
    opts = parse_args(args)
    points = load_baseline(opts.baseline)

    ensure_output_path(opts.output, opts.overwrite)

    tau = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)
    models_solver = Main.Models.NLsolveGapSolver(
        method=:trust_region,
        jacobian=:finite,
        xtol=1e-10,
        ftol=1e-10,
    )

    breaches = String[]

    open(opts.output, "w") do io
        println(io, "# baseline=" * opts.baseline)
        println(io, @sprintf("# rtol=%.6e, atol=%.6e", opts.rtol, opts.atol))
        write_header(io)

        for pt in points
            res_legacy, res_models = run_one(pt.T, pt.mu, pt.xi, tau, models_solver)

            legacy = res_legacy.transport
            models = res_models.transport

            compare_metric!(breaches, pt, :legacy, :eta, legacy.eta, pt.eta, opts.rtol, opts.atol)
            compare_metric!(breaches, pt, :legacy, :sigma, legacy.sigma, pt.sigma, opts.rtol, opts.atol)
            compare_metric!(breaches, pt, :legacy, :zeta, legacy.zeta, pt.zeta, opts.rtol, opts.atol)

            compare_metric!(breaches, pt, :models, :eta, models.eta, pt.eta, opts.rtol, opts.atol)
            compare_metric!(breaches, pt, :models, :sigma, models.sigma, pt.sigma, opts.rtol, opts.atol)
            compare_metric!(breaches, pt, :models, :zeta, models.zeta, pt.zeta, opts.rtol, opts.atol)

            write_row(io, pt, :legacy, legacy.eta, legacy.sigma, legacy.zeta, opts.atol, opts.rtol)
            write_row(io, pt, :models, models.eta, models.sigma, models.zeta, opts.atol, opts.rtol)
        end
    end

    println("transport fixedpoint regression finished")
    println("baseline: " * opts.baseline)
    println("output:   " * opts.output)
    println(@sprintf("points: %d", length(points)))
    println(@sprintf("checks: %d", 6 * length(points)))

    if !isempty(breaches)
        println()
        println("[FAIL] Baseline threshold breaches detected:")
        for msg in breaches
            println("  - " * msg)
        end
        error(@sprintf("transport fixedpoint regression failed: %d breaches", length(breaches)))
    end

    println("[PASS] All points are within tolerance.")
end

main(ARGS)
