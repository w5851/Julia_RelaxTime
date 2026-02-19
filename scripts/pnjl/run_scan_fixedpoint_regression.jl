#!/usr/bin/env julia

"""
PNJL 扫描固定点回归：读取 baseline，对 T-μ / T-ρ 各 1 组固定点执行回归校验。

示例：
  julia --project=. scripts/pnjl/run_scan_fixedpoint_regression.jl
  julia --project=. scripts/pnjl/run_scan_fixedpoint_regression.jl --rtol 8e-2 --atol 1e-6
"""

using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_BASELINE = joinpath(PROJECT_ROOT, "tests", "baselines", "pnjl", "baseline_pnjl_scan_fixedpoints_v1.csv")
const DEFAULT_OUTPUT = joinpath(PROJECT_ROOT, "data", "outputs", "results", "pnjl", "regression", "scan_fixedpoint_regression_latest.csv")

include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))
using .PNJL: run_tmu_scan, run_trho_scan

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
            println("Usage: julia --project=. scripts/pnjl/run_scan_fixedpoint_regression.jl [--baseline <csv>] [--output <csv>] [--rtol <x>] [--atol <x>] [--no-overwrite]")
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end

    return RegressionOptions(baseline, output, rtol, atol, overwrite)
end

function ensure_output_path(path::String, overwrite::Bool)
    mkpath(dirname(path))
    if isfile(path) && !overwrite
        error("output already exists: $path")
    end
end

function load_baseline(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    lines = readlines(path)
    isempty(lines) && error("baseline CSV is empty: $path")

    rows = NamedTuple[]
    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 10 || error("invalid baseline row: $line")
        push!(rows, (
            kind=strip(cols[1]),
            T=parse(Float64, strip(cols[2])),
            mu=parse(Float64, strip(cols[3])),
            rho=parse(Float64, strip(cols[4])),
            xi=parse(Float64, strip(cols[5])),
            pressure=parse(Float64, strip(cols[6])),
            rho_out=parse(Float64, strip(cols[7])),
            mu_avg=parse(Float64, strip(cols[8])),
            entropy=parse(Float64, strip(cols[9])),
            energy=parse(Float64, strip(cols[10])),
        ))
    end
    isempty(rows) && error("baseline has no data rows")
    return rows
end

function _parse_single_row(path::String)
    lines = readlines(path)
    length(lines) >= 2 || error("invalid csv output: $path")
    header = split(lines[1], ',')
    cols = split(lines[2], ',')
    return Dict(header[i] => cols[i] for i in eachindex(header))
end

function _f(row::AbstractDict, key::String)
    haskey(row, key) || error("missing column '$key'")
    return parse(Float64, row[key])
end

function compare_metric!(breaches::Vector{String}, label::String, metric::String, actual::Float64, expected::Float64, rtol::Float64, atol::Float64)
    if !isapprox(actual, expected; rtol=rtol, atol=atol)
        push!(breaches,
            @sprintf("%s | %s | actual=%.16e expected=%.16e abs_diff=%.3e",
                label, metric, actual, expected, abs(actual - expected)))
    end
end

function write_header(io)
    println(io, "kind,T,mu,rho,xi,pressure_actual,pressure_baseline,entropy_actual,entropy_baseline,energy_actual,energy_baseline,rho_out_actual,rho_out_baseline,mu_avg_actual,mu_avg_baseline,pass")
end

function write_row(io, row, actual; pass::Bool)
    @printf(io,
        "%s,%.6f,%.6f,%.6f,%.6f,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%s\n",
        row.kind, row.T, row.mu, row.rho, row.xi,
        actual.pressure, row.pressure,
        actual.entropy, row.entropy,
        actual.energy, row.energy,
        actual.rho_out, row.rho_out,
        actual.mu_avg, row.mu_avg,
        pass ? "true" : "false"
    )
end

function run_row(row)
    tmp_dir = mktempdir()
    out = joinpath(tmp_dir, "single.csv")

    if row.kind == "tmu"
        run_tmu_scan(
            T_values=[row.T],
            mu_values=[row.mu],
            xi_values=[row.xi],
            output_path=out,
            overwrite=true,
            resume=false,
            use_phase_aware=false,
            thermo_backend=:legacy,
            p_num=8,
            t_num=4,
            iterations=80,
        )
        r = _parse_single_row(out)
        return (pressure=_f(r, "pressure_fm4"), entropy=_f(r, "entropy_fm3"), energy=_f(r, "energy_fm4"), rho_out=_f(r, "rho"), mu_avg=NaN)
    elseif row.kind == "trho"
        run_trho_scan(
            T_values=[row.T],
            rho_values=[row.rho],
            xi_values=[row.xi],
            output_path=out,
            overwrite=true,
            resume=false,
            reverse_rho=false,
            thermo_backend=:legacy,
            p_num=8,
            t_num=4,
            iterations=80,
        )
        r = _parse_single_row(out)
        return (pressure=_f(r, "pressure_fm4"), entropy=_f(r, "entropy_fm3"), energy=_f(r, "energy_fm4"), rho_out=NaN, mu_avg=_f(r, "mu_avg_MeV"))
    else
        error("unsupported kind: $(row.kind)")
    end
end

function main(args::Vector{String})
    opts = parse_args(args)
    ensure_output_path(opts.output, opts.overwrite)
    baseline = load_baseline(opts.baseline)
    breaches = String[]

    open(opts.output, "w") do io
        println(io, "# baseline=" * opts.baseline)
        println(io, @sprintf("# rtol=%.6e, atol=%.6e", opts.rtol, opts.atol))
        write_header(io)

        for row in baseline
            actual = run_row(row)
            label = @sprintf("%s(T=%.3f,mu=%.3f,rho=%.3f,xi=%.3f)", row.kind, row.T, row.mu, row.rho, row.xi)

            compare_metric!(breaches, label, "pressure", actual.pressure, row.pressure, opts.rtol, opts.atol)
            compare_metric!(breaches, label, "entropy", actual.entropy, row.entropy, opts.rtol, opts.atol)
            compare_metric!(breaches, label, "energy", actual.energy, row.energy, opts.rtol, opts.atol)
            if row.kind == "tmu" && isfinite(row.rho_out)
                compare_metric!(breaches, label, "rho_out", actual.rho_out, row.rho_out, opts.rtol, opts.atol)
            elseif row.kind == "trho" && isfinite(row.mu_avg)
                compare_metric!(breaches, label, "mu_avg", actual.mu_avg, row.mu_avg, opts.rtol, opts.atol)
            end

            pass = true
            if row.kind == "tmu"
                pass &= isapprox(actual.rho_out, row.rho_out; rtol=opts.rtol, atol=opts.atol)
            else
                pass &= isapprox(actual.mu_avg, row.mu_avg; rtol=opts.rtol, atol=opts.atol)
            end
            pass &= isapprox(actual.pressure, row.pressure; rtol=opts.rtol, atol=opts.atol)
            pass &= isapprox(actual.entropy, row.entropy; rtol=opts.rtol, atol=opts.atol)
            pass &= isapprox(actual.energy, row.energy; rtol=opts.rtol, atol=opts.atol)

            write_row(io, row, actual; pass=pass)
        end
    end

    println("pnjl scan fixedpoint regression finished")
    println("baseline: " * opts.baseline)
    println("output:   " * opts.output)

    if !isempty(breaches)
        println("[FAIL] baseline threshold breaches detected:")
        for msg in breaches
            println("  - " * msg)
        end
        error(@sprintf("pnjl scan fixedpoint regression failed: %d breaches", length(breaches)))
    end

    println("[PASS] all fixedpoints are within tolerance")
end

main(ARGS)
