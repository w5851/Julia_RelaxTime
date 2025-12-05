#!/usr/bin/env julia

"""
High-density T-ρ scan helper.

Adds extra TrhoScan samples near a target temperature window and appends the
results to the main CSV for CEP/Maxwell analysis.
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL.TrhoScan

struct DenseScanOptions
    output::String
    xi_values::Vector{Float64}
    tmin::Float64
    tmax::Float64
    tstep::Float64
    rho_max::Float64
    coarse_step::Float64
    medium_switch::Float64
    medium_step::Float64
    fine_switch::Float64
    fine_step::Float64
    ultra_switch::Float64
    ultra_step::Float64
    overwrite::Bool
    resume::Bool
    p_num::Int
    t_num::Int
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :output => joinpath("data", "outputs", "results", "pnjl", "trho_scan.csv"),
        :xi_values => Float64[0.0],
        :tmin => 130.0,
        :tmax => 132.0,
        :tstep => 0.2,
        :rho_max => 3.0,
        :coarse_step => 0.02,
        :medium_switch => 0.8,
        :medium_step => 0.01,
        :fine_switch => 0.3,
        :fine_step => 0.005,
        :ultra_switch => 0.15,
        :ultra_step => 0.0025,
        :overwrite => false,
        :resume => true,
        :p_num => 24,
        :t_num => 8,
    )
    i = 1
    while i <= length(args)
        arg = args[i]
        function require_value()
            i == length(args) && error("missing value for $arg")
            val = args[i + 1]
            i += 1
            return val
        end
        if arg == "--output"
            opts[:output] = require_value()
        elseif arg == "--xi"
            val = parse(Float64, require_value())
            if opts[:xi_values] == Float64[0.0]
                opts[:xi_values] = Float64[]
            end
            push!(opts[:xi_values], val)
        elseif arg == "--xi-list"
            raw = split(require_value(), ',')
            vals = Float64[
                parse(Float64, strip(v)) for v in raw if !isempty(strip(v))
            ]
            opts[:xi_values] = vals
        elseif arg == "--tmin"
            opts[:tmin] = parse(Float64, require_value())
        elseif arg == "--tmax"
            opts[:tmax] = parse(Float64, require_value())
        elseif arg == "--tstep"
            opts[:tstep] = parse(Float64, require_value())
        elseif arg == "--rho-max"
            opts[:rho_max] = parse(Float64, require_value())
        elseif arg == "--coarse-step"
            opts[:coarse_step] = parse(Float64, require_value())
        elseif arg == "--medium-switch"
            opts[:medium_switch] = parse(Float64, require_value())
        elseif arg == "--medium-step"
            opts[:medium_step] = parse(Float64, require_value())
        elseif arg == "--fine-switch"
            opts[:fine_switch] = parse(Float64, require_value())
        elseif arg == "--fine-step"
            opts[:fine_step] = parse(Float64, require_value())
        elseif arg == "--ultra-switch"
            opts[:ultra_switch] = parse(Float64, require_value())
        elseif arg == "--ultra-step"
            opts[:ultra_step] = parse(Float64, require_value())
        elseif arg == "--overwrite"
            opts[:overwrite] = true
        elseif arg == "--no-resume"
            opts[:resume] = false
        elseif arg == "--p-num"
            opts[:p_num] = parse(Int, require_value())
        elseif arg == "--t-num"
            opts[:t_num] = parse(Int, require_value())
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    xi_vals = opts[:xi_values]
    if isempty(xi_vals)
        xi_vals = Float64[0.0]
    end
    xi_vals = unique(sort(xi_vals))
    opts[:tstep] > 0 || error("tstep must be positive")
    opts[:rho_max] > 0 || error("rho_max must be positive")
    opts[:coarse_step] > 0 || error("coarse_step must be positive")
    opts[:medium_step] > 0 || error("medium_step must be positive")
    opts[:fine_step] > 0 || error("fine_step must be positive")
    opts[:ultra_step] > 0 || error("ultra_step must be positive")
    return DenseScanOptions(
        String(opts[:output]),
        Float64.(xi_vals),
        Float64(opts[:tmin]),
        Float64(opts[:tmax]),
        Float64(opts[:tstep]),
        Float64(opts[:rho_max]),
        Float64(opts[:coarse_step]),
        Float64(opts[:medium_switch]),
        Float64(opts[:medium_step]),
        Float64(opts[:fine_switch]),
        Float64(opts[:fine_step]),
        Float64(opts[:ultra_switch]),
        Float64(opts[:ultra_step]),
        Bool(opts[:overwrite]),
        Bool(opts[:resume]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
    )
end

function print_usage()
    println("Usage: julia scripts/pnjl/run_dense_trho_scan.jl [options]\n")
    println("Options:")
    println("  --output <path>             Destination CSV (default data/outputs/results/pnjl/trho_scan.csv)")
    println("  --xi <value>                Append an additional ξ value (default 0.0)")
    println("  --xi-list v1,v2,...         Replace ξ list with comma-separated values")
    println("  --tmin/--tmax <MeV>         Temperature window (default 130–132 MeV)")
    println("  --tstep <MeV>               Temperature increment (default 0.2 MeV)")
    println("  --rho-max <value>           Maximum normalized density (default 3.0)")
    println("  --coarse-step <value>       Base Δρ for the full range (default 0.02)")
    println("  --medium-switch <value>     ρ where medium step kicks in (default 0.8)")
    println("  --medium-step <value>       Medium region Δρ (default 0.01)")
    println("  --fine-switch <value>       ρ where fine step kicks in (default 0.3)")
    println("  --fine-step <value>         Fine region Δρ (default 0.005)")
    println("  --ultra-switch <value>      ρ where ultra-fine step kicks in (default 0.15)")
    println("  --ultra-step <value>        Ultra-fine Δρ (default 0.0025)")
    println("  --overwrite                 Replace destination CSV instead of appending")
    println("  --no-resume                 Disable resume/skip logic")
    println("  --p-num <int>               Momentum grid size (pass-through)")
    println("  --t-num <int>               Theta grid size (pass-through)")
    println("  -h, --help                  Show this help text")
end

function build_dense_rho_values(opts::DenseScanOptions)
    return TrhoScan.build_default_rho_grid(
        rho_max=opts.rho_max,
        coarse_step=opts.coarse_step,
        medium_switch=opts.medium_switch,
        medium_step=opts.medium_step,
        fine_switch=opts.fine_switch,
        fine_step=opts.fine_step,
        ultra_fine_switch=opts.ultra_switch,
        ultra_fine_step=opts.ultra_step,
    )
end

function run_dense_scan(opts::DenseScanOptions)
    T_values = collect(range(opts.tmin; stop=opts.tmax, step=opts.tstep))
    isempty(T_values) && error("temperature range is empty")
    rho_values = build_dense_rho_values(opts)
    println("Running dense scan for T in [$(opts.tmin), $(opts.tmax)] with $(length(T_values)) slices")
    println("ρ grid points: $(length(rho_values)) (Δρ adjusted for low-density region)")
    stats = TrhoScan.run_trho_scan(
        ; T_values=T_values,
          rho_values=rho_values,
          xi_values=opts.xi_values,
          output_path=opts.output,
          overwrite=opts.overwrite,
          resume=opts.resume,
          p_num=opts.p_num,
          t_num=opts.t_num,
    )
    println("Dense scan finished: total=$(stats.total) success=$(stats.success) failure=$(stats.failure) skipped=$(stats.skipped)")
    println("Output written to $(stats.output)")
end

function main()
    opts = parse_args(copy(ARGS))
    run_dense_scan(opts)
end

main()
