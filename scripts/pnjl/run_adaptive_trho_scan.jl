#!/usr/bin/env julia

"""
自适应 T-ρ 扫描：读取既有 TrhoScan CSV，检测斜率趋零的区段，自动补充 ρ 采样。
"""

using CSV

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL.TrhoScan
using .PNJL.AdaptiveRhoRefinement
using .PNJL.SeedCache: DEFAULT_SEED_PATH

struct AdaptiveCLIOptions
    source::String
    output::String
    xi::Float64
    xi_tol::Float64
    slope_tol::Float64
    min_gap::Float64
    max_points::Int
    passes::Int
    tmin::Union{Nothing, Float64}
    tmax::Union{Nothing, Float64}
    resume::Bool
    p_num::Int
    t_num::Int
    seed_path::String
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :source => joinpath("data", "outputs", "results", "pnjl", "trho_scan.csv"),
        :output => joinpath("data", "outputs", "results", "pnjl", "trho_scan.csv"),
        :xi => 0.0,
        :xi_tol => 1e-6,
        :slope_tol => 5.0,
        :min_gap => 0.003,
        :max_points => 48,
        :passes => 1,
        :tmin => nothing,
        :tmax => nothing,
        :resume => true,
        :p_num => 24,
        :t_num => 8,
        :seed_path => DEFAULT_SEED_PATH,
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
        if arg == "--source"
            opts[:source] = require_value()
        elseif arg == "--output"
            opts[:output] = require_value()
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--xi-tol"
            opts[:xi_tol] = parse(Float64, require_value())
        elseif arg == "--slope-tol"
            opts[:slope_tol] = parse(Float64, require_value())
        elseif arg == "--min-gap"
            opts[:min_gap] = parse(Float64, require_value())
        elseif arg == "--max-points"
            opts[:max_points] = parse(Int, require_value())
        elseif arg == "--passes"
            opts[:passes] = parse(Int, require_value())
        elseif arg == "--tmin"
            opts[:tmin] = parse(Float64, require_value())
        elseif arg == "--tmax"
            opts[:tmax] = parse(Float64, require_value())
        elseif arg == "--no-resume"
            opts[:resume] = false
        elseif arg == "--p-num"
            opts[:p_num] = parse(Int, require_value())
        elseif arg == "--t-num"
            opts[:t_num] = parse(Int, require_value())
        elseif arg == "--seed-path"
            opts[:seed_path] = require_value()
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return AdaptiveCLIOptions(
        String(opts[:source]),
        String(opts[:output]),
        Float64(opts[:xi]),
        Float64(opts[:xi_tol]),
        Float64(opts[:slope_tol]),
        Float64(opts[:min_gap]),
        Int(opts[:max_points]),
        Int(opts[:passes]),
        opts[:tmin],
        opts[:tmax],
        Bool(opts[:resume]),
        Int(opts[:p_num]),
        Int(opts[:t_num]),
        String(opts[:seed_path]),
    )
end

function print_usage()
    println("Usage: julia scripts/pnjl/run_adaptive_trho_scan.jl [options]\n")
    println("Options:")
    println("  --source <path>           Existing TrhoScan CSV (default data/.../trho_scan.csv)")
    println("  --output <path>           Output CSV (default同 source)")
    println("  --xi <value>              目标 ξ (default 0.0)")
    println("  --xi-tol <value>          ξ 过滤容差 (default 1e-6)")
    println("  --slope-tol <value>       触发加密的 |Δμ/Δρ| 阈值 (MeV)")
    println("  --min-gap <value>         仅在大于该间距的区段补点")
    println("  --max-points <int>        每条曲线最多新增的 ρ 数量")
    println("  --passes <int>            重复迭代次数 (默认 1)")
    println("  --tmin/--tmax <MeV>       限制温度窗口")
    println("  --no-resume               禁用 skip 逻辑，强制重算")
    println("  --p-num / --t-num <int>   透传给 TrhoScan 的积分节点")
    println("  --seed-path <path>        自定义连续种子文件")
end

function _maybe_float(row, key::Symbol)
    val = try
        row[key]
    catch
        return nothing
    end
    val isa Missing && return nothing
    return Float64(val)
end

function load_curves(path::AbstractString; xi::Float64, tol::Float64, tmin, tmax)
    isfile(path) || error("Trho CSV not found: $path")
    grouped = Dict{Float64, Vector{Tuple{Float64, Float64}}}()
    for row in CSV.File(path)
        T_val = _maybe_float(row, :T_MeV)
        T_val === nothing && continue
        if tmin !== nothing && T_val < tmin
            continue
        end
        if tmax !== nothing && T_val > tmax
            continue
        end
        xi_val = _maybe_float(row, :xi)
        xi_val === nothing && continue
        abs(xi_val - xi) <= tol || continue
        rho_val = _maybe_float(row, :rho)
        mu_val = _maybe_float(row, :mu_avg_MeV)
        if mu_val === nothing
            mu_val = _maybe_float(row, :mu_MeV)
        end
        (rho_val === nothing || mu_val === nothing) && continue
        bucket = get!(grouped, T_val) do
            Vector{Tuple{Float64, Float64}}()
        end
        push!(bucket, (rho_val, mu_val))
    end
    return grouped
end

function plan_refinement(curves::Dict{Float64, Vector{Tuple{Float64, Float64}}}, config::AdaptiveRhoConfig)
    plan = Dict{Float64, Vector{Float64}}()
    for (T, samples) in curves
        length(samples) < 2 && continue
        rho_vals = Float64[s[1] for s in samples]
        mu_vals = Float64[s[2] for s in samples]
        suggestions = suggest_refinement_points(rho_vals, mu_vals; config=config)
        if isempty(suggestions)
            continue
        end
        existing = Set(round.(rho_vals; digits=config.digits))
        filtered = Float64[]
        for val in suggestions
            rounded = round(val; digits=config.digits)
            rounded in existing && continue
            push!(filtered, rounded)
        end
        isempty(filtered) && continue
        plan[T] = filtered
    end
    return plan
end

function run_adaptive_scan(opts::AdaptiveCLIOptions)
    config = AdaptiveRhoConfig(; slope_tol=opts.slope_tol, min_gap=opts.min_gap, max_points=opts.max_points)
    source_path = opts.source
    for pass in 1:opts.passes
        curves = load_curves(source_path; xi=opts.xi, tol=opts.xi_tol, tmin=opts.tmin, tmax=opts.tmax)
        plan = plan_refinement(curves, config)
        if isempty(plan)
            println("[pass $(pass)] no refinement targets found, stopping.")
            break
        end
        println("[pass $(pass)] refining $(length(plan)) temperature slices")
        for T in sort(collect(keys(plan)))
            rho_values = plan[T]
            println("  T=$(T) -> $(length(rho_values)) new ρ samples")
            stats = TrhoScan.run_trho_scan(
                ; T_values=[T],
                  rho_values=rho_values,
                  xi_values=[opts.xi],
                  output_path=opts.output,
                  seed_path=opts.seed_path,
                  overwrite=false,
                  resume=opts.resume,
                  p_num=opts.p_num,
                  t_num=opts.t_num,
            )
            println("    stats: total=$(stats.total) success=$(stats.success) skipped=$(stats.skipped)")
        end
        source_path = opts.output
    end
end

function main()
    opts = parse_args(copy(ARGS))
    run_adaptive_scan(opts)
end

main()
