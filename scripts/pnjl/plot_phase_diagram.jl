#!/usr/bin/env julia

"""
PNJL 相图绘制脚本（步骤 7）。

读取 `TrhoScan` CSV，使用步骤 4-6 的分析模块（S 形检测、CEP 搜索、Maxwell 等面积）
推导共存线与自旋odal，并输出 `T-μ` 相图。
"""

using CSV
using Plots

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL.PhaseTransition: detect_s_shape
using .PNJL.CEPFinder: build_curves, find_cep, CEPResult
using .PNJL.MaxwellRhoMu: build_phase_boundary, MaxwellResult
using .PNJL.TrhoScan

struct CLIOptions
    trho_csv::String
    xi::Float64
    xi_tol::Float64
    output::String
    min_samples::Int
    area_tol::Float64
    candidate_steps::Int
    max_iter::Int
    tmin::Union{Nothing, Float64}
    tmax::Union{Nothing, Float64}
    dump_csv::Union{Nothing, String}
    processed_dir::String
    processed_fig_dir::String
    force_update_cep::Bool
end

function parse_args(args::Vector{String})
    opts = Dict{Symbol, Any}(
        :trho_csv => joinpath("data", "outputs", "results", "pnjl", "trho_scan.csv"),
        :xi => 0.0,
        :xi_tol => 1e-6,
        :output => joinpath("data", "outputs", "figures", "pnjl", "pnjl_phase_diagram.png"),
        :min_samples => 12,
        :area_tol => 1e-4,
        :candidate_steps => 64,
        :max_iter => 60,
        :tmin => nothing,
        :tmax => nothing,
        :dump_csv => nothing,
        :processed_dir => joinpath("data", "processed", "results", "pnjl"),
        :processed_fig_dir => joinpath("data", "processed", "figures", "pnjl"),
        :force_update_cep => false,
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
        if arg == "--trho-csv"
            opts[:trho_csv] = require_value()
        elseif arg == "--xi"
            opts[:xi] = parse(Float64, require_value())
        elseif arg == "--xi-tol"
            opts[:xi_tol] = parse(Float64, require_value())
        elseif arg == "--output"
            opts[:output] = require_value()
        elseif arg == "--min-samples"
            opts[:min_samples] = parse(Int, require_value())
        elseif arg == "--area-tol"
            opts[:area_tol] = parse(Float64, require_value())
        elseif arg == "--candidate-steps"
            opts[:candidate_steps] = parse(Int, require_value())
        elseif arg == "--max-iter"
            opts[:max_iter] = parse(Int, require_value())
        elseif arg == "--tmin"
            opts[:tmin] = parse(Float64, require_value())
        elseif arg == "--tmax"
            opts[:tmax] = parse(Float64, require_value())
        elseif arg == "--dump-boundary"
            opts[:dump_csv] = require_value()
        elseif arg == "--processed-dir"
            opts[:processed_dir] = require_value()
        elseif arg == "--processed-figures-dir"
            opts[:processed_fig_dir] = require_value()
        elseif arg == "--force-update-cep"
            opts[:force_update_cep] = true
        elseif arg in ("-h", "--help")
            print_usage()
            exit(0)
        else
            error("unknown option: $arg")
        end
        i += 1
    end
    return CLIOptions(
        String(opts[:trho_csv]),
        Float64(opts[:xi]),
        Float64(opts[:xi_tol]),
        String(opts[:output]),
        Int(opts[:min_samples]),
        Float64(opts[:area_tol]),
        Int(opts[:candidate_steps]),
        Int(opts[:max_iter]),
        opts[:tmin],
        opts[:tmax],
        opts[:dump_csv],
        String(opts[:processed_dir]),
        String(opts[:processed_fig_dir]),
        Bool(opts[:force_update_cep]),
    )
end

function print_usage()
    println("Usage: julia scripts/pnjl/plot_phase_diagram.jl [options]\n")
    println("Options:")
    println("  --trho-csv <path>         Path to Trho scan CSV (default data/outputs/results/pnjl/trho_scan.csv)")
    println("  --xi <value>              ξ filter (default 0.0)")
    println("  --xi-tol <value>          Allowed deviation for ξ when filtering (default 1e-6)")
    println("  --output <path>           Output image path (default data/outputs/figures/pnjl/pnjl_phase_diagram.png)")
    println("  --min-samples <int>       Minimum samples per curve for Maxwell (default 12)")
    println("  --area-tol <float>        Area tolerance for Maxwell equal-area (default 1e-4)")
    println("  --candidate-steps <int>   Initial μ scan steps (default 64)")
    println("  --max-iter <int>          Bisection iterations (default 60)")
    println("  --tmin/--tmax <float>     Optional temperature window (MeV)")
    println("  --dump-boundary <path>    Optional CSV dump for coexistence points")
    println("  --processed-dir <path>    Directory for intermediate CSV outputs (default data/processed/results/pnjl)")
    println("  --processed-figures-dir <path>  Directory for intermediate figures (default data/processed/figures/pnjl)")
    println("  --force-update-cep        Recompute CEP even if cached for the same ξ")
    println("  -h, --help                Show this help message")
end

function load_trho_curves(path::AbstractString; xi::Float64=0.0, tol::Float64=1e-6)
    isfile(path) || error("Trho CSV not found: $path")
    grouped = Dict{Float64, Vector{Tuple{Float64, Float64}}}()
    for row in CSV.File(path)
        T_val = try
            Float64(row[:T_MeV])
        catch
            continue
        end
        xi_val = _maybe_float(row, :xi)
        isnothing(xi_val) && continue
        abs(xi_val - xi) <= tol || continue
        if !_is_converged(row)
            continue
        end
        rho_val = _maybe_float(row, :rho)
        mu_val = _maybe_float(row, :mu_avg_MeV)
        if mu_val === nothing
            mu_val = _maybe_float(row, :mu_MeV)
        end
        if rho_val === nothing || mu_val === nothing
            continue
        end
        bucket = get!(grouped, T_val) do
            Vector{Tuple{Float64, Float64}}()
        end
        push!(bucket, (mu_val, rho_val))
    end
    return grouped
end

_maybe_float(row, key::Symbol) = try
    val = row[key]
    val isa Missing && return nothing
    return Float64(val)
catch
    return nothing
end

function _is_converged(row)
    val = try row[:converged] catch; missing end
    if val isa Missing
        return true
    elseif val isa Bool
        return val
    elseif val isa AbstractString
        str = lowercase(strip(val))
        return str == "true"
    else
        return Bool(val)
    end
end

function filter_curves(curves::Dict{Float64, T}, tmin, tmax) where {T}
    if tmin === nothing && tmax === nothing
        return curves
    end
    filtered = Dict{Float64, T}()
    for (Tval, data) in curves
        if (tmin === nothing || Tval >= tmin) && (tmax === nothing || Tval <= tmax)
            filtered[Tval] = data
        end
    end
    return filtered
end

function extract_phase_boundary(results::Dict{Float64, MaxwellResult})
    T_vals = Float64[]
    mu_vals = Float64[]
    rho_gas = Float64[]
    rho_liquid = Float64[]
    for T in sort(collect(keys(results)))
        res = results[T]
        if res.converged && res.mu_coex_MeV !== nothing && res.rho_gas !== nothing && res.rho_liquid !== nothing
            push!(T_vals, T)
            push!(mu_vals, res.mu_coex_MeV)
            push!(rho_gas, res.rho_gas)
            push!(rho_liquid, res.rho_liquid)
        end
    end
    return (T=T_vals, mu=mu_vals, rho_g=rho_gas, rho_l=rho_liquid)
end

function extract_spinodals(curves)
    T_vals = Float64[]
    mu_low = Float64[]
    mu_high = Float64[]
    for T in sort(collect(keys(curves)))
        mu_arr, rho_arr = curves[T]
        result = detect_s_shape(mu_arr, rho_arr)
        if result.has_s_shape && result.mu_spinodal_low !== nothing && result.mu_spinodal_high !== nothing
            push!(T_vals, T)
            push!(mu_low, result.mu_spinodal_low)
            push!(mu_high, result.mu_spinodal_high)
        end
    end
    return (T=T_vals, low=mu_low, high=mu_high)
end

function dump_boundary_csv(path::AbstractString, boundary)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "T_MeV,mu_coex_MeV,rho_gas,rho_liquid")
        for i in eachindex(boundary.T)
            println(io, join((boundary.T[i], boundary.mu[i], boundary.rho_g[i], boundary.rho_l[i]), ','))
        end
    end
end

function plot_phase_diagram(opts::CLIOptions)
    grouped = load_trho_curves(opts.trho_csv; xi=opts.xi, tol=opts.xi_tol)
    isempty(grouped) && error("no curves found for ξ=$(opts.xi)")
    curves = build_curves(grouped)
    curves = filter_curves(curves, opts.tmin, opts.tmax)
    isempty(curves) && error("no curves remain after temperature filter")

    curve_fetcher = T -> _fetch_curve_via_trho!(curves, opts, T)
    cep = _get_or_compute_cep(opts.processed_dir, curves;
        xi=opts.xi, xi_tol=opts.xi_tol, temp_tol=0.01,
        force_update=opts.force_update_cep, curve_fetcher=curve_fetcher)

    curves = filter_curves(curves, opts.tmin, opts.tmax)
    isempty(curves) && error("no curves remain after CEP refinement")

    boundary_results = build_phase_boundary(curves;
        min_samples=opts.min_samples,
        tol_area=opts.area_tol,
        candidate_steps=opts.candidate_steps,
        max_iter=opts.max_iter,
    )
    boundary = extract_phase_boundary(boundary_results)
    spinodals = extract_spinodals(curves)
    dump_processed_data(opts.processed_dir, curves, boundary_results, boundary, spinodals, cep; xi=opts.xi)
    plot_processed_figures(opts.processed_fig_dir, curves, boundary, spinodals, boundary_results, cep)

    if opts.dump_csv !== nothing
        dump_boundary_csv(opts.dump_csv, boundary)
        @info "Saved coexistence table" opts.dump_csv
    end

    mkpath(dirname(opts.output))
    default(; legend=:topright, framestyle=:box, size=(900, 600))
    plt = plot(; xlabel="μ (MeV)", ylabel="T (MeV)", title="PNJL Phase Diagram", grid=true)
    if !isempty(spinodals.T)
        plot!(plt, spinodals.low, spinodals.T; label="Lower spinodal", linestyle=:dash, color=:gray)
        plot!(plt, spinodals.high, spinodals.T; label="Upper spinodal", linestyle=:dash, color=:gray)
    end
    if !isempty(boundary.T)
        plot!(plt, boundary.mu, boundary.T; label="Coexistence line", color=:red, linewidth=2)
    end
    if cep.has_cep && cep.T_cep_MeV !== nothing && cep.mu_cep_MeV !== nothing
        scatter!(plt, [cep.mu_cep_MeV], [cep.T_cep_MeV]; label="CEP", color=:blue, marker=:star5, markersize=10)
    end
    savefig(plt, opts.output)
    println("Saved phase diagram to $(opts.output)")
end

function dump_processed_data(dir::AbstractString, curves, boundary_results, boundary, spinodals, cep; xi::Float64=0.0)
    mkpath(dir)
    _dump_curves(joinpath(dir, "curves.csv"), curves, xi)
    _dump_maxwell(joinpath(dir, "maxwell_results.csv"), boundary_results, xi)
    _dump_maxwell_inputs(joinpath(dir, "maxwell_inputs.csv"), curves, xi)
    _dump_spinodals(joinpath(dir, "spinodals.csv"), spinodals, xi)
end

function plot_processed_figures(dir::AbstractString, curves, boundary, spinodals, boundary_results, cep)
    mkpath(dir)
    _plot_curves(joinpath(dir, "curves.png"), curves)
    _plot_spinodals(joinpath(dir, "spinodals.png"), spinodals)
    _plot_maxwell(joinpath(dir, "maxwell.png"), boundary, boundary_results, cep)
end

function _plot_curves(path, curves)
    isempty(curves) && return
    plt = plot(; xlabel="ρ/ρ₀", ylabel="μ (MeV)", title="PNJL Curves", legend=:outerright, size=(900, 600))
    colors = get_color_palette(:auto, length(curves))
    for (idx, T) in enumerate(sort(collect(keys(curves))))
        mu_vec, rho_vec = curves[T]
        color = colors[mod1(idx, length(colors))]
        plot!(plt, rho_vec, mu_vec; label="T=$(T) MeV", color=color)
    end
    savefig(plt, path)
end

function _plot_spinodals(path, spinodals)
    isempty(spinodals.T) && return
    plt = plot(; xlabel="μ (MeV)", ylabel="T (MeV)", title="Spinodals", legend=:bottomright)
    plot!(plt, spinodals.low, spinodals.T; label="Lower", linestyle=:dash, color=:gray)
    plot!(plt, spinodals.high, spinodals.T; label="Upper", linestyle=:dashdot, color=:gray)
    savefig(plt, path)
end

function _plot_maxwell(path, boundary, boundary_results, cep)
    plt = plot(; xlabel="μ (MeV)", ylabel="T (MeV)", title="Maxwell Status", legend=:topright)
    successes = Float64[]
    mu_success = Float64[]
    failures = Float64[]
    failure_mu = Float64[]
    for T in sort(collect(keys(boundary_results)))
        res = boundary_results[T]
        if res.converged && res.mu_coex_MeV !== nothing
            push!(successes, T)
            push!(mu_success, res.mu_coex_MeV)
        else
            reason = get(res.details, :reason, "")
            if res.mu_coex_MeV !== nothing
                push!(failures, T)
                push!(failure_mu, res.mu_coex_MeV)
            elseif reason == "no_sign_change"
                push!(failures, T)
                push!(failure_mu, NaN)
            end
        end
    end
    if !isempty(successes)
        scatter!(plt, mu_success, successes; label="Coexistence", color=:red, marker=:circle)
    end
    if !isempty(failures)
        scatter!(plt, failure_mu, failures; label="Failed", color=:gray, marker=:x)
    end
    if cep.has_cep && cep.mu_cep_MeV !== nothing && cep.T_cep_MeV !== nothing
        scatter!(plt, [cep.mu_cep_MeV], [cep.T_cep_MeV]; label="CEP", color=:blue, marker=:star5, markersize=10)
    end
    savefig(plt, path)
end

function _dump_curves(path, curves, xi::Float64)
    open(path, "w") do io
        println(io, "T_MeV,mu_MeV,rho,xi")
        for T in sort(collect(keys(curves)))
            mu_vec, rho_vec = curves[T]
            n = min(length(mu_vec), length(rho_vec))
            for i in 1:n
                println(io, join((T, mu_vec[i], rho_vec[i], xi), ','))
            end
        end
    end
end

function _dump_maxwell(path, boundary_results::Dict{Float64, MaxwellResult}, xi::Float64)
    open(path, "w") do io
        println(io, "T_MeV,xi,status,mu_coex_MeV,rho_gas,rho_liquid,area_residual,iterations,reason")
        for T in sort(collect(keys(boundary_results)))
            res = boundary_results[T]
            status = res.converged ? "success" : "failure"
            reason = get(res.details, :reason, "")
            values = (
                T,
                xi,
                status,
                something(res.mu_coex_MeV, ""),
                something(res.rho_gas, ""),
                something(res.rho_liquid, ""),
                something(res.area_residual, ""),
                res.iterations,
                reason,
            )
            println(io, join(values, ','))
        end
    end
end

function _dump_maxwell_inputs(path, curves, xi::Float64)
    open(path, "w") do io
        println(io, "T_MeV,xi,sample_index,rho,mu_MeV,has_s_shape,mu_spinodal_low,mu_spinodal_high")
        for T in sort(collect(keys(curves)))
            mu_vec, rho_vec = curves[T]
            rho_sorted, mu_sorted = _sorted_curve(mu_vec, rho_vec)
            hint = detect_s_shape(mu_vec, rho_vec)
            for idx in eachindex(rho_sorted)
                println(io, join((T, xi, idx, rho_sorted[idx], mu_sorted[idx], hint.has_s_shape,
                    something(hint.mu_spinodal_low, ""), something(hint.mu_spinodal_high, "")), ','))
            end
        end
    end
end

function _sorted_curve(mu_vals, rho_vals)
    n = min(length(mu_vals), length(rho_vals))
    pairs = Vector{Tuple{Float64, Float64}}()
    sizehint!(pairs, n)
    for i in 1:n
        mu = Float64(mu_vals[i])
        rho = Float64(rho_vals[i])
        (isfinite(mu) && isfinite(rho)) || continue
        push!(pairs, (rho, mu))
    end
    sort!(pairs; by=first)
    return Float64[first(p) for p in pairs], Float64[last(p) for p in pairs]
end

function _dump_spinodals(path, spinodals, xi::Float64)
    open(path, "w") do io
        println(io, "T_MeV,mu_low_MeV,mu_high_MeV,xi")
        n = length(spinodals.T)
        for i in 1:n
            println(io, join((spinodals.T[i], spinodals.low[i], spinodals.high[i], xi), ','))
        end
    end
end

function _get_or_compute_cep(processed_dir::AbstractString, curves;
        xi::Float64=0.0, xi_tol::Float64=1e-6, temp_tol::Real=0.01,
        force_update::Bool=false, curve_fetcher=nothing)
    cep_path = joinpath(processed_dir, "cep.csv")
    if !force_update && isfile(cep_path)
        for row in CSV.File(cep_path)
            xi_val = try Float64(row[:xi]) catch; continue end
            if abs(xi_val - xi) <= xi_tol
                has = try String(row[:has_cep]) catch; "false" end
                if lowercase(has) in ("true", "t", "1")
                    mu = try Float64(row[:mu_cep_MeV]) catch; nothing end
                    T = try Float64(row[:T_cep_MeV]) catch; nothing end
                    return CEPResult(true, T, mu, nothing, nothing, Dict(:cached => true))
                else
                    return CEPResult(false, nothing, nothing, nothing, nothing, Dict(:cached => true))
                end
            end
        end
    end
    cep = find_cep(curves; temp_tol=temp_tol, curve_fetcher=curve_fetcher)
    _write_cep_row(cep_path, cep, xi; force=force_update, xi_tol=xi_tol)
    return cep
end

function _write_cep_row(path::AbstractString, cep, xi; force::Bool, xi_tol::Float64)
    mkpath(dirname(path))
    if force && isfile(path)
        retained = Vector{NTuple{4, String}}()
        for row in CSV.File(path)
            xi_val = try Float64(row[:xi]) catch; continue end
            if abs(xi_val - xi) <= xi_tol
                continue
            end
            push!(retained, (
                string(row[:has_cep]),
                string(row[:mu_cep_MeV]),
                string(row[:T_cep_MeV]),
                string(row[:xi]),
            ))
        end
        open(path, "w") do io
            println(io, "has_cep,mu_cep_MeV,T_cep_MeV,xi")
            for r in retained
                println(io, join(r, ','))
            end
            println(io, join(_cep_row_values(cep, xi), ','))
        end
        return
    end
    exists = isfile(path)
    open(path, exists ? "a" : "w") do io
        if !exists
            println(io, "has_cep,mu_cep_MeV,T_cep_MeV,xi")
        end
        println(io, join(_cep_row_values(cep, xi), ','))
    end
end

function _cep_row_values(cep, xi)
    if cep.has_cep && cep.mu_cep_MeV !== nothing && cep.T_cep_MeV !== nothing
        return ("true", string(cep.mu_cep_MeV), string(cep.T_cep_MeV), string(xi))
    else
        return ("false", "", "", string(xi))
    end
end

function _fetch_curve_via_trho!(curves, opts::CLIOptions, T_target::Float64)
    tmp = tempname() * ".csv"
    try
        TrhoScan.run_trho_scan(
            ; T_values=[T_target],
              rho_values=TrhoScan.DEFAULT_RHO_VALUES,
              xi_values=[opts.xi],
              output_path=tmp,
              overwrite=true,
              resume=false,
        )
        grouped = load_trho_curves(tmp; xi=opts.xi, tol=opts.xi_tol)
        match_key = _nearest_temperature_key(grouped, T_target)
        match_key === nothing && return nothing
        curve_dict = build_curves(Dict(match_key => grouped[match_key]))
        new_curve = curve_dict[match_key]
        _append_trho_csv!(tmp, opts.trho_csv)
        curves[T_target] = new_curve
        return new_curve
    finally
        isfile(tmp) && rm(tmp; force=true)
    end
end

function _nearest_temperature_key(curves::Dict{Float64, T}, target::Float64; tol::Float64=1e-4) where {T}
    haskey(curves, target) && return target
    rounded = round(target; digits=6)
    haskey(curves, rounded) && return rounded
    best = nothing
    best_diff = tol
    for key in keys(curves)
        diff = abs(key - target)
        if diff <= best_diff
            best = key
            best_diff = diff
        end
    end
    return best
end

function _append_trho_csv!(src::AbstractString, dest::AbstractString)
    mkpath(dirname(dest))
    if !isfile(dest)
        cp(src, dest; force=true)
        return
    end
    first_line = true
    open(dest, "a") do out
        for line in eachline(src)
            if first_line
                first_line = false
                continue
            end
            isempty(strip(line)) && continue
            println(out, line)
        end
    end
end
function main()
    opts = parse_args(copy(ARGS))
    plot_phase_diagram(opts)
end

main()
