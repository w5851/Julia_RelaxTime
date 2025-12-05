module CEPFinder

using ..PhaseTransition: SShapeResult, detect_s_shape
using ..TrhoScan: run_trho_scan, DEFAULT_RHO_VALUES, DEFAULT_OUTPUT_PATH
using ..SeedCache: DEFAULT_SEED_PATH
using ..AdaptiveRhoRefinement: AdaptiveRhoConfig, suggest_refinement_points, merge_rho_values

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_PROCESSED_CURVES_PATH = normpath(joinpath(PROJECT_ROOT, "data", "processed", "results", "pnjl", "curves.csv"))

export CEPResult, find_cep, build_curves

struct CEPResult
    has_cep::Bool
    T_cep_MeV::Union{Nothing, Float64}
    mu_cep_MeV::Union{Nothing, Float64}
    rho_cep::Union{Nothing, Float64}
    bracket::Union{Nothing, Tuple{Float64, Float64}}
    details::Dict{Symbol, Any}
end

CEPResult() = CEPResult(false, nothing, nothing, nothing, nothing, Dict{Symbol, Any}())

"""Convert grouped curve data (Dict{T => Vector{Tuple{μ,ρ}}}) into sorted arrays."""
function build_curves(grouped::Dict{Float64, Vector{Tuple{Float64, Float64}}})
    curves = Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}()
    for (T, samples) in grouped
        if length(samples) < 3
            continue
        end
        mu_vals = Float64[first(s) for s in samples]
        rho_vals = Float64[last(s) for s in samples]
        push!(curves, T => (mu_vals, rho_vals))
    end
    return curves
end

"""
    find_cep(curves; min_points=6, eps=0.0, temp_tol=0.01)

Given a dictionary `curves` keyed by temperature (MeV) with `(μ, ρ)` vectors,
attempt to locate the critical end point (CEP) by identifying the last temperature
with an S-shaped curve and the first one without.
"""
function find_cep(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}};
    min_points::Int=6, eps::Real=0.0, temp_tol::Real=0.01,
    curve_fetcher::Union{Nothing, Function}=nothing, auto_fetch=:default)
    fetch_cb = curve_fetcher
    if fetch_cb === nothing
        if auto_fetch === nothing
            fetch_cb = nothing
        else
            config = auto_fetch === :default ? _default_auto_fetch_config() : auto_fetch
            fetch_cb = _build_auto_fetcher(config, curves)
        end
    end
    temperatures = sort(collect(keys(curves)))
    isempty(temperatures) && return CEPResult()

    last_with = nothing
    first_without = nothing
    shape_cache = Dict{Float64, SShapeResult}()
    for T in temperatures
        mu_vals, rho_vals = curves[T]
        result = detect_s_shape(mu_vals, rho_vals; eps=eps, min_points=min_points)
        shape_cache[T] = result
        if result.has_s_shape
            last_with = (T, result)
        elseif last_with !== nothing
            first_without = (T, result)
            break
        end
    end

    if last_with === nothing || first_without === nothing
        return CEPResult(false, nothing, nothing, nothing, nothing, Dict(:reason => "no_bracket"))
    end

    T_low, res_low = last_with
    T_high, _ = first_without
    orig_T_low = T_low
    orig_T_high = T_high

    # Continuous temperature bisection. At each step, evaluate S-shape at
    # the midpoint temperature by either taking the exact curve (if present)
    # or linearly interpolating the μ values between the two neighboring
    # temperature curves (assumes same ρ sampling across temperatures).
    while (T_high - T_low) > temp_tol
        T_mid = 0.5 * (T_low + T_high)
        mu_mid, rho_mid = _curve_at_T!(curves, T_mid, fetch_cb)
        mu_mid === nothing && break
        res_mid = detect_s_shape(mu_mid, rho_mid; eps=eps, min_points=min_points)
        shape_cache[T_mid] = res_mid
        if res_mid.has_s_shape
            T_low = T_mid
            res_low = res_mid
        else
            T_high = T_mid
        end
    end

    T_cep = T_low + (T_high - T_low) / 2
    mu_cep = (something(res_low.mu_spinodal_low, NaN) + something(res_low.mu_spinodal_high, NaN)) / 2
    rho_cep = (something(res_low.rho_spinodal_low, NaN) + something(res_low.rho_spinodal_high, NaN)) / 2

    details = Dict(
        :bracket_low => T_low,
        :bracket_high => T_high,
        :spinodal_low => res_low.mu_spinodal_low,
        :spinodal_high => res_low.mu_spinodal_high,
        :temp_tol => temp_tol,
    )

    # For backward compatibility with tests that expect the original
    # discrete-temperature bracket, return `bracket` as the initial
    # (last_with, first_without) pair while recording the refined
    # bracket in details.
    return CEPResult(true, T_cep, mu_cep, rho_cep, (orig_T_low, orig_T_high), details)
end

function _curve_at_T!(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}},
        T_target::Float64, curve_fetcher)
    if haskey(curves, T_target)
        return curves[T_target]
    end
    curve_fetcher === nothing && return nothing, nothing
    fetched = curve_fetcher(T_target)
    fetched === nothing && return nothing, nothing
    mu_vals, rho_vals = fetched
    curves[T_target] = (mu_vals, rho_vals)
    return curves[T_target]
end

function _curve_at_T(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}}, T_target::Float64)
    # direct lookup
    if haskey(curves, T_target)
        mu, rho = curves[T_target]
        return mu, rho
    end
    # find nearest lower and upper available temperatures
    temps = sort(collect(keys(curves)))
    lower = nothing
    upper = nothing
    for t in temps
        if t < T_target
            lower = t
        elseif t > T_target
            upper = t
            break
        end
    end
    if lower === nothing || upper === nothing
        return nothing, nothing
    end
    mu_lo, rho_lo = curves[lower]
    mu_hi, rho_hi = curves[upper]
    # Require same rho sampling to interpolate; otherwise abort
    length(rho_lo) == length(rho_hi) || return nothing, nothing
    # linear interpolation in temperature for μ at each sample index
    α = (T_target - lower) / (upper - lower)
    mu_mid = similar(mu_lo)
    for i in eachindex(mu_lo)
        mu_mid[i] = mu_lo[i] + α * (mu_hi[i] - mu_lo[i])
    end
    return mu_mid, rho_lo
end

_default_auto_fetch_config() = (
    xi = 0.0,
    rho_values = DEFAULT_RHO_VALUES,
    seed_path = DEFAULT_SEED_PATH,
    output_path = DEFAULT_OUTPUT_PATH,
    processed_path = DEFAULT_PROCESSED_CURVES_PATH,
    adaptive = true,
    adaptive_config = AdaptiveRhoConfig(),
    adaptive_window = 5.0,
    rho_digits = 6,
)

function _build_auto_fetcher(config, curves)
    config isa NamedTuple || config isa AbstractDict || error("auto_fetch must be a NamedTuple or Dict")
    xi = Float64(_get_config(config, :xi, 0.0))
    rho_values = _normalize_vector(_get_config(config, :rho_values, DEFAULT_RHO_VALUES))
    seed_path = String(_get_config(config, :seed_path, DEFAULT_SEED_PATH))
    p_num = Int(_get_config(config, :p_num, 24))
    t_num = Int(_get_config(config, :t_num, 8))
    xi_tol = Float64(_get_config(config, :xi_tol, 1e-6))
    output_path = _normalize_path(_get_config(config, :output_path, DEFAULT_OUTPUT_PATH))
    processed_path = _normalize_path(_get_config(config, :processed_path, DEFAULT_PROCESSED_CURVES_PATH))
    solver_kwargs = _get_config(config, :solver_kwargs, NamedTuple())
    solver_kwargs isa NamedTuple || error("auto_fetch solver_kwargs must be a NamedTuple")
    runner = _get_config(config, :runner, run_trho_scan)
    loader = _get_config(config, :loader, _extract_curve_from_csv)
    generator = _get_config(config, :generator, nothing)
    adaptive = Bool(_get_config(config, :adaptive, true))
    rho_digits = Int(_get_config(config, :rho_digits, 6))
    adaptive_cfg = _normalize_adaptive_config(_get_config(config, :adaptive_config, AdaptiveRhoConfig()))
    adaptive_window = _get_config(config, :adaptive_window, 5.0)
    return T_target -> _auto_fetch_curve(Float64(T_target), curves, (;
        xi=xi,
        rho_values=rho_values,
        seed_path=seed_path,
        p_num=p_num,
        t_num=t_num,
        xi_tol=xi_tol,
        output_path=output_path,
        processed_path=processed_path,
        solver_kwargs=solver_kwargs,
        runner=runner,
        loader=loader,
        generator=generator,
        adaptive=adaptive,
        rho_digits=rho_digits,
        adaptive_config=adaptive_cfg,
        adaptive_window=isnothing(adaptive_window) ? nothing : Float64(adaptive_window),
    ))
end

function _auto_fetch_curve(T_target::Float64, curves, cfg)
    generator = cfg.generator
    if generator !== nothing
        return _normalize_curve_tuple(generator(T_target))
    end
    rho_values = cfg.rho_values
    if cfg.adaptive
        extras = _adaptive_rho_candidates(curves, T_target, cfg)
        if !isempty(extras)
            rho_values = merge_rho_values(rho_values, extras; digits=cfg.rho_digits)
        end
    end
    scratch = tempname() * ".csv"
    runner_kwargs = (; T_values=[T_target], rho_values=rho_values, xi_values=[cfg.xi],
        output_path=scratch, seed_path=cfg.seed_path, overwrite=true, resume=false,
        p_num=cfg.p_num, t_num=cfg.t_num)
    cfg.runner(; runner_kwargs..., cfg.solver_kwargs...)
    curve = cfg.loader(scratch; xi=cfg.xi, tol=cfg.xi_tol)
    normalized = _normalize_curve_tuple(curve)
    if normalized !== nothing
        cfg.output_path !== nothing && _append_trho_csv!(scratch, cfg.output_path)
        cfg.processed_path !== nothing && _append_processed_curve!(cfg.processed_path, T_target, normalized, cfg.xi)
    end
    isfile(scratch) && rm(scratch; force=true)
    return normalized
end

function _extract_curve_from_csv(path::AbstractString; xi::Float64, tol::Float64=1e-6)
    mu_vals = Float64[]
    rho_vals = Float64[]
    isfile(path) || return nothing
    open(path, "r") do io
        first_line = true
        for line in eachline(io)
            if first_line
                first_line = false
                continue
            end
            isempty(strip(line)) && continue
            cols = split(line, ','; limit=19)
            length(cols) < 18 && continue
            xi_val = _maybe_parse_float(cols[3])
            rho_val = _maybe_parse_float(cols[2])
            mu_val = _maybe_parse_float(cols[7])
            xi_val === nothing && continue
            rho_val === nothing && continue
            mu_val === nothing && continue
            (isfinite(xi_val) && isfinite(rho_val) && isfinite(mu_val)) || continue
            abs(xi_val - xi) <= tol || continue
            converged = lowercase(strip(cols[18]))
            if !(converged in ("true", "t", "1"))
                continue
            end
            push!(mu_vals, mu_val)
            push!(rho_vals, rho_val)
        end
    end
    isempty(mu_vals) && return nothing
    return (mu_vals, rho_vals)
end

function _adaptive_rho_candidates(curves, T_target, cfg)
    additions = Float64[]
    window = cfg.adaptive_window
    for (T, (mu_vals, rho_vals)) in curves
        if window !== nothing && abs(T - T_target) > window
            continue
        end
        extra = suggest_refinement_points(rho_vals, mu_vals; config=cfg.adaptive_config)
        isempty(extra) && continue
        append!(additions, extra)
    end
    return additions
end

function _append_trho_csv!(src::AbstractString, dest::AbstractString)
    dest === nothing && return
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

function _append_processed_curve!(path::AbstractString, T::Float64, curve, xi::Float64)
    path === nothing && return
    mkpath(dirname(path))
    mu_vals, rho_vals = curve
    exists = isfile(path)
    open(path, exists ? "a" : "w") do io
        if !exists
            println(io, "T_MeV,mu_MeV,rho,xi")
        end
        for i in eachindex(mu_vals)
            println(io, join((T, mu_vals[i], rho_vals[i], xi), ','))
        end
    end
end

function _normalize_curve_tuple(data)
    data === nothing && return nothing
    (data isa Tuple && length(data) == 2) || error("curve fetcher must return (mu_vals, rho_vals)")
    mu_vals = Float64.(collect(data[1]))
    rho_vals = Float64.(collect(data[2]))
    length(mu_vals) == length(rho_vals) || error("curve fetcher returned mismatched lengths")
    return (mu_vals, rho_vals)
end

function _normalize_vector(values)
    values === nothing && return DEFAULT_RHO_VALUES
    return Float64.(collect(values))
end

function _normalize_adaptive_config(value)
    if value isa AdaptiveRhoConfig
        return value
    elseif value isa NamedTuple
        return AdaptiveRhoConfig(; value...)
    elseif value isa AbstractDict
        pairs = Pair{Symbol, Any}[]
        for (k, v) in value
            push!(pairs, Symbol(k) => v)
        end
        return AdaptiveRhoConfig(; pairs...)
    else
        return AdaptiveRhoConfig()
    end
end

function _normalize_path(value)
    value === nothing && return nothing
    str = String(value)
    isempty(str) && return nothing
    return normpath(str)
end

function _maybe_parse_float(value)
    try
        return parse(Float64, strip(value))
    catch
        return nothing
    end
end

function _get_config(config, key::Symbol, default)
    if config isa NamedTuple
        return hasproperty(config, key) ? getfield(config, key) : default
    elseif config isa AbstractDict
        return get(config, key, default)
    else
        return default
    end
end

end # module
