#!/usr/bin/env julia

"""Run one isolated PNJL CEP narrow-window pilot job.

The script intentionally owns only the pilot sampler.  It does not write
`data/reference/pnjl` or invoke the production phase pipeline.  Each job is
one `(xi, method)` pair; the workflow supplies the immutable calculation SHA.
"""

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using CSV
using Dates
using JSON3
using SHA
using Statistics

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const HBARC_MEV_FM = 197.3269804
const METHODS = (:c2_dense_baseline, :rho_support_cascade, :high_resolution_oracle)
const CANONICAL_CEP_PATH = joinpath(PROJECT_ROOT, "data", "reference", "pnjl", "cep.csv")

# JSON (unlike CSV) does not permit IEEE NaN/Inf literals.  A pilot job is
# allowed to finish without finding a CEP bracket; those diagnostic fields are
# represented as `null` in JSON rather than turning an otherwise valid solver
# run into a failed Actions job.  Keep this conversion local to the artifact
# boundary so the numerical and CSV paths continue to expose their values.
_json_safe(x::Nothing) = nothing
_json_safe(x::Bool) = x
_json_safe(x::Integer) = x
_json_safe(x::AbstractFloat) = isfinite(x) ? x : nothing
_json_safe(x::AbstractString) = String(x)
_json_safe(x::Symbol) = String(x)
_json_safe(x::NamedTuple) = Dict(string(name) => _json_safe(getproperty(x, name)) for name in keys(x))
_json_safe(x::AbstractDict) = Dict(string(key) => _json_safe(value) for (key, value) in x)
_json_safe(x::AbstractArray) = [_json_safe(value) for value in x]
_json_safe(x::Tuple) = [_json_safe(value) for value in x]
_json_safe(x) = x

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
using Main.Models

if !isdefined(Main, :PNJLCriticalityFeasibility)
    include(joinpath(@__DIR__, "pnjl_criticality_feasibility.jl"))
end
const Criticality = Main.PNJLCriticalityFeasibility

Base.@kwdef struct PilotConfig
    xi::Float64 = 0.0
    method::Symbol = :c2_dense_baseline
    output_dir::String = joinpath(PROJECT_ROOT, "data", "outputs", "pnjl_cep_narrow_pilot")
    tag::String = "cep_narrow_pilot_v1"
    calculation_sha::String = ""
    p_num::Int = 24
    t_num::Int = 8
    thermo_quadrature_policy::Symbol = :rs_reduced_adaptive
    thermo_quadrature_rtol::Float64 = 1e-8
    thermo_quadrature_atol::Float64 = 1e-10
    thermo_quadrature_maxevals::Int = 10^7
    iterations::Int = 80
    rho_min::Float64 = 0.0
    rho_max::Float64 = 4.0
    baseline_rho_step::Float64 = 0.00625
    oracle_rho_step::Float64 = 0.00625
    oracle_refine_step::Float64 = 0.003125
    initial_window_MeV::Float64 = 4.0
    max_window_MeV::Float64 = 8.0
    initial_T_step_MeV::Float64 = 1.0
    bracket_tol_MeV::Float64 = 0.125
    targeted_max_points::Int = 12
end

mutable struct PilotMemo
    model::Models.AbstractQCDModel
    config::PilotConfig
    telemetry::Models.SolverWorkTelemetry
    point_cache::Dict{Tuple{Float64, Float64, Float64}, NamedTuple}
    slice_cache::Dict{Float64, NamedTuple}
    previous_seeds::Dict{Float64, Vector{Float64}}
    point_requests::Int
    cache_hits::Int
    unique_solves::Int
    targeted_additions::Int
    scan_retries::Int
end

function PilotMemo(config::PilotConfig)
    return PilotMemo(
        Models.create_model(:PNJL),
        config,
        Models.SolverWorkTelemetry(),
        Dict{Tuple{Float64, Float64, Float64}, NamedTuple}(),
        Dict{Float64, NamedTuple}(),
        Dict{Float64, Vector{Float64}}(),
        0,
        0,
        0,
        0,
        0,
    )
end

@inline function _parse_cli(args)
    values = Dict{String, String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        if startswith(arg, "--") && occursin("=", arg)
            key, value = split(arg[3:end], "="; limit=2)
            values[key] = value
        elseif startswith(arg, "--")
            key = arg[3:end]
            i == length(args) && throw(ArgumentError("missing value for --$(key)"))
            i += 1
            values[key] = args[i]
        else
            throw(ArgumentError("unexpected argument: $(arg)"))
        end
        i += 1
    end
    return values
end

@inline _get_float(values, key, default) = haskey(values, key) ? parse(Float64, values[key]) : default
@inline _get_int(values, key, default) = haskey(values, key) ? parse(Int, values[key]) : default
@inline _get_string(values, key, default) = get(values, key, default)

function parse_config(args)
    values = _parse_cli(args)
    xi = _get_float(values, "xi", 0.0)
    method = Symbol(_get_string(values, "method", "c2_dense_baseline"))
    method in METHODS || throw(ArgumentError("method must be one of $(METHODS), got $(method)"))
    output_dir = abspath(_get_string(values, "output-dir", joinpath(PROJECT_ROOT, "data", "outputs", "pnjl_cep_narrow_pilot")))
    calculation_sha = _get_string(values, "calculation-sha", "")
    isempty(calculation_sha) || occursin(r"^[0-9a-fA-F]{40}$", calculation_sha) ||
        throw(ArgumentError("calculation-sha must be an immutable 40-character SHA"))
    return PilotConfig(
        xi=xi,
        method=method,
        output_dir=output_dir,
        tag=_get_string(values, "tag", "cep_narrow_pilot_v1"),
        calculation_sha=calculation_sha,
        p_num=_get_int(values, "p-num", 24),
        t_num=_get_int(values, "t-num", 8),
        thermo_quadrature_policy=Symbol(_get_string(values, "thermo-quadrature-policy", "rs_reduced_adaptive")),
        thermo_quadrature_rtol=_get_float(values, "thermo-quadrature-rtol", 1e-8),
        thermo_quadrature_atol=_get_float(values, "thermo-quadrature-atol", 1e-10),
        thermo_quadrature_maxevals=_get_int(values, "thermo-quadrature-maxevals", 10^7),
        iterations=_get_int(values, "iterations", 80),
    )
end

function _canonical_cep(xi::Float64)
    isfile(CANONICAL_CEP_PATH) || throw(ArgumentError("canonical CEP file missing: $(CANONICAL_CEP_PATH)"))
    rows = collect(CSV.File(CANONICAL_CEP_PATH))
    candidates = [row for row in rows if hasproperty(row, :xi) && abs(Float64(row.xi) - xi) <= 1e-9]
    isempty(candidates) && throw(ArgumentError("canonical CEP has no xi=$(xi) row"))
    row = first(candidates)
    return (T=Float64(row.T_CEP_MeV), muq=Float64(row.muq_CEP_MeV))
end

@inline function _rho_grid(lo::Float64, hi::Float64, step::Float64)
    values = collect(lo:step:hi)
    if isempty(values) || values[end] < hi - 64eps(max(abs(lo), abs(hi), abs(step), 1.0))
        push!(values, hi)
    else
        values[end] = hi
    end
    return values
end

function _default_seed(rho::Float64)
    if rho < 0.5
        return copy(Models.HADRON_SEED_8)
    elseif rho < 2.0
        return copy(Models.MEDIUM_SEED_8)
    end
    return copy(Models.HIGH_DENSITY_SEED_8)
end

function _seed_for(memo::PilotMemo, rho::Float64)
    haskey(memo.previous_seeds, rho) && return copy(memo.previous_seeds[rho])
    return _default_seed(rho)
end

function _solve_point!(memo::PilotMemo, T_MeV::Float64, rho::Float64; targeted::Bool=false)
    config = memo.config
    key = (T_MeV, config.xi, rho)
    memo.point_requests += 1
    if haskey(memo.point_cache, key)
        memo.cache_hits += 1
        return memo.point_cache[key]
    end

    seed = _seed_for(memo, rho)
    T_fm = T_MeV / HBARC_MEV_FM
    result = Models.solve(
        memo.model,
        Models.FixedRho(rho),
        T_fm;
        seed_guess=seed,
        xi=config.xi,
        p_num=config.p_num,
        t_num=config.t_num,
        thermo_quadrature_policy=config.thermo_quadrature_policy,
        thermo_quadrature_rtol=config.thermo_quadrature_rtol,
        thermo_quadrature_atol=config.thermo_quadrature_atol,
        thermo_quadrature_maxevals=config.thermo_quadrature_maxevals,
        residual_norm_max=1e-6,
        iterations=config.iterations,
        semantic_mode=:ground_state,
        nlsolve_method=:newton,
        trust_region_fallback=true,
        work_telemetry=memo.telemetry,
    )
    converged = Models.solver_result_is_success(result; residual_norm_max=1e-6)
    finite = all(isfinite, (result.pressure, result.residual_norm)) &&
             all(isfinite, result.solution) && all(isfinite, result.mu_vec)
    row = (
        T_MeV=T_MeV,
        xi=config.xi,
        rho=rho,
        muq_MeV=finite ? HBARC_MEV_FM * mean(result.mu_vec) : NaN,
        pressure_fm4=Float64(result.pressure),
        residual_norm=Float64(result.residual_norm),
        iterations=Int(result.iterations),
        converged=Bool(converged && finite),
        finite=Bool(finite),
        solution=Float64.(result.solution),
    )
    memo.point_cache[key] = row
    memo.unique_solves += 1
    targeted && (memo.targeted_additions += 1)
    if row.converged
        memo.previous_seeds[rho] = copy(row.solution)
    end
    if !row.converged
        Models.record_scan_retry!(memo.telemetry)
        memo.scan_retries += 1
        retry_seed = copy(Models.HIGH_DENSITY_SEED_8)
        retry = Models.solve(
            memo.model,
            Models.FixedRho(rho),
            T_fm;
            seed_guess=retry_seed,
            xi=config.xi,
            p_num=config.p_num,
            t_num=config.t_num,
            thermo_quadrature_policy=config.thermo_quadrature_policy,
            thermo_quadrature_rtol=config.thermo_quadrature_rtol,
            thermo_quadrature_atol=config.thermo_quadrature_atol,
            thermo_quadrature_maxevals=config.thermo_quadrature_maxevals,
            residual_norm_max=1e-6,
            iterations=config.iterations,
            semantic_mode=:ground_state,
            nlsolve_method=:trust_region,
            trust_region_fallback=false,
            work_telemetry=memo.telemetry,
        )
        retry_finite = all(isfinite, (retry.pressure, retry.residual_norm)) && all(isfinite, retry.solution)
        if Models.solver_result_is_success(retry; residual_norm_max=1e-6) && retry_finite
            row = merge(row, (
                muq_MeV=HBARC_MEV_FM * mean(retry.mu_vec),
                pressure_fm4=Float64(retry.pressure),
                residual_norm=Float64(retry.residual_norm),
                iterations=Int(retry.iterations),
                converged=true,
                finite=true,
                solution=Float64.(retry.solution),
            ))
            memo.point_cache[key] = row
            memo.previous_seeds[rho] = copy(row.solution)
        end
    end
    return row
end

@inline function _curve_status(rho::Vector{Float64}, mu::Vector{Float64})
    length(rho) == length(mu) || return (status=:failed, mu_transition=NaN, area=NaN, s_shape=false, reason="length_mismatch")
    all(isfinite, mu) || return (status=:failed, mu_transition=NaN, area=NaN, s_shape=false, reason="nonfinite_mu")
    sres = Models.detect_s_shape(mu, rho)
    sres.has_s_shape || return (status=:monotone, mu_transition=NaN, area=NaN, s_shape=false, reason="no_s_shape")
    mres = Models.maxwell_construction(mu, rho; spinodal_hint=sres, tol_area=1e-4)
    if mres.converged && mres.mu_transition !== nothing
        return (status=:resolved_s_shape, mu_transition=Float64(mres.mu_transition), area=Float64(something(mres.area_residual, NaN)), s_shape=true, reason="maxwell")
    end
    return (status=:near_critical, mu_transition=NaN, area=NaN, s_shape=true, reason="maxwell_unresolved")
end

function _baseline_slice!(memo::PilotMemo, T_MeV::Float64; step::Float64=memo.config.baseline_rho_step)
    rho = _rho_grid(memo.config.rho_min, memo.config.rho_max, step)
    points = [_solve_point!(memo, T_MeV, value) for value in rho]
    mu = [point.muq_MeV for point in points]
    status = _curve_status(rho, mu)
    return merge(status, (T_MeV=T_MeV, rho=rho, mu=mu, points=points, targeted_points=0, method_note="dense secant + Maxwell"))
end

function _cascade_slice!(memo::PilotMemo, T_MeV::Float64, prior)
    rho = _rho_grid(memo.config.rho_min, memo.config.rho_max, 0.05)
    points = [_solve_point!(memo, T_MeV, value) for value in rho]
    mu = [point.muq_MeV for point in points]
    targeted_values = Float64[]
    sampler = value -> begin
        push!(targeted_values, Float64(value))
        return _solve_point!(memo, T_MeV, Float64(value); targeted=true).muq_MeV
    end
    assessment = Criticality.analyze_curve_cascade(
        rho,
        mu;
        sample_mu=sampler,
        prior=prior,
        config=Criticality.RhoSupportConfig(max_extra_points=memo.config.targeted_max_points),
    )
    # Re-run the physical Maxwell check on the union of coarse and targeted
    # points. This keeps the cascade classifier and Maxwell input auditable.
    all_rho = sort(unique(vcat(rho, targeted_values)))
    all_points = [_solve_point!(memo, T_MeV, value) for value in all_rho]
    all_mu = [point.muq_MeV for point in all_points]
    status = _curve_status(all_rho, all_mu)
    if assessment.status == :near_critical && status.status == :monotone
        status = merge(status, (status=:near_critical, s_shape=false, reason="cascade_near_critical"))
    end
    return merge(status, (
        T_MeV=T_MeV,
        rho=all_rho,
        mu=all_mu,
        points=all_points,
        targeted_points=length(unique(targeted_values)),
        method_note="rho-support cascade + re-Maxwell",
        cascade_status=assessment.status,
        spinodal_rho_center=assessment.spinodal_rho_center,
        spinodal_rho_gap=assessment.spinodal_rho_gap,
    ))
end

function _slice!(memo::PilotMemo, T_MeV::Float64, prior)
    haskey(memo.slice_cache, T_MeV) && return memo.slice_cache[T_MeV]
    started = time_ns()
    before_requests = memo.point_requests
    before_unique = memo.unique_solves
    before_hits = memo.cache_hits
    before_targeted = memo.targeted_additions
    before_retries = memo.scan_retries
    before_telemetry = Models.solver_work_snapshot(memo.telemetry)
    curve = if memo.config.method == :rho_support_cascade
        _cascade_slice!(memo, T_MeV, prior)
    else
        step = memo.config.method == :high_resolution_oracle ? memo.config.oracle_rho_step : memo.config.baseline_rho_step
        _baseline_slice!(memo, T_MeV; step=step)
    end
    after_telemetry = Models.solver_work_snapshot(memo.telemetry)
    curve = merge(curve, (
        wall_seconds=(time_ns() - started) / 1e9,
        requests=memo.point_requests - before_requests,
        unique_solves=memo.unique_solves - before_unique,
        cache_hits=memo.cache_hits - before_hits,
        targeted_additions=memo.targeted_additions - before_targeted,
        residual_calls=(after_telemetry.nlsolve_f_calls + after_telemetry.postprocess_residual_calls) -
            (before_telemetry.nlsolve_f_calls + before_telemetry.postprocess_residual_calls),
        jacobian_calls=after_telemetry.nlsolve_g_calls - before_telemetry.nlsolve_g_calls,
        newton_iterations=after_telemetry.newton_iterations - before_telemetry.newton_iterations,
        trust_region_iterations=after_telemetry.trust_region_iterations - before_telemetry.trust_region_iterations,
        fallback_count=after_telemetry.root_fallbacks - before_telemetry.root_fallbacks,
        governed_rescue_count=after_telemetry.governed_rescues - before_telemetry.governed_rescues,
        retry_count=memo.scan_retries - before_retries,
    ))
    memo.slice_cache[T_MeV] = curve
    return curve
end

function _valid_low_side(curve)
    return curve.status == :resolved_s_shape || curve.status == :near_critical
end

function _status_bracket(slice_cache)
    temps = sort(collect(keys(slice_cache)))
    low = [T for T in temps if _valid_low_side(slice_cache[T])]
    high = [
        T for T in temps if !_valid_low_side(slice_cache[T]) &&
            any(_valid_low_side(slice_cache[τ]) for τ in temps if τ < T)
    ]
    isempty(low) && return nothing
    isempty(high) && return nothing
    T_low = maximum(low)
    higher = filter(>(T_low), high)
    isempty(higher) && return nothing
    return (T_low=T_low, T_high=minimum(higher))
end

function _find_bracket!(memo::PilotMemo, center::Float64)
    prior = nothing
    window = memo.config.initial_window_MeV
    while window <= memo.config.max_window_MeV + 1e-9
        initial = collect((center - window):memo.config.initial_T_step_MeV:(center + window))
        for T in initial
            _slice!(memo, Float64(T), prior)
            curve = memo.slice_cache[Float64(T)]
            if memo.config.method == :rho_support_cascade &&
               curve.spinodal_rho_center !== nothing && curve.spinodal_rho_gap !== nothing
                prior = Criticality.RhoSupportPrior(Float64(T), Float64(curve.spinodal_rho_center), Float64(curve.spinodal_rho_gap))
            end
        end
        bracket = _status_bracket(memo.slice_cache)
        bracket !== nothing && return bracket
        window += 4.0
    end
    return nothing
end

function _refine_bracket!(memo::PilotMemo, bracket)
    bracket === nothing && return nothing
    low = Float64(bracket.T_low)
    high = Float64(bracket.T_high)
    while high - low > memo.config.bracket_tol_MeV + 1e-10
        mid = 0.5 * (low + high)
        curve = _slice!(memo, mid, nothing)
        if _valid_low_side(curve)
            low = mid
        else
            high = mid
        end
    end
    return (T_low=low, T_high=high, T_CEP=0.5 * (low + high))
end

function _run_job(config::PilotConfig)
    started = time_ns()
    canonical = _canonical_cep(config.xi)
    memo = PilotMemo(config)
    bracket = _refine_bracket!(memo, _find_bracket!(memo, canonical.T))
    cep = bracket === nothing ?
        (found=false, T_CEP_MeV=NaN, muq_CEP_MeV=NaN, T_bracket_low_MeV=NaN, T_bracket_high_MeV=NaN, bracket_width_T_MeV=NaN) :
        begin
            low_curve = memo.slice_cache[bracket.T_low]
            muq = low_curve.mu_transition
            (found=true, T_CEP_MeV=bracket.T_CEP, muq_CEP_MeV=Float64(muq), T_bracket_low_MeV=bracket.T_low, T_bracket_high_MeV=bracket.T_high, bracket_width_T_MeV=bracket.T_high - bracket.T_low)
        end

    oracle_refine = nothing
    if config.method == :high_resolution_oracle && bracket !== nothing
        coarse_low = memo.slice_cache[bracket.T_low]
        coarse_high = memo.slice_cache[bracket.T_high]
        low = _baseline_slice!(memo, bracket.T_low; step=config.oracle_refine_step)
        high = _baseline_slice!(memo, bracket.T_high; step=config.oracle_refine_step)
        oracle_refine = (
            coarse_low_status=String(coarse_low.status),
            coarse_high_status=String(coarse_high.status),
            low_status=String(low.status),
            high_status=String(high.status),
            low_mu_transition=low.mu_transition,
            high_mu_transition=high.mu_transition,
            stable=(coarse_low.status == low.status && coarse_high.status == high.status),
        )
    end

    rows = NamedTuple[]
    for T in sort(collect(keys(memo.slice_cache)))
        curve = memo.slice_cache[T]
        push!(rows, (
            xi=config.xi,
            method=String(config.method),
            T_MeV=T,
            requests=Int(curve.requests),
            unique_solve=Int(curve.unique_solves),
            targeted_additions=Int(curve.targeted_additions),
            cache_hits=Int(curve.cache_hits),
            residual_calls=Int(curve.residual_calls),
            jacobian_calls=Int(curve.jacobian_calls),
            newton_iterations=Int(curve.newton_iterations),
            trust_region_iterations=Int(curve.trust_region_iterations),
            fallback_count=Int(curve.fallback_count),
            governed_rescue_count=Int(curve.governed_rescue_count),
            retry_count=Int(curve.retry_count),
            failure_count=count(point -> !point.converged, curve.points),
            status=String(curve.status),
            mu_transition_MeV=Float64(curve.mu_transition),
            wall_seconds=Float64(curve.wall_seconds),
        ))
    end

    total_seconds = (time_ns() - started) / 1e9
    snapshot = Models.solver_work_snapshot(memo.telemetry)
    costs = [(
        xi=config.xi,
        method=String(config.method),
        calculation_sha=config.calculation_sha,
        equilibrium_requests=snapshot.equilibrium_requests,
        fixedrho_requests=snapshot.fixedrho_requests,
        unique_solves=memo.unique_solves,
        requested_point_calls=memo.point_requests,
        uncached_equivalent_requests=memo.point_requests,
        cache_hits=memo.cache_hits,
        targeted_additions=memo.targeted_additions,
        residual_calls=snapshot.nlsolve_f_calls + snapshot.postprocess_residual_calls,
        jacobian_calls=snapshot.nlsolve_g_calls,
        newton_iterations=snapshot.newton_iterations,
        trust_region_iterations=snapshot.trust_region_iterations,
        nonconverged_attempts=snapshot.nonconverged_attempts,
        fallback_count=snapshot.root_fallbacks,
        governed_rescue_count=snapshot.governed_rescues,
        retry_count=memo.scan_retries,
        exception_count=snapshot.exceptions,
        runner_seconds=total_seconds,
    )]
    accuracy = [(
        xi=config.xi,
        method=String(config.method),
        canonical_T_CEP_MeV=canonical.T,
        canonical_muq_CEP_MeV=canonical.muq,
        estimated_T_CEP_MeV=Float64(cep.T_CEP_MeV),
        estimated_muq_CEP_MeV=Float64(cep.muq_CEP_MeV),
        delta_T_MeV=Float64(cep.T_CEP_MeV) - canonical.T,
        delta_muq_MeV=Float64(cep.muq_CEP_MeV) - canonical.muq,
        bracket_low_MeV=Float64(cep.T_bracket_low_MeV),
        bracket_high_MeV=Float64(cep.T_bracket_high_MeV),
        bracket_width_MeV=Float64(cep.bracket_width_T_MeV),
        found=Bool(cep.found),
        oracle_refine=oracle_refine === nothing ? "" : JSON3.write(_json_safe(oracle_refine)),
        oracle_refine_stable=oracle_refine === nothing ? false : Bool(oracle_refine.stable),
    )]

    mkpath(config.output_dir)
    CSV.write(joinpath(config.output_dir, "slice_metrics.csv"), rows)
    CSV.write(joinpath(config.output_dir, "method_costs.csv"), costs)
    CSV.write(joinpath(config.output_dir, "cep_accuracy.csv"), accuracy)
    summary = Dict(
        "schema_version" => "cep_narrow_pilot_v1",
        "xi" => config.xi,
        "method" => String(config.method),
        "tag" => config.tag,
        "calculation_sha" => config.calculation_sha,
        "canonical_cep_path" => relpath(CANONICAL_CEP_PATH, PROJECT_ROOT),
        "canonical_cep_sha256" => bytes2hex(sha256(read(CANONICAL_CEP_PATH))),
        "parameters" => Dict(string(field) => getfield(config, field) for field in fieldnames(PilotConfig)),
        "cep" => Dict(string(field) => getproperty(cep, field) for field in propertynames(cep)),
        "telemetry" => Dict(string(field) => getproperty(snapshot, field) for field in propertynames(snapshot)),
        "cache_hits" => memo.cache_hits,
        "unique_solves" => memo.unique_solves,
        "targeted_additions" => memo.targeted_additions,
        "scan_retries" => memo.scan_retries,
        "runner_seconds" => total_seconds,
        "finite_and_converged" => all(point -> point.converged && point.finite, values(memo.point_cache)),
        "slice_count" => length(memo.slice_cache),
    )
    open(joinpath(config.output_dir, "job_summary.json"), "w") do io
        write(io, JSON3.write(_json_safe(summary)))
        write(io, '\n')
    end
    println(JSON3.write(_json_safe(summary)))
    return summary
end

function main(args=ARGS)
    config = parse_config(args)
    _run_job(config)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
