#!/usr/bin/env julia

"""Run one PNJL CEP narrow-window pilot v2 job.

The v2 runner is intentionally diagnostic-only.  It reuses the v1 pilot's
request-scoped solver telemetry and exact point memoization, but performs two
rho evidence levels and reports the three-state CEP contract.  It never
writes `data/reference/pnjl` or invokes the production phase pipeline.
"""

using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..")))

using CSV
using Dates
using JSON3
using SHA
using Statistics

const PROJECT_ROOT_V2 = normpath(joinpath(@__DIR__, "..", ".."))
const HBARC_MEV_FM_V2 = 197.3269804
const V2_CANONICAL_CEP_PATH = joinpath(
    PROJECT_ROOT_V2, "data", "reference", "pnjl", "issue130_phase_reference_v2",
    "accepted", "tables", "cep_boundary_accepted_phase_map_v1.csv",
)
const V2_STATUS_VALUES = (:confirmed_first_order, :confirmed_monotone, :ambiguous_near_critical)

# The v1 script supplies Models, Criticality, SolverWorkTelemetry, the stable
# solver wrapper, JSON NaN conversion and rho-grid helpers.  Its main guard is
# false when this file is the entrypoint, so including it has no side effects.
if !isdefined(Main, :PilotConfig)
    include(joinpath(@__DIR__, "pnjl_cep_narrow_pilot.jl"))
end

Base.@kwdef struct PilotV2Config
    xi::Float64 = 0.0
    method::Symbol = :rho_support_cascade
    stage::Symbol = :validation
    output_dir::String = joinpath(PROJECT_ROOT_V2, "data", "outputs", "pnjl_cep_narrow_pilot_v2")
    tag::String = "cep_narrow_pilot_v2"
    calculation_sha::String = ""
    window_file::String = ""
    p_num::Int = 24
    t_num::Int = 8
    thermo_quadrature_policy::Symbol = :rs_reduced_adaptive
    thermo_quadrature_rtol::Float64 = 1e-8
    thermo_quadrature_atol::Float64 = 1e-10
    thermo_quadrature_maxevals::Int = 10^7
    iterations::Int = 80
    rho_min::Float64 = 0.0
    rho_max::Float64 = 4.0
    cascade_coarse_step::Float64 = 0.05
    cascade_fine_step::Float64 = 0.025
    dense_coarse_step::Float64 = 0.0125
    dense_fine_step::Float64 = 0.00625
    oracle_coarse_step::Float64 = 0.00625
    oracle_fine_step::Float64 = 0.003125
    initial_window_MeV::Float64 = 4.0
    low_extension_MeV::Float64 = 8.0
    regular_high_extension_MeV::Float64 = 8.0
    xi05_high_extension_MeV::Float64 = 32.0
    expansion_step_MeV::Float64 = 4.0
    initial_T_step_MeV::Float64 = 1.0
    temperature_resolution_target_MeV::Float64 = 0.125
    max_bisect_iter::Int = 20
    unknown_budget::Int = 5
    targeted_max_points::Int = 12
    rho_position_tol_MeV::Float64 = 0.025
    rho_density_tol::Float64 = 0.0025
    rho_maxwell_area_tol::Float64 = 5e-5
    area_tol_good::Float64 = 1e-4
    area_tol_bad::Float64 = 5e-4
end

mutable struct PilotV2Memo
    base::PilotMemo
    config::PilotV2Config
    point_rows::Dict{Tuple{Float64, Float64, Int}, NamedTuple}
    slice_cache::Dict{Float64, NamedTuple}
    targeted_by_slice::Dict{Float64, Set{Float64}}
end

function _v2_base_config(config::PilotV2Config)
    method = config.method == :rho_support_cascade ? :rho_support_cascade :
        config.method == :high_resolution_oracle ? :high_resolution_oracle : :c2_dense_baseline
    return PilotConfig(
        xi=config.xi,
        method=method,
        output_dir=config.output_dir,
        tag=config.tag,
        calculation_sha=config.calculation_sha,
        p_num=config.p_num,
        t_num=config.t_num,
        thermo_quadrature_policy=config.thermo_quadrature_policy,
        thermo_quadrature_rtol=config.thermo_quadrature_rtol,
        thermo_quadrature_atol=config.thermo_quadrature_atol,
        thermo_quadrature_maxevals=config.thermo_quadrature_maxevals,
        iterations=config.iterations,
        rho_min=config.rho_min,
        rho_max=config.rho_max,
        targeted_max_points=config.targeted_max_points,
    )
end

function PilotV2Memo(config::PilotV2Config)
    return PilotV2Memo(
        PilotMemo(_v2_base_config(config)),
        config,
        Dict{Tuple{Float64, Float64, Int}, NamedTuple}(),
        Dict{Float64, NamedTuple}(),
        Dict{Float64, Set{Float64}}(),
    )
end

@inline _v2_symbol(value) = value isa Symbol ? value : Symbol(value)

function _v2_validate_config(config::PilotV2Config)
    config.method in (:c2_dense_baseline, :rho_support_cascade, :high_resolution_oracle) ||
        throw(ArgumentError("method must be c2_dense_baseline, rho_support_cascade or high_resolution_oracle"))
    config.stage in (:discovery, :validation) || throw(ArgumentError("stage must be discovery or validation"))
    isfinite(config.xi) || throw(ArgumentError("xi must be finite"))
    config.temperature_resolution_target_MeV > 0 || throw(ArgumentError("temperature resolution target must be positive"))
    config.max_bisect_iter > 0 || throw(ArgumentError("max_bisect_iter must be positive"))
    config.unknown_budget >= 0 || throw(ArgumentError("unknown_budget must be nonnegative"))
    config.targeted_max_points >= 0 || throw(ArgumentError("targeted_max_points must be nonnegative"))
    config.rho_position_tol_MeV > 0 && config.rho_density_tol > 0 && config.rho_maxwell_area_tol > 0 ||
        throw(ArgumentError("rho geometry tolerances must be positive"))
    isempty(config.calculation_sha) || occursin(r"^[0-9a-fA-F]{40}$", config.calculation_sha) ||
        throw(ArgumentError("calculation-sha must be an immutable 40-character SHA"))
    return config
end

function _v2_read_canonical(xi::Float64)
    isfile(V2_CANONICAL_CEP_PATH) || throw(ArgumentError("canonical CEP file missing: $(V2_CANONICAL_CEP_PATH)"))
    rows = collect(CSV.File(V2_CANONICAL_CEP_PATH))
    candidates = [row for row in rows if hasproperty(row, :xi) && abs(Float64(row.xi) - xi) <= 1e-9]
    isempty(candidates) && throw(ArgumentError("canonical CEP has no xi=$(xi) row"))
    row = first(candidates)
    if hasproperty(row, :T_CEP_MeV)
        return (T=Float64(row.T_CEP_MeV), muq=Float64(row.muq_CEP_MeV))
    end
    return (T=Float64(row.T_midpoint_MeV), muq=Float64(row.mu_CEP_proxy_MeV))
end

function _v2_json_read(path::AbstractString)
    isempty(path) && return nothing
    isfile(path) || throw(ArgumentError("window file does not exist: $(path)"))
    return JSON3.read(read(path, String))
end

function _v2_float(value, default=NaN)
    value === nothing && return default
    value isa AbstractString && isempty(value) && return default
    return Float64(value)
end

function _v2_window(config::PilotV2Config, canonical)
    raw = _v2_json_read(config.window_file)
    raw === nothing && return (
        T_min=canonical.T - config.low_extension_MeV,
        T_max=canonical.T + (abs(config.xi - 0.5) <= 1e-9 ? config.xi05_high_extension_MeV : config.regular_high_extension_MeV),
        required_anchors=Float64[],
        source="local_default",
    )
    rows = hasproperty(raw, :windows) ? raw.windows : raw
    selected = nothing
    for row in rows
        abs(Float64(row.xi) - config.xi) <= 1e-9 && (selected = row; break)
    end
    selected === nothing && throw(ArgumentError("window file has no xi=$(config.xi) row"))
    anchors = hasproperty(selected, :required_T_anchors) ? Float64.(collect(selected.required_T_anchors)) : Float64[]
    return (
        T_min=Float64(selected.T_min_MeV),
        T_max=Float64(selected.T_max_MeV),
        required_anchors=sort(unique(anchors)),
        source="frozen_validation_window",
    )
end

function _v2_prior(memo::PilotV2Memo, T::Float64)
    candidates = NamedTuple[]
    for record in values(memo.slice_cache)
        center = record.spinodal_rho_center
        gap = record.spinodal_rho_gap
        center === nothing && continue
        gap === nothing && continue
        isfinite(Float64(center)) && isfinite(Float64(gap)) && Float64(gap) > 0 || continue
        push!(candidates, record)
    end
    isempty(candidates) && return nothing
    previous = argmin(record -> abs(Float64(record.T_MeV) - T), candidates)
    return Criticality.RhoSupportPrior(
        Float64(previous.T_MeV),
        Float64(previous.spinodal_rho_center),
        Float64(previous.spinodal_rho_gap),
    )
end

function _v2_point!(memo::PilotV2Memo, T::Float64, rho::Float64, level::Int; targeted::Bool=false)
    set = get!(memo.targeted_by_slice, T, Set{Float64}())
    key = (T, rho, level)
    before_cache = haskey(memo.base.point_cache, (T, memo.config.xi, rho))
    point = _solve_point!(memo.base, T, rho; targeted=targeted)
    role = targeted ? "targeted" : "base"
    if haskey(memo.point_rows, key)
        old = memo.point_rows[key]
        previous_role = old.sampling_role == "targeted" ? "targeted" : role
        memo.point_rows[key] = merge(old, (sampling_role=previous_role,))
    else
        memo.point_rows[key] = (
            xi=memo.config.xi,
            method=String(memo.config.method),
            T_MeV=T,
            rho_level=level,
            rho=rho,
            muq_MeV=Float64(point.muq_MeV),
            pressure_fm4=Float64(point.pressure_fm4),
            residual_norm=Float64(point.residual_norm),
            iterations=Int(point.iterations),
            converged=Bool(point.converged),
            finite=Bool(point.finite),
            sampling_role=role,
            cache_reused=before_cache,
            solution_json=JSON3.write(_json_safe(point.solution)),
        )
    end
    targeted && push!(set, rho)
    return point
end

function _v2_level_result(memo::PilotV2Memo, T::Float64, rho::Vector{Float64}, level::Int, targeted_values::Vector{Float64})
    points = [_v2_point!(memo, T, value, level) for value in rho]
    mu = [Float64(point.muq_MeV) for point in points]
    all_converged = all(point -> point.converged && point.finite, points)
    if !all_converged
        return (
            status=:failed,
            mu_transition=nothing,
            area_residual=nothing,
            sres=Models.SShapeResult(),
            rho_hadron=nothing,
            rho_quark=nothing,
            reason="solver_or_curve_failure",
            curve_rho=rho,
            curve_mu=mu,
            points=points,
            s_shape_status=:unknown,
            maxwell_status=:not_run,
            all_converged=false,
            targeted_values=targeted_values,
            level=level,
        )
    end

    cres = Models._classify_s_curve(
        mu,
        rho;
        area_tol_good=memo.config.area_tol_good,
        area_tol_bad=memo.config.area_tol_bad,
    )
    maxwell = cres.sres.has_s_shape ? Models.maxwell_construction(mu, rho; spinodal_hint=cres.sres) : nothing
    valid_maxwell = cres.status == :valid && maxwell !== nothing && maxwell.converged && maxwell.mu_transition !== nothing
    return (
        status=Symbol(cres.status),
        mu_transition=cres.mu_transition,
        area_residual=cres.area_residual,
        sres=cres.sres,
        rho_hadron=valid_maxwell ? maxwell.rho_hadron : nothing,
        rho_quark=valid_maxwell ? maxwell.rho_quark : nothing,
        reason=String(cres.reason),
        curve_rho=rho,
        curve_mu=mu,
        points=points,
        s_shape_status=cres.sres.has_s_shape ?
            (cres.status == :weak_s_shape ? :weak : :present) : :none,
        maxwell_status=valid_maxwell ? :converged :
            (cres.sres.has_s_shape ? :failed : :not_applicable),
        all_converged=true,
        targeted_values=targeted_values,
        level=level,
    )
end

function _v2_target_config(remaining::Int)
    remaining < 5 && return nothing
    count = min(9, remaining)
    iseven(count) && (count -= 1)
    count < 5 && return nothing
    return Criticality.RhoSupportConfig(target_point_count=count, max_extra_points=count)
end

function _v2_cascade_level!(memo::PilotV2Memo, T::Float64, level::Int, step::Float64, prior)
    rho = _rho_grid(memo.config.rho_min, memo.config.rho_max, step)
    targeted_values = Float64[]
    total_targeted = length(get!(memo.targeted_by_slice, T, Set{Float64}()))
    remaining = memo.config.targeted_max_points - total_targeted
    cascade_cfg = _v2_target_config(remaining)
    sampler = nothing
    if cascade_cfg !== nothing
        sampler = value -> begin
            value = Float64(value)
            push!(targeted_values, value)
            targeted_set = get!(memo.targeted_by_slice, T, Set{Float64}())
            is_new = !(value in targeted_set) && length(targeted_set) < memo.config.targeted_max_points
            _v2_point!(memo, T, value, level; targeted=is_new).muq_MeV
        end
    end
    base_points = [_v2_point!(memo, T, value, level) for value in rho]
    base_mu = [Float64(point.muq_MeV) for point in base_points]
    assessment = if sampler === nothing
        nothing
    else
        Criticality.analyze_curve_cascade(
            rho,
            base_mu;
            sample_mu=sampler,
            prior=prior,
            config=cascade_cfg,
        )
    end
    all_rho = sort(unique(vcat(rho, targeted_values)))
    result = _v2_level_result(memo, T, all_rho, level, unique(targeted_values))
    return merge(result, (
        cascade_status=assessment === nothing ? :not_run : assessment.status,
        cascade_reason=assessment === nothing ? "target_budget_unavailable" : assessment.reason,
        spinodal_rho_center=assessment === nothing ? nothing : assessment.spinodal_rho_center,
        spinodal_rho_gap=assessment === nothing ? nothing : assessment.spinodal_rho_gap,
        base_point_count=length(rho),
    ))
end

function _v2_dense_level!(memo::PilotV2Memo, T::Float64, level::Int, step::Float64)
    rho = _rho_grid(memo.config.rho_min, memo.config.rho_max, step)
    result = _v2_level_result(memo, T, rho, level, Float64[])
    return merge(result, (
        cascade_status=:not_applicable,
        cascade_reason="dense_grid",
        spinodal_rho_center=nothing,
        spinodal_rho_gap=nothing,
        base_point_count=length(rho),
    ))
end

function _v2_slice!(memo::PilotV2Memo, T::Float64)
    haskey(memo.slice_cache, T) && return memo.slice_cache[T]
    started = time_ns()
    before_requests = memo.base.point_requests
    before_unique = memo.base.unique_solves
    before_hits = memo.base.cache_hits
    before_targeted = memo.base.targeted_additions
    before_retries = memo.base.scan_retries
    before_telemetry = Models.solver_work_snapshot(memo.base.telemetry)
    prior = memo.config.method == :rho_support_cascade ? _v2_prior(memo, T) : nothing
    levels = if memo.config.method == :rho_support_cascade
        (
            _v2_cascade_level!(memo, T, 0, memo.config.cascade_coarse_step, prior),
            _v2_cascade_level!(memo, T, 1, memo.config.cascade_fine_step, prior),
        )
    elseif memo.config.method == :high_resolution_oracle
        (
            _v2_dense_level!(memo, T, 0, memo.config.oracle_coarse_step),
            _v2_dense_level!(memo, T, 1, memo.config.oracle_fine_step),
        )
    else
        (
            _v2_dense_level!(memo, T, 0, memo.config.dense_coarse_step),
            _v2_dense_level!(memo, T, 1, memo.config.dense_fine_step),
        )
    end
    coarse, fine = levels
    geometry = if coarse.status == :invalid && fine.status == :invalid &&
                     coarse.reason == "no_s_shape" && fine.reason == "no_s_shape"
        Models.PhaseGeometryError(
            comparable=true,
            converged=true,
            position_MeV=0.0,
            density=0.0,
            maxwell_area=0.0,
            response_rtol=0.0,
            reason="stable_no_s_shape",
        )
    else
        Models._compare_phase_geometry(
            coarse,
            fine,
            Models.PhaseGeometryTolerances(
                position_MeV=memo.config.rho_position_tol_MeV,
                density=memo.config.rho_density_tol,
                maxwell_area=memo.config.rho_maxwell_area_tol,
            ),
        )
    end
    status = if coarse.status == :valid && fine.status == :valid && geometry.converged
        :confirmed_first_order
    elseif coarse.status == :invalid && fine.status == :invalid &&
           coarse.reason == "no_s_shape" && fine.reason == "no_s_shape" && geometry.converged
        :confirmed_monotone
    else
        :ambiguous_near_critical
    end
    candidate_state(level) = level.status == :valid ? :first_order_candidate :
        (level.status == :invalid && level.reason == "no_s_shape" ? :monotone_candidate : :ambiguous_candidate)
    after_telemetry = Models.solver_work_snapshot(memo.base.telemetry)
    snapshot_diff = (
        requests=memo.base.point_requests - before_requests,
        unique_solves=memo.base.unique_solves - before_unique,
        cache_hits=memo.base.cache_hits - before_hits,
        targeted_additions=memo.base.targeted_additions - before_targeted,
        residual_calls=(after_telemetry.nlsolve_f_calls + after_telemetry.postprocess_residual_calls) -
            (before_telemetry.nlsolve_f_calls + before_telemetry.postprocess_residual_calls),
        jacobian_calls=after_telemetry.nlsolve_g_calls - before_telemetry.nlsolve_g_calls,
        newton_iterations=after_telemetry.newton_iterations - before_telemetry.newton_iterations,
        trust_region_iterations=after_telemetry.trust_region_iterations - before_telemetry.trust_region_iterations,
        fallback_count=after_telemetry.root_fallbacks - before_telemetry.root_fallbacks,
        governed_rescue_count=after_telemetry.governed_rescues - before_telemetry.governed_rescues,
        retry_count=memo.base.scan_retries - before_retries,
    )
    fine_failed = count(point -> !(point.converged && point.finite), fine.points)
    record = (
        xi=memo.config.xi,
        method=String(memo.config.method),
        T_MeV=T,
        coarse_rho_step=memo.config.method == :rho_support_cascade ? memo.config.cascade_coarse_step :
            memo.config.method == :high_resolution_oracle ? memo.config.oracle_coarse_step : memo.config.dense_coarse_step,
        fine_rho_step=memo.config.method == :rho_support_cascade ? memo.config.cascade_fine_step :
            memo.config.method == :high_resolution_oracle ? memo.config.oracle_fine_step : memo.config.dense_fine_step,
        coarse_status=String(coarse.status),
        fine_status=String(fine.status),
        coarse_candidate_status=String(candidate_state(coarse)),
        fine_candidate_status=String(candidate_state(fine)),
        coarse_s_shape_status=String(coarse.s_shape_status),
        fine_s_shape_status=String(fine.s_shape_status),
        coarse_maxwell_status=String(coarse.maxwell_status),
        fine_maxwell_status=String(fine.maxwell_status),
        coarse_reason=String(coarse.reason),
        fine_reason=String(fine.reason),
        result_status=String(status),
        geometry_converged=Bool(geometry.converged),
        geometry_position_error_MeV=Float64(geometry.position_MeV),
        geometry_density_error=Float64(geometry.density),
        geometry_maxwell_area=Float64(geometry.maxwell_area),
        targeted_additions=snapshot_diff.targeted_additions,
        failure_count=fine_failed,
        requests=snapshot_diff.requests,
        unique_solves=snapshot_diff.unique_solves,
        cache_hits=snapshot_diff.cache_hits,
        residual_calls=snapshot_diff.residual_calls,
        jacobian_calls=snapshot_diff.jacobian_calls,
        newton_iterations=snapshot_diff.newton_iterations,
        trust_region_iterations=snapshot_diff.trust_region_iterations,
        fallback_count=snapshot_diff.fallback_count,
        governed_rescue_count=snapshot_diff.governed_rescue_count,
        retry_count=snapshot_diff.retry_count,
        wall_seconds=(time_ns() - started) / 1e9,
        mu_last_first_order_MeV=(status == :confirmed_first_order ? Float64(something(fine.mu_transition, NaN)) : NaN),
        spinodal_rho_center=fine.spinodal_rho_center,
        spinodal_rho_gap=fine.spinodal_rho_gap,
        coarse_level=coarse,
        fine_level=fine,
    )
    memo.slice_cache[T] = record
    return record
end

function _v2_required_temperatures(config::PilotV2Config, canonical, window)
    isempty(window.required_anchors) || return window.required_anchors
    low = canonical.T - config.initial_window_MeV
    high = canonical.T + config.initial_window_MeV
    values = collect(low:config.initial_T_step_MeV:high)
    push!(values, canonical.T - config.low_extension_MeV)
    cursor = high
    while cursor < window.T_max - 1e-9
        cursor = min(cursor + config.expansion_step_MeV, window.T_max)
        push!(values, cursor)
    end
    return sort(unique(Float64[clamp(value, window.T_min, window.T_max) for value in values]))
end

function _v2_initial_temperatures(config::PilotV2Config, canonical, window)
    values = _v2_required_temperatures(config, canonical, window)
    low = canonical.T - config.initial_window_MeV
    high = canonical.T + config.initial_window_MeV
    for value in low:config.initial_T_step_MeV:high
        push!(values, Float64(value))
    end
    return sort(unique(Float64[clamp(value, window.T_min, window.T_max) for value in values]))
end

@inline function _v2_in_window(T::Float64, window)
    return !(T < window.T_min - 1e-8 || T > window.T_max + 1e-8)
end

function _v2_evaluate!(memo::PilotV2Memo, T::Float64, window)
    if !_v2_in_window(T, window)
        return nothing
    end
    return _v2_slice!(memo, T)
end

function _v2_unknown_count(memo::PilotV2Memo)
    return count(record ->
        record.coarse_status in ("failed", "unknown") || record.fine_status in ("failed", "unknown"),
        values(memo.slice_cache))
end

function _v2_frontier_pair(memo::PilotV2Memo, state::Symbol)
    records = collect(values(memo.slice_cache))
    if state == :first_order
        lows = [record for record in records if record.result_status == "confirmed_first_order"]
        isempty(lows) && return nothing
        low = argmax(record -> record.T_MeV, lows)
        higher = [record for record in records if record.T_MeV > low.T_MeV && record.result_status != "confirmed_first_order"]
        high = isempty(higher) ? nothing : argmin(record -> record.T_MeV, higher)
        return low.T_MeV, high === nothing ? nothing : high.T_MeV
    end
    highs = [record for record in records if record.result_status == "confirmed_monotone"]
    isempty(highs) && return nothing
    high = argmin(record -> record.T_MeV, highs)
    lower = [record for record in records if record.T_MeV < high.T_MeV && record.result_status != "confirmed_monotone"]
    low = isempty(lower) ? nothing : argmax(record -> record.T_MeV, lower)
    return low === nothing ? nothing : (low.T_MeV, high.T_MeV)
end

function _v2_search!(memo::PilotV2Memo, canonical, window)
    for T in _v2_initial_temperatures(memo.config, canonical, window)
        _v2_evaluate!(memo, T, window)
    end
    # Ensure the low side is represented when a canonical-centered initial
    # probe does not contain a confirmed first-order slice.
    if _v2_frontier_pair(memo, :first_order) === nothing
        _v2_evaluate!(memo, clamp(canonical.T - memo.config.low_extension_MeV, window.T_min, window.T_max), window)
    end
    for _ in 1:memo.config.max_bisect_iter
        _v2_unknown_count(memo) > memo.config.unknown_budget && break
        changed = false
        first_pair = _v2_frontier_pair(memo, :first_order)
        if first_pair !== nothing
            low, high = first_pair
            high === nothing && (high = window.T_max)
            if high - low > memo.config.temperature_resolution_target_MeV + 1e-10
                _v2_evaluate!(memo, 0.5 * (low + high), window)
                changed = true
            end
        end
        monotone_pair = _v2_frontier_pair(memo, :monotone)
        if monotone_pair !== nothing
            low, high = monotone_pair
            if high - low > memo.config.temperature_resolution_target_MeV + 1e-10
                _v2_evaluate!(memo, 0.5 * (low + high), window)
                changed = true
            end
        end
        changed || break
    end
    return memo
end

function _v2_candidate_endpoints(memo::PilotV2Memo, level::Symbol)
    records = sort(collect(values(memo.slice_cache)); by=record -> record.T_MeV)
    state_name = level == :coarse ? :coarse_candidate_status : :fine_candidate_status
    first_order = [record for record in records if getproperty(record, state_name) == "first_order_candidate"]
    monotone = [record for record in records if getproperty(record, state_name) == "monotone_candidate"]
    last_first = isempty(first_order) ? nothing : argmax(record -> record.T_MeV, first_order)
    first_mono = if last_first === nothing
        isempty(monotone) ? nothing : argmin(record -> record.T_MeV, monotone)
    else
        candidates = [record for record in monotone if record.T_MeV > last_first.T_MeV]
        isempty(candidates) ? nothing : argmin(record -> record.T_MeV, candidates)
    end
    endpoint_level = level == :coarse ? :coarse_level : :fine_level
    return (
        T_last_first_order_MeV=last_first === nothing ? NaN : Float64(last_first.T_MeV),
        muq_last_first_order_MeV=last_first === nothing ? NaN : Float64(something(getproperty(last_first, endpoint_level).mu_transition, NaN)),
        T_first_monotone_MeV=first_mono === nothing ? NaN : Float64(first_mono.T_MeV),
    )
end

function _v2_endpoint_record(memo::PilotV2Memo)
    records = sort(collect(values(memo.slice_cache)); by=record -> record.T_MeV)
    first_order = [record for record in records if record.result_status == "confirmed_first_order"]
    monotone = [record for record in records if record.result_status == "confirmed_monotone"]
    last_first = isempty(first_order) ? nothing : argmax(record -> record.T_MeV, first_order)
    first_mono = if last_first === nothing
        isempty(monotone) ? nothing : argmin(record -> record.T_MeV, monotone)
    else
        candidates = [record for record in monotone if record.T_MeV > last_first.T_MeV]
        isempty(candidates) ? nothing : argmin(record -> record.T_MeV, candidates)
    end
    result_status = if last_first !== nothing && first_mono !== nothing
        :ambiguous
    elseif last_first === nothing && first_mono !== nothing
        :not_found
    else
        :ambiguous
    end
    low_T = last_first === nothing ? NaN : Float64(last_first.T_MeV)
    high_T = first_mono === nothing ? NaN : Float64(first_mono.T_MeV)
    low_mu = last_first === nothing ? NaN : Float64(last_first.mu_last_first_order_MeV)
    coarse_endpoints = _v2_candidate_endpoints(memo, :coarse)
    fine_endpoints = _v2_candidate_endpoints(memo, :fine)
    return (
        result_status=result_status,
        found=false,
        T_cep_MeV=NaN,
        muq_CEP_MeV=NaN,
        T_last_first_order_MeV=low_T,
        muq_last_first_order_MeV=low_mu,
        T_first_monotone_MeV=high_T,
        ambiguity_width_T_MeV=(isfinite(low_T) && isfinite(high_T) ? high_T - low_T : NaN),
        temperature_resolution_target_MeV=memo.config.temperature_resolution_target_MeV,
        first_order_count=length(first_order),
        monotone_count=length(monotone),
        evaluated_temperature_count=length(records),
        coarse_T_last_first_order_MeV=coarse_endpoints.T_last_first_order_MeV,
        coarse_muq_last_first_order_MeV=coarse_endpoints.muq_last_first_order_MeV,
        coarse_T_first_monotone_MeV=coarse_endpoints.T_first_monotone_MeV,
        fine_T_last_first_order_MeV=fine_endpoints.T_last_first_order_MeV,
        fine_muq_last_first_order_MeV=fine_endpoints.muq_last_first_order_MeV,
        fine_T_first_monotone_MeV=fine_endpoints.T_first_monotone_MeV,
    )
end

@inline _v2_csv_value(_, value) = value === nothing ? missing : value

function _v2_write_outputs(memo::PilotV2Memo, canonical, window, endpoints, started_ns::UInt64)
    config = memo.config
    mkpath(config.output_dir)
    point_rows = collect(values(memo.point_rows))
    sort!(point_rows; by=row -> (row.T_MeV, row.rho_level, row.rho))
    slice_rows = NamedTuple[]
    for record in sort(collect(values(memo.slice_cache)); by=record -> record.T_MeV)
        push!(slice_rows, (; (key => value for (key, value) in pairs(record) if key ∉ (:coarse_level, :fine_level))...))
    end
    snapshot = Models.solver_work_snapshot(memo.base.telemetry)
    total_seconds = (time_ns() - started_ns) / 1e9
    costs = [(
        xi=config.xi,
        method=String(config.method),
        calculation_sha=config.calculation_sha,
        equilibrium_requests=snapshot.equilibrium_requests,
        fixedrho_requests=snapshot.fixedrho_requests,
        unique_solves=memo.base.unique_solves,
        requested_point_calls=memo.base.point_requests,
        uncached_equivalent_requests=memo.base.point_requests,
        cache_hits=memo.base.cache_hits,
        targeted_additions=memo.base.targeted_additions,
        residual_calls=snapshot.nlsolve_f_calls + snapshot.postprocess_residual_calls,
        jacobian_calls=snapshot.nlsolve_g_calls,
        newton_iterations=snapshot.newton_iterations,
        trust_region_iterations=snapshot.trust_region_iterations,
        nonconverged_attempts=snapshot.nonconverged_attempts,
        fallback_count=snapshot.root_fallbacks,
        governed_rescue_count=snapshot.governed_rescues,
        retry_count=memo.base.scan_retries,
        exception_count=snapshot.exceptions,
        runner_seconds=total_seconds,
    )]
    accuracy = [(
        xi=config.xi,
        method=String(config.method),
        canonical_T_CEP_MeV=canonical.T,
        canonical_muq_CEP_MeV=canonical.muq,
        result_status=String(endpoints.result_status),
        found=Bool(endpoints.found),
        T_last_first_order_MeV=endpoints.T_last_first_order_MeV,
        muq_last_first_order_MeV=endpoints.muq_last_first_order_MeV,
        T_first_monotone_MeV=endpoints.T_first_monotone_MeV,
        ambiguity_width_T_MeV=endpoints.ambiguity_width_T_MeV,
        temperature_resolution_target_MeV=endpoints.temperature_resolution_target_MeV,
        delta_T_last_first_order_MeV=isfinite(endpoints.T_last_first_order_MeV) ? endpoints.T_last_first_order_MeV - canonical.T : NaN,
        delta_muq_last_first_order_MeV=isfinite(endpoints.muq_last_first_order_MeV) ? endpoints.muq_last_first_order_MeV - canonical.muq : NaN,
        evaluated_temperature_count=endpoints.evaluated_temperature_count,
        first_order_count=endpoints.first_order_count,
        monotone_count=endpoints.monotone_count,
        coarse_T_last_first_order_MeV=endpoints.coarse_T_last_first_order_MeV,
        coarse_muq_last_first_order_MeV=endpoints.coarse_muq_last_first_order_MeV,
        coarse_T_first_monotone_MeV=endpoints.coarse_T_first_monotone_MeV,
        fine_T_last_first_order_MeV=endpoints.fine_T_last_first_order_MeV,
        fine_muq_last_first_order_MeV=endpoints.fine_muq_last_first_order_MeV,
        fine_T_first_monotone_MeV=endpoints.fine_T_first_monotone_MeV,
    )]
    CSV.write(joinpath(config.output_dir, "curve_points.csv"), point_rows; transform=_v2_csv_value)
    CSV.write(joinpath(config.output_dir, "slice_metrics.csv"), slice_rows; transform=_v2_csv_value)
    CSV.write(joinpath(config.output_dir, "method_costs.csv"), costs; transform=_v2_csv_value)
    CSV.write(joinpath(config.output_dir, "cep_accuracy.csv"), accuracy; transform=_v2_csv_value)
    summary = Dict(
        "schema_version" => "cep_narrow_pilot_v2",
        "xi" => config.xi,
        "method" => String(config.method),
        "stage" => String(config.stage),
        "tag" => config.tag,
        "calculation_sha" => config.calculation_sha,
        "canonical_cep_path" => relpath(V2_CANONICAL_CEP_PATH, PROJECT_ROOT_V2),
        "canonical_cep_sha256" => bytes2hex(sha256(read(V2_CANONICAL_CEP_PATH))),
        "window" => Dict("T_min_MeV" => window.T_min, "T_max_MeV" => window.T_max, "source" => window.source),
        "parameters" => Dict(string(field) => getfield(config, field) for field in fieldnames(PilotV2Config)),
        "cep" => Dict(string(field) => getproperty(endpoints, field) for field in propertynames(endpoints)),
        "telemetry" => Dict(string(field) => getproperty(snapshot, field) for field in propertynames(snapshot)),
        "cache_hits" => memo.base.cache_hits,
        "unique_solves" => memo.base.unique_solves,
        "targeted_additions" => memo.base.targeted_additions,
        "scan_retries" => memo.base.scan_retries,
        "runner_seconds" => total_seconds,
        "finite_and_converged_final" => all(record -> record.fine_level.all_converged, values(memo.slice_cache)),
        "slice_count" => length(memo.slice_cache),
    )
    open(joinpath(config.output_dir, "job_summary.json"), "w") do io
        write(io, JSON3.write(_json_safe(summary)))
        write(io, '\n')
    end
    if config.stage == :discovery
        proposal = Dict(
            "schema_version" => "cep_narrow_pilot_v2_window_proposal",
            "xi" => config.xi,
            "canonical_T_CEP_MeV" => canonical.T,
            "canonical_muq_CEP_MeV" => canonical.muq,
            "T_min_MeV" => window.T_min,
            "T_max_MeV" => window.T_max,
            "required_T_anchors" => window.required_anchors,
            "cep" => Dict(string(field) => getproperty(endpoints, field) for field in propertynames(endpoints)),
            "calculation_sha" => config.calculation_sha,
        )
        open(joinpath(config.output_dir, "window_proposal.json"), "w") do io
            write(io, JSON3.write(_json_safe(proposal)))
            write(io, '\n')
        end
    end
    return summary
end

function _v2_parse_config(args)
    values = _parse_cli(args)
    method = Symbol(_get_string(values, "method", "rho_support_cascade"))
    stage = Symbol(_get_string(values, "stage", "validation"))
    window_file = _get_string(values, "window-file", "")
    config = PilotV2Config(
        xi=_get_float(values, "xi", 0.0),
        method=method,
        stage=stage,
        output_dir=abspath(_get_string(values, "output-dir", joinpath(PROJECT_ROOT_V2, "data", "outputs", "pnjl_cep_narrow_pilot_v2"))),
        tag=_get_string(values, "tag", "cep_narrow_pilot_v2"),
        calculation_sha=_get_string(values, "calculation-sha", ""),
        window_file=isempty(window_file) ? "" : abspath(window_file),
    )
    return _v2_validate_config(config)
end

function _v2_run(config::PilotV2Config)
    started = time_ns()
    canonical = _v2_read_canonical(config.xi)
    window = _v2_window(config, canonical)
    memo = PilotV2Memo(config)
    _v2_search!(memo, canonical, window)
    endpoints = _v2_endpoint_record(memo)
    summary = _v2_write_outputs(memo, canonical, window, endpoints, started)
    println(JSON3.write(_json_safe(summary)))
    return summary
end

function main(args=ARGS)
    _v2_run(_v2_parse_config(args))
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
