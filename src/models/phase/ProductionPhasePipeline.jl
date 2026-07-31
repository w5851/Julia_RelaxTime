function _production_temperature_grid(T_start::Float64, T_end::Float64, dT::Float64)
    dT > 0 || error("dT must be positive")
    T_end >= T_start || error("T_end must be >= T_start")
    temps = collect(T_start:dT:T_end)
    if isempty(temps) || temps[end] < T_end
        push!(temps, T_end)
    end
    return unique(Float64.(temps))
end

function _production_eval_filename(T::Float64, level::Int, pass::Int)
    t_token = replace(@sprintf("%.6f", T), "." => "p", "-" => "m")
    return "prod_eval_T$(t_token)_L$(level)_P$(pass).csv"
end

function _append_scan_csv!(aggregate_csv::String, source_csv::String)
    isfile(source_csv) || return
    lines = readlines(source_csv)
    isempty(lines) && return
    if !isfile(aggregate_csv)
        mkpath(dirname(aggregate_csv))
        open(aggregate_csv, "w") do io
            for line in lines
                println(io, line)
            end
        end
        return
    end
    open(aggregate_csv, "a") do io
        for line in lines[2:end]
            println(io, line)
        end
    end
end

function _run_single_temperature_production_scan(
        aggregate_csv::String,
        eval_dir::String,
        T_mid::Float64,
        rho_eval::Vector{Float64},
        xi::Float64,
        reverse_rho::Bool,
        seed_policy::Symbol,
        solver_backend::Symbol,
        model_kind::Symbol,
        p_num::Int,
        t_num::Int,
        thermo_quadrature_kwargs::NamedTuple,
        iterations::Int,
        level::Int,
        pass::Int,
        work_telemetry::Union{Nothing, SolverWorkTelemetry}=nothing)
    out_csv = joinpath(eval_dir, _production_eval_filename(T_mid, level, pass))
    stats = run_trho_scan(;
        T_values=[Float64(T_mid)],
        rho_values=rho_eval,
        xi_values=[Float64(xi)],
        output_path=out_csv,
        overwrite=true,
        resume=false,
        reverse_rho=reverse_rho,
        seed_policy=seed_policy,
        solver_backend=solver_backend,
        model_kind=model_kind,
        p_num=p_num,
        t_num=t_num,
        thermo_quadrature_kwargs...,
        iterations=iterations,
        work_telemetry=work_telemetry,
    )
    local_curves = load_curves_from_trho_csv(out_csv; xi=Float64(xi), min_points=3)
    if isempty(local_curves)
        return (curve=nothing, stats=stats, out_csv=out_csv)
    end
    local_key = first(sort(collect(keys(local_curves))))
    return (curve=local_curves[local_key], stats=stats, out_csv=out_csv)
end

function _rho_support_fine_grid(base_rho_grid::Vector{Float64}, fine_step::Float64)
    first_rho = first(base_rho_grid)
    last_rho = last(base_rho_grid)
    fine_step > 0 || throw(ArgumentError("rho_support_fine_step must be positive"))
    span = last_rho - first_rho
    count = span / fine_step
    isapprox(count, round(count); atol=1e-9, rtol=0.0) || throw(ArgumentError(
        "rho_support_fine_step must divide the rho span exactly",
    ))
    return Float64.(collect(range(first_rho, last_rho; length=Int(round(count)) + 1)))
end

function _validate_rho_support_grid(base_rho_grid::Vector{Float64}, fine_step::Float64)
    length(base_rho_grid) >= 2 || throw(ArgumentError("rho_support_cascade requires at least two rho points"))
    all(diff(base_rho_grid) .> 0) || throw(ArgumentError(
        "rho_support_cascade requires a strictly increasing coarse rho_grid",
    ))
    coarse_step = first(diff(base_rho_grid))
    all(isapprox(step, coarse_step; atol=1e-9, rtol=0.0) for step in diff(base_rho_grid)) ||
        throw(ArgumentError("rho_support_cascade requires a uniformly spaced coarse rho_grid"))
    fine_step < coarse_step || throw(ArgumentError(
        "rho_support_fine_step must be strictly smaller than the coarse rho step",
    ))
    ratio = coarse_step / fine_step
    isapprox(ratio, round(ratio); atol=1e-9, rtol=0.0) || throw(ArgumentError(
        "rho_support_fine_step must divide every coarse rho interval exactly",
    ))
    _rho_support_fine_grid(base_rho_grid, fine_step)
    return nothing
end

function _production_session_rows(session, T::Float64, xi::Float64)
    rows = NamedTuple[]
    for ((row_T, row_xi, _), row) in session.cache
        row_T == T && row_xi == xi && push!(rows, row)
    end
    sort!(rows; by=row -> row.rho)
    return rows
end

function _existing_session_keys(path::String)
    keys = Set{Tuple{Float64, Float64, Float64}}()
    isfile(path) || return keys
    for (line_number, line) in enumerate(readlines(path))
        line_number == 1 && continue
        fields = split(line, ','; limit=4)
        length(fields) >= 3 || continue
        try
            T = parse(Float64, fields[1])
            rho = parse(Float64, fields[2])
            xi = parse(Float64, fields[3])
            push!(keys, (T, xi, rho))
        catch
            # A malformed historical row is left for the normal validator;
            # it must not prevent the opt-in session from starting.
        end
    end
    return keys
end

function _production_session_curve(session, T::Float64, xi::Float64)
    rows = _production_session_rows(session, T, xi)
    isempty(rows) && return nothing, rows
    all(row -> row.converged && row.finite, rows) || return nothing, rows
    return ([row.muq_MeV for row in rows], [row.rho for row in rows]), rows
end

function _write_session_row!(io, written::Set{Tuple{Float64, Float64, Float64}}, row,
        T::Float64, xi::Float64, p_num::Int, t_num::Int, thermo_quadrature_kwargs::NamedTuple)
    key = (T, xi, row.rho)
    key in written && return false
    TrhoScan._write_row(io, T, row.rho, xi, row.result, row.message;
        p_num=p_num, t_num=t_num, thermo_quadrature_kwargs...)
    push!(written, key)
    return true
end

function _production_session_scan!(session, io, written, T::Float64, xi::Float64,
        rho_values::Vector{Float64}, reverse_rho::Bool, p_num::Int, t_num::Int,
        thermo_quadrature_kwargs::NamedTuple)
    order = reverse_rho ? reverse(rho_values) : rho_values
    for rho in order
        row = TrhoScan.rho_point!(session, T, xi, rho)
        _write_session_row!(io, written, row, T, xi, p_num, t_num, thermo_quadrature_kwargs)
    end
    return _production_session_curve(session, T, xi)
end

function _production_classify_temperature_cascade(
        aggregate_csv::String,
        eval_dir::String,
        T_mid::Float64,
        base_rho_grid::Vector{Float64},
        xi::Float64,
        reverse_rho::Bool,
        p_num::Int,
        t_num::Int,
        thermo_quadrature_kwargs::NamedTuple,
        cfg::ProductionPipelineConfig,
        session;
        prior::Union{Nothing, RhoSupportRefinement.RhoSupportPrior}=nothing)
    token = replace(@sprintf("%.6f", T_mid), "." => "p", "-" => "m")
    out_csv = joinpath(eval_dir, "prod_eval_T$(token)_cascade.csv")
    io_mode = isfile(out_csv) ? "a" : "w"
    written = Set{Tuple{Float64, Float64, Float64}}()
    if io_mode == "a"
        union!(written, _existing_session_keys(out_csv))
        for row in _production_session_rows(session, T_mid, xi)
            push!(written, (T_mid, xi, row.rho))
        end
    end

    open(out_csv, io_mode) do io
        io_mode == "w" && println(io, TrhoScan.HEADER)
        level_results = NamedTuple[]
        targeted_total = 0
        targeted_seen = Set{Float64}()
        rho_levels = (copy(base_rho_grid), _rho_support_fine_grid(base_rho_grid, cfg.rho_support_fine_step))
        before = TrhoScan.rho_session_snapshot(session)

        for (level, rho_grid) in zip(0:1, rho_levels)
            targeted_values = Float64[]
            remaining = max(cfg.rho_support_targeted_cap - targeted_total, 0)
            target_count = min(cfg.rho_support_config.target_point_count, remaining)
            iseven(target_count) && (target_count -= 1)
            support_cfg = if target_count >= 5
                RhoSupportRefinement.RhoSupportConfig(
                    support_slope_tol=cfg.rho_support_config.support_slope_tol,
                    positive_slope_margin=cfg.rho_support_config.positive_slope_margin,
                    negative_slope_margin=cfg.rho_support_config.negative_slope_margin,
                    minimum_negative_secant_run=cfg.rho_support_config.minimum_negative_secant_run,
                    target_point_count=target_count,
                    max_extra_points=max(target_count, min(cfg.rho_support_config.max_extra_points, remaining)),
                    support_expansion_gaps=cfg.rho_support_config.support_expansion_gaps,
                    local_fit_rmse_tol=cfg.rho_support_config.local_fit_rmse_tol,
                    near_critical_slope_tol=cfg.rho_support_config.near_critical_slope_tol,
                )
            else
                nothing
            end

            curve, _ = _production_session_scan!(session, io, written, T_mid, xi,
                rho_grid, reverse_rho, p_num, t_num, thermo_quadrature_kwargs)
            sample_mu = if support_cfg === nothing || curve === nothing
                nothing
            else
                value -> begin
                    value = Float64(value)
                    push!(targeted_values, value)
                    row = TrhoScan.rho_point!(session, T_mid, xi, value; targeted=true)
                    _write_session_row!(io, written, row, T_mid, xi, p_num, t_num, thermo_quadrature_kwargs)
                    row.muq_MeV
                end
            end
            assessment = if sample_mu === nothing
                nothing
            else
                RhoSupportRefinement.analyze_rho_support_cascade(
                    curve[2], curve[1]; sample_mu=sample_mu,
                    prior=prior,
                    config=support_cfg,
                )
            end
            for value in unique(targeted_values)
                value in targeted_seen && continue
                push!(targeted_seen, value)
                targeted_total += 1
            end
            final_curve, _ = _production_session_curve(session, T_mid, xi)

            current = if final_curve === nothing
                (
                    status=:invalid, mu_transition=nothing, area_residual=nothing,
                    sres=SShapeResult(), rho_hadron=nothing, rho_quark=nothing,
                    reason="solver_or_curve_failure", level=level, curve=nothing,
                    out_csv=out_csv,
                )
            else
                mu_vals, rho_vals = final_curve
                cres = _classify_s_curve(mu_vals, rho_vals;
                    area_tol_good=cfg.area_tol_good, area_tol_bad=cfg.area_tol_bad)
                maxwell = maxwell_construction(mu_vals, rho_vals; spinodal_hint=cres.sres)
                (
                    status=cres.status,
                    mu_transition=cres.mu_transition,
                    area_residual=cres.area_residual,
                    sres=cres.sres,
                    rho_hadron=maxwell.converged ? maxwell.rho_hadron : nothing,
                    rho_quark=maxwell.converged ? maxwell.rho_quark : nothing,
                    reason=String(cres.reason), level=level, curve=final_curve,
                    out_csv=out_csv,
                )
            end
            push!(level_results, merge(current, (
                cascade_status=assessment === nothing ? :not_run : assessment.status,
                cascade_support_center=assessment === nothing ? nothing : assessment.spinodal_rho_center,
                cascade_support_gap=assessment === nothing ? nothing : assessment.spinodal_rho_gap,
                cascade_support_low=assessment === nothing ? nothing : assessment.support_low,
                cascade_support_high=assessment === nothing ? nothing : assessment.support_high,
                cascade_support_origin=assessment === nothing ? :none : assessment.support_origin,
                targeted_count=targeted_total,
            )))
        end

        coarse, fine = level_results
        geometry = if coarse.status == :invalid && fine.status == :invalid &&
                         coarse.reason == "no_s_shape" && fine.reason == "no_s_shape"
            PhaseGeometryError(comparable=true, converged=true, position_MeV=0.0,
                density=0.0, maxwell_area=0.0, response_rtol=0.0, reason="stable_no_s_shape")
        else
            _compare_phase_geometry(coarse, fine, PhaseGeometryTolerances(
                position_MeV=cfg.rho_position_tol_MeV,
                density=cfg.rho_density_tol,
                maxwell_area=cfg.rho_maxwell_area_tol))
        end
        final = merge(fine, (
            geometry_converged=geometry.converged,
            position_error_MeV=geometry.position_MeV,
            density_error=geometry.density,
            maxwell_area_gate=geometry.maxwell_area,
        ))
        semantic_status = if coarse.status == :valid && fine.status == :valid && geometry.converged
            :confirmed_first_order
        elseif coarse.status == :invalid && fine.status == :invalid &&
               coarse.reason == "no_s_shape" && fine.reason == "no_s_shape" && geometry.converged
            :confirmed_monotone
        else
            :ambiguous_near_critical
        end
        after = TrhoScan.rho_session_snapshot(session)
        unique_delta = after.unique_solves - before.unique_solves
        failure_delta = after.failed_points - before.failed_points
        _append_scan_csv!(aggregate_csv, out_csv)
        return merge(final, (
            raw_status=final.status,
            slice_status=semantic_status,
            coarse_status=coarse.status,
            fine_status=fine.status,
            coarse_reason=coarse.reason,
            fine_reason=fine.reason,
            stats=(total=unique_delta, success=unique_delta - failure_delta, failure=failure_delta),
            cache_stats=(
                point_requests=after.point_requests - before.point_requests,
                cache_hits=after.cache_hits - before.cache_hits,
                unique_solves=unique_delta,
                targeted_additions=after.targeted_additions - before.targeted_additions,
            ),
            rho_convergence_records=[(
                axis="rho", level=1, left=0.0, right=1.0, midpoint=1.0,
                position_error_MeV=isfinite(geometry.position_MeV) ? geometry.position_MeV : nothing,
                density_error=isfinite(geometry.density) ? geometry.density : nothing,
                maxwell_area=isfinite(geometry.maxwell_area) ? geometry.maxwell_area : nothing,
                response_rtol=isfinite(geometry.response_rtol) ? geometry.response_rtol : nothing,
                converged=geometry.converged, reason=geometry.reason,
            )],
            cascade_targeted_count=targeted_total,
            cascade_coarse_point_count=length(base_rho_grid),
            cascade_fine_point_count=length(rho_levels[2]),
            cascade_support_low=get(fine, :cascade_support_low, nothing),
            cascade_support_high=get(fine, :cascade_support_high, nothing),
            cascade_support_origin=get(fine, :cascade_support_origin, :none),
        ))
    end
end

function _production_classify_temperature_memoized_uniform(
        aggregate_csv::String,
        eval_dir::String,
        T_mid::Float64,
        base_rho_grid::Vector{Float64},
        xi::Float64,
        reverse_rho::Bool,
        p_num::Int,
        t_num::Int,
        thermo_quadrature_kwargs::NamedTuple,
        cfg::ProductionPipelineConfig,
        session;
        level_start::Int=0,
        level_stop::Int=cfg.max_refine_level_rho,
        output_suffix::AbstractString="memoized")
    0 <= level_start <= level_stop || throw(ArgumentError(
        "memoized rho refinement level range is invalid: $(level_start):$(level_stop)",
    ))
    token = replace(@sprintf("%.6f", T_mid), "." => "p", "-" => "m")
    out_csv = joinpath(eval_dir, "prod_eval_T$(token)_$(output_suffix).csv")
    io_mode = isfile(out_csv) ? "a" : "w"
    written = Set{Tuple{Float64, Float64, Float64}}()
    if io_mode == "a"
        union!(written, _existing_session_keys(out_csv))
        foreach(row -> push!(written, (T_mid, xi, row.rho)), _production_session_rows(session, T_mid, xi))
    end
    open(out_csv, io_mode) do io
        io_mode == "w" && println(io, TrhoScan.HEADER)
        before = TrhoScan.rho_session_snapshot(session)
        level_results = NamedTuple[]
        previous = nothing
        for level in level_start:level_stop
            rho_grid = level == 0 ? copy(base_rho_grid) : _refine_rho_grid(base_rho_grid, level)
            _production_session_scan!(session, io, written, T_mid, xi, rho_grid, reverse_rho,
                p_num, t_num, thermo_quadrature_kwargs)
            curve, _ = _production_session_curve(session, T_mid, xi)
            current = if curve === nothing
                (status=:invalid, mu_transition=nothing, area_residual=nothing,
                 sres=SShapeResult(), rho_hadron=nothing, rho_quark=nothing,
                 reason="solver_or_curve_failure", level=level, curve=nothing, out_csv=out_csv)
            else
                mu_vals, rho_vals = curve
                cres = _classify_s_curve(mu_vals, rho_vals;
                    area_tol_good=cfg.area_tol_good, area_tol_bad=cfg.area_tol_bad)
                maxwell = maxwell_construction(mu_vals, rho_vals; spinodal_hint=cres.sres)
                (status=cres.status, mu_transition=cres.mu_transition,
                 area_residual=cres.area_residual, sres=cres.sres,
                 rho_hadron=maxwell.converged ? maxwell.rho_hadron : nothing,
                 rho_quark=maxwell.converged ? maxwell.rho_quark : nothing,
                 reason=String(cres.reason), level=level, curve=curve, out_csv=out_csv)
            end
            geometry = previous === nothing ? PhaseGeometryError(reason="coarse_level_only") :
                _compare_phase_geometry(previous, current, PhaseGeometryTolerances(
                    position_MeV=cfg.rho_position_tol_MeV,
                    density=cfg.rho_density_tol,
                    maxwell_area=cfg.rho_maxwell_area_tol))
            result = merge(current, (
                geometry_converged=geometry.converged,
                position_error_MeV=geometry.position_MeV,
                density_error=geometry.density,
                maxwell_area_gate=geometry.maxwell_area,
            ))
            push!(level_results, result)
            previous = current
            geometry.converged && level > level_start && break
        end
        coarse = first(level_results)
        fine = last(level_results)
        geometry = get(fine, :geometry_converged, false)
        semantic_status = if coarse.status == :valid && fine.status == :valid && geometry
            :confirmed_first_order
        elseif coarse.status == :invalid && fine.status == :invalid &&
               coarse.reason == "no_s_shape" && fine.reason == "no_s_shape" && geometry
            :confirmed_monotone
        else
            :ambiguous_near_critical
        end
        after = TrhoScan.rho_session_snapshot(session)
        _append_scan_csv!(aggregate_csv, out_csv)
        return merge(fine, (
            raw_status=fine.status,
            slice_status=semantic_status,
            coarse_status=coarse.status,
            fine_status=fine.status,
            coarse_reason=coarse.reason,
            fine_reason=fine.reason,
            stats=(total=after.unique_solves - before.unique_solves,
                success=(after.unique_solves - before.unique_solves) - (after.failed_points - before.failed_points),
                failure=after.failed_points - before.failed_points),
            cache_stats=(
                point_requests=after.point_requests - before.point_requests,
                cache_hits=after.cache_hits - before.cache_hits,
                unique_solves=after.unique_solves - before.unique_solves,
                targeted_additions=after.targeted_additions - before.targeted_additions,
            ),
            rho_convergence_records=[(
                axis="rho", level=level_start == 0 ? 1 : Int(fine.level),
                left=level_start == 0 ? 0.0 : Float64(max(level_start, fine.level - 1)),
                right=level_start == 0 ? 1.0 : Float64(fine.level),
                midpoint=level_start == 0 ? 1.0 : Float64(fine.level),
                position_error_MeV=isfinite(fine.position_error_MeV) ? fine.position_error_MeV : nothing,
                density_error=isfinite(fine.density_error) ? fine.density_error : nothing,
                maxwell_area=isfinite(fine.maxwell_area_gate) ? fine.maxwell_area_gate : nothing,
                response_rtol=nothing, converged=geometry,
                reason=geometry ? "ok" : "rho_geometry_not_converged",
            )],
        ))
    end
end

function _production_session_curve_for_grid(
        session,
        T::Float64,
        xi::Float64,
        rho_grid::AbstractVector)
    rows = NamedTuple[]
    for rho in sort(unique(Float64.(rho_grid)))
        key = (T, xi, rho)
        haskey(session.cache, key) || return nothing, rows
        push!(rows, session.cache[key])
    end
    isempty(rows) && return nothing, rows
    all(row -> row.converged && row.finite, rows) || return nothing, rows
    return ([row.muq_MeV for row in rows], [row.rho for row in rows]), rows
end

function _classify_production_curve(curve, cfg::ProductionPipelineConfig,
        level::Int, out_csv::String)
    curve === nothing && return (
        status=:invalid,
        mu_transition=nothing,
        area_residual=nothing,
        sres=SShapeResult(),
        rho_hadron=nothing,
        rho_quark=nothing,
        reason="solver_or_curve_failure",
        level=level,
        curve=nothing,
        out_csv=out_csv,
    )
    mu_vals, rho_vals = curve
    cres = _classify_s_curve(mu_vals, rho_vals;
        area_tol_good=cfg.area_tol_good, area_tol_bad=cfg.area_tol_bad)
    maxwell = maxwell_construction(mu_vals, rho_vals; spinodal_hint=cres.sres)
    return (
        status=cres.status,
        mu_transition=cres.mu_transition,
        area_residual=cres.area_residual,
        sres=cres.sres,
        rho_hadron=maxwell.converged ? maxwell.rho_hadron : nothing,
        rho_quark=maxwell.converged ? maxwell.rho_quark : nothing,
        reason=String(cres.reason),
        level=level,
        curve=curve,
        out_csv=out_csv,
    )
end

function _hybrid_support_grid(
        base_rho_grid::Vector{Float64},
        stage_a,
        stage_b;
        local_step::Float64=0.003125,
        padding::Float64=0.025)
    local_step > 0 || throw(ArgumentError("hybrid local rho step must be positive"))
    padding >= 0 || throw(ArgumentError("hybrid rho support padding must be nonnegative"))
    first_rho, last_rho = first(base_rho_grid), last(base_rho_grid)
    support_values = Float64[]
    support_sources = Symbol[]
    for (label, result) in ((:cascade, stage_a), (:dense, stage_b))
        for (field, source) in ((:cascade_support_low, :cascade_support),
                (:cascade_support_high, :cascade_support),
                (:rho_hadron, :coexistence_density),
                (:rho_quark, :coexistence_density))
            value = get(result, field, nothing)
            value === nothing && continue
            value = Float64(value)
            isfinite(value) && first_rho <= value <= last_rho || continue
            push!(support_values, value)
            push!(support_sources, source)
        end
        sres = get(result, :sres, SShapeResult())
        for (field, source) in ((:rho_spinodal_hadron, :spinodal_density),
                (:rho_spinodal_quark, :spinodal_density))
            value = getproperty(sres, field)
            value === nothing && continue
            value = Float64(value)
            isfinite(value) && first_rho <= value <= last_rho || continue
            push!(support_values, value)
            push!(support_sources, source)
        end
    end
    isempty(support_values) && return nothing
    low_raw = max(first_rho, minimum(support_values) - padding)
    high_raw = min(last_rho, maximum(support_values) + padding)
    high_raw > low_raw || return nothing

    # Align both edges to the local oracle grid.  Rounding outwards is
    # deliberate: no Stage-C point may lie outside the declared support.
    origin = first_rho
    low_index = floor(Int, (low_raw - origin) / local_step + 1e-10)
    high_index = ceil(Int, (high_raw - origin) / local_step - 1e-10)
    low = max(first_rho, origin + low_index * local_step)
    high = min(last_rho, origin + high_index * local_step)
    high > low || return nothing
    grid = Float64.(collect(low:local_step:high))
    isempty(grid) && return nothing
    return (
        low=low,
        high=high,
        grid=grid,
        local_step=local_step,
        padding=padding,
        source=sort!(unique(support_sources)),
    )
end

function _hybrid_stats(before, after)
    unique_delta = after.unique_solves - before.unique_solves
    failure_delta = after.failed_points - before.failed_points
    return (
        stats=(total=unique_delta, success=unique_delta - failure_delta, failure=failure_delta),
        cache_stats=(
            point_requests=after.point_requests - before.point_requests,
            cache_hits=after.cache_hits - before.cache_hits,
            unique_solves=unique_delta,
            targeted_additions=after.targeted_additions - before.targeted_additions,
        ),
    )
end

function _production_classify_temperature_hybrid(
        aggregate_csv::String,
        eval_dir::String,
        T_mid::Float64,
        base_rho_grid::Vector{Float64},
        xi::Float64,
        reverse_rho::Bool,
        p_num::Int,
        t_num::Int,
        thermo_quadrature_kwargs::NamedTuple,
        cfg::ProductionPipelineConfig,
        session;
        prior::Union{Nothing, RhoSupportRefinement.RhoSupportPrior}=nothing)
    before = TrhoScan.rho_session_snapshot(session)
    stage_a = _production_classify_temperature_cascade(
        aggregate_csv, eval_dir, T_mid, base_rho_grid, xi, reverse_rho,
        p_num, t_num, thermo_quadrature_kwargs, cfg, session; prior=prior)

    common = (
        stage_a_status=stage_a.slice_status,
        stage_b_status=:not_run,
        stage_c_status=:not_run,
        stage_used=:stage_a,
        hybrid_upgrade_reason="stage_a_certificate",
        hybrid_support_low=nothing,
        hybrid_support_high=nothing,
        hybrid_support_source=Symbol[],
        hybrid_verification_point_count=0,
    )

    # Stage A is allowed to certify monotonicity because it already contains
    # both cascade rho layers.  A first-order result always continues to
    # Stage B; this prevents the cheaper support heuristic from publishing a
    # first-order certificate on its own.
    if stage_a.slice_status == :confirmed_monotone
        after = TrhoScan.rho_session_snapshot(session)
        return merge(stage_a, common, _hybrid_stats(before, after))
    end

    stage_b = _production_classify_temperature_memoized_uniform(
        aggregate_csv, eval_dir, T_mid, base_rho_grid, xi, reverse_rho,
        p_num, t_num, thermo_quadrature_kwargs, cfg, session;
        level_start=2, level_stop=3, output_suffix="hybrid_stage_b")
    stage_a_first = stage_a.slice_status == :confirmed_first_order && stage_a.raw_status == :valid
    # Stage C is specifically allowed to repair a Stage-B geometry mismatch;
    # it therefore needs a valid Maxwell curve, not a pre-existing Stage-B
    # semantic certificate.
    stage_b_first = stage_b.raw_status == :valid
    cross_ab = _compare_phase_geometry(stage_a, stage_b, PhaseGeometryTolerances(
        position_MeV=cfg.rho_position_tol_MeV,
        density=cfg.rho_density_tol,
        maxwell_area=cfg.rho_maxwell_area_tol,
    ))

    if stage_a_first && stage_b_first && cross_ab.converged
        after = TrhoScan.rho_session_snapshot(session)
        return merge(stage_b, (
            slice_status=:confirmed_first_order,
            geometry_converged=true,
            position_error_MeV=cross_ab.position_MeV,
            density_error=cross_ab.density,
            maxwell_area_gate=cross_ab.maxwell_area,
            stage_b_status=stage_b.slice_status,
            stage_c_status=:not_run,
            stage_used=:stage_b,
            hybrid_upgrade_reason="stage_a_stage_b_first_order_geometry_pass",
            hybrid_support_low=nothing,
            hybrid_support_high=nothing,
            hybrid_support_source=Symbol[],
            hybrid_verification_point_count=0,
        ), _hybrid_stats(before, after))
    end

    support = _hybrid_support_grid(base_rho_grid, stage_a, stage_b)
    if support === nothing
        after = TrhoScan.rho_session_snapshot(session)
        return merge(stage_b, (
            slice_status=:ambiguous_near_critical,
            geometry_converged=false,
            stage_a_status=stage_a.slice_status,
            stage_b_status=stage_b.slice_status,
            stage_c_status=:not_run,
            stage_used=:stage_b,
            hybrid_upgrade_reason="no_reliable_support_for_local_oracle",
            hybrid_support_low=nothing,
            hybrid_support_high=nothing,
            hybrid_support_source=Symbol[],
            hybrid_verification_point_count=0,
        ), _hybrid_stats(before, after))
    end

    token = replace(@sprintf("%.6f", T_mid), "." => "p", "-" => "m")
    out_csv = joinpath(eval_dir, "prod_eval_T$(token)_hybrid_stage_c.csv")
    io_mode = isfile(out_csv) ? "a" : "w"
    written = Set{Tuple{Float64, Float64, Float64}}()
    io_mode == "a" && union!(written, _existing_session_keys(out_csv))
    # Stage C reuses all Stage-B points inside the declared support and adds
    # only the local 0.003125 grid.  Build the union before scanning so
    # floating-point representations of shared nested points are materialized
    # under the exact keys used by the verification curve.
    stage_b_grid = _refine_rho_grid(base_rho_grid, 3)
    verification_grid = sort!(unique(vcat(
        [rho for rho in stage_b_grid if support.low <= rho <= support.high],
        support.grid,
    )))
    open(out_csv, io_mode) do io
        io_mode == "w" && println(io, TrhoScan.HEADER)
        _production_session_scan!(session, io, written, T_mid, xi, verification_grid,
            reverse_rho, p_num, t_num, thermo_quadrature_kwargs)
    end
    # It never falls back to a global oracle: every verification point is
    # either a Stage-B point within support or part of support.grid.
    curve_c, _ = _production_session_curve_for_grid(session, T_mid, xi, verification_grid)
    stage_c = _classify_production_curve(curve_c, cfg, 4, out_csv)
    geometry_bc = _compare_phase_geometry(stage_b, stage_c, PhaseGeometryTolerances(
        position_MeV=cfg.rho_position_tol_MeV,
        density=cfg.rho_density_tol,
        maxwell_area=cfg.rho_maxwell_area_tol,
    ))
    stage_c_valid = stage_c.status == :valid && stage_b_first && geometry_bc.converged
    _append_scan_csv!(aggregate_csv, out_csv)
    after = TrhoScan.rho_session_snapshot(session)
    stats = _hybrid_stats(before, after)
    return merge(stage_c, stats, (
        raw_status=stage_c.status,
        slice_status=stage_c_valid ? :confirmed_first_order : :ambiguous_near_critical,
        coarse_status=stage_a.raw_status,
        fine_status=stage_c.status,
        coarse_reason=stage_a.reason,
        fine_reason=stage_c.reason,
        geometry_converged=geometry_bc.converged && stage_c_valid,
        position_error_MeV=geometry_bc.position_MeV,
        density_error=geometry_bc.density,
        maxwell_area_gate=geometry_bc.maxwell_area,
        stage_a_status=stage_a.slice_status,
        stage_b_status=stage_b.slice_status,
        stage_c_status=stage_c_valid ? :confirmed_first_order : stage_c.status,
        stage_used=:stage_c,
        hybrid_upgrade_reason=stage_c_valid ?
            "local_stage_c_first_order_geometry_pass" :
            "stage_c_local_oracle_did_not_close_geometry",
        hybrid_support_low=support.low,
        hybrid_support_high=support.high,
        hybrid_support_source=support.source,
        hybrid_verification_point_count=length(verification_grid),
        rho_convergence_records=vcat(
            collect(get(stage_a, :rho_convergence_records, NamedTuple[])),
            collect(get(stage_b, :rho_convergence_records, NamedTuple[])),
            [(
                axis="rho", level=4, left=3.0, right=4.0, midpoint=4.0,
                position_error_MeV=isfinite(geometry_bc.position_MeV) ? geometry_bc.position_MeV : nothing,
                density_error=isfinite(geometry_bc.density) ? geometry_bc.density : nothing,
                maxwell_area=isfinite(geometry_bc.maxwell_area) ? geometry_bc.maxwell_area : nothing,
                response_rtol=nothing, converged=stage_c_valid,
                reason=stage_c_valid ? "ok" : "hybrid_stage_c_not_converged",
            )],
        ),
    ))
end

function _production_classify_temperature(
        aggregate_csv::String,
        eval_dir::String,
        T_mid::Float64,
        base_rho_grid::Vector{Float64},
        xi::Float64,
        reverse_rho::Bool,
        seed_policy::Symbol,
        solver_backend::Symbol,
        model_kind::Symbol,
        p_num::Int,
        t_num::Int,
        thermo_quadrature_kwargs::NamedTuple,
        iterations::Int,
        cfg::ProductionPipelineConfig;
        rho_session=nothing,
        rho_prior::Union{Nothing, RhoSupportRefinement.RhoSupportPrior}=nothing,
        work_telemetry::Union{Nothing, SolverWorkTelemetry}=nothing)
    if cfg.rho_refinement_policy === :rho_support_hybrid
        return _production_classify_temperature_hybrid(
            aggregate_csv, eval_dir, T_mid, base_rho_grid, xi, reverse_rho,
            p_num, t_num, thermo_quadrature_kwargs, cfg, rho_session;
            prior=rho_prior,
        )
    elseif cfg.rho_refinement_policy === :rho_support_cascade
        return _production_classify_temperature_cascade(
            aggregate_csv, eval_dir, T_mid, base_rho_grid, xi, reverse_rho,
            p_num, t_num, thermo_quadrature_kwargs, cfg, rho_session;
            prior=rho_prior,
        )
    elseif rho_session !== nothing
        return _production_classify_temperature_memoized_uniform(
            aggregate_csv, eval_dir, T_mid, base_rho_grid, xi, reverse_rho,
            p_num, t_num, thermo_quadrature_kwargs, cfg, rho_session,
        )
    end
    adaptive_cfg = AdaptiveRhoRefinement.AdaptiveRhoConfig(
        slope_tol=cfg.adaptive_slope_tol,
        min_gap=cfg.adaptive_min_gap,
        max_points=cfg.adaptive_max_points,
        digits=cfg.adaptive_digits,
    )

    rho_eval = copy(base_rho_grid)
    last_result = nothing
    previous_result = nothing
    convergence_records = NamedTuple[]
    pass = 1
    scan_total = 0
    scan_success = 0
    scan_failure = 0
    geometry_tol = PhaseGeometryTolerances(
        position_MeV=cfg.rho_position_tol_MeV,
        density=cfg.rho_density_tol,
        maxwell_area=cfg.rho_maxwell_area_tol,
    )

    for level in 0:cfg.max_refine_level_rho
        scan = _run_single_temperature_production_scan(
            aggregate_csv,
            eval_dir,
            T_mid,
            rho_eval,
            xi,
            reverse_rho,
            seed_policy,
            solver_backend,
            model_kind,
            p_num,
            t_num,
            thermo_quadrature_kwargs,
            iterations,
            level,
            pass,
            work_telemetry,
        )
        pass += 1
        scan_total += scan.stats.total
        scan_success += scan.stats.success
        scan_failure += scan.stats.failure
        if scan.curve === nothing
            last_result = (
                status=:invalid,
                mu_transition=nothing,
                area_residual=nothing,
                sres=SShapeResult(),
                rho_hadron=nothing,
                rho_quark=nothing,
                reason="no_curve",
                level=level,
                curve=nothing,
                geometry_converged=false,
                position_error_MeV=Inf,
                density_error=Inf,
                maxwell_area_gate=Inf,
                out_csv=scan.out_csv,
            )
            break
        end

        mu_vals, rho_vals = scan.curve
        cres = _classify_s_curve(
            mu_vals,
            rho_vals;
            area_tol_good=cfg.area_tol_good,
            area_tol_bad=cfg.area_tol_bad,
        )
        maxwell = maxwell_construction(mu_vals, rho_vals; spinodal_hint=cres.sres)
        current_result = (
            status=cres.status,
            mu_transition=cres.mu_transition,
            area_residual=cres.area_residual,
            sres=cres.sres,
            rho_hadron=maxwell.converged ? maxwell.rho_hadron : nothing,
            rho_quark=maxwell.converged ? maxwell.rho_quark : nothing,
            reason=String(cres.reason),
            level=level,
            curve=scan.curve,
            out_csv=scan.out_csv,
        )

        geometry_error = previous_result === nothing ? PhaseGeometryError(reason="coarse_level_only") :
            _compare_phase_geometry(previous_result, current_result, geometry_tol)
        if previous_result !== nothing
            push!(convergence_records, (
                axis="rho",
                level=level,
                left=Float64(level - 1),
                right=Float64(level),
                midpoint=Float64(level),
                position_error_MeV=isfinite(geometry_error.position_MeV) ? geometry_error.position_MeV : nothing,
                density_error=isfinite(geometry_error.density) ? geometry_error.density : nothing,
                maxwell_area=isfinite(geometry_error.maxwell_area) ? geometry_error.maxwell_area : nothing,
                response_rtol=isfinite(geometry_error.response_rtol) ? geometry_error.response_rtol : nothing,
                converged=geometry_error.converged,
                reason=geometry_error.reason,
            ))
        end

        last_result = merge(current_result, (
            geometry_converged=geometry_error.converged,
            position_error_MeV=geometry_error.position_MeV,
            density_error=geometry_error.density,
            maxwell_area_gate=geometry_error.maxwell_area,
        ))

        if cfg.rho_geometry_convergence
            geometry_error.converged && break
            if level == cfg.max_refine_level_rho
                if last_result.status == :valid
                    last_result = merge(last_result, (
                        status=:unknown,
                        reason="rho_geometry_not_converged:$(geometry_error.reason)",
                    ))
                end
                break
            end
            previous_result = current_result
            rho_eval = _refine_rho_grid(base_rho_grid, level + 1)
            continue
        end

        if cres.status != :unknown || !cfg.adaptive_rho || level == cfg.max_refine_level_rho
            break
        end

        additions = AdaptiveRhoRefinement.suggest_refinement_points(rho_vals, mu_vals; config=adaptive_cfg)
        isempty(additions) && break
        previous_result = current_result
        rho_eval = AdaptiveRhoRefinement.merge_rho_values(rho_eval, additions; digits=adaptive_cfg.digits)
    end

    last_result === nothing && error("production rho refinement produced no classification at T=$(T_mid)")
    _append_scan_csv!(aggregate_csv, last_result.out_csv)
    # `rho_geometry_convergence=false` is an explicit legacy/diagnostic
    # opt-out.  In that mode a valid Maxwell slice may remain a confirmed
    # first-order candidate for backwards-compatible boundary consumers, but
    # a monotone certificate is never manufactured from the single rho layer.
    # The default production path keeps the strict two-layer geometry gate.
    geometry_gate_satisfied = !cfg.rho_geometry_convergence || last_result.geometry_converged
    semantic_status = if last_result.status == :valid && geometry_gate_satisfied
        :confirmed_first_order
    elseif last_result.status == :invalid &&
           last_result.reason == "no_s_shape" &&
           last_result.geometry_converged &&
           any(record -> record.reason == "stable_no_s_shape", convergence_records)
        :confirmed_monotone
    else
        :ambiguous_near_critical
    end
    return merge(last_result, (
        raw_status=last_result.status,
        slice_status=semantic_status,
        stats=(total=scan_total, success=scan_success, failure=scan_failure),
        rho_convergence_records=convergence_records,
    ))
end

function _push_production_boundary!(rows::Vector{NamedTuple}, T::Float64, res)
    push!(rows, (
        T_MeV=Float64(T),
        mu_transition_MeV=Float64(something(res.mu_transition, NaN)),
        rho_hadron=Float64(something(res.rho_hadron, NaN)),
        rho_quark=Float64(something(res.rho_quark, NaN)),
        area_residual=Float64(something(res.area_residual, NaN)),
        converged=true,
    ))
end

function _push_production_spinodal!(rows::Vector{NamedTuple}, T::Float64, res)
    push!(rows, (
        T_MeV=Float64(T),
        mu_spinodal_hadron_MeV=Float64(something(res.sres.mu_spinodal_hadron, NaN)),
        mu_spinodal_quark_MeV=Float64(something(res.sres.mu_spinodal_quark, NaN)),
        rho_spinodal_hadron=Float64(something(res.sres.rho_spinodal_hadron, NaN)),
        rho_spinodal_quark=Float64(something(res.sres.rho_spinodal_quark, NaN)),
    ))
end

function _materialize_sweep_result(records::AbstractVector{<:NamedTuple}, first_point_fallback::Bool, fallback_start_T_MeV::Float64,
        unknown_count::Int, forced_invalid_count::Int)
    return FirstOrderSweepResult(
        temperatures_MeV=[Float64(r.T_MeV) for r in records],
        statuses=[Symbol(r.status) for r in records],
        mu_transitions_MeV=[Float64(r.mu_transition_MeV) for r in records],
        area_residuals=[Float64(r.area_residual) for r in records],
        reasons=[String(r.reason) for r in records],
        first_point_fallback=first_point_fallback,
        fallback_start_T_MeV=fallback_start_T_MeV,
        unknown_count=unknown_count,
        forced_invalid_count=forced_invalid_count,
    )
end

function _find_production_bracket(records::AbstractVector{<:NamedTuple})
    prev_first_order = nothing
    for record in records
        status = if record.status == :valid
            :confirmed_first_order
        elseif record.status == :invalid
            :confirmed_monotone
        else
            record.status
        end
        if status == :confirmed_first_order
            prev_first_order = record
        elseif status == :confirmed_monotone && prev_first_order !== nothing
            return (T_low=prev_first_order.T_MeV, mu_low=prev_first_order.mu_transition_MeV, T_high=record.T_MeV)
        end
    end
    return nothing
end

function _refine_production_cep_frontiers(
        bracket,
        evaluate_temperature::Function,
        cfg::ProductionPipelineConfig;
        unknown_count_fn::Function=()->0)
    bracket === nothing && return nothing
    low = (T=Float64(bracket.T_low), mu=Float64(bracket.mu_low))
    high = (T=Float64(bracket.T_high),)
    low_search_hi = high
    high_search_lo = low
    budget_exhausted = false

    for _ in 1:cfg.cep_max_bisect_iter
        if unknown_count_fn() > cfg.unknown_budget
            budget_exhausted = true
            break
        end
        changed = false
        if low_search_hi.T - low.T > cfg.temperature_resolution_target_MeV
            T_mid = 0.5 * (low.T + low_search_hi.T)
            res = evaluate_temperature(T_mid)
            if res.slice_status == :confirmed_first_order
                low = (T=T_mid, mu=Float64(something(res.mu_transition, low.mu)))
            else
                low_search_hi = (T=T_mid,)
                res.slice_status == :confirmed_monotone && (high = (T=T_mid,))
            end
            changed = true
        end
        if unknown_count_fn() > cfg.unknown_budget
            budget_exhausted = true
            break
        end
        if high.T - high_search_lo.T > cfg.temperature_resolution_target_MeV
            T_mid = 0.5 * (high_search_lo.T + high.T)
            res = evaluate_temperature(T_mid)
            if res.slice_status == :confirmed_monotone
                high = (T=T_mid,)
            else
                high_search_lo = (T=T_mid,)
                res.slice_status == :confirmed_first_order &&
                    (low = (T=T_mid, mu=Float64(something(res.mu_transition, low.mu))))
            end
            changed = true
        end
        changed || break
        low_search_hi.T - low.T <= cfg.temperature_resolution_target_MeV && high.T - high_search_lo.T <= cfg.temperature_resolution_target_MeV && break
    end
    return (low=low, high=high, budget_exhausted=budget_exhausted, unknown_count=unknown_count_fn())
end

function _adaptive_production_temperature_refinement!(
        temperatures::Vector{Float64},
        evaluate_temperature::Function,
        cfg::ProductionPipelineConfig;
        stop_refinement::Function=()->false)
    cfg.adaptive_temperature || return sort(unique(copy(temperatures))), NamedTuple[]
    cfg.temperature_max_refine_level > 0 || return sort(unique(copy(temperatures))), NamedTuple[]

    tol = PhaseGeometryTolerances(
        position_MeV=cfg.temperature_position_tol_MeV,
        density=cfg.temperature_density_tol,
        maxwell_area=cfg.temperature_maxwell_area_tol,
    )
    resolved = sort(unique(copy(temperatures)))
    intervals = Tuple{Float64, Float64}[]
    for i in 1:(length(resolved) - 1)
        left = evaluate_temperature(resolved[i])
        right = evaluate_temperature(resolved[i + 1])
        if left.status == :valid && right.status == :valid
            push!(intervals, (resolved[i], resolved[i + 1]))
        end
    end

    records = NamedTuple[]
    for level in 1:cfg.temperature_max_refine_level
        stop_refinement() && break
        isempty(intervals) && break
        midpoints = sort(unique(Float64[0.5 * (left + right) for (left, right) in intervals]))
        for midpoint in midpoints
            evaluate_temperature(Float64(midpoint))
        end
        stop_refinement() && break

        next_intervals = Tuple{Float64, Float64}[]
        for (left_T, right_T) in intervals
            midpoint_T = 0.5 * (left_T + right_T)
            left = evaluate_temperature(left_T)
            midpoint = evaluate_temperature(midpoint_T)
            right = evaluate_temperature(right_T)
            error = _phase_geometry_midpoint_error(left, midpoint, right, tol)
            push!(records, (
                axis="temperature",
                level=level,
                left=left_T,
                right=right_T,
                midpoint=midpoint_T,
                position_error_MeV=isfinite(error.position_MeV) ? error.position_MeV : nothing,
                density_error=isfinite(error.density) ? error.density : nothing,
                maxwell_area=isfinite(error.maxwell_area) ? error.maxwell_area : nothing,
                response_rtol=isfinite(error.response_rtol) ? error.response_rtol : nothing,
                converged=error.converged,
                reason=error.reason,
            ))
            push!(resolved, midpoint_T)
            if !error.converged && midpoint.status == :valid
                push!(next_intervals, (left_T, midpoint_T))
                push!(next_intervals, (midpoint_T, right_T))
            end
        end
        sort!(resolved)
        unique!(resolved)
        intervals = next_intervals
    end
    return resolved, records
end

function _production_crossover_temperature_bounds(
        cfg::ProductionPipelineConfig,
        cep::CEPResult,
        boundary::AbstractVector{<:NamedTuple})
    T_max_mev = isfinite(cfg.crossover_T_max_MeV) ? cfg.crossover_T_max_MeV : cfg.T_end
    T_min_mev = cfg.T_start

    if !isempty(boundary)
        boundary_top = maximum(row.T_MeV for row in boundary if isfinite(row.T_MeV))
        T_min_mev = max(cfg.T_start, boundary_top)
    elseif cep.found && isfinite(cep.T_cep_MeV)
        T_min_mev = max(cfg.T_start, cep.T_cep_MeV - something(cep.uncertainty_T_MeV, 0.0))
    end

    return (T_min_MeV=T_min_mev, T_max_MeV=T_max_mev)
end

function run_production_phase_pipeline(model_kind::Symbol=:PNJL;
        T_start::Float64,
        T_end::Float64,
        dT::Float64=5.0,
        rho_grid::AbstractVector=[0.2],
        xi::Real=0.0,
        output_dir::Union{Nothing, String}=nothing,
        profile=:default,
        run_id::Union{Nothing, String}=nothing,
        policy::Symbol=:processed_first,
        solver_backend::Symbol=:models,
        reverse_rho::Bool=true,
        seed_policy::Symbol=:hybrid_continuity,
        p_num::Int=12,
        t_num::Int=4,
        thermo_quadrature_policy::Symbol=:tensor_gauss,
        thermo_quadrature_rtol::Float64=1e-8,
        thermo_quadrature_atol::Float64=1e-10,
        thermo_quadrature_maxevals::Int=10^7,
        iterations::Int=80,
        compute_crossover::Bool=false,
        crossover_method::Symbol=:peak,
        crossover_variable::Symbol=:phi_u,
        crossover_n_mu::Int=12,
        crossover_mu0_only::Bool=false,
        crossover_T_max_MeV::Float64=NaN,
        cep_tol::Float64=0.1,
        temperature_resolution_target_MeV::Float64=NaN,
        cep_max_bisect_iter::Int=20,
        area_tol_good::Float64=1e-4,
        area_tol_bad::Float64=5e-4,
        unknown_budget::Int=5,
        promote_reference::Bool=false,
        promotion_gate_options::NamedTuple=(;),
        cep_adaptive_rho::Bool=true,
        cep_adaptive_slope_tol::Float64=5.0,
        cep_adaptive_min_gap::Float64=0.002,
        cep_adaptive_max_points::Int=32,
        cep_adaptive_digits::Int=6,
        cep_max_refine_level_rho::Int=2,
        rho_geometry_convergence::Bool=true,
        rho_position_tol_MeV::Float64=0.05,
        rho_density_tol::Float64=0.005,
        rho_maxwell_area_tol::Float64=1e-4,
        rho_refinement_policy::Symbol=:uniform_nested,
        rho_support_fine_step::Float64=0.025,
        rho_support_targeted_cap::Int=12,
        rho_support_config::RhoSupportRefinement.RhoSupportConfig=RhoSupportRefinement.RhoSupportConfig(),
        work_telemetry::Union{Nothing, SolverWorkTelemetry}=nothing,
        memoize_uniform::Bool=false,
        adaptive_temperature::Bool=false,
        temperature_max_refine_level::Int=2,
        temperature_position_tol_MeV::Float64=0.10,
        temperature_density_tol::Float64=0.01,
        temperature_maxwell_area_tol::Float64=1e-4)

    thermo_quadrature_kwargs = _phase_thermo_quadrature_kwargs(
        thermo_quadrature_policy,
        thermo_quadrature_rtol,
        thermo_quadrature_atol,
        thermo_quadrature_maxevals,
    )
    model_kind === :PNJL && PNJLIntegrals.validate_rs_anisotropy(xi)
    if thermo_quadrature_policy === :rs_reduced_adaptive && model_kind !== :PNJL
        throw(ArgumentError("rs_reduced_adaptive thermal quadrature is supported only for model_kind=:PNJL"))
    end
    T_start > 0.0 || throw(ArgumentError(
        "production phase pipeline requires T_start > 0 MeV; strict T=0 PNJL five-field solves are Polyakov-degenerate",
    ))
    T_end >= T_start || throw(ArgumentError("T_end must be >= T_start, got $(T_end) < $(T_start)"))
    resolution_target = isfinite(temperature_resolution_target_MeV) ? temperature_resolution_target_MeV : cep_tol
    resolution_target > 0 || throw(ArgumentError("temperature_resolution_target_MeV must be positive, got $(resolution_target)"))
    cep_max_bisect_iter > 0 || throw(ArgumentError("cep_max_bisect_iter must be positive"))
    unknown_budget >= 0 || throw(ArgumentError("unknown_budget must be nonnegative, got $(unknown_budget)"))
    cep_max_refine_level_rho >= 0 || throw(ArgumentError("cep_max_refine_level_rho must be nonnegative"))
    rho_refinement_policy in (:uniform_nested, :rho_support_cascade, :rho_support_hybrid) || throw(ArgumentError(
        "rho_refinement_policy must be :uniform_nested, :rho_support_cascade, or :rho_support_hybrid, got $(rho_refinement_policy)",
    ))
    if rho_refinement_policy in (:rho_support_cascade, :rho_support_hybrid)
        rho_support_targeted_cap >= rho_support_config.target_point_count || throw(ArgumentError(
            "rho_support_targeted_cap must cover rho_support_config.target_point_count",
        ))
    end
    rho_support_fine_step > 0 || throw(ArgumentError("rho_support_fine_step must be positive"))
    if rho_refinement_policy in (:rho_support_cascade, :rho_support_hybrid)
        model_kind === :PNJL || throw(ArgumentError("$(rho_refinement_policy) is currently supported only for model_kind=:PNJL"))
        rho_geometry_convergence || throw(ArgumentError("$(rho_refinement_policy) requires rho_geometry_convergence=true"))
        required_level = rho_refinement_policy === :rho_support_hybrid ? 4 : 1
        cep_max_refine_level_rho == required_level || throw(ArgumentError(
            "$(rho_refinement_policy) requires cep_max_refine_level_rho=$(required_level)",
        ))
        if rho_refinement_policy === :rho_support_hybrid
            rho_support_targeted_cap <= 12 || throw(ArgumentError(
                "rho_support_hybrid Stage A targeted cap must be <= 12",
            ))
        end
    end
    if rho_geometry_convergence && cep_max_refine_level_rho < 1
        throw(ArgumentError("rho geometry convergence requires cep_max_refine_level_rho >= 1"))
    end
    temperature_max_refine_level >= 0 || throw(ArgumentError("temperature_max_refine_level must be nonnegative"))
    if isfinite(crossover_T_max_MeV)
        crossover_T_max_MeV >= T_start || throw(ArgumentError(
            "crossover_T_max_MeV must be >= T_start, got $(crossover_T_max_MeV)",
        ))
    elseif !isnan(crossover_T_max_MeV)
        throw(ArgumentError("crossover_T_max_MeV must be finite or NaN (use NaN to inherit T_end)"))
    end

    target = resolve_phase_output_target(model_kind; profile=profile, run_id=run_id, policy=policy)
    run_dir = isnothing(output_dir) ? target.run_dir : output_dir
    effective_run_id = isnothing(run_id) ? target.run_id : run_id
    mkpath(run_dir)

    cfg = ProductionPipelineConfig(
        T_start=Float64(T_start),
        T_end=Float64(T_end),
        dT_initial=Float64(dT),
        temperature_resolution_target_MeV=Float64(resolution_target),
        cep_tol_MeV=Float64(resolution_target),
        cep_max_bisect_iter=Int(cep_max_bisect_iter),
        area_tol_good=Float64(area_tol_good),
        area_tol_bad=Float64(area_tol_bad),
        unknown_budget=Int(unknown_budget),
        max_refine_level_rho=Int(cep_max_refine_level_rho),
        adaptive_rho=Bool(cep_adaptive_rho),
        adaptive_slope_tol=Float64(cep_adaptive_slope_tol),
        adaptive_min_gap=Float64(cep_adaptive_min_gap),
        adaptive_max_points=Int(cep_adaptive_max_points),
        adaptive_digits=Int(cep_adaptive_digits),
        rho_geometry_convergence=Bool(rho_geometry_convergence),
        rho_position_tol_MeV=Float64(rho_position_tol_MeV),
        rho_density_tol=Float64(rho_density_tol),
        rho_maxwell_area_tol=Float64(rho_maxwell_area_tol),
        rho_refinement_policy=rho_refinement_policy,
        rho_support_fine_step=Float64(rho_support_fine_step),
        rho_support_targeted_cap=Int(rho_support_targeted_cap),
        rho_support_config=rho_support_config,
        adaptive_temperature=Bool(adaptive_temperature),
        temperature_max_refine_level=Int(temperature_max_refine_level),
        temperature_position_tol_MeV=Float64(temperature_position_tol_MeV),
        temperature_density_tol=Float64(temperature_density_tol),
        temperature_maxwell_area_tol=Float64(temperature_maxwell_area_tol),
        crossover_T_max_MeV=Float64(crossover_T_max_MeV),
    )
    _validate_phase_geometry_tolerances(PhaseGeometryTolerances(
        position_MeV=cfg.rho_position_tol_MeV,
        density=cfg.rho_density_tol,
        maxwell_area=cfg.rho_maxwell_area_tol,
    ))
    _validate_phase_geometry_tolerances(PhaseGeometryTolerances(
        position_MeV=cfg.temperature_position_tol_MeV,
        density=cfg.temperature_density_tol,
        maxwell_area=cfg.temperature_maxwell_area_tol,
    ))

    temps = _production_temperature_grid(cfg.T_start, cfg.T_end, cfg.dT_initial)
    rho_base = collect(Float64.(rho_grid))
    if cfg.rho_refinement_policy in (:rho_support_cascade, :rho_support_hybrid)
        sort!(rho_base)
        unique!(rho_base)
        length(rho_base) >= 2 || throw(ArgumentError("rho_grid must contain at least two unique values"))
        _validate_rho_support_grid(rho_base, cfg.rho_support_fine_step)
    end
    eval_dir = joinpath(run_dir, "production_eval")
    aggregate_csv = joinpath(run_dir, "trho_scan.csv")
    mkpath(eval_dir)

    rho_session = (cfg.rho_refinement_policy in (:rho_support_cascade, :rho_support_hybrid) || memoize_uniform) ?
        TrhoScan.new_rho_point_session(
            model_kind=model_kind,
            reverse_rho=reverse_rho,
            seed_policy=seed_policy,
            solver_backend=solver_backend,
            p_num=p_num,
            t_num=t_num,
            iterations=iterations,
            thermo_quadrature_policy=thermo_quadrature_policy,
            thermo_quadrature_rtol=thermo_quadrature_rtol,
            thermo_quadrature_atol=thermo_quadrature_atol,
            thermo_quadrature_maxevals=thermo_quadrature_maxevals,
            telemetry=work_telemetry,
        ) : nothing
    eval_cache = Dict{Float64, NamedTuple}()
    function rho_support_prior(T::Float64)
        candidates = NamedTuple[]
        for (cached_T, record) in eval_cache
            center = get(record, :cascade_support_center, nothing)
            gap = get(record, :cascade_support_gap, nothing)
            center === nothing || gap === nothing || continue
            isfinite(Float64(center)) && isfinite(Float64(gap)) && Float64(gap) > 0 || continue
            push!(candidates, (distance=abs(cached_T - T), T=cached_T, center=Float64(center), gap=Float64(gap)))
        end
        isempty(candidates) && return nothing
        sort!(candidates; by=item -> (item.distance, item.T))
        selected = first(candidates)
        return RhoSupportRefinement.RhoSupportPrior(selected.T, selected.center, selected.gap)
    end
    function evaluate_temperature(T::Float64)
        key = round(Float64(T); digits=8)
        if haskey(eval_cache, key)
            return eval_cache[key]
        end
        value = _production_classify_temperature(
            aggregate_csv,
            eval_dir,
            Float64(T),
            rho_base,
            Float64(xi),
            reverse_rho,
            seed_policy,
            solver_backend,
            model_kind,
            p_num,
            t_num,
            thermo_quadrature_kwargs,
            iterations,
            cfg,
            rho_session=rho_session,
            rho_prior=cfg.rho_refinement_policy in (:rho_support_cascade, :rho_support_hybrid) ? rho_support_prior(Float64(T)) : nothing,
            work_telemetry=work_telemetry,
        )
        eval_cache[key] = value
        return value
    end

    unknown_count_fn() = count(res -> get(res, :status, :unknown) == :unknown, values(eval_cache))
    unknown_budget_exhausted() = unknown_count_fn() > cfg.unknown_budget

    first_point_fallback = false
    fallback_start_T_MeV = NaN
    start_idx = 1
    first_eval = evaluate_temperature(temps[1])
    if first_eval.reason == "no_curve"
        for idx in 2:length(temps)
            probe = evaluate_temperature(temps[idx])
            if probe.reason != "no_curve"
                first_point_fallback = true
                fallback_start_T_MeV = temps[idx]
                start_idx = idx
                break
            end
        end
    end

    sweep_order = vcat(temps[start_idx:end], reverse(temps[1:(start_idx - 1)]))
    for T in sweep_order
        evaluate_temperature(Float64(T))
    end
    resolved_temps, temperature_convergence_records =
        _adaptive_production_temperature_refinement!(temps, evaluate_temperature, cfg;
            stop_refinement=unknown_budget_exhausted)

    sweep_records = NamedTuple[]
    boundary = NamedTuple[]
    spinodal = NamedTuple[]
    unknown_count = 0
    forced_invalid_count = 0

    for T in sort(resolved_temps)
        res = evaluate_temperature(Float64(T))
        status = get(res, :slice_status, :ambiguous_near_critical)
        raw_status = res.status
        reason = res.reason
        if raw_status == :unknown
            unknown_count += 1
        end

        if status == :confirmed_first_order && raw_status == :valid
            _push_production_boundary!(boundary, T, res)
            _push_production_spinodal!(spinodal, T, res)
        end

        push!(sweep_records, (
            T_MeV=Float64(T),
            status=status,
            raw_status=raw_status,
            coarse_status=Symbol(get(res, :coarse_status, raw_status)),
            fine_status=Symbol(get(res, :fine_status, raw_status)),
            coarse_reason=String(get(res, :coarse_reason, reason)),
            fine_reason=String(get(res, :fine_reason, reason)),
            cascade_status=Symbol(get(res, :cascade_status, :not_run)),
            cascade_targeted_count=Int(get(res, :cascade_targeted_count, 0)),
            stage_used=Symbol(get(res, :stage_used, :uniform_nested)),
            hybrid_upgrade_reason=String(get(res, :hybrid_upgrade_reason, "not_applicable")),
            hybrid_stage_a_status=Symbol(get(res, :stage_a_status, :not_run)),
            hybrid_stage_b_status=Symbol(get(res, :stage_b_status, :not_run)),
            hybrid_stage_c_status=Symbol(get(res, :stage_c_status, :not_run)),
            hybrid_support_low=Float64(something(get(res, :hybrid_support_low, nothing), NaN)),
            hybrid_support_high=Float64(something(get(res, :hybrid_support_high, nothing), NaN)),
            hybrid_support_source=String.(get(res, :hybrid_support_source, Symbol[])),
            hybrid_verification_point_count=Int(get(res, :hybrid_verification_point_count, 0)),
            geometry_converged=Bool(get(res, :geometry_converged, false)),
            position_error_MeV=Float64(get(res, :position_error_MeV, Inf)),
            density_error=Float64(get(res, :density_error, Inf)),
            maxwell_area_gate=Float64(get(res, :maxwell_area_gate, Inf)),
            area_residual=Float64(something(res.area_residual, NaN)),
            rho_hadron=Float64(something(res.rho_hadron, NaN)),
            rho_quark=Float64(something(res.rho_quark, NaN)),
            mu_spinodal_hadron_MeV=Float64(something(res.sres.mu_spinodal_hadron, NaN)),
            mu_spinodal_quark_MeV=Float64(something(res.sres.mu_spinodal_quark, NaN)),
            rho_spinodal_hadron=Float64(something(res.sres.rho_spinodal_hadron, NaN)),
            rho_spinodal_quark=Float64(something(res.sres.rho_spinodal_quark, NaN)),
            stats_total=Int(get(res.stats, :total, 0)),
            stats_success=Int(get(res.stats, :success, 0)),
            stats_failure=Int(get(res.stats, :failure, 0)),
            cache_stats=get(res, :cache_stats, (point_requests=0, cache_hits=0, unique_solves=0, targeted_additions=0)),
            mu_transition_MeV=Float64(something(res.mu_transition, NaN)),
            reason=String(reason),
            level=Int(res.level),
        ))
    end

    sort!(sweep_records; by=r -> r.T_MeV)
    sort!(boundary; by=r -> r.T_MeV)
    sort!(spinodal; by=r -> r.T_MeV)
    bracket = _find_production_bracket(sweep_records)
    frontiers = _refine_production_cep_frontiers(
        bracket,
        evaluate_temperature,
        cfg;
        unknown_count_fn=unknown_count_fn,
    )
    frontier_budget_exhausted = frontiers !== nothing && frontiers.budget_exhausted
    unknown_count = unknown_count_fn()
    sweep_result = _materialize_sweep_result(
        sweep_records,
        first_point_fallback,
        fallback_start_T_MeV,
        unknown_count,
        forced_invalid_count,
    )
    cep = if frontiers === nothing
        has_first_order = any(row -> row.status == :confirmed_first_order, sweep_records)
        has_monotone = any(row -> row.status == :confirmed_monotone, sweep_records)
        last_first_order = findlast(row -> row.status == :confirmed_first_order, sweep_records)
        low_T = last_first_order === nothing ? NaN : Float64(sweep_records[last_first_order].T_MeV)
        low_mu = last_first_order === nothing ? NaN : Float64(sweep_records[last_first_order].mu_transition_MeV)
        base_reason = if has_first_order
            "confirmed_first_order_without_monotone_anchor"
        elseif has_monotone
            "no_confirmed_first_order"
        else
            "no_confirmed_phase_anchor"
        end
        CEPResult(
            result_status=has_first_order ? :ambiguous : (has_monotone ? :not_found : :ambiguous),
            T_last_first_order_MeV=low_T,
            mu_last_first_order_MeV=low_mu,
            eval_count=length(eval_cache),
            unknown_count=unknown_count,
            reason=unknown_budget_exhausted() ? "unknown_budget_exceeded:$(base_reason)" : base_reason,
            method=:production_three_state_frontier,
            temperature_resolution_target_MeV=cfg.temperature_resolution_target_MeV,
        )
    else
        low = frontiers.low
        high = frontiers.high
        CEPResult(
            result_status=:ambiguous,
            T_bracket_low_MeV=low.T,
            T_bracket_high_MeV=high.T,
            bracket_width_T_MeV=high.T - low.T,
            T_last_first_order_MeV=low.T,
            mu_last_first_order_MeV=low.mu,
            T_first_monotone_MeV=high.T,
            ambiguity_width_T_MeV=high.T - low.T,
            eval_count=length(eval_cache),
            unknown_count=unknown_count,
            reason=(frontier_budget_exhausted || unknown_budget_exhausted()) ?
                "unknown_budget_exceeded:ambiguous_interval_between_confirmed_first_order_and_monotone" :
                "ambiguous_interval_between_confirmed_first_order_and_monotone",
            method=:production_three_state_frontier,
            temperature_resolution_target_MeV=cfg.temperature_resolution_target_MeV,
        )
    end

    crossover = if compute_crossover
        mu_max = if cep.found && isfinite(cep.mu_cep_MeV)
            cep.mu_cep_MeV
        elseif !isempty(boundary)
            maximum(row.mu_transition_MeV for row in boundary if isfinite(row.mu_transition_MeV))
        else
            300.0
        end
        crossover_bounds = _production_crossover_temperature_bounds(cfg, cep, boundary)
        T_min_mev = crossover_bounds.T_min_MeV
        T_max_mev = crossover_bounds.T_max_MeV
        if T_min_mev < T_max_mev
            build_crossover_line(;
                mu_max_MeV=mu_max,
                T_min_MeV=T_min_mev,
                T_max_MeV=T_max_mev,
                xi=Float64(xi),
                n_mu=crossover_n_mu,
                mu0_only=crossover_mu0_only,
                method=crossover_method,
                variable=crossover_variable,
                model_kind=model_kind,
                solver_backend=solver_backend,
                p_num=p_num,
                t_num=t_num,
                thermo_quadrature_kwargs...,
            )
        else
            NamedTuple[]
        end
    else
        NamedTuple[]
    end

    scan_total = sum((res.stats.total for res in values(eval_cache)); init=0)
    scan_success = sum((res.stats.success for res in values(eval_cache)); init=0)
    scan_failure = sum((res.stats.failure for res in values(eval_cache)); init=0)
    grid_convergence_records = NamedTuple[]
    for record in temperature_convergence_records
        push!(grid_convergence_records, merge((T_MeV=Float64(record.midpoint), xi=Float64(xi)), record))
    end
    for T in sort(collect(keys(eval_cache)))
        res = eval_cache[T]
        for record in res.rho_convergence_records
            push!(grid_convergence_records, merge((T_MeV=Float64(T), xi=Float64(xi)), record))
        end
    end
    rho_unconverged_count = count(record -> record.axis == "rho" && !record.converged, grid_convergence_records)
    temperature_unconverged_count = count(record -> record.axis == "temperature" && !record.converged, grid_convergence_records)

    config_snapshot = Dict(
        "model_kind" => String(model_kind),
        "profile" => String(profile),
        "mode" => "production",
        "xi" => Float64(xi),
        "T_start" => cfg.T_start,
        "T_end" => cfg.T_end,
        "dT_initial" => cfg.dT_initial,
        "rho_grid" => rho_base,
        "solver_backend" => String(solver_backend),
        "p_num" => p_num,
        "t_num" => t_num,
        "iterations" => iterations,
        "thermo_quadrature_policy" => String(thermo_quadrature_policy),
        "thermo_quadrature_rtol" => thermo_quadrature_rtol,
        "thermo_quadrature_atol" => thermo_quadrature_atol,
        "thermo_quadrature_maxevals" => thermo_quadrature_maxevals,
        "unknown_budget" => cfg.unknown_budget,
        "temperature_resolution_target_MeV" => cfg.temperature_resolution_target_MeV,
        "cep_tol_MeV" => cfg.cep_tol_MeV,
        "max_refine_level_rho" => cfg.max_refine_level_rho,
        "adaptive_rho" => cfg.adaptive_rho,
        "rho_geometry_convergence" => cfg.rho_geometry_convergence,
        "rho_position_tol_MeV" => cfg.rho_position_tol_MeV,
        "rho_density_tol" => cfg.rho_density_tol,
        "rho_maxwell_area_tol" => cfg.rho_maxwell_area_tol,
        "rho_refinement_policy" => String(cfg.rho_refinement_policy),
        "rho_support_fine_step" => cfg.rho_support_fine_step,
        "rho_support_targeted_cap" => cfg.rho_support_targeted_cap,
        "rho_support_config" => Dict(
            "support_slope_tol" => cfg.rho_support_config.support_slope_tol,
            "positive_slope_margin" => cfg.rho_support_config.positive_slope_margin,
            "negative_slope_margin" => cfg.rho_support_config.negative_slope_margin,
            "minimum_negative_secant_run" => cfg.rho_support_config.minimum_negative_secant_run,
            "target_point_count" => cfg.rho_support_config.target_point_count,
            "max_extra_points" => cfg.rho_support_config.max_extra_points,
            "support_expansion_gaps" => cfg.rho_support_config.support_expansion_gaps,
            "local_fit_rmse_tol" => cfg.rho_support_config.local_fit_rmse_tol,
            "near_critical_slope_tol" => cfg.rho_support_config.near_critical_slope_tol,
        ),
        "adaptive_temperature" => cfg.adaptive_temperature,
        "temperature_max_refine_level" => cfg.temperature_max_refine_level,
        "temperature_position_tol_MeV" => cfg.temperature_position_tol_MeV,
        "temperature_density_tol" => cfg.temperature_density_tol,
        "temperature_maxwell_area_tol" => cfg.temperature_maxwell_area_tol,
        "crossover_mu0_only" => crossover_mu0_only,
        "crossover_T_max_MeV" => (isfinite(cfg.crossover_T_max_MeV) ? cfg.crossover_T_max_MeV : cfg.T_end),
    )
    hash_options = (
        mode=:production,
        profile=profile,
        xi=xi,
        T_start=cfg.T_start,
        T_end=cfg.T_end,
        dT=cfg.dT_initial,
        rho_grid=join(rho_base, ","),
        solver_backend=solver_backend,
        p_num=p_num,
        t_num=t_num,
        iterations=iterations,
        thermo_quadrature_policy=thermo_quadrature_policy,
        thermo_quadrature_rtol=thermo_quadrature_rtol,
        thermo_quadrature_atol=thermo_quadrature_atol,
        thermo_quadrature_maxevals=thermo_quadrature_maxevals,
        temperature_resolution_target_MeV=cfg.temperature_resolution_target_MeV,
        unknown_budget=cfg.unknown_budget,
        rho_geometry_convergence=cfg.rho_geometry_convergence,
        rho_position_tol_MeV=cfg.rho_position_tol_MeV,
        rho_density_tol=cfg.rho_density_tol,
        rho_maxwell_area_tol=cfg.rho_maxwell_area_tol,
        rho_refinement_policy=cfg.rho_refinement_policy,
        rho_support_fine_step=cfg.rho_support_fine_step,
        rho_support_targeted_cap=cfg.rho_support_targeted_cap,
        rho_support_config=join(string(getproperty(cfg.rho_support_config, field)) for field in fieldnames(typeof(cfg.rho_support_config))),
        adaptive_temperature=cfg.adaptive_temperature,
        temperature_max_refine_level=cfg.temperature_max_refine_level,
        temperature_position_tol_MeV=cfg.temperature_position_tol_MeV,
        temperature_density_tol=cfg.temperature_density_tol,
        temperature_maxwell_area_tol=cfg.temperature_maxwell_area_tol,
        crossover_mu0_only=crossover_mu0_only,
        crossover_T_max_MeV=(isfinite(cfg.crossover_T_max_MeV) ? cfg.crossover_T_max_MeV : cfg.T_end),
    )
    if cfg.rho_refinement_policy === :rho_support_hybrid
        config_snapshot["rho_hybrid_local_step"] = 0.003125
        config_snapshot["rho_hybrid_support_padding"] = 0.025
        config_snapshot["rho_hybrid_stage_chain"] = ["stage_a_cascade", "stage_b_dense", "stage_c_local_oracle"]
        hash_options = merge(hash_options, (rho_hybrid_local_step=0.003125, rho_hybrid_support_padding=0.025))
    end
    config_snapshot["config_hash"] = _config_hash(model_kind; hash_options...)

    session_snapshot = rho_session === nothing ? nothing : TrhoScan.rho_session_snapshot(rho_session)
    diagnostics = Dict{String, Any}(
        "mode" => "production",
        "scan_total" => scan_total,
        "scan_success" => scan_success,
        "scan_failure" => scan_failure,
        "curve_count" => length(eval_cache),
        "boundary_count" => length(boundary),
        "spinodal_count" => length(spinodal),
        "crossover_count" => length(crossover),
        "cep_method" => String(cep.method),
        "cep_result_status" => String(cep.result_status),
        "cep_eval_count" => cep.eval_count,
        "cep_unknown_count" => cep.unknown_count,
        "cep_temperature_resolution_target_MeV" => (isfinite(cep.temperature_resolution_target_MeV) ? cep.temperature_resolution_target_MeV : nothing),
        "cep_T_last_first_order_MeV" => (isfinite(cep.T_last_first_order_MeV) ? cep.T_last_first_order_MeV : nothing),
        "cep_T_first_monotone_MeV" => (isfinite(cep.T_first_monotone_MeV) ? cep.T_first_monotone_MeV : nothing),
        "unknown_budget" => cfg.unknown_budget,
        "unknown_budget_observed" => unknown_count,
        "unknown_budget_exhausted" => unknown_budget_exhausted() || frontier_budget_exhausted,
        "first_point_fallback" => sweep_result.first_point_fallback,
        "fallback_start_T_MeV" => (isfinite(sweep_result.fallback_start_T_MeV) ? sweep_result.fallback_start_T_MeV : nothing),
        "forced_invalid_count" => forced_invalid_count,
        "sweep_temperatures_MeV" => sweep_result.temperatures_MeV,
        "sweep_statuses" => String.(sweep_result.statuses),
        "sweep_reasons" => sweep_result.reasons,
        "sweep_records" => sweep_records,
        "rho_unconverged_count" => rho_unconverged_count,
        "temperature_unconverged_count" => temperature_unconverged_count,
        "grid_convergence_records" => grid_convergence_records,
        "rho_refinement_policy" => String(cfg.rho_refinement_policy),
        "rho_support_fine_step" => cfg.rho_support_fine_step,
        "rho_support_targeted_cap" => cfg.rho_support_targeted_cap,
        "rho_support_cache" => session_snapshot === nothing ? nothing : Dict(
            "point_requests" => session_snapshot.point_requests,
            "cache_hits" => session_snapshot.cache_hits,
            "unique_solves" => session_snapshot.unique_solves,
            "targeted_additions" => session_snapshot.targeted_additions,
            "failed_points" => session_snapshot.failed_points,
        ),
    )
    if cfg.rho_refinement_policy === :rho_support_hybrid
        diagnostics["rho_hybrid_local_step"] = 0.003125
        diagnostics["rho_hybrid_support_padding"] = 0.025
        diagnostics["rho_hybrid_stage_counts"] = Dict(
            String(stage) => count(row -> row.stage_used == stage, sweep_records)
            for stage in unique(row.stage_used for row in sweep_records)
        )
    end

    base_result = PhasePipelineResult(
        model_kind=model_kind,
        model_variant="default",
        xi=Float64(xi),
        run_id=effective_run_id,
        cep=cep,
        first_order_boundary=boundary,
        spinodal=spinodal,
        crossover_line=crossover,
        diagnostics=diagnostics,
        config_snapshot=config_snapshot,
    )

    artifact_paths = build_phase_artifacts(base_result; output_dir=run_dir)

    promotion_status = if promote_reference
        promote_phase_artifacts(run_dir;
            reference_root=target.reference_root,
            gate_options=merge((; profile=String(profile)), promotion_gate_options),
            write_reference=true)
    else
        PromotionResult()
    end

    return PhasePipelineResult(
        model_kind=base_result.model_kind,
        model_variant=base_result.model_variant,
        xi=base_result.xi,
        run_id=base_result.run_id,
        cep=base_result.cep,
        first_order_boundary=base_result.first_order_boundary,
        spinodal=base_result.spinodal,
        crossover_line=base_result.crossover_line,
        diagnostics=base_result.diagnostics,
        config_snapshot=base_result.config_snapshot,
        artifact_paths=artifact_paths,
        promotion_status=promotion_status,
    )
end
