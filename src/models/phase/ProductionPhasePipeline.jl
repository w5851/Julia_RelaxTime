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
        pass::Int)
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
    )
    _append_scan_csv!(aggregate_csv, out_csv)
    local_curves = load_curves_from_trho_csv(out_csv; xi=Float64(xi), min_points=3)
    if isempty(local_curves)
        return (curve=nothing, stats=stats, out_csv=out_csv)
    end
    local_key = first(sort(collect(keys(local_curves))))
    return (curve=local_curves[local_key], stats=stats, out_csv=out_csv)
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
        cfg::ProductionPipelineConfig)
    adaptive_cfg = AdaptiveRhoRefinement.AdaptiveRhoConfig(
        slope_tol=cfg.adaptive_slope_tol,
        min_gap=cfg.adaptive_min_gap,
        max_points=cfg.adaptive_max_points,
        digits=cfg.adaptive_digits,
    )

    rho_eval = copy(base_rho_grid)
    last_result = nothing
    pass = 1

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
        )
        pass += 1
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
                stats=scan.stats,
                curve=nothing,
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
        last_result = (
            status=cres.status,
            mu_transition=cres.mu_transition,
            area_residual=cres.area_residual,
            sres=cres.sres,
            rho_hadron=maxwell.converged ? maxwell.rho_hadron : nothing,
            rho_quark=maxwell.converged ? maxwell.rho_quark : nothing,
            reason=String(cres.reason),
            level=level,
            stats=scan.stats,
            curve=scan.curve,
        )

        if cres.status != :unknown || !cfg.adaptive_rho || level == cfg.max_refine_level_rho
            break
        end

        additions = AdaptiveRhoRefinement.suggest_refinement_points(rho_vals, mu_vals; config=adaptive_cfg)
        isempty(additions) && break
        rho_eval = AdaptiveRhoRefinement.merge_rho_values(rho_eval, additions; digits=adaptive_cfg.digits)
    end

    return last_result
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
    prev_valid = nothing
    for record in records
        if record.status == :valid
            prev_valid = record
        elseif record.status == :invalid && prev_valid !== nothing
            return (T_low=prev_valid.T_MeV, mu_low=prev_valid.mu_transition_MeV, T_high=record.T_MeV)
        end
    end
    return nothing
end

function _production_crossover_temperature_bounds(
        cfg::ProductionPipelineConfig,
        cep::CEPResult,
        boundary::AbstractVector{<:NamedTuple})
    T_max_mev = min(cfg.T_end, 220.0)
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
        cep_tol::Float64=0.5,
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
        cep_max_refine_level_rho::Int=2)

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

    target = resolve_phase_output_target(model_kind; profile=profile, run_id=run_id, policy=policy)
    run_dir = isnothing(output_dir) ? target.run_dir : output_dir
    effective_run_id = isnothing(run_id) ? target.run_id : run_id
    mkpath(run_dir)

    cfg = ProductionPipelineConfig(
        T_start=Float64(T_start),
        T_end=Float64(T_end),
        dT_initial=Float64(dT),
        cep_tol_MeV=Float64(cep_tol),
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
    )

    temps = _production_temperature_grid(cfg.T_start, cfg.T_end, cfg.dT_initial)
    rho_base = collect(Float64.(rho_grid))
    eval_dir = joinpath(run_dir, "production_eval")
    aggregate_csv = joinpath(run_dir, "trho_scan.csv")
    mkpath(eval_dir)

    eval_cache = Dict{Float64, NamedTuple}()
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
        )
        eval_cache[key] = value
        return value
    end

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
    seen_t = Set{Float64}()
    sweep_records = NamedTuple[]
    boundary = NamedTuple[]
    spinodal = NamedTuple[]
    unknown_count = 0
    forced_invalid_count = 0
    scan_total = 0
    scan_success = 0
    scan_failure = 0

    for T in sweep_order
        T_key = round(Float64(T); digits=8)
        T_key in seen_t && continue
        push!(seen_t, T_key)
        res = evaluate_temperature(T)
        scan_total += res.stats.total
        scan_success += res.stats.success
        scan_failure += res.stats.failure

        status = res.status
        reason = res.reason
        if status == :unknown
            unknown_count += 1
            if unknown_count > cfg.unknown_budget
                status = :invalid
                reason = "unknown_budget_exceeded:$(res.reason)"
                forced_invalid_count += 1
            end
        end

        if status == :valid
            _push_production_boundary!(boundary, T, res)
            _push_production_spinodal!(spinodal, T, res)
        end

        push!(sweep_records, (
            T_MeV=Float64(T),
            status=status,
            mu_transition_MeV=Float64(something(res.mu_transition, NaN)),
            area_residual=Float64(something(res.area_residual, NaN)),
            reason=String(reason),
            level=Int(res.level),
        ))
    end

    sort!(sweep_records; by=r -> r.T_MeV)
    sort!(boundary; by=r -> r.T_MeV)
    sort!(spinodal; by=r -> r.T_MeV)
    sweep_result = _materialize_sweep_result(
        sweep_records,
        first_point_fallback,
        fallback_start_T_MeV,
        unknown_count,
        forced_invalid_count,
    )

    bracket = _find_production_bracket(sweep_records)
    cep = CEPResult(method=:production_no_valid_invalid_transition)
    if bracket !== nothing
        T_low = Float64(bracket.T_low)
        T_high = Float64(bracket.T_high)
        last_mu = Float64(bracket.mu_low)
        bisect_count = 0
        refine_unknown_count = 0

        while (T_high - T_low) > cfg.cep_tol_MeV && bisect_count < cfg.cep_max_bisect_iter
            T_mid = 0.5 * (T_low + T_high)
            res = evaluate_temperature(T_mid)
            scan_total += res.stats.total
            scan_success += res.stats.success
            scan_failure += res.stats.failure

            status = res.status
            if status == :unknown
                refine_unknown_count += 1
                if unknown_count + refine_unknown_count > cfg.unknown_budget
                    status = :invalid
                    forced_invalid_count += 1
                else
                    status = :invalid
                end
            end

            if status == :valid
                T_low = T_mid
                if res.mu_transition !== nothing && isfinite(res.mu_transition)
                    last_mu = Float64(res.mu_transition)
                end
            else
                T_high = T_mid
            end
            bisect_count += 1
        end

        cep = CEPResult(
            found=true,
            T_cep_MeV=0.5 * (T_low + T_high),
            mu_cep_MeV=last_mu,
            uncertainty_T_MeV=T_high - T_low,
            eval_count=length(eval_cache),
            unknown_count=unknown_count + refine_unknown_count,
            reason=forced_invalid_count > 0 ? "unknown_budget_forced_invalid_present" : nothing,
            method=:production_bisect_last_valid_maxwell,
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
                method=crossover_method,
                variable=crossover_variable,
                model_kind=model_kind,
                solver_backend=solver_backend,
                thermo_quadrature_kwargs...,
            )
        else
            NamedTuple[]
        end
    else
        NamedTuple[]
    end

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
        "thermo_quadrature_policy" => String(thermo_quadrature_policy),
        "thermo_quadrature_rtol" => thermo_quadrature_rtol,
        "thermo_quadrature_atol" => thermo_quadrature_atol,
        "thermo_quadrature_maxevals" => thermo_quadrature_maxevals,
        "unknown_budget" => cfg.unknown_budget,
        "cep_tol_MeV" => cfg.cep_tol_MeV,
        "max_refine_level_rho" => cfg.max_refine_level_rho,
        "adaptive_rho" => cfg.adaptive_rho,
    )
    config_snapshot["config_hash"] = _config_hash(model_kind;
        mode=:production,
        profile=profile,
        xi=xi,
        T_start=cfg.T_start,
        T_end=cfg.T_end,
        dT=cfg.dT_initial,
        rho_grid=join(rho_base, ","),
        solver_backend=solver_backend,
        thermo_quadrature_policy=thermo_quadrature_policy,
        thermo_quadrature_rtol=thermo_quadrature_rtol,
        thermo_quadrature_atol=thermo_quadrature_atol,
        thermo_quadrature_maxevals=thermo_quadrature_maxevals)

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
        "cep_eval_count" => cep.eval_count,
        "cep_unknown_count" => cep.unknown_count,
        "unknown_budget" => cfg.unknown_budget,
        "first_point_fallback" => sweep_result.first_point_fallback,
        "fallback_start_T_MeV" => (isfinite(sweep_result.fallback_start_T_MeV) ? sweep_result.fallback_start_T_MeV : nothing),
        "forced_invalid_count" => sweep_result.forced_invalid_count,
        "sweep_temperatures_MeV" => sweep_result.temperatures_MeV,
        "sweep_statuses" => String.(sweep_result.statuses),
        "sweep_reasons" => sweep_result.reasons,
    )

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
