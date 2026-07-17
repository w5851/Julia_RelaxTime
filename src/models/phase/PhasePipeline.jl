using SHA
using Printf

function _config_hash(model_kind::Symbol; kwargs...)
    payload = string(model_kind, "|", join(sort(collect(string(k) * "=" * string(v) for (k, v) in kwargs)), ";"))
    return bytes2hex(sha1(payload))
end

@inline function _phase_thermo_quadrature_kwargs(
        policy::Symbol,
        rtol::Float64,
        atol::Float64,
        maxevals::Int)
    PNJLIntegrals.validate_thermal_quadrature_policy(policy)
    PNJLIntegrals.validate_thermal_quadrature_controls(rtol, atol, maxevals)
    return (
        thermo_quadrature_policy=policy,
        thermo_quadrature_rtol=rtol,
        thermo_quadrature_atol=atol,
        thermo_quadrature_maxevals=maxevals,
    )
end

function _refine_rho_grid(rho_grid::Vector{Float64}, level::Int)
    level <= 0 && return sort(unique(copy(rho_grid)))
    refined = sort(unique(copy(rho_grid)))
    for _ in 1:level
        if length(refined) < 2
            return refined
        end
        next_grid = Float64[]
        sizehint!(next_grid, 2 * length(refined) - 1)
        for i in 1:(length(refined) - 1)
            left = refined[i]
            right = refined[i + 1]
            push!(next_grid, left)
            push!(next_grid, 0.5 * (left + right))
        end
        push!(next_grid, refined[end])
        refined = next_grid
    end
    return refined
end

function _cep_eval_filename(T::Float64, level::Int)
    t_token = replace(@sprintf("%.6f", T), "." => "p", "-" => "m")
    return "cep_eval_T$(t_token)_L$(level).csv"
end

function _infer_temperature_step(T_grid)
    temps = sort(unique(Float64.(T_grid)))
    if length(temps) < 2
        return 5.0
    end
    diffs = diff(temps)
    positive = filter(>(0.0), diffs)
    return isempty(positive) ? 5.0 : minimum(positive)
end

function _run_single_temperature_scan(
        out_csv::String,
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
        iterations::Int)
    run_trho_scan(;
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
    local_curves = load_curves_from_trho_csv(out_csv; xi=Float64(xi), min_points=3)
    isempty(local_curves) && return nothing
    local_key = first(sort(collect(keys(local_curves))))
    return local_curves[local_key]
end

function _collect_first_order_boundary(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}})
    rows = NamedTuple[]
    for T in sort(collect(keys(curves)))
        mu_vals, rho_vals = curves[T]
        sres = detect_s_shape(mu_vals, rho_vals)
        sres.has_s_shape || continue

        mres = maxwell_construction(mu_vals, rho_vals; spinodal_hint=sres)
        mres.converged || continue

        mu_transition = something(mres.mu_transition, NaN)
        rho_hadron = something(mres.rho_hadron, NaN)
        rho_quark = something(mres.rho_quark, NaN)
        area_residual = something(mres.area_residual, NaN)

        push!(rows, (
            T_MeV=Float64(T),
            mu_transition_MeV=Float64(mu_transition),
            rho_hadron=Float64(rho_hadron),
            rho_quark=Float64(rho_quark),
            area_residual=Float64(area_residual),
            converged=true,
        ))
    end
    return rows
end

function _collect_spinodal(curves::Dict{Float64, Tuple{Vector{Float64}, Vector{Float64}}})
    rows = NamedTuple[]
    for T in sort(collect(keys(curves)))
        mu_vals, rho_vals = curves[T]
        sres = detect_s_shape(mu_vals, rho_vals)
        sres.has_s_shape || continue
        push!(rows, (
            T_MeV=Float64(T),
            mu_spinodal_hadron_MeV=Float64(something(sres.mu_spinodal_hadron, NaN)),
            mu_spinodal_quark_MeV=Float64(something(sres.mu_spinodal_quark, NaN)),
            rho_spinodal_hadron=Float64(something(sres.rho_spinodal_hadron, NaN)),
            rho_spinodal_quark=Float64(something(sres.rho_spinodal_quark, NaN)),
        ))
    end
    return rows
end

function _run_phase_pipeline_core(model_kind::Symbol=:PNJL;
        mode::Symbol=:production,
        T_grid::AbstractVector=[150.0],
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
        cep_strategy::Symbol=:interpolate,
        cep_interpolate_use_direct_eval::Bool=false,
        cep_tol::Float64=0.01,
        cep_max_bisect_iter::Int=20,
        cep_area_tol_good::Float64=1e-4,
        cep_area_tol_bad::Float64=5e-4,
        cep_max_refine_level::Int=2,
        cep_adaptive_rho::Bool=true,
        cep_adaptive_slope_tol::Float64=5.0,
        cep_adaptive_min_gap::Float64=0.002,
        cep_adaptive_max_points::Int=32,
        cep_adaptive_digits::Int=6,
        cep_direct_bracket_mode::Symbol=:directional,
        cep_direct_start::Symbol=:low,
        cep_direct_initial_step::Float64=NaN,
        cep_direct_expand_factor::Float64=2.0,
        cep_direct_max_expand_steps::Int=8,
        cep_direct_fallback_scan::Bool=true,
        promote_reference::Bool=false,
        promotion_gate_options::NamedTuple=(;),
        T_start::Float64=NaN,
        T_end::Float64=NaN,
        dT::Float64=NaN,
        unknown_budget::Int=5)

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
    minimum(Float64.(T_grid)) > 0.0 || throw(ArgumentError(
        "phase pipeline requires T_grid > 0 MeV; strict T=0 PNJL five-field solves are Polyakov-degenerate",
    ))

    if mode == :production
        T_start_eff = isfinite(T_start) ? T_start : Float64(minimum(T_grid))
        T_end_eff = isfinite(T_end) ? T_end : Float64(maximum(T_grid))
        dT_eff = isfinite(dT) ? dT : _infer_temperature_step(T_grid)
        return run_production_phase_pipeline(model_kind;
            T_start=T_start_eff,
            T_end=T_end_eff,
            dT=dT_eff,
            rho_grid=rho_grid,
            xi=xi,
            output_dir=output_dir,
            profile=profile,
            run_id=run_id,
            policy=policy,
            solver_backend=solver_backend,
            reverse_rho=reverse_rho,
            seed_policy=seed_policy,
            p_num=p_num,
            t_num=t_num,
            thermo_quadrature_kwargs...,
            iterations=iterations,
            compute_crossover=compute_crossover,
            crossover_method=crossover_method,
            crossover_variable=crossover_variable,
            crossover_n_mu=crossover_n_mu,
            cep_tol=isfinite(cep_tol) && cep_tol > 0 ? cep_tol : 0.5,
            cep_max_bisect_iter=cep_max_bisect_iter,
            area_tol_good=cep_area_tol_good,
            area_tol_bad=cep_area_tol_bad,
            unknown_budget=unknown_budget,
            promote_reference=promote_reference,
            promotion_gate_options=promotion_gate_options,
        )
    elseif mode != :research
        throw(ArgumentError("Unknown phase pipeline mode: $mode. Use :production or :research"))
    end

    # Research mode: original pipeline logic
    target = resolve_phase_output_target(model_kind; profile=profile, run_id=run_id, policy=policy)
    run_dir = isnothing(output_dir) ? target.run_dir : output_dir
    effective_run_id = isnothing(run_id) ? target.run_id : run_id

    mkpath(run_dir)

    scan_csv = joinpath(run_dir, "trho_scan.csv")
    stats = run_trho_scan(;
        T_values=collect(Float64.(T_grid)),
        rho_values=collect(Float64.(rho_grid)),
        xi_values=[Float64(xi)],
        output_path=scan_csv,
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

    curves = load_curves_from_trho_csv(scan_csv; xi=Float64(xi), min_points=3)
    boundary = _collect_first_order_boundary(curves)
    spinodal = _collect_spinodal(curves)

    cep_evaluator = if cep_strategy == :direct
        base_rho_grid = collect(Float64.(rho_grid))
        eval_dir = joinpath(run_dir, "cep_direct")
        mkpath(eval_dir)
        adaptive_cfg = AdaptiveRhoRefinement.AdaptiveRhoConfig(
            slope_tol=cep_adaptive_slope_tol,
            min_gap=cep_adaptive_min_gap,
            max_points=cep_adaptive_max_points,
            digits=cep_adaptive_digits,
        )
        function (T_mid::Float64, level::Int)
            rho_eval = _refine_rho_grid(base_rho_grid, level)
            out_csv = joinpath(eval_dir, _cep_eval_filename(T_mid, level))

            curve = _run_single_temperature_scan(
                out_csv,
                T_mid,
                rho_eval,
                Float64(xi),
                reverse_rho,
                seed_policy,
                solver_backend,
                model_kind,
                p_num,
                t_num,
                thermo_quadrature_kwargs,
                iterations,
            )
            curve === nothing && return nothing

            if cep_adaptive_rho && level > 0
                for adapt_round in 1:level
                    mu_vals, rho_vals = curve
                    additions = AdaptiveRhoRefinement.suggest_refinement_points(rho_vals, mu_vals; config=adaptive_cfg)
                    isempty(additions) && break
                    rho_eval = AdaptiveRhoRefinement.merge_rho_values(rho_eval, additions; digits=adaptive_cfg.digits)
                    out_csv = joinpath(eval_dir, _cep_eval_filename(T_mid, level) * "_A$(adapt_round).csv")
                    curve = _run_single_temperature_scan(
                        out_csv,
                        T_mid,
                        rho_eval,
                        Float64(xi),
                        reverse_rho,
                        seed_policy,
                        solver_backend,
                        model_kind,
                        p_num,
                        t_num,
                        thermo_quadrature_kwargs,
                        iterations,
                    )
                    curve === nothing && return nothing
                end
            end

            return curve
        end
    elseif cep_strategy == :interpolate && cep_interpolate_use_direct_eval
        base_rho_grid = collect(Float64.(rho_grid))
        eval_dir = joinpath(run_dir, "cep_interpolate_direct")
        mkpath(eval_dir)
        function (T_mid::Float64, level::Int)
            rho_eval = level <= 0 ? base_rho_grid : _refine_rho_grid(base_rho_grid, level)
            out_csv = joinpath(eval_dir, _cep_eval_filename(T_mid, level))
            return _run_single_temperature_scan(
                out_csv,
                T_mid,
                rho_eval,
                Float64(xi),
                reverse_rho,
                seed_policy,
                solver_backend,
                model_kind,
                p_num,
                t_num,
                thermo_quadrature_kwargs,
                iterations,
            )
        end
    else
        nothing
    end

    cep = find_cep(curves;
        strategy=cep_strategy,
        evaluate_at_T=cep_evaluator,
        tol=cep_tol,
        max_bisect_iter=cep_max_bisect_iter,
        area_tol_good=cep_area_tol_good,
        area_tol_bad=cep_area_tol_bad,
        max_refine_level=cep_max_refine_level,
        direct_bracket_mode=cep_direct_bracket_mode,
        direct_start=cep_direct_start,
        direct_initial_step=cep_direct_initial_step,
        direct_expand_factor=cep_direct_expand_factor,
        direct_max_expand_steps=cep_direct_max_expand_steps,
        direct_fallback_scan=cep_direct_fallback_scan)

    crossover = if compute_crossover
        mu_max = if cep.found && isfinite(cep.mu_cep_MeV)
            cep.mu_cep_MeV
        elseif !isempty(boundary)
            maximum(row.mu_transition_MeV for row in boundary if isfinite(row.mu_transition_MeV))
        else
            300.0
        end
        T_min_mev = minimum(Float64.(T_grid))
        T_max_mev = min(maximum(Float64.(T_grid)), 220.0)
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

    config_snapshot = Dict(
        "model_kind" => String(model_kind),
        "profile" => String(profile),
        "mode" => "research",
        "xi" => Float64(xi),
        "T_grid" => collect(Float64.(T_grid)),
        "rho_grid" => collect(Float64.(rho_grid)),
        "solver_backend" => String(solver_backend),
        "thermo_quadrature_policy" => String(thermo_quadrature_policy),
        "thermo_quadrature_rtol" => thermo_quadrature_rtol,
        "thermo_quadrature_atol" => thermo_quadrature_atol,
        "thermo_quadrature_maxevals" => thermo_quadrature_maxevals,
        "cep_strategy" => String(cep_strategy),
        "cep_interpolate_use_direct_eval" => cep_interpolate_use_direct_eval,
        "cep_adaptive_rho" => cep_adaptive_rho,
        "cep_adaptive_slope_tol" => cep_adaptive_slope_tol,
        "cep_adaptive_min_gap" => cep_adaptive_min_gap,
        "cep_adaptive_max_points" => cep_adaptive_max_points,
        "cep_direct_bracket_mode" => String(cep_direct_bracket_mode),
        "cep_direct_start" => String(cep_direct_start),
    )
    config_snapshot["config_hash"] = _config_hash(model_kind;
        mode=:research, profile=profile, xi=xi, T_grid=join(T_grid, ","), rho_grid=join(rho_grid, ","), solver_backend=solver_backend,
        thermo_quadrature_policy=thermo_quadrature_policy, thermo_quadrature_rtol=thermo_quadrature_rtol,
        thermo_quadrature_atol=thermo_quadrature_atol, thermo_quadrature_maxevals=thermo_quadrature_maxevals)

    diagnostics = Dict{String, Any}(
        "scan_total" => getproperty(stats, :total),
        "scan_success" => getproperty(stats, :success),
        "scan_failure" => getproperty(stats, :failure),
        "curve_count" => length(curves),
        "boundary_count" => length(boundary),
        "spinodal_count" => length(spinodal),
        "crossover_count" => length(crossover),
        "cep_method" => String(cep.method),
        "cep_eval_count" => cep.eval_count,
        "cep_unknown_count" => cep.unknown_count,
        "cep_adaptive_rho" => cep_adaptive_rho,
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
