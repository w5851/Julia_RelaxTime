function _pm_normalize_mev_grid(values, name::String)
    values isa AbstractVector || throw(ArgumentError("$name must be an AbstractVector"))
    normalized = collect(Float64, values)
    isempty(normalized) && throw(ArgumentError("$name must not be empty"))
    all(isfinite, normalized) || throw(ArgumentError("$name must contain only finite values"))
    return normalized
end

const PM_FATAL_ENDPOINT_CAUSES = (:max_iter, :nan_guard, :nonconvergence, :hard_constraint_failed)

@inline function _pm_solver_attempt_origin(result)
    if hasproperty(result, :governed_attempt_origin)
        return Symbol(getproperty(result, :governed_attempt_origin))
    elseif hasproperty(result, :fixedrho_joint_attempt_origin)
        return Symbol(getproperty(result, :fixedrho_joint_attempt_origin))
    elseif hasproperty(result, :entropy_attempt_origin)
        return Symbol(getproperty(result, :entropy_attempt_origin))
    elseif hasproperty(result, :sigma_attempt_origin)
        return Symbol(getproperty(result, :sigma_attempt_origin))
    elseif hasproperty(result, :asym_attempt_origin)
        return Symbol(getproperty(result, :asym_attempt_origin))
    end
    return :fallback
end

@inline function _pm_extract_solver_diagnostic(result; seed_source::Union{Symbol,Nothing}=nothing)
    _diag_view(diag) = if diag isa NamedTuple
        diag
    elseif applicable(to_namedtuple, diag)
        to_namedtuple(diag)
    else
        (
            attempt_origin=(hasproperty(diag, :attempt_origin) ? Symbol(getproperty(diag, :attempt_origin)) : :fallback),
            seed_source=(hasproperty(diag, :seed_source) ? getproperty(diag, :seed_source) : nothing),
            hard_constraint_ok=(hasproperty(diag, :hard_constraint_ok) ? getproperty(diag, :hard_constraint_ok) : nothing),
            failed_constraints=(hasproperty(diag, :failed_constraints) ? Symbol.(getproperty(diag, :failed_constraints)) : Symbol[]),
            selection_reason=(hasproperty(diag, :selection_reason) ? Symbol(getproperty(diag, :selection_reason)) : :none),
            endpoint_cause=(hasproperty(diag, :endpoint_cause) ? getproperty(diag, :endpoint_cause) : :nonconvergence),
            continuity_distance=(hasproperty(diag, :continuity_distance) ? getproperty(diag, :continuity_distance) : nothing),
        )
    end

    if hasproperty(result, :diagnostic)
        diag = _diag_view(getproperty(result, :diagnostic))
        if haskey(diag, :seed_source)
            return diag
        end
        return merge(diag, (seed_source=seed_source,))
    end

    endpoint_cause = if Bool(getproperty(result, :converged))
        :converged
    elseif hasproperty(result, :failed_constraints) && !isempty(getproperty(result, :failed_constraints))
        :hard_constraint_failed
    else
        :nonconvergence
    end

    hard_ok = hasproperty(result, :hard_constraint_ok) ? Bool(getproperty(result, :hard_constraint_ok)) : nothing
    failed = hasproperty(result, :failed_constraints) ? Symbol.(getproperty(result, :failed_constraints)) : Symbol[]
    continuity = hasproperty(result, :continuity_distance) ? Float64(getproperty(result, :continuity_distance)) : nothing

    return (
        attempt_origin=_pm_solver_attempt_origin(result),
        seed_source=seed_source,
        hard_constraint_ok=hard_ok,
        failed_constraints=failed,
        endpoint_cause=endpoint_cause,
        continuity_distance=continuity,
    )
end

@inline function _pm_infer_phase_status(diag)::Symbol
    diag = diag isa NamedTuple ? diag : (applicable(to_namedtuple, diag) ? to_namedtuple(diag) : diag)
    hard_ok = get(diag, :hard_constraint_ok, nothing)
    failed = Symbol.(get(diag, :failed_constraints, Symbol[]))
    endpoint = get(diag, :endpoint_cause, nothing)

    if hard_ok === true && isempty(failed) && !(endpoint in PM_FATAL_ENDPOINT_CAUSES)
        return :valid
    elseif hard_ok === false || !isempty(failed)
        return :invalid
    end
    return :unknown
end

function _pm_interpolate_transition_mu(rows)
    length(rows) >= 2 || return nothing
    left = rows[1]
    right = rows[2]
    Δl = Float64(left.delta_omega)
    Δr = Float64(right.delta_omega)
    μl = Float64(left.mu_MeV)
    μr = Float64(right.mu_MeV)
    Δl == Δr && return 0.5 * (μl + μr)
    return μl - Δl * (μr - μl) / (Δr - Δl)
end

function _pm_refine_transition_bracket(rows; mu_refine_step::Float64=0.01)
    length(rows) >= 2 || return rows
    left = rows[1]
    right = rows[2]
    μl = Float64(left.mu_MeV)
    μr = Float64(right.mu_MeV)
    μr > μl || return rows

    refined_grid = collect(μl:mu_refine_step:μr)
    last(refined_grid) < μr && push!(refined_grid, μr)
    Δl = Float64(left.delta_omega)
    Δr = Float64(right.delta_omega)
    Δpl = Float64(left.delta_pressure)
    Δpr = Float64(right.delta_pressure)

    return [(
        mu_MeV=μ,
        delta_omega=Δl + (Δr - Δl) * ((μ - μl) / (μr - μl)),
        delta_pressure=Δpl + (Δpr - Δpl) * ((μ - μl) / (μr - μl)),
        hadron_exists=left.hadron_exists && right.hadron_exists,
        quark_exists=left.quark_exists && right.quark_exists,
        hadron_status=get(left, :hadron_status, :accepted),
        quark_status=get(left, :quark_status, :accepted),
    ) for μ in refined_grid]
end

function _pm_has_bistable_window(rows)
    return any(getproperty(row, :hadron_exists) && getproperty(row, :quark_exists) for row in rows)
end

function _pm_extract_endpoints(rows)
    overlap = [row for row in rows if row.hadron_exists && row.quark_exists]
    isempty(overlap) && return (
        hadron_endpoint_mu_MeV=nothing,
        quark_endpoint_mu_MeV=nothing,
        bistable_window_width_MeV=0.0,
        branch_disappears_first=:none,
    )

    hadron_endpoint = maximum(row.mu_MeV for row in overlap)
    quark_endpoint = minimum(row.mu_MeV for row in overlap)
    branch_disappears_first = :none
    overlap_end = maximum(row.mu_MeV for row in overlap)
    for row in rows
        if row.mu_MeV > overlap_end
            if !row.hadron_exists && row.quark_exists
                branch_disappears_first = :hadron
                break
            elseif row.hadron_exists && !row.quark_exists
                branch_disappears_first = :quark
                break
            end
        end
    end
    return (
        hadron_endpoint_mu_MeV=Float64(hadron_endpoint),
        quark_endpoint_mu_MeV=Float64(quark_endpoint),
        bistable_window_width_MeV=Float64(hadron_endpoint - quark_endpoint),
        branch_disappears_first=branch_disappears_first,
    )
end

function _pm_compare_with_maxwell(mu_transition_pm, mu_transition_maxwell; comparison_mu_tol::Float64=0.05)
    if mu_transition_pm === nothing && mu_transition_maxwell === nothing
        return (comparison_status=:neither, delta_mu_pm_minus_maxwell_MeV=nothing)
    elseif mu_transition_pm !== nothing && mu_transition_maxwell === nothing
        return (comparison_status=:pm_only, delta_mu_pm_minus_maxwell_MeV=nothing)
    elseif mu_transition_pm === nothing && mu_transition_maxwell !== nothing
        return (comparison_status=:maxwell_only, delta_mu_pm_minus_maxwell_MeV=nothing)
    end

    delta = Float64(mu_transition_pm) - Float64(mu_transition_maxwell)
    status = abs(delta) <= comparison_mu_tol ? :both : :neither
    return (comparison_status=status, delta_mu_pm_minus_maxwell_MeV=delta)
end

function _pm_pressure_crosscheck(rows)
    length(rows) >= 2 || return (consistent=false,)
    left = rows[1]
    right = rows[2]
    consistent = signbit(Float64(left.delta_omega)) != signbit(Float64(right.delta_omega)) &&
                 signbit(Float64(left.delta_pressure)) != signbit(Float64(right.delta_pressure))
    return (consistent=consistent,)
end

function _pm_accept_branch_point(result; residual_accept_tol::Float64=1e-6)
    finite_thermo = isfinite(Float64(result.omega)) && isfinite(Float64(result.pressure)) &&
                    isfinite(Float64(result.rho_norm)) && isfinite(Float64(result.residual_norm))
    if !finite_thermo
        return (accepted=false, branch_status=:invalid_thermo)
    end
    if !Bool(result.converged) || Float64(result.residual_norm) > residual_accept_tol
        return (accepted=false, branch_status=:nonconverged)
    end
    return (accepted=true, branch_status=:accepted)
end

function _pm_check_branch_continuity(prev_x, curr_x, prev_rho::Real, curr_rho::Real;
        continuity_x_tol::Float64=0.25,
        continuity_rho_tol::Float64=0.15)
    xdist = sqrt(sum((Float64(a) - Float64(b))^2 for (a, b) in zip(prev_x, curr_x)))
    rhodist = abs(Float64(curr_rho) - Float64(prev_rho))
    if xdist <= continuity_x_tol && rhodist <= continuity_rho_tol
        return (continuity_ok=true, branch_status=:accepted, endpoint_cause=nothing)
    end
    return (continuity_ok=false, branch_status=:branch_jump, endpoint_cause=:branch_jump)
end

function _pm_solve_single_branch(T_MeV::Real, mu_MeV::Real, branch::Symbol, seed_state::AbstractVector{<:Real};
        xi::Real=0.0,
        solver_backend::Symbol=:auto,
        seed_source::Union{Symbol,Nothing}=nothing,
        p_num::Int=24,
        t_num::Int=8,
        residual_accept_tol::Float64=1e-6)
    branch in (:hadron, :quark) || throw(ArgumentError("unsupported branch: $branch"))

    T_fm = Float64(T_MeV) / Main.Constants_PNJL.ħc_MeV_fm
    mu_fm = Float64(mu_MeV) / Main.Constants_PNJL.ħc_MeV_fm

    solver_backend === :models || throw(ArgumentError("solver_backend=:legacy has been removed from PM phase diagnostic path; use solver_backend=:models"))
    result = begin
        model = create_model(:PNJL)
        solve_constraint(
            model,
            FixedMu(),
            T_fm;
            μ_fm=mu_fm,
            seed_guess=Float64.(seed_state),
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            residual_norm_max=residual_accept_tol,
            diagnostic_level=:summary,
        )
    end

    accept = _pm_accept_branch_point(result; residual_accept_tol=residual_accept_tol)
    diagnostic = _pm_extract_solver_diagnostic(result; seed_source=seed_source)
    diag_status = _pm_infer_phase_status(diagnostic)
    accepted = accept.accepted && diag_status == :valid
    branch_status = if accepted
        accept.branch_status
    elseif diag_status == :invalid
        :nonconverged
    else
        accept.branch_status
    end
    return (
        branch=branch,
        T_MeV=Float64(T_MeV),
        mu_MeV=Float64(mu_MeV),
        x_state=Float64.(collect(getproperty(result, :x_state))),
        pressure=Float64(getproperty(result, :pressure)),
        omega=Float64(getproperty(result, :omega)),
        rho_norm=Float64(getproperty(result, :rho_norm)),
        residual_norm=Float64(getproperty(result, :residual_norm)),
        converged=Bool(getproperty(result, :converged)),
        branch_status=branch_status,
        accepted=accepted,
        diagnostic_status=diag_status,
        diagnostic=diagnostic,
        raw_result=result,
    )
end

function _pm_branch_rows_for_temperature(T_MeV::Real, mu_grid, seed_pair::PMSeedPair;
        xi::Real=0.0,
        solver_backend::Symbol=:auto,
        p_num::Int=24,
        t_num::Int=8,
        residual_accept_tol::Float64=1e-6,
        continuity_x_tol::Float64=0.25,
        continuity_rho_tol::Float64=0.15)
    rows = NamedTuple[]
    hadron_seed = copy(seed_pair.hadron_seed0)
    quark_seed = copy(seed_pair.quark_seed0)
    prev_branch_row = Dict{Symbol, Any}()

    for mu_MeV in mu_grid
        for branch in (:hadron, :quark)
            seed_state = branch === :hadron ? hadron_seed : quark_seed
            seed_source = haskey(prev_branch_row, branch) ? :previous_same_branch : :seed0

            solve_row = _pm_solve_single_branch(T_MeV, mu_MeV, branch, seed_state;
                xi=xi,
                solver_backend=solver_backend,
                seed_source=seed_source,
                p_num=p_num,
                t_num=t_num,
                residual_accept_tol=residual_accept_tol)

            continuity = if haskey(prev_branch_row, branch) && solve_row.accepted
                prev = prev_branch_row[branch]
                _pm_check_branch_continuity(
                    prev.x_state,
                    solve_row.x_state,
                    prev.rho_norm,
                    solve_row.rho_norm;
                    continuity_x_tol=continuity_x_tol,
                    continuity_rho_tol=continuity_rho_tol,
                )
            else
                (continuity_ok=true, branch_status=solve_row.branch_status, endpoint_cause=nothing)
            end

            accepted = solve_row.accepted && continuity.continuity_ok
            branch_status = accepted ? solve_row.branch_status : continuity.branch_status
            endpoint_cause = if accepted
                nothing
            elseif solve_row.branch_status == :nonconverged
                get(solve_row.diagnostic, :endpoint_cause, :nonconvergence)
            elseif solve_row.branch_status == :invalid_thermo
                :physical_loss_candidate
            else
                continuity.endpoint_cause
            end

            row = (
                T_MeV=Float64(T_MeV),
                mu_MeV=Float64(mu_MeV),
                branch=branch,
                branch_status=branch_status,
                seed_source=seed_source,
                continuity_ok=Bool(continuity.continuity_ok),
                converged=Bool(solve_row.converged),
                residual_norm=Float64(solve_row.residual_norm),
                rho_norm=Float64(solve_row.rho_norm),
                pressure=Float64(solve_row.pressure),
                omega=Float64(solve_row.omega),
                endpoint_cause=endpoint_cause,
                accepted=accepted,
                x_state=solve_row.x_state,
            )
            push!(rows, row)

            if accepted
                prev_branch_row[branch] = row
                if branch === :hadron
                    hadron_seed = copy(solve_row.x_state)
                else
                    quark_seed = copy(solve_row.x_state)
                end
            end
        end
    end

    return rows
end

function _pm_temperature_summary(T_MeV::Real, rows;
        comparison_mu_tol::Float64=0.05,
        residual_accept_tol::Float64=1e-6,
        xi::Union{Nothing, Float64}=nothing,
        solver_backend::Union{Nothing, Symbol}=nothing,
        p_num::Union{Nothing, Int}=nothing,
        t_num::Union{Nothing, Int}=nothing,
        continuity_x_tol::Float64=0.25,
        continuity_rho_tol::Float64=0.15)
    mu_values = sort(unique(Float64(row.mu_MeV) for row in rows))
    compare_rows = NamedTuple[]

    for mu_MeV in mu_values
        hadron_row = findfirst(row -> row.mu_MeV == mu_MeV && row.branch == :hadron, rows)
        quark_row = findfirst(row -> row.mu_MeV == mu_MeV && row.branch == :quark, rows)
        hadron = isnothing(hadron_row) ? nothing : rows[hadron_row]
        quark = isnothing(quark_row) ? nothing : rows[quark_row]

        if hadron === nothing || quark === nothing
            continue
        end

        push!(compare_rows, (
            mu_MeV=mu_MeV,
            hadron_exists=Bool(hadron.accepted),
            quark_exists=Bool(quark.accepted),
            hadron_status=hadron.branch_status,
            quark_status=quark.branch_status,
            delta_omega=Float64(hadron.omega - quark.omega),
            delta_pressure=Float64(hadron.pressure - quark.pressure),
        ))
    end

    overlap_rows = [row for row in compare_rows if row.hadron_exists && row.quark_exists]
    transition_rows = NamedTuple[]
    if length(overlap_rows) >= 2
        for i in 1:(length(overlap_rows) - 1)
            left = overlap_rows[i]
            right = overlap_rows[i + 1]
            if Float64(left.delta_omega) == 0.0
                push!(transition_rows, left)
                break
            elseif Float64(left.delta_omega) * Float64(right.delta_omega) <= 0.0
                push!(transition_rows, left)
                push!(transition_rows, right)
                break
            end
        end
    end

    refined_transition_rows = length(transition_rows) == 2 ?
        _pm_refine_transition_bracket(transition_rows) : transition_rows

    mu_transition_pm = if length(refined_transition_rows) == 1
        Float64(transition_rows[1].mu_MeV)
    else
        _pm_interpolate_transition_mu(refined_transition_rows)
    end
    endpoints = _pm_extract_endpoints(compare_rows)
    maxwell = if isnothing(xi) || isnothing(solver_backend) || isnothing(p_num) || isnothing(t_num)
        (mu_transition_maxwell_MeV=nothing, comparison_status=:neither, delta_mu_pm_minus_maxwell_MeV=nothing)
    else
        _pm_maxwell_reference_from_rows(rows;
            T_MeV=T_MeV,
            xi=xi,
            solver_backend=solver_backend,
            p_num=p_num,
            t_num=t_num,
        )
    end
    comparison = _pm_compare_with_maxwell(mu_transition_pm, maxwell.mu_transition_maxwell_MeV; comparison_mu_tol=comparison_mu_tol)

    return (
        T_MeV=Float64(T_MeV),
        mu_transition_pm_MeV=mu_transition_pm,
        mu_transition_maxwell_MeV=maxwell.mu_transition_maxwell_MeV,
        delta_mu_pm_minus_maxwell_MeV=comparison.delta_mu_pm_minus_maxwell_MeV,
        hadron_endpoint_mu_MeV=endpoints.hadron_endpoint_mu_MeV,
        quark_endpoint_mu_MeV=endpoints.quark_endpoint_mu_MeV,
        bistable_window_width_MeV=endpoints.bistable_window_width_MeV,
        branch_disappears_first=endpoints.branch_disappears_first,
        comparison_status=comparison.comparison_status,
        comparison_mu_tol_MeV=Float64(comparison_mu_tol),
        residual_accept_tol=Float64(residual_accept_tol),
        continuity_x_tol=Float64(continuity_x_tol),
        continuity_rho_tol=Float64(continuity_rho_tol),
    ), compare_rows
end

function _pm_maxwell_reference_from_rows(rows;
        T_MeV::Real,
        xi::Real,
        solver_backend::Symbol,
        p_num::Int,
        t_num::Int)
    filtered = [row for row in rows if
        get(row, :accepted, false) &&
        get(row, :T_MeV, T_MeV) == T_MeV &&
        get(row, :xi, xi) == xi &&
        get(row, :solver_backend, solver_backend) == solver_backend &&
        get(row, :p_num, p_num) == p_num &&
        get(row, :t_num, t_num) == t_num]

    if length(filtered) < 6
        return (mu_transition_maxwell_MeV=nothing, comparison_status=:neither, delta_mu_pm_minus_maxwell_MeV=nothing)
    end

    sorted = sort(filtered; by=row -> Float64(row.rho_norm))
    mu_vals = Float64[row.mu_MeV for row in sorted]
    rho_vals = Float64[row.rho_norm for row in sorted]
    sres = detect_s_shape(mu_vals, rho_vals; min_points=5)
    sres.has_s_shape || return (mu_transition_maxwell_MeV=nothing, comparison_status=:neither, delta_mu_pm_minus_maxwell_MeV=nothing)
    mres = maxwell_construction(mu_vals, rho_vals; min_samples=6, spinodal_hint=sres)
    mres.converged || return (mu_transition_maxwell_MeV=nothing, comparison_status=:neither, delta_mu_pm_minus_maxwell_MeV=nothing)

    return (
        mu_transition_maxwell_MeV=Float64(mres.mu_transition),
        comparison_status=:maxwell_only,
        delta_mu_pm_minus_maxwell_MeV=nothing,
    )
end

function analyze_pm_branch_competition(;
        T_values,
        mu_grid,
        xi::Real=0.0,
        solver_backend::Symbol=:auto,
        p_num::Int=24,
        t_num::Int=8,
        seed_pair=nothing,
        output_dir::String,
        comparison_mu_tol::Float64=0.05,
        residual_accept_tol::Float64=1e-6,
        continuity_x_tol::Float64=0.25,
        continuity_rho_tol::Float64=0.15)
    T_values_vec = _pm_normalize_mev_grid(T_values, "T_values")
    mu_grid_vec = _pm_normalize_mev_grid(mu_grid, "mu_grid")
    issorted(mu_grid_vec) || throw(ArgumentError("mu_grid must be sorted in ascending order"))

    branch_rows = NamedTuple[]
    temperature_summaries = NamedTuple[]
    comparison_rows = NamedTuple[]

    for T_MeV in T_values_vec
        local_seed_pair = isnothing(seed_pair) ? derive_pm_seed_pair(T_MeV, mu_grid_vec;
            xi=xi,
            solver_backend=solver_backend,
            p_num=p_num,
            t_num=t_num,
            residual_accept_tol=residual_accept_tol) : normalize_pm_seed_pair(seed_pair)

        rows = _pm_branch_rows_for_temperature(T_MeV, mu_grid_vec, local_seed_pair;
            xi=xi,
            solver_backend=solver_backend,
            p_num=p_num,
            t_num=t_num,
            residual_accept_tol=residual_accept_tol,
            continuity_x_tol=continuity_x_tol,
            continuity_rho_tol=continuity_rho_tol)
        rows_with_context = [merge(row, (
            xi=Float64(xi),
            solver_backend=solver_backend,
            p_num=Int(p_num),
            t_num=Int(t_num),
        )) for row in rows]
        append!(branch_rows, rows_with_context)

        summary, _ = _pm_temperature_summary(T_MeV, rows_with_context;
            comparison_mu_tol=comparison_mu_tol,
            residual_accept_tol=residual_accept_tol,
            xi=Float64(xi),
            solver_backend=solver_backend,
            p_num=Int(p_num),
            t_num=Int(t_num),
            continuity_x_tol=continuity_x_tol,
            continuity_rho_tol=continuity_rho_tol)
        push!(temperature_summaries, summary)
        push!(comparison_rows, summary)
    end

    artifact_paths = _pm_write_artifacts(output_dir, branch_rows, temperature_summaries, comparison_rows)
    return (
        branch_rows=branch_rows,
        temperature_summaries=temperature_summaries,
        comparison_rows=comparison_rows,
        artifacts=artifact_paths,
    )
end
