module GapTransportScanPhaseEquilibrium

using StaticArrays

mutable struct LocalPhaseTracker
    boundary_data::Union{Nothing, NamedTuple{(:T_values, :mu_values, :T_CEP, :mu_CEP, :xi), Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64, Float64}}}
    previous_solution::Union{Nothing, Vector{Float64}}
    previous_phase::Symbol
    hadron_seed::Vector{Float64}
    quark_seed::Vector{Float64}
end

function load_phase_boundary_data(xi::Float64; boundary_path::String=Main.DEFAULT_PHASE_BOUNDARY_PATH, cep_path::String=Main.DEFAULT_PHASE_CEP_PATH)
    T_CEP = NaN
    mu_CEP = NaN
    if isfile(cep_path)
        for line in eachline(cep_path)
            startswith(line, "xi") && continue
            parts = split(line, ',')
            length(parts) >= 3 || continue
            xi_val = tryparse(Float64, parts[1])
            xi_val === nothing && continue
            abs(xi_val - xi) > 1e-6 && continue
            T_CEP = tryparse(Float64, parts[2])
            mu_CEP = tryparse(Float64, parts[3])
            break
        end
    end

    T_values = Float64[]
    mu_values = Float64[]
    if isfile(boundary_path)
        for line in eachline(boundary_path)
            startswith(line, "xi") && continue
            parts = split(line, ',')
            length(parts) >= 3 || continue
            xi_val = tryparse(Float64, parts[1])
            xi_val === nothing && continue
            abs(xi_val - xi) > 1e-6 && continue
            T_val = tryparse(Float64, parts[2])
            mu_val = tryparse(Float64, parts[3])
            (T_val === nothing || mu_val === nothing) && continue
            push!(T_values, T_val)
            push!(mu_values, mu_val)
        end
    end

    if !isempty(T_values)
        order = sortperm(T_values)
        T_values = T_values[order]
        mu_values = mu_values[order]
    end

    return (T_values=T_values, mu_values=mu_values, T_CEP=T_CEP, mu_CEP=mu_CEP, xi=Float64(xi))
end

function interpolate_boundary_mu_c(data, T_mev::Float64)
    if !isnan(data.T_CEP) && T_mev > data.T_CEP
        return NaN
    end
    isempty(data.T_values) && return NaN
    Ts = data.T_values
    mus = data.mu_values
    if T_mev <= Ts[1]
        return mus[1]
    elseif T_mev >= Ts[end]
        return mus[end]
    end
    for i in 1:(length(Ts)-1)
        if Ts[i] <= T_mev <= Ts[i+1]
            t = (T_mev - Ts[i]) / (Ts[i+1] - Ts[i])
            return mus[i] + t * (mus[i+1] - mus[i])
        end
    end
    return NaN
end

function current_phase_hint(tracker::LocalPhaseTracker, T_mev::Float64, muq_mev::Float64)
    data = tracker.boundary_data
    data === nothing && return :unknown
    if !isnan(data.T_CEP) && T_mev > data.T_CEP
        return :crossover
    end
    mu_c_mev = interpolate_boundary_mu_c(data, T_mev)
    isnan(mu_c_mev) && return :unknown
    return muq_mev < mu_c_mev ? :hadron : :quark
end

@inline function local_is_phase_transition(prev_phase::Symbol, current_phase::Symbol)
    return (prev_phase == :hadron && current_phase == :quark) ||
           (prev_phase == :quark && current_phase == :hadron)
end

function tracker_seed(tracker::LocalPhaseTracker, T_fm::Float64, muq_fm::Float64)
    mode = Main.Models.FixedMu()
    if tracker.previous_solution !== nothing
        if length(tracker.previous_solution) == Main.Models.state_dim(mode)
            return copy(tracker.previous_solution)
        elseif length(tracker.previous_solution) >= 5
            return Main.Models.extend_seed(tracker.previous_solution, mode)
        end
    end
    hint = Main.Models.auto_phase_hint(T_fm, muq_fm)
    base = (hint === :quark) ? tracker.quark_seed : tracker.hadron_seed
    return Main.Models.extend_seed(base, mode)
end

function available_phase_boundary_xis()
    if Main.PHASE_BOUNDARY_XI_CACHE[] !== nothing
        return Main.PHASE_BOUNDARY_XI_CACHE[]
    end
    xis = Float64[]
    if isfile(Main.DEFAULT_PHASE_BOUNDARY_PATH)
        for line in eachline(Main.DEFAULT_PHASE_BOUNDARY_PATH)
            startswith(line, "xi") && continue
            parts = split(line, ',')
            length(parts) >= 1 || continue
            xi_val = tryparse(Float64, parts[1])
            xi_val === nothing && continue
            push!(xis, xi_val)
        end
    end
    Main.PHASE_BOUNDARY_XI_CACHE[] = unique(sort(xis))
    return Main.PHASE_BOUNDARY_XI_CACHE[]
end

function nearest_phase_boundary_xi(xi::Float64)
    xis = available_phase_boundary_xis()
    isempty(xis) && return nothing
    distances = abs.(xis .- xi)
    return xis[argmin(distances)]
end

function available_phase_crossover_xis()
    if Main.PHASE_CROSSOVER_XI_CACHE[] !== nothing
        return Main.PHASE_CROSSOVER_XI_CACHE[]
    end
    xis = Float64[]
    if isfile(Main.DEFAULT_PHASE_CROSSOVER_PATH)
        for line in eachline(Main.DEFAULT_PHASE_CROSSOVER_PATH)
            startswith(line, "xi") && continue
            parts = split(line, ',')
            length(parts) >= 1 || continue
            xi_val = tryparse(Float64, parts[1])
            xi_val === nothing && continue
            push!(xis, xi_val)
        end
    end
    Main.PHASE_CROSSOVER_XI_CACHE[] = unique(sort(xis))
    return Main.PHASE_CROSSOVER_XI_CACHE[]
end

function nearest_phase_crossover_xi(xi::Float64)
    xis = available_phase_crossover_xis()
    isempty(xis) && return nothing
    distances = abs.(xis .- xi)
    return xis[argmin(distances)]
end

function load_crossover_reference(xi::Float64)
    isfile(Main.DEFAULT_PHASE_CROSSOVER_PATH) || return nothing, nothing

    header = nothing
    rows = NamedTuple{(:mu_MeV, :T_crossover_MeV), Tuple{Float64, Float64}}[]
    xi_used = xi
    exact_match_found = false
    nearest_xi = nearest_phase_crossover_xi(xi)

    open(Main.DEFAULT_PHASE_CROSSOVER_PATH, "r") do io
        for raw_line in eachline(io)
            line = strip(raw_line)
            isempty(line) && continue
            if header === nothing
                header = split(line, ',')
                continue
            end

            parts = split(line, ',')
            length(parts) == length(header) || continue
            idx_xi = findfirst(==("xi"), header)
            idx_mu = findfirst(==("mu_MeV"), header)
            idx_T = findfirst(==("T_crossover_MeV"), header)
            if idx_T === nothing
                idx_T = findfirst(==("T_crossover_chiral_MeV"), header)
            end
            (idx_xi === nothing || idx_mu === nothing || idx_T === nothing) && return nothing, nothing

            row_xi = tryparse(Float64, parts[idx_xi])
            row_mu = tryparse(Float64, parts[idx_mu])
            row_T = tryparse(Float64, parts[idx_T])
            (row_xi === nothing || row_mu === nothing || row_T === nothing) && continue

            if abs(row_xi - xi) <= 1e-6
                exact_match_found = true
                push!(rows, (mu_MeV=Float64(row_mu), T_crossover_MeV=Float64(row_T)))
                xi_used = Float64(row_xi)
            elseif !exact_match_found && nearest_xi !== nothing && abs(row_xi - nearest_xi) <= 1e-6
                push!(rows, (mu_MeV=Float64(row_mu), T_crossover_MeV=Float64(row_T)))
                xi_used = Float64(row_xi)
            end
        end
    end

    filtered = filter(row -> isfinite(row.mu_MeV) && isfinite(row.T_crossover_MeV), rows)
    isempty(filtered) && return nothing, nothing
    sort!(filtered, by=row -> row.mu_MeV)
    return filtered, xi_used
end

function interpolate_crossover_temperature(xi::Float64, muq_mev::Float64)
    data, xi_used = load_crossover_reference(xi)
    data === nothing && return NaN, xi_used
    length(data) == 1 && return data[1].T_crossover_MeV, xi_used

    if muq_mev <= data[1].mu_MeV
        return data[1].T_crossover_MeV, xi_used
    elseif muq_mev >= data[end].mu_MeV
        return data[end].T_crossover_MeV, xi_used
    end

    for i in 1:(length(data) - 1)
        left = data[i]
        right = data[i + 1]
        if left.mu_MeV <= muq_mev <= right.mu_MeV
            weight = (muq_mev - left.mu_MeV) / (right.mu_MeV - left.mu_MeV)
            return left.T_crossover_MeV + weight * (right.T_crossover_MeV - left.T_crossover_MeV), xi_used
        end
    end

    return NaN, xi_used
end

function tracker_phase(tracker::LocalPhaseTracker, T_mev::Float64, muq_mev::Float64, xi::Float64)
    boundary_phase = current_phase_hint(tracker, T_mev, muq_mev)
    if boundary_phase in (:hadron, :quark)
        return boundary_phase
    end

    Tc_mev, _ = interpolate_crossover_temperature(xi, muq_mev)
    if isfinite(Tc_mev)
        if abs(T_mev - Tc_mev) <= 2.0
            return :crossover
        elseif T_mev < Tc_mev
            return :hadron
        else
            return :quark
        end
    end

    if boundary_phase == :crossover
        return :quark
    end
    return boundary_phase
end

function is_phase_transition(prev_phase::Symbol, current_phase::Symbol)
    return local_is_phase_transition(prev_phase, current_phase)
end

function describe_seed_source(tracker, current_phase::Symbol)
    if tracker.previous_solution === nothing
        return "phase_aware_default_$(String(current_phase))"
    elseif is_phase_transition(tracker.previous_phase, current_phase)
        return "phase_aware_phase_switch_$(String(tracker.previous_phase))_to_$(String(current_phase))"
    else
        return "phase_aware_continuity_$(String(current_phase))"
    end
end

function phase_structure(tracker, T_mev::Float64, muq_mev::Float64, xi::Float64)
    data = tracker.boundary_data
    if data !== nothing && !isempty(data.T_values) && !isnan(data.T_CEP) && T_mev <= data.T_CEP
        return :first_order_possible
    end

    Tc_mev, _ = interpolate_crossover_temperature(xi, muq_mev)
    if isfinite(Tc_mev)
        return T_mev <= Tc_mev ? :crossover_possible : :no_transition
    end
    return :unknown
end

@inline function _seed_state_5(seed_state)
    return Float64.(seed_state[1:min(5, length(seed_state))])
end

function _normalize_equilibrium_result(raw; solver_backend::Symbol=:models)
    Bool(raw.converged) || return nothing
    x_state = SVector{5}(Tuple(Float64.(raw.x_state[1:5])))
    mu_vec = SVector{3}(Tuple(Float64.(raw.mu_vec[1:3])))
    masses = SVector{3}(Tuple(Float64.(raw.masses[1:3])))
    (all(isfinite, masses) && all(>(0.0), masses)) || return nothing
    return (
        converged=true,
        x_state=x_state,
        mu_vec=mu_vec,
        masses=masses,
        iterations=Int(raw.iterations),
        residual_norm=Float64(raw.residual_norm),
        solver_backend=solver_backend,
        omega=Float64(raw.omega),
    )
end

function _solve_fixedmu_via_models_solve(T_fm::Float64, muq_fm::Float64, xi::Float64, seed_state, opts)
    return Main.Models.solve(
        Main.PNJL_MODEL,
        Main.Models.FixedMu(),
        T_fm,
        muq_fm;
        seed_guess=_seed_state_5(seed_state),
        xi=xi,
        p_num=opts.p_num,
        t_num=opts.t_num,
        residual_norm_max=1e-4,
    )
end

function _solve_fixedmu_via_models_constraint(T_fm::Float64, muq_fm::Float64, xi::Float64, seed_state, opts)
    return Main.Models.solve_constraint(
        Main.PNJL_MODEL,
        Main.Models.FixedMu(),
        T_fm;
        μ_fm=muq_fm,
        seed_guess=_seed_state_5(seed_state),
        xi=xi,
        p_num=opts.p_num,
        t_num=opts.t_num,
        residual_norm_max=1e-4,
        physicality_check=Main.Models.is_physical_solution,
    )
end

function solve_models_equilibrium(T_fm::Float64, muq_fm::Float64, xi::Float64, seed_state, opts)
    models_err = nothing
    try
        raw = _solve_fixedmu_via_models_solve(T_fm, muq_fm, xi, seed_state, opts)
        eq = _normalize_equilibrium_result(raw; solver_backend=:models)
        eq !== nothing && return eq
    catch err
        models_err = err
    end

    try
        raw = _solve_fixedmu_via_models_constraint(T_fm, muq_fm, xi, seed_state, opts)
        return _normalize_equilibrium_result(raw; solver_backend=:models)
    catch
        if models_err !== nothing
            rethrow(models_err)
        end
        rethrow()
    end
end

function solve_with_multiseed_governance(T_fm::Float64, muq_fm::Float64, xi::Float64, opts)
    result = try
        Main.Models.solve_multi(
            Main.PNJL_MODEL,
            Main.Models.FixedMu(),
            T_fm,
            muq_fm;
            xi=xi,
            p_num=opts.p_num,
            t_num=opts.t_num,
            residual_norm_max=1e-4,
            evaluate_all_attempts=true,
        )
    catch
        nothing
    end
    result === nothing && return nothing
    return _normalize_equilibrium_result(result; solver_backend=:models_multi)
end

function classify_phase_from_solution(eq)
    m_u = eq.masses[1]
    Phi = eq.x_state[4]
    if m_u < 0.8 || Phi > 0.1
        return :quark
    else
        return :hadron
    end
end

function solve_equilibrium_with_diagnostics(T_mev::Float64, muB_mev::Float64, xi::Float64, opts;
    previous_solution=nothing,
    previous_phase::Symbol=:unknown,
)
    T_fm = T_mev / Main.ħc_MeV_fm
    muq_mev = muB_mev / 3.0
    muq_fm = muq_mev / Main.ħc_MeV_fm

    tracker, boundary_xi_used = build_phase_tracker(xi, previous_solution, previous_phase)
    phase_prev = tracker.previous_phase
    phase_curr_hint = tracker_phase(tracker, T_mev, muq_mev, xi)
    structure = phase_structure(tracker, T_mev, muq_mev, xi)
    seed_source = describe_seed_source(tracker, phase_curr_hint)
    seed_state = tracker_seed(tracker, T_fm, muq_fm)

    eq = nothing
    phase_curr = phase_curr_hint
    models_err = nothing

    if tracker.previous_solution !== nothing && is_phase_transition(phase_prev, phase_curr_hint) && phase_curr_hint in (:hadron, :quark)
        eq_multi = solve_with_multiseed_governance(T_fm, muq_fm, xi, opts)
        if eq_multi !== nothing
            eq = eq_multi
            phase_curr = if phase_curr_hint in (:hadron, :quark)
                phase_curr_hint
            else
                classify_phase_from_solution(eq)
            end
            seed_source = "phase_aware_multiseed_governance_$(String(phase_curr))"
        end
    end

    if eq === nothing
        try
            eq = solve_models_equilibrium(T_fm, muq_fm, xi, seed_state, opts)
        catch err
            models_err = err
        end
    end

    if eq === nothing
        if models_err !== nothing
            rethrow(models_err)
        end
        throw(ArgumentError("models equilibrium solver returned no valid candidate"))
    end

    next_solution = collect(Float64, eq.x_state)
    next_phase = phase_curr
    return eq, next_solution, next_phase, (
        seed_source=seed_source,
        phase_prev=phase_prev,
        phase_curr=next_phase,
        phase_structure=structure,
        phase_boundary_xi_used=boundary_xi_used,
    )
end

export LocalPhaseTracker
export build_phase_tracker
export tracker_phase
export phase_structure
export describe_seed_source
export classify_phase_from_solution
export solve_models_equilibrium
export solve_equilibrium_with_diagnostics
export is_phase_transition

end # module GapTransportScanPhaseEquilibrium
