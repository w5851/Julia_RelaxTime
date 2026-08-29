module GapTransportScanPhaseEquilibrium

using StaticArrays

mutable struct LocalPhaseTracker
    boundary_data::Union{Nothing, NamedTuple{(:T_values, :mu_values, :T_CEP, :mu_CEP, :muq_CEP, :muB_CEP, :xi), Tuple{Vector{Float64}, Vector{Float64}, Float64, Float64, Float64, Float64, Float64}}}
    previous_solution::Union{Nothing, Vector{Float64}}
    previous_phase::Symbol
    hadron_seed::Vector{Float64}
    quark_seed::Vector{Float64}
end

function _column_index(header::AbstractVector{<:AbstractString}, names::Tuple{Vararg{String}})
    for name in names
        idx = findfirst(==(name), header)
        idx !== nothing && return idx
    end
    return nothing
end

function load_phase_boundary_data(xi::Float64;
    boundary_path::String=Main.DEFAULT_PHASE_BOUNDARY_PATH,
    cep_path::String=Main.DEFAULT_PHASE_CEP_PATH,
    phase_reference=nothing,
    phase_reference_mode::Symbol=:runtime,
)
    if phase_reference !== nothing
        phase_reference_mode in (:runtime, :strict, :diagnostic) ||
            throw(ArgumentError("phase_reference_mode must be :runtime, :strict, or :diagnostic"))
        return Main.PhaseReferenceAdapter.boundary_data(
            phase_reference,
            xi;
            require_certified=phase_reference_mode in (:runtime, :strict),
        )
    end
    T_CEP = NaN
    muq_CEP = NaN
    muB_CEP = NaN
    if isfile(cep_path)
        header = SubString{String}[]
        for line in eachline(cep_path)
            if startswith(line, "xi")
                header = split(strip(line), ',')
                continue
            end
            parts = split(line, ',')
            length(parts) >= 3 || continue
            xi_val = tryparse(Float64, parts[1])
            xi_val === nothing && continue
            abs(xi_val - xi) > 1e-6 && continue

            idx_T = isempty(header) ? 2 : something(_column_index(header, ("T_CEP_MeV",)), 2)
            idx_muq = isempty(header) ? 3 : _column_index(header, ("muq_CEP_MeV", "mu_CEP_MeV"))
            idx_muB = isempty(header) ? nothing : _column_index(header, ("muB_CEP_MeV",))

            T_val = idx_T <= length(parts) ? tryparse(Float64, parts[idx_T]) : nothing
            muq_val = idx_muq !== nothing && idx_muq <= length(parts) ? tryparse(Float64, parts[idx_muq]) : nothing
            muB_val = idx_muB !== nothing && idx_muB <= length(parts) ? tryparse(Float64, parts[idx_muB]) : nothing

            T_CEP = T_val === nothing ? NaN : T_val
            if muq_val !== nothing
                muq_CEP = muq_val
                muB_CEP = muB_val === nothing ? 3.0 * muq_CEP : muB_val
            elseif muB_val !== nothing
                muB_CEP = muB_val
                muq_CEP = muB_CEP / 3.0
            end
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

    return (
        T_values=T_values,
        mu_values=mu_values,
        T_CEP=T_CEP,
        mu_CEP=muq_CEP, # compatibility alias: phase-reference mu is mu_q, not mu_B
        muq_CEP=muq_CEP,
        muB_CEP=muB_CEP,
        xi=Float64(xi),
    )
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

function available_phase_boundary_xis(; phase_reference=nothing)
    phase_reference !== nothing && return Main.PhaseReferenceAdapter.available_xi(phase_reference, :boundary)
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

function nearest_phase_boundary_xi(xi::Float64; phase_reference=nothing)
    xis = available_phase_boundary_xis(; phase_reference=phase_reference)
    isempty(xis) && return nothing
    distances = abs.(xis .- xi)
    return xis[argmin(distances)]
end

function available_phase_crossover_xis(; phase_reference=nothing)
    phase_reference !== nothing && return Main.PhaseReferenceAdapter.available_xi(phase_reference, :crossover)
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

function nearest_phase_crossover_xi(xi::Float64; phase_reference=nothing)
    xis = available_phase_crossover_xis(; phase_reference=phase_reference)
    isempty(xis) && return nothing
    distances = abs.(xis .- xi)
    return xis[argmin(distances)]
end

function build_phase_tracker(xi::Float64, previous_solution=nothing, previous_phase::Symbol=:unknown;
    phase_reference=nothing,
    phase_reference_mode::Symbol=:runtime,
)
    boundary_xi = nearest_phase_boundary_xi(xi; phase_reference=phase_reference)
    boundary_data = boundary_xi === nothing ? nothing : load_phase_boundary_data(
        boundary_xi;
        phase_reference=phase_reference,
        phase_reference_mode=phase_reference_mode,
    )

    hadron_seed = Float64.(Main.Models.HADRON_SEED_5)
    quark_seed = Float64.(Main.Models.QUARK_SEED_5)
    normalized_previous = previous_solution === nothing ? nothing : collect(Float64, previous_solution)

    if normalized_previous !== nothing
        seed5 = _seed_state_5(normalized_previous)
        if previous_phase == :hadron || previous_phase == :crossover
            hadron_seed = copy(seed5)
        elseif previous_phase == :quark
            quark_seed = copy(seed5)
        end
    end

    tracker = LocalPhaseTracker(
        boundary_data,
        normalized_previous,
        previous_phase,
        hadron_seed,
        quark_seed,
    )
    return tracker, boundary_xi
end

function load_crossover_reference(xi::Float64; phase_reference=nothing, phase_reference_mode::Symbol=:runtime)
    if phase_reference !== nothing
        phase_reference_mode in (:runtime, :strict, :diagnostic) ||
            throw(ArgumentError("phase_reference_mode must be :runtime, :strict, or :diagnostic"))
        rows = Main.PhaseReferenceAdapter.crossover_rows(
            phase_reference,
            xi;
            require_certified=phase_reference_mode in (:runtime, :strict),
        )
        isempty(rows) && return nothing, nothing
        return [(mu_MeV=row.muq_MeV, T_crossover_MeV=row.T_MeV) for row in rows],
            first(rows).xi
    end
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

function interpolate_crossover_temperature(xi::Float64, muq_mev::Float64;
    phase_reference=nothing,
    phase_reference_mode::Symbol=:runtime,
)
    data, xi_used = load_crossover_reference(
        xi;
        phase_reference=phase_reference,
        phase_reference_mode=phase_reference_mode,
    )
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

function tracker_phase(tracker::LocalPhaseTracker, T_mev::Float64, muq_mev::Float64, xi::Float64;
    phase_reference=nothing,
    phase_reference_mode::Symbol=:runtime,
)
    boundary_phase = current_phase_hint(tracker, T_mev, muq_mev)
    if boundary_phase in (:hadron, :quark)
        return boundary_phase
    end

    Tc_mev, _ = interpolate_crossover_temperature(
        xi,
        muq_mev;
        phase_reference=phase_reference,
        phase_reference_mode=phase_reference_mode,
    )
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

function phase_structure(tracker, T_mev::Float64, muq_mev::Float64, xi::Float64;
    phase_reference=nothing,
    phase_reference_mode::Symbol=:runtime,
)
    data = tracker.boundary_data
    if data !== nothing && !isempty(data.T_values) && !isnan(data.T_CEP) && T_mev <= data.T_CEP
        return :first_order_possible
    end

    Tc_mev, _ = interpolate_crossover_temperature(
        xi,
        muq_mev;
        phase_reference=phase_reference,
        phase_reference_mode=phase_reference_mode,
    )
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

function solve_with_multiseed_governance(T_fm::Float64, muq_fm::Float64, xi::Float64, opts;
    continuation_seed=nothing,
)
    seeds = Vector{Vector{Float64}}()
    if continuation_seed !== nothing
        push!(seeds, _seed_state_5(continuation_seed))
    end
    push!(seeds, Float64.(Main.Models.HADRON_SEED_5))
    push!(seeds, Float64.(Main.Models.QUARK_SEED_5))
    result = try
        Main.Models.solve_multi(
            Main.PNJL_MODEL,
            Main.Models.FixedMu(),
            T_fm,
            muq_fm;
            seeds=seeds,
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
    return m_u < 0.8 ? :quark : :hadron
end

@inline function _relative_delta(a::Real, b::Real)
    return abs(Float64(a) - Float64(b)) /
        max(abs(Float64(a)), abs(Float64(b)), eps(Float64))
end

function solve_two_branch_candidates(
    T_fm::Float64,
    muq_fm::Float64,
    xi::Float64,
    opts;
    branch_mass_rel_tol::Real=1e-3,
    hadron_seed=Main.Models.HADRON_SEED_5,
    quark_seed=Main.Models.QUARK_SEED_5,
)
    candidate_h = solve_models_equilibrium(T_fm, muq_fm, xi, hadron_seed, opts)
    candidate_q = solve_models_equilibrium(T_fm, muq_fm, xi, quark_seed, opts)
    (candidate_h === nothing || candidate_q === nothing) &&
        throw(ArgumentError("failed to obtain both PNJL phase-branch candidates"))

    high, low = candidate_h.masses[1] >= candidate_q.masses[1] ?
        (candidate_h, candidate_q) :
        (candidate_q, candidate_h)
    distinct = _relative_delta(high.masses[1], low.masses[1]) > Float64(branch_mass_rel_tol)
    delta_omega = Float64(high.omega) - Float64(low.omega)
    stable = delta_omega <= 0.0 ? high : low
    stable_phase = delta_omega <= 0.0 ? :hadron : :quark
    return (
        distinct=distinct,
        high_mass_hadron=high,
        low_mass_quark=low,
        delta_omega_high_minus_low=delta_omega,
        stable=stable,
        stable_phase=stable_phase,
    )
end

function direct_coexistence_anchor(
    muB_MeV::Real,
    reference_T_MeV::Real,
    opts;
    xi::Real=0.0,
    scan_half_width_MeV::Real=6.0,
    scan_step_MeV::Real=0.25,
    bisection_steps::Int=24,
)
    scan_half_width_MeV > 0 || throw(ArgumentError("scan_half_width_MeV must be positive"))
    scan_step_MeV > 0 || throw(ArgumentError("scan_step_MeV must be positive"))
    bisection_steps > 0 || throw(ArgumentError("bisection_steps must be positive"))

    muq_fm = (Float64(muB_MeV) / 3.0) / Main.ħc_MeV_fm
    xi0 = Float64(xi)
    previous = nothing
    lower = nothing
    upper = nothing
    offsets = collect(-Float64(scan_half_width_MeV):Float64(scan_step_MeV):Float64(scan_half_width_MeV))
    for offset in offsets
        T_MeV = Float64(reference_T_MeV) + offset
        candidates = solve_two_branch_candidates(
            T_MeV / Main.ħc_MeV_fm,
            muq_fm,
            xi0,
            opts,
        )
        candidates.distinct || continue
        current = (
            T_MeV=T_MeV,
            delta_omega=candidates.delta_omega_high_minus_low,
        )
        if previous !== nothing && previous.delta_omega * current.delta_omega <= 0.0
            lower = previous
            upper = current
            break
        end
        previous = current
    end
    (lower === nothing || upper === nothing) &&
        throw(ArgumentError("failed to bracket the direct PNJL two-branch coexistence temperature near $(reference_T_MeV) MeV"))

    for _ in 1:bisection_steps
        midpoint_T_MeV = 0.5 * (lower.T_MeV + upper.T_MeV)
        candidates = solve_two_branch_candidates(
            midpoint_T_MeV / Main.ħc_MeV_fm,
            muq_fm,
            xi0,
            opts,
        )
        candidates.distinct || throw(ArgumentError("direct coexistence bisection lost the two-root branch bracket"))
        midpoint = (
            T_MeV=midpoint_T_MeV,
            delta_omega=candidates.delta_omega_high_minus_low,
        )
        if lower.delta_omega * midpoint.delta_omega <= 0.0
            upper = midpoint
        else
            lower = midpoint
        end
    end

    return (
        T_lower_MeV=lower.T_MeV,
        T_upper_MeV=upper.T_MeV,
        T_mid_MeV=0.5 * (lower.T_MeV + upper.T_MeV),
        bracket_width_MeV=upper.T_MeV - lower.T_MeV,
        delta_omega_lower=lower.delta_omega,
        delta_omega_upper=upper.delta_omega,
        reference_T_MeV=Float64(reference_T_MeV),
        reference_xi=xi0,
        p_num=Int(opts.p_num),
        t_num=Int(opts.t_num),
        method=:direct_two_branch_equal_omega_bisection,
    )
end

function _default_coexistence_convergence_numerics(opts)
    p_num = Int(opts.p_num)
    t_num = Int(opts.t_num)
    return (
        p_num=p_num < 24 ? 24 : p_num + 8,
        t_num=t_num < 8 ? 8 : t_num + 2,
    )
end

function certify_coexistence_side_points(
    anchor,
    muB_MeV::Real,
    opts;
    delta_xi_candidates=(0.003, 0.005, 0.007, 0.01, 0.015, 0.02),
    convergence_p_num::Int=_default_coexistence_convergence_numerics(opts).p_num,
    convergence_t_num::Int=_default_coexistence_convergence_numerics(opts).t_num,
    anchor_convergence_tol_MeV::Real=0.1,
)
    base_config = (p_num=Int(opts.p_num), t_num=Int(opts.t_num))
    convergence_config = (p_num=convergence_p_num, t_num=convergence_t_num)
    if convergence_config.p_num < base_config.p_num ||
       convergence_config.t_num < base_config.t_num ||
       convergence_config == base_config
        throw(ArgumentError("coexistence certification requires an independent, non-lower thermodynamic node configuration: base=$(base_config), convergence=$(convergence_config)"))
    end
    node_configs = [
        base_config,
        convergence_config,
    ]
    muq_fm = (Float64(muB_MeV) / 3.0) / Main.ħc_MeV_fm
    convergence_cfg = last(node_configs)
    convergence_anchor = if convergence_cfg.p_num == anchor.p_num && convergence_cfg.t_num == anchor.t_num
        anchor
    else
        direct_coexistence_anchor(
            muB_MeV,
            anchor.T_mid_MeV,
            convergence_cfg;
            xi=anchor.reference_xi,
        )
    end
    anchor_convergence_delta_MeV = convergence_anchor.T_mid_MeV - anchor.T_mid_MeV
    anchor_convergence_certified = abs(anchor_convergence_delta_MeV) <= Float64(anchor_convergence_tol_MeV)
    anchor_convergence_certified || throw(ArgumentError("direct coexistence anchor did not converge across thermodynamic node configurations: delta_T_MeV=$(anchor_convergence_delta_MeV), tolerance=$(anchor_convergence_tol_MeV)"))

    for delta_raw in delta_xi_candidates
        delta = Float64(delta_raw)
        delta > 0 || continue
        evidence = NamedTuple[]
        certified = true
        for cfg in node_configs
            numerics = (p_num=cfg.p_num, t_num=cfg.t_num)
            for T_MeV in (anchor.T_lower_MeV, anchor.T_upper_MeV)
                minus = solve_two_branch_candidates(
                    T_MeV / Main.ħc_MeV_fm,
                    muq_fm,
                    -delta,
                    numerics,
                )
                plus = solve_two_branch_candidates(
                    T_MeV / Main.ħc_MeV_fm,
                    muq_fm,
                    delta,
                    numerics,
                )
                minus_ok = minus.distinct && minus.delta_omega_high_minus_low > 0.0
                plus_ok = plus.distinct && plus.delta_omega_high_minus_low < 0.0
                certified &= minus_ok && plus_ok
                push!(evidence, (
                    p_num=cfg.p_num,
                    t_num=cfg.t_num,
                    T_MeV=T_MeV,
                    delta_xi=delta,
                    minus_delta_omega=minus.delta_omega_high_minus_low,
                    plus_delta_omega=plus.delta_omega_high_minus_low,
                    minus_certified=minus_ok,
                    plus_certified=plus_ok,
                ))
            end
        end
        certified && return (
            certified=true,
            delta_xi=delta,
            minus_xi=-delta,
            plus_xi=delta,
            node_configs=node_configs,
            evidence=evidence,
            convergence_anchor=convergence_anchor,
            anchor_convergence_delta_MeV=anchor_convergence_delta_MeV,
            anchor_convergence_tol_MeV=Float64(anchor_convergence_tol_MeV),
            anchor_convergence_certified=anchor_convergence_certified,
            method=:adaptive_two_sided_phase_certification,
        )
    end
    throw(ArgumentError("failed to certify quark/hadron coexistence-side xi points across thermodynamic node configurations"))
end

function solve_equilibrium_with_diagnostics(T_mev::Float64, muB_mev::Float64, xi::Float64, opts;
    previous_solution=nothing,
    previous_phase::Symbol=:unknown,
    phase_reference=nothing,
    phase_reference_mode::Symbol=:runtime,
)
    T_fm = T_mev / Main.ħc_MeV_fm
    muq_mev = muB_mev / 3.0
    muq_fm = muq_mev / Main.ħc_MeV_fm

    tracker, boundary_xi_used = build_phase_tracker(xi, previous_solution, previous_phase;
        phase_reference=phase_reference,
        phase_reference_mode=phase_reference_mode,
    )
    phase_prev = tracker.previous_phase
    phase_curr_hint = tracker_phase(tracker, T_mev, muq_mev, xi;
        phase_reference=phase_reference,
        phase_reference_mode=phase_reference_mode,
    )
    structure = phase_structure(tracker, T_mev, muq_mev, xi;
        phase_reference=phase_reference,
        phase_reference_mode=phase_reference_mode,
    )
    seed_source = describe_seed_source(tracker, phase_curr_hint)
    seed_state = tracker_seed(tracker, T_fm, muq_fm)

    eq = nothing
    phase_curr = phase_curr_hint
    models_err = nothing

    stable_multiseed_required = structure === :first_order_possible
    hinted_switch_multiseed = tracker.previous_solution !== nothing &&
        is_phase_transition(phase_prev, phase_curr_hint) &&
        phase_curr_hint in (:hadron, :quark)
    selection_policy = :continuation_seed

    if stable_multiseed_required || hinted_switch_multiseed
        eq_multi = solve_with_multiseed_governance(
            T_fm,
            muq_fm,
            xi,
            opts;
            continuation_seed=tracker.previous_solution,
        )
        if eq_multi !== nothing
            eq = eq_multi
            phase_curr = classify_phase_from_solution(eq)
            selection_policy = :pressure_max_under_constraints
            seed_source = stable_multiseed_required ?
                "first_order_stable_multiseed_$(String(phase_curr))" :
                "phase_aware_multiseed_governance_$(String(phase_curr))"
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
    next_phase = phase_curr in (:hadron, :quark) ? phase_curr : classify_phase_from_solution(eq)
    return eq, next_solution, next_phase, (
        seed_source=seed_source,
        phase_prev=phase_prev,
        phase_curr=next_phase,
        phase_curr_hint=phase_curr_hint,
        phase_structure=structure,
        phase_boundary_xi_used=boundary_xi_used,
        equilibrium_selection_policy=selection_policy,
    )
end

export LocalPhaseTracker
export build_phase_tracker
export tracker_phase
export phase_structure
export describe_seed_source
export classify_phase_from_solution
export solve_two_branch_candidates
export direct_coexistence_anchor
export certify_coexistence_side_points
export solve_models_equilibrium
export solve_equilibrium_with_diagnostics
export is_phase_transition

end # module GapTransportScanPhaseEquilibrium
