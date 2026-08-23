using LinearAlgebra: eigvals, norm, Symmetric
using NLsolve
using ForwardDiff

"""One distinct stationary-point candidate found by the magnetic solver."""
struct MagneticGapCandidate
    seed_index::Int
    seed::SVector{5, Float64}
    x_state::SVector{5, Float64}
    omega::Float64
    residual::SVector{5, Float64}
    residual_norm::Float64
    converged::Bool
    physical::Bool
    method::Symbol
    iterations::Int
    n_max::Int
    stability::Symbol
    branch_label::Symbol
end

"""Branch-aware result for a magnetic five-field equilibrium solve.

`candidates` contains distinct converged roots. The ordinary `solve_gap` API
returns `state`, while branch/stability workflows can inspect every candidate.
Stability classification is opt-in because a finite-difference Hessian is much
more expensive than the root solve itself.
"""
struct MagneticGapResult
    state::Union{Nothing, MeanFieldState{Float64}}
    selected_index::Int
    candidates::Vector{MagneticGapCandidate}
    converged::Bool
    attempt_count::Int
    failed_attempts::Int
    stability_classified::Bool
end

"""ForwardDiff stationarity residual for the magnetic five-field Omega.

The Landau cutoff is discrete, so callers must provide a fixed `n_max` for
the differentiable kernel. `solve_magnetic_gap` resolves that cutoff once per
seed before entering NLsolve and reuses it for both the primary and fallback
attempts.
"""
function gap_residual(model::PNJLMagneticModel, x, T_fm, mu_vec; kwargs...)
    return magnetic_gap_residual(model, x, T_fm, mu_vec; kwargs...)
end

@inline function _magnetic_solver_validate(model::PNJLMagneticModel, T_fm, mu_vec, xi)
    T_value = Float64(T_fm)
    isfinite(T_value) && T_value > 0.0 || throw(ArgumentError(
        "magnetic five-field gap solve requires finite T_fm > 0, got $(T_fm)",
    ))
    xi_value = Float64(xi)
    isfinite(xi_value) && abs(xi_value) <= 1e-14 || throw(ArgumentError(
        "PNJLMagneticModel currently supports only xi=0, got xi=$(xi)",
    ))
    μ = normalize_mu_vec(mu_vec)
    all(isfinite, μ) || throw(ArgumentError("magnetic mu_vec must be finite, got $(mu_vec)"))
    _magnetic_thermodynamics_module().validate_magnetic_eB(model.magnetic.eB_fm2)
    return μ
end

@inline function _magnetic_seed_vec(seed)
    if seed isa MeanFieldState
        return SVector{5, Float64}(Tuple(state_vector(seed)))
    elseif seed isa AbstractVector
        length(seed) >= 5 || throw(ArgumentError("magnetic seed must have length >= 5, got $(length(seed))"))
        return SVector{5, Float64}(Float64(seed[1]), Float64(seed[2]), Float64(seed[3]), Float64(seed[4]), Float64(seed[5]))
    elseif seed isa NamedTuple
        return _magnetic_seed_vec(MeanFieldState(seed))
    end
    throw(ArgumentError("unsupported magnetic seed type $(typeof(seed))"))
end

function _magnetic_seed_pool(T_fm, mu_vec; initial_guess=nothing, seed_candidates=nothing, include_default_seeds::Bool=true)
    pool = SVector{5, Float64}[]
    if initial_guess !== nothing
        push!(pool, _magnetic_seed_vec(initial_guess))
    end
    if seed_candidates !== nothing
        for seed in seed_candidates
            push!(pool, _magnetic_seed_vec(seed))
        end
    end

    # Reuse the project-wide PNJL seed catalog, but always append it even when
    # the caller supplied a primary seed: magnetic branch selection must not
    # collapse to one initial condition.
    if include_default_seeds && isdefined(@__MODULE__, :seed_catalog) && isdefined(@__MODULE__, :FixedMu)
        μ_seed = sum(Float64.(Tuple(mu_vec))) / 3.0
        for seed in seed_catalog(FixedMu(), [Float64(T_fm), μ_seed])
            push!(pool, _magnetic_seed_vec(seed))
        end
    elseif include_default_seeds
        push!(pool, SVector{5, Float64}(-1.84329, -1.84329, -2.22701, 1.0e-5, 4.0e-5))
        push!(pool, SVector{5, Float64}(-0.73192, -0.73192, -1.79539, 0.60532, 0.60532))
    end

    unique_pool = SVector{5, Float64}[]
    for seed in pool
        all(isfinite, seed) || continue
        any(existing -> norm(seed - existing) <= 1e-12, unique_pool) || push!(unique_pool, seed)
    end
    isempty(unique_pool) && throw(ArgumentError("magnetic solver received no finite seed"))
    return unique_pool
end

@inline function _magnetic_state_vec(x)
    if x isa MeanFieldState
        return state_vector(x)
    elseif x isa AbstractVector
        length(x) >= 5 || throw(ArgumentError("magnetic state must have length >= 5, got $(length(x))"))
        Tx = promote_type(typeof(x[1]), typeof(x[2]), typeof(x[3]), typeof(x[4]), typeof(x[5]))
        return SVector{5, Tx}(Tx(x[1]), Tx(x[2]), Tx(x[3]), Tx(x[4]), Tx(x[5]))
    end
    throw(ArgumentError("unsupported magnetic state type $(typeof(x))"))
end

"""Construct the magnetic five-field stationarity residual with ForwardDiff."""
function magnetic_gap_residual(
    model::PNJLMagneticModel,
    x,
    T_fm,
    mu_vec;
    xi::Real=0.0,
    p_num::Union{Nothing, Int}=nothing,
    t_num::Int=8,
    pz_max::Union{Nothing, Real}=nothing,
    n_max::Union{Nothing, Int}=nothing,
    cutoff_N::Union{Nothing, Int}=nothing,
)
    n_max_value = n_max === nothing ? model.magnetic.n_max : n_max
    n_max_value !== nothing || throw(ArgumentError(
        "magnetic AD residual requires a fixed n_max; solve_magnetic_gap resolves it per seed",
    ))
    n_max_value >= 0 || throw(ArgumentError("magnetic AD residual requires n_max >= 0, got $(n_max_value)"))
    μ = _magnetic_solver_validate(model, T_fm, mu_vec, xi)
    x0 = _magnetic_state_vec(x)
    thermo = _magnetic_thermodynamics_module()
    omega_fn = y -> thermo.calculate_magnetic_omega_components(
        y,
        μ,
        T_fm,
        model.magnetic;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        pz_max=pz_max,
        n_max=n_max_value,
        cutoff_N=cutoff_N,
    ).omega
    gradient = ForwardDiff.gradient(omega_fn, x0)
    return SVector{5, eltype(gradient)}(Tuple(gradient))
end

"""Compatibility alias for the historical explicit AD residual name."""
@inline magnetic_gap_residual_autodiff(args...; kwargs...) = magnetic_gap_residual(args...; kwargs...)

function _resolve_magnetic_attempt_nmax(
    model::PNJLMagneticModel,
    seed::SVector{5, Float64},
    T_fm,
    μ;
    xi::Real,
    p_num::Int,
    t_num::Int,
    pz_max,
    n_max,
    cutoff_N,
)
    explicit = n_max === nothing ? model.magnetic.n_max : n_max
    explicit !== nothing && return Int(explicit)

    thermo = _magnetic_thermodynamics_module()
    components = thermo.calculate_magnetic_omega_components(
        seed,
        μ,
        T_fm,
        model.magnetic;
        xi=xi,
        p_num=p_num,
        t_num=t_num,
        pz_max=pz_max,
        n_max=nothing,
        cutoff_N=cutoff_N,
    )
    resolved = Int(components.n_max)
    resolved >= 0 || throw(ArgumentError("resolved magnetic n_max must be >= 0, got $(resolved)"))
    return resolved
end

function _magnetic_attempt(
    model::PNJLMagneticModel,
    seed::SVector{5, Float64},
    T_fm,
    μ,
    method::Symbol;
    xtol::Real,
    ftol::Real,
    iterations::Int,
    residual_norm_max::Real,
    xi::Real,
    p_num::Int,
    t_num::Int,
    pz_max,
    n_max,
    cutoff_N,
)
    residual! = (F, x) -> begin
        r = magnetic_gap_residual(model, x, T_fm, μ;
            xi=xi, p_num=p_num, t_num=t_num, pz_max=pz_max, n_max=n_max,
            cutoff_N=cutoff_N)
        @inbounds for i in 1:5
            F[i] = r[i]
        end
        nothing
    end
    try
        result = nlsolve(
            residual!,
            collect(seed);
            autodiff=:forward,
            method=method,
            xtol=Float64(xtol),
            ftol=Float64(ftol),
            iterations=iterations,
        )
        x_state = SVector{5, Float64}(Tuple(Float64.(result.zero)))
        residual = magnetic_gap_residual(model, x_state, T_fm, μ;
            xi=xi, p_num=p_num, t_num=t_num, pz_max=pz_max, n_max=n_max,
            cutoff_N=cutoff_N)
        residual_norm = norm(residual)
        comp = _magnetic_thermodynamics_module().calculate_magnetic_omega_components(
            x_state, μ, T_fm, model.magnetic;
            xi=xi, p_num=p_num, t_num=t_num, pz_max=pz_max, n_max=n_max, cutoff_N=cutoff_N,
        )
        physical = all(m -> isfinite(m) && m > 0.0, comp.masses) &&
            all(v -> isfinite(v) && 0.0 <= v <= 1.0, x_state[4:5]) &&
            isfinite(comp.omega)
        converged = Bool(result.f_converged) && isfinite(residual_norm) &&
            residual_norm <= Float64(residual_norm_max) && physical
        iterations_used = hasproperty(result, :iterations) ? Int(result.iterations) : 0
        return (
            seed=seed,
            x_state=x_state,
            omega=Float64(comp.omega),
            residual=SVector{5, Float64}(residual),
            residual_norm=Float64(residual_norm),
            converged=converged,
            physical=physical,
            method=method,
            iterations=iterations_used,
            n_max=Int(comp.n_max),
        )
    catch err
        err isa InterruptException && rethrow()
        return (
            seed=seed,
            x_state=seed,
            omega=Inf,
            residual=SVector{5, Float64}(fill(NaN, 5)),
            residual_norm=Inf,
            converged=false,
            physical=false,
            method=method,
            iterations=0,
            n_max=0,
        )
    end
end

function _magnetic_hessian(
    model::PNJLMagneticModel,
    x_state,
    T_fm,
    μ;
    finite_difference_step::Real,
    xi::Real,
    p_num::Int,
    t_num::Int,
    pz_max,
    n_max,
    cutoff_N,
)
    x0 = _magnetic_seed_vec(x_state)
    H = Matrix{Float64}(undef, 5, 5)
    @inbounds for j in 1:5
        h = 0.5 * Float64(finite_difference_step) * max(1.0, abs(x0[j]))
        xp = MVector{5, Float64}(x0)
        xm = MVector{5, Float64}(x0)
        xp[j] += h
        xm[j] -= h
        rp = magnetic_gap_residual(model, xp, T_fm, μ;
            xi=xi, p_num=p_num, t_num=t_num, pz_max=pz_max, n_max=n_max,
            cutoff_N=cutoff_N)
        rm = magnetic_gap_residual(model, xm, T_fm, μ;
            xi=xi, p_num=p_num, t_num=t_num, pz_max=pz_max, n_max=n_max,
            cutoff_N=cutoff_N)
        H[:, j] = (rp - rm) / (2 * h)
    end
    return 0.5 .* (H + H')
end

function _magnetic_stability(
    model::PNJLMagneticModel,
    candidate,
    T_fm,
    μ;
    classify_stability::Bool,
    finite_difference_step::Real,
    xi::Real,
    p_num::Int,
    t_num::Int,
    pz_max,
    n_max,
    cutoff_N,
)
    classify_stability || return :not_evaluated
    try
        H = _magnetic_hessian(model, candidate.x_state, T_fm, μ;
            finite_difference_step=finite_difference_step, xi=xi, p_num=p_num,
            t_num=t_num, pz_max=pz_max, n_max=n_max, cutoff_N=cutoff_N)
        all(isfinite, H) || return :unknown
        λ = eigvals(Symmetric(H))
        scale = max(1.0, norm(H))
        tol = 1e-7 * scale
        all(real(v) > tol for v in λ) && return :local_minimum
        return :saddle_or_maximum
    catch
        return :unknown
    end
end

function solve_magnetic_gap(
    model::PNJLMagneticModel,
    T_fm,
    mu_vec;
    xi::Real=0.0,
    p_num::Union{Nothing, Int}=nothing,
    t_num::Int=8,
    pz_max::Union{Nothing, Real}=nothing,
    n_max::Union{Nothing, Int}=nothing,
    cutoff_N::Union{Nothing, Int}=nothing,
    initial_guess=nothing,
    seed_candidates=nothing,
    solver=nothing,
    method::Symbol=:trust_region,
    fallback_method::Union{Nothing, Symbol}=:newton,
    xtol::Real=1e-9,
    ftol::Real=1e-9,
    iterations::Int=600,
    residual_norm_max::Real=1e-6,
    finite_difference_step::Real=1e-5,
    root_merge_tol::Real=1e-5,
    classify_stability::Bool=false,
    include_default_seeds::Bool=true,
)
    μ = _magnetic_solver_validate(model, T_fm, mu_vec, xi)
    p_num_eff = p_num === nothing ? model.magnetic.p_num : Int(p_num)
    p_num_eff >= 4 || throw(ArgumentError("magnetic p_num must be >= 4, got $(p_num_eff)"))
    method_eff = solver === nothing ? method : Symbol(getproperty(solver, :method))
    xtol_eff = solver === nothing ? xtol : Float64(getproperty(solver, :xtol))
    ftol_eff = solver === nothing ? ftol : Float64(getproperty(solver, :ftol))

    seeds = _magnetic_seed_pool(T_fm, μ;
        initial_guess=initial_guess,
        seed_candidates=seed_candidates,
        include_default_seeds=include_default_seeds,
    )
    attempts = NamedTuple[]
    for seed in seeds
        n_max_seed = _resolve_magnetic_attempt_nmax(model, seed, T_fm, μ;
            xi=xi, p_num=p_num_eff, t_num=t_num, pz_max=pz_max, n_max=n_max, cutoff_N=cutoff_N)
        primary = _magnetic_attempt(model, seed, T_fm, μ, method_eff;
            xtol=xtol_eff, ftol=ftol_eff, iterations=iterations,
            residual_norm_max=residual_norm_max,
            xi=xi, p_num=p_num_eff, t_num=t_num, pz_max=pz_max, n_max=n_max_seed, cutoff_N=cutoff_N)
        push!(attempts, primary)
        if !primary.converged && fallback_method !== nothing && fallback_method != method_eff
            fallback = _magnetic_attempt(model, seed, T_fm, μ, fallback_method;
                xtol=xtol_eff, ftol=ftol_eff, iterations=iterations,
                residual_norm_max=residual_norm_max,
                xi=xi, p_num=p_num_eff, t_num=t_num, pz_max=pz_max, n_max=n_max_seed, cutoff_N=cutoff_N)
            push!(attempts, fallback)
        end
    end

    roots = NamedTuple[]
    for attempt in attempts
        attempt.converged || continue
        duplicate = findfirst(existing -> norm(attempt.x_state - existing.x_state) <= Float64(root_merge_tol), roots)
        if duplicate === nothing
            push!(roots, attempt)
        elseif attempt.residual_norm < roots[duplicate].residual_norm
            roots[duplicate] = attempt
        end
    end
    isempty(roots) && throw(ErrorException(
        "magnetic five-field solve did not converge: $(length(attempts)) attempts at T=$(T_fm), eB=$(model.magnetic.eB_fm2)",
    ))

    sort!(roots; by=attempt -> attempt.omega)
    candidates = MagneticGapCandidate[]
    stable_idx = 0
    for (idx, root) in enumerate(roots)
        stability = _magnetic_stability(model, root, T_fm, μ;
            classify_stability=classify_stability, finite_difference_step=finite_difference_step,
            xi=xi, p_num=p_num_eff, t_num=t_num, pz_max=pz_max, n_max=root.n_max, cutoff_N=cutoff_N)
        if stability == :local_minimum && stable_idx == 0
            stable_idx = idx
        end
        label = if !classify_stability
            idx == 1 ? :equilibrium : :candidate
        elseif stability == :local_minimum
            idx == stable_idx ? :stable : :metastable
        elseif stability == :saddle_or_maximum
            :unstable
        else
            :unclassified
        end
        seed_index = findfirst(seed -> norm(seed - root.seed) <= 1e-12, seeds)
        push!(candidates, MagneticGapCandidate(
            Int(seed_index === nothing ? idx : seed_index), root.seed, root.x_state,
            root.omega, root.residual, root.residual_norm, root.converged, root.physical,
            root.method, root.iterations, root.n_max, stability, label,
        ))
    end
    # Hessian classification is diagnostic metadata, not a default PNJL
    # admissibility filter. Keep the ordinary convenience selection (the
    # lowest-Omega converged candidate) even when every candidate is labelled
    # `saddle_or_maximum`; branch-aware callers can inspect all labels and
    # apply an explicit policy of their own.
    selected_idx = 1
    selected = candidates[selected_idx]
    selected_state = MeanFieldState(
        SVector{3, Float64}(selected.x_state[1], selected.x_state[2], selected.x_state[3]);
        Phi=selected.x_state[4], PhiBar=selected.x_state[5],
    )
    return MagneticGapResult(
        selected_state, selected_idx, candidates, true, length(attempts),
        count(attempt -> !attempt.converged, attempts), classify_stability,
    )
end

function solve_gap(model::PNJLMagneticModel, T_fm, mu_vec; kwargs...)
    result = solve_magnetic_gap(model, T_fm, mu_vec; kwargs...)
    result.state === nothing && throw(ErrorException("magnetic gap solve returned no equilibrium state"))
    return result.state
end
