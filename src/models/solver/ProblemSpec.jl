"""
    ProblemSpec

约束问题契约：统一描述 mode、维度以及求解/后处理回调。
"""
struct ProblemSpec{M,R,F,C,U,P,H,S}
    mode::M
    x_dim::Int
    theta_dim::Int
    residual!::R
    forward_solve::F
    conditions::C
    unpack_solution::U
    postprocess::P
    hard_rules::H
    selector::S
end

@inline function ProblemSpec(
    mode::ConstraintMode;
    x_dim::Int=state_dim(mode),
    theta_dim::Int=param_dim(mode),
    residual! = _unimplemented_residual,
    forward_solve=_unimplemented_forward_solve,
    conditions=_unimplemented_conditions,
    unpack_solution=identity,
    postprocess=identity,
    hard_rules=default_hard_constraint_rules(),
    selector=select_pressure_max_candidate,
)
    x_dim > 0 || throw(ArgumentError("x_dim must be positive, got $(x_dim)"))
    theta_dim > 0 || throw(ArgumentError("theta_dim must be positive, got $(theta_dim)"))
    return ProblemSpec(mode, x_dim, theta_dim, residual!, forward_solve, conditions, unpack_solution, postprocess, hard_rules, selector)
end

@inline _unimplemented_residual(F, x, theta, cfg, mode) = throw(ArgumentError("residual! is not configured for $(typeof(mode))"))
@inline _unimplemented_forward_solve(theta, cfg, mode) = throw(ArgumentError("forward_solve is not configured for $(typeof(mode))"))
@inline _unimplemented_conditions(theta, x, meta, cfg, mode) = throw(ArgumentError("conditions is not configured for $(typeof(mode))"))

function _fixedrho_problem_spec_forward_solve(model::AbstractQCDModel, mode::FixedRho, T_fm::Real; fwd_kwargs...)
    kwargs = Dict{Symbol,Any}(pairs(fwd_kwargs))

    seed_guess = get(kwargs, :seed_guess, nothing)
    seed_guess === nothing && throw(ArgumentError("seed_guess is required for ProblemSpec FixedRho forward_solve"))

    seed_pool = if haskey(kwargs, :seed_candidates)
        [Float64.(s) for s in kwargs[:seed_candidates]]
    else
        _build_default_seed_candidates(seed_guess)
    end

    physicality_check = get(kwargs, :physicality_check, ((_, _) -> true))
    hard_constraints = get(kwargs, :hard_constraints, default_hard_constraint_rules(; physicality_check=physicality_check))

    candidates = NamedTuple[]
    for (seed_index, seed) in pairs(seed_pool)
        local raw
        try
            local_kwargs = Dict{Symbol,Any}(kwargs)
            local_kwargs[:seed_guess] = seed
            delete!(local_kwargs, :seed_candidates)
            delete!(local_kwargs, :hard_constraints)

            solved = _solve_constraint_fixedrho(model, T_fm, mode.rho_target; pairs(local_kwargs)...)
            raw = (; solved..., residual_norm_max=get(local_kwargs, :residual_norm_max, 1e-6))
            ok, failed = evaluate_hard_constraints(raw, hard_constraints)
            push!(candidates, (; raw..., hard_constraint_ok=ok, failed_constraints=failed, converged=ok, seed_index=Int(seed_index)))
        catch
            raw = (
                converged=false,
                solution=Float64[],
                x_state=zeros(5),
                mu_vec=zeros(3),
                omega=NaN,
                pressure=-Inf,
                rho_norm=NaN,
                entropy=NaN,
                energy=NaN,
                masses=zeros(3),
                iterations=0,
                residual_norm=Inf,
                residual_norm_max=Float64(get(kwargs, :residual_norm_max, 1e-6)),
            )
            push!(candidates, (; raw..., hard_constraint_ok=false, failed_constraints=Symbol[:solver_failed], seed_index=Int(seed_index)))
        end
    end

    selected = select_pressure_max_candidate(candidates)
    s = selected.selected_candidate
    return (
        converged=Bool(s.converged),
        solution=Vector{Float64}(s.solution),
        x_state=s.x_state,
        mu_vec=s.mu_vec,
        omega=s.omega,
        pressure=s.pressure,
        rho_norm=s.rho_norm,
        entropy=s.entropy,
        energy=s.energy,
        masses=s.masses,
        iterations=s.iterations,
        residual_norm=s.residual_norm,
        hard_constraint_ok=s.hard_constraint_ok,
        failed_constraints=s.failed_constraints,
        selection_reason=selected.selection_reason,
        selected_index=selected.selected_index,
        candidate_count=length(candidates),
    )
end

@inline function build_problem_spec(mode::ConstraintMode; kwargs...)
    if isempty(kwargs)
        components = build_constraint_components(mode)
        forward_solve = if mode isa FixedRho
            (model, T_fm; fwd_kwargs...) -> _fixedrho_problem_spec_forward_solve(model, mode, T_fm; fwd_kwargs...)
        else
            (model, T_fm; fwd_kwargs...) -> solve_constraint(model, mode, T_fm; fwd_kwargs...)
        end
        return ProblemSpec(
            mode;
            x_dim=constraint_total_dim(components),
            conditions=params -> build_conditions(mode, params),
            forward_solve=forward_solve,
        )
    end
    return ProblemSpec(mode; kwargs...)
end
