function _solve_constraint_fixedasymrho(
    model::AbstractQCDModel,
    T_fm::Real,
    rho_target::Real,
    ud_ratio_target::Real,
    s_target::Real;
    seed_guess::AbstractVector,
    mu0::AbstractVector=Float64[0.2, 0.2, 0.2],
    solver_primary::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:forward, xtol=1e-9, ftol=1e-9),
    solver_secondary::Union{Nothing, AbstractGapSolver}=nothing,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    nlsolve_method::Symbol=:trust_region,
    rho0::Real,
    enforce_physicality::Bool=false,
    physicality_check::Function=((_, _) -> true),
    nlsolve_kwargs...,
)
    _ = solver_primary, solver_secondary, nlsolve_method, enforce_physicality

    mode = FixedAsymmetricRho(Float64(rho_target), Float64(ud_ratio_target), Float64(s_target))
    model_kind = _model_kind_for_shared_core(model)
    joint_params = GapParams(Float64(T_fm), cached_nodes(p_num, t_num), Float64(xi);
        p_num=p_num,
        t_num=t_num,
        model_kind=model_kind,
    )
    joint_residual! = build_residual!(mode, joint_params)

    candidate_seeds = Vector{Vector{Float64}}()
    if length(seed_guess) >= 8
        push!(candidate_seeds, Float64.(seed_guess[1:8]))
    else
        push!(candidate_seeds, extend_seed(Float64.(seed_guess), mode))
    end
    if length(seed_guess) >= 5 && length(mu0) >= 3
        push!(candidate_seeds, Float64[Float64.(seed_guess[1:5])..., Float64.(mu0[1:3])...])
    end
    for seed in seed_catalog(mode, [T_fm])
        if length(seed) >= 8
            push!(candidate_seeds, Float64.(seed[1:8]))
        else
            push!(candidate_seeds, extend_seed(Float64.(seed), mode))
        end
    end

    uniq = Dict{String, Vector{Float64}}()
    for s in candidate_seeds
        key = join(round.(s; digits=10), ",")
        haskey(uniq, key) || (uniq[key] = s)
    end

    nls_kwargs = Dict{Symbol,Any}(pairs(nlsolve_kwargs))
    joint_iterations = Int(get(nls_kwargs, :iterations, 600))
    best_joint = nothing

    for x0 in values(uniq)
        joint = try
            nlsolve(
                joint_residual!,
                x0;
                autodiff=:forward,
                method=:trust_region,
                xtol=1e-9,
                ftol=1e-9,
                iterations=max(joint_iterations, 300),
            )
        catch err
            err isa InterruptException && rethrow()
            nothing
        end
        joint === nothing && continue

        joint_solution = Float64.(joint.zero)
        joint_thermo = compute_thermo_from_solution(
            model,
            joint_solution,
            T_fm;
            xi=xi,
            p_num=p_num,
            t_num=t_num,
            rho0_scale=rho0,
            state_n=5,
            mu_n=3,
        )

        rho_u, rho_d, rho_s = joint_thermo.rho_vec[1], joint_thermo.rho_vec[2], joint_thermo.rho_vec[3]
        nB_joint = sum(joint_thermo.rho_vec) / (3.0 * rho0)
        ud_ratio_joint = rho_u - rho_d * ud_ratio_target

        joint_residual_norm = _compose_mode_residual_norm(
            model,
            joint_thermo.x_state,
            joint_thermo.mu_vec,
            T_fm,
            (nB_joint, rho_target),
            (ud_ratio_joint, 0.0),
            (rho_s, s_target);
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )

        joint_thermo_finite = isfinite(joint_thermo.omega) && isfinite(joint_thermo.pressure) &&
            isfinite(joint_thermo.rho_norm) && isfinite(joint_thermo.entropy) && isfinite(joint_thermo.energy)
        joint_phys = physicality_check(joint_thermo.x_state, joint_thermo.masses) &&
            joint_thermo_finite && _is_mass_positive(joint_thermo.masses) && _is_phi_in_range(joint_thermo.x_state)
        joint_converged = Bool(joint.f_converged) && joint_phys && isfinite(joint_residual_norm) &&
            joint_residual_norm <= max(Float64(residual_norm_max), 1e-3)

        joint_result = _build_mode_result_from_outer_state(
            joint_thermo.x_state,
            joint_thermo.mu_vec,
            joint_thermo.pressure,
            joint_thermo.rho_norm,
            joint_thermo.entropy,
            joint_thermo.energy,
            joint_thermo.masses,
            joint_residual_norm;
            iterations=joint.iterations,
            converged=joint_converged,
            legacy_fallback_used=false,
        )

        if best_joint === nothing ||
           (joint_result.converged != best_joint.converged && joint_result.converged) ||
           (joint_result.converged == best_joint.converged && joint_result.residual_norm < best_joint.residual_norm)
            best_joint = joint_result
        end
    end

    if best_joint !== nothing
        return best_joint
    end

    return _build_mode_failure_candidate(
        residual_norm_max=Float64(residual_norm_max),
    )
end
