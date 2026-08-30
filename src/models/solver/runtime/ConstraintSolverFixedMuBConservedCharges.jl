"""
Joint 8-dimensional executor for `FixedMuBConservedCharges`.

The unknown vector remains `[mean_fields(5), mu_u, mu_d, mu_s]`. This keeps the
new mode on the same unified non-fixedmu solve semantics as the existing
constrained modes; `mu_Q` and `mu_S` are derived, not nested outer-loop state.
"""
function _solve_constraint_fixedmub_conserved(
    model::AbstractQCDModel,
    T_fm::Real,
    mode;
    seed_guess::AbstractVector,
    mu0::AbstractVector=let flavor = flavor_mu_from_bqs(mode.muB_fm, 0.0, mode.muB_fm / 3)
        Float64[flavor.mu_u, flavor.mu_d, flavor.mu_s]
    end,
    solver_primary::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:forward, xtol=1e-9, ftol=1e-9),
    solver_secondary::Union{Nothing,AbstractGapSolver}=nothing,
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
    _ = solver_primary, solver_secondary, enforce_physicality
    isfinite(rho0) && rho0 > 0 || throw(ArgumentError("rho0 must be finite and positive, got $(rho0)"))
    nlsolve_method in (:newton, :trust_region) || throw(ArgumentError(
        "nlsolve_method must be :newton or :trust_region, got $(nlsolve_method)",
    ))

    model_kind = _model_kind_for_shared_core(model)
    joint_params = GapParams(
        Float64(T_fm),
        cached_nodes(p_num, t_num),
        Float64(xi);
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
        push!(candidate_seeds, length(seed) >= 8 ? Float64.(seed[1:8]) : extend_seed(Float64.(seed), mode))
    end

    unique_seeds = Dict{String,Vector{Float64}}()
    for seed in candidate_seeds
        key = join(round.(seed; digits=10), ",")
        haskey(unique_seeds, key) || (unique_seeds[key] = seed)
    end

    nls_kwargs = Dict{Symbol,Any}(pairs(nlsolve_kwargs))
    joint_iterations = Int(get(nls_kwargs, :iterations, 600))
    best_joint = nothing

    for x0 in values(unique_seeds)
        joint = try
            nlsolve(
                joint_residual!,
                x0;
                autodiff=:forward,
                method=nlsolve_method,
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
        conserved_mu = conserved_mu_from_flavor(joint_thermo.mu_vec...)
        conserved_rho = conserved_densities_from_flavor(joint_thermo.rho_vec)
        residual_norm = _compose_mode_residual_norm(
            model,
            joint_thermo.x_state,
            joint_thermo.mu_vec,
            T_fm,
            (conserved_mu.mu_B, mode.muB_fm),
            (
                conserved_rho.rho_Q / rho0,
                mode.charge_to_baryon_ratio * conserved_rho.rho_B / rho0,
            ),
            (
                conserved_rho.rho_S / rho0,
                mode.strangeness_density_target / rho0,
            );
            xi=xi,
            p_num=p_num,
            t_num=t_num,
        )

        thermo_finite = (
            isfinite(joint_thermo.omega) && isfinite(joint_thermo.pressure) &&
            isfinite(joint_thermo.rho_norm) && isfinite(joint_thermo.entropy) &&
            isfinite(joint_thermo.energy)
        )
        physical = (
            physicality_check(joint_thermo.x_state, joint_thermo.masses) &&
            thermo_finite && _is_mass_positive(joint_thermo.masses) &&
            _is_phi_in_range(joint_thermo.x_state)
        )
        converged = (
            Bool(joint.f_converged) && physical && isfinite(residual_norm) &&
            residual_norm <= Float64(residual_norm_max)
        )

        candidate = _build_mode_result_from_outer_state(
            joint_thermo.x_state,
            joint_thermo.mu_vec,
            joint_thermo.pressure,
            joint_thermo.rho_norm,
            joint_thermo.entropy,
            joint_thermo.energy,
            joint_thermo.masses,
            residual_norm;
            iterations=joint.iterations,
            converged=converged,
            legacy_fallback_used=false,
        )

        if best_joint === nothing ||
           (candidate.converged != best_joint.converged && candidate.converged) ||
           (candidate.converged == best_joint.converged && candidate.residual_norm < best_joint.residual_norm)
            best_joint = candidate
        end
    end

    best_joint !== nothing && return best_joint
    return _build_mode_failure_candidate(residual_norm_max=Float64(residual_norm_max))
end
