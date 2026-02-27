"""constraint_solver.jl

Models 侧约束求解（阶段 B 内核下沉）：
- 先以 FixedEntropy 为起点，将单变量 μ 外层约束求解从 PNJL 兼容层下沉到 models。
"""

using StaticArrays
using NLsolve
using ForwardDiff
using Statistics: mean

@inline function default_mu0_from_seed(seed)::Float64
    if length(seed) >= 8
        return mean(Float64.(seed[6:8]))
    elseif length(seed) >= 6
        return Float64(seed[6])
    end
    return 0.2
end

@inline function default_muvec0_from_seed(seed)::Vector{Float64}
    if length(seed) >= 8
        return Float64.(seed[6:8])
    end
    μ0 = default_mu0_from_seed(seed)
    return [μ0, μ0, μ0]
end

function solve_fixedmu_constraint(
    model::AbstractQCDModel,
    T_fm::Real,
    μ_fm::Real;
    seed_guess::AbstractVector,
    solver::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:finite, xtol=1e-9, ftol=1e-9),
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    physicality_check::Function=((_, _) -> true),
)
    st = solve_gap(
        model,
        T_fm,
        μ_fm;
        solver_backend=:models,
        solver=solver,
        initial_guess=Float64.(seed_guess),
        residual_norm_max=residual_norm_max,
        xi=xi,
        p_num=p_num,
        t_num=t_num,
    )

    x_state = SVector{5}(Tuple(state_vector(st)))
    mu_vec = normalize_mu_vec(μ_fm)

    pressure = -omega(model, x_state, T_fm, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
    pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
    rho_vec = ForwardDiff.gradient(pressure_mu, mu_vec)
    rho_norm = sum(rho_vec) / 3.0
    pressure_T = τ -> -omega(model, x_state, τ, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
    entropy = ForwardDiff.derivative(pressure_T, T_fm)
    energy = -pressure + sum(mu_vec .* rho_vec) + T_fm * entropy
    omega_val = -pressure
    masses = calculate_mass_vec(model, SVector{3}(x_state[1], x_state[2], x_state[3]))

    residual_vec = gap_residual(model, st, T_fm, mu_vec; xi=xi, p_num=p_num, t_num=t_num)
    residual_norm = sqrt(sum(abs2, residual_vec))

    thermo_finite = isfinite(omega_val) && isfinite(pressure) && isfinite(rho_norm) && isfinite(entropy) && isfinite(energy)
    converged = physicality_check(x_state, masses) && thermo_finite && isfinite(residual_norm) && residual_norm <= Float64(residual_norm_max)

    return (
        converged=converged,
        solution=Vector{Float64}(x_state),
        x_state=x_state,
        mu_vec=SVector{3}(Tuple(mu_vec)),
        omega=omega_val,
        pressure=pressure,
        rho_norm=rho_norm,
        entropy=entropy,
        energy=energy,
        masses=masses,
        iterations=0,
        residual_norm=residual_norm,
    )
end

function solve_fixedrho_constraint(
    model::AbstractQCDModel,
    T_fm::Real,
    rho_target::Real;
    seed_guess::AbstractVector,
    mu0::Real=0.2,
    solver_primary::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:finite, xtol=1e-9, ftol=1e-9),
    solver_secondary::Union{Nothing, AbstractGapSolver}=nothing,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    nlsolve_method::Symbol=:newton,
    physicality_check::Function=((_, _) -> true),
    nlsolve_kwargs...,
)
    st_ref = Ref{Any}(nothing)
    x_state_ref = Ref(SVector{5, Float64}(fill(NaN, 5)))
    mu_vec_ref = Ref(SVector{3, Float64}(fill(NaN, 3)))
    pressure_ref = Ref(NaN)
    rho_norm_ref = Ref(NaN)
    entropy_ref = Ref(NaN)
    energy_ref = Ref(NaN)
    masses_ref = Ref(SVector{3, Float64}(fill(NaN, 3)))

    residual_fn! = (F, x) -> begin
        μ = Float64(x[1])
        μ_vec = normalize_mu_vec(μ)

        initial_guess = if st_ref[] === nothing
            Float64.(seed_guess)
        else
            Vector{Float64}(state_vector(st_ref[]))
        end

        st = try
            solve_gap(
                model,
                T_fm,
                μ;
                solver_backend=:models,
                solver=solver_primary,
                initial_guess=initial_guess,
                residual_norm_max=residual_norm_max,
                xi=xi,
                p_num=p_num,
                t_num=t_num,
            )
        catch
            if solver_secondary === nothing
                nothing
            else
                try
                    solve_gap(
                        model,
                        T_fm,
                        μ;
                        solver_backend=:models,
                        solver=solver_secondary,
                        initial_guess=initial_guess,
                        residual_norm_max=residual_norm_max,
                        xi=xi,
                        p_num=p_num,
                        t_num=t_num,
                    )
                catch
                    nothing
                end
            end
        end

        if st === nothing
            F[1] = 1e6
            return nothing
        end

        x_state = SVector{5}(Tuple(state_vector(st)))
        pressure = -omega(model, x_state, T_fm, μ_vec; p_num=p_num, t_num=t_num, xi=xi)

        pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
        rho_vec = ForwardDiff.gradient(pressure_mu, μ_vec)
        rho_norm = sum(rho_vec) / 3.0

        pressure_T = τ -> -omega(model, x_state, τ, μ_vec; p_num=p_num, t_num=t_num, xi=xi)
        entropy = ForwardDiff.derivative(pressure_T, T_fm)

        energy = -pressure + sum(μ_vec .* rho_vec) + T_fm * entropy
        masses = calculate_mass_vec(model, SVector{3}(x_state[1], x_state[2], x_state[3]))

        st_ref[] = st
        x_state_ref[] = x_state
        mu_vec_ref[] = SVector{3}(Tuple(μ_vec))
        pressure_ref[] = pressure
        rho_norm_ref[] = rho_norm
        entropy_ref[] = entropy
        energy_ref[] = energy
        masses_ref[] = masses

        F[1] = rho_norm - rho_target
        return nothing
    end

    res = try
        nlsolve(
            residual_fn!,
            [Float64(mu0)];
            autodiff=:finite,
            method=nlsolve_method,
            xtol=1e-9,
            ftol=1e-9,
            nlsolve_kwargs...,
        )
    catch
        if nlsolve_method === :trust_region
            rethrow()
        end
        nlsolve(
            residual_fn!,
            [Float64(mu0)];
            autodiff=:finite,
            method=:trust_region,
            xtol=1e-9,
            ftol=1e-9,
            nlsolve_kwargs...,
        )
    end

    if st_ref[] === nothing
        tmpF = zeros(1)
        residual_fn!(tmpF, res.zero)
    end

    gap_vec = gap_residual(model, st_ref[], T_fm, mu_vec_ref[]; xi=xi, p_num=p_num, t_num=t_num)
    gap_norm = sqrt(sum(abs2, gap_vec))
    rho_residual = abs(rho_norm_ref[] - rho_target)
    residual_norm = max(gap_norm, rho_residual)

    omega_val = -pressure_ref[]
    thermo_finite = isfinite(omega_val) && isfinite(pressure_ref[]) && isfinite(rho_norm_ref[]) && isfinite(entropy_ref[]) && isfinite(energy_ref[])
    phys = physicality_check(x_state_ref[], masses_ref[]) && thermo_finite
    converged = res.f_converged && phys && isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-4)

    return (
        converged=converged,
        solution=Float64[
            x_state_ref[][1], x_state_ref[][2], x_state_ref[][3], x_state_ref[][4], x_state_ref[][5],
            mu_vec_ref[][1], mu_vec_ref[][2], mu_vec_ref[][3],
        ],
        x_state=x_state_ref[],
        mu_vec=mu_vec_ref[],
        omega=omega_val,
        pressure=pressure_ref[],
        rho_norm=rho_norm_ref[],
        entropy=entropy_ref[],
        energy=energy_ref[],
        masses=masses_ref[],
        iterations=res.iterations,
        residual_norm=residual_norm,
    )
end

function solve_fixedentropy_constraint(
    model::AbstractQCDModel,
    T_fm::Real,
    s_target::Real;
    seed_guess::AbstractVector,
    mu0::Real=0.2,
    solver_primary::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:finite, xtol=1e-9, ftol=1e-9),
    solver_secondary::Union{Nothing, AbstractGapSolver}=nothing,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    rho0::Real,
    physicality_check::Function=((_, _) -> true),
    mass_positive_constraint::Bool=true,
    nlsolve_kwargs...,
)
    st_ref = Ref{Any}(nothing)
    x_state_ref = Ref(SVector{5, Float64}(fill(NaN, 5)))
    mu_vec_ref = Ref(SVector{3, Float64}(fill(NaN, 3)))
    pressure_ref = Ref(NaN)
    rho_norm_ref = Ref(NaN)
    entropy_ref = Ref(NaN)
    energy_ref = Ref(NaN)
    masses_ref = Ref(SVector{3, Float64}(fill(NaN, 3)))

    residual_fn! = (F, x) -> begin
        μ = Float64(x[1])
        μ_vec = normalize_mu_vec(μ)

        initial_guess = if st_ref[] === nothing
            Float64.(seed_guess)
        else
            Vector{Float64}(state_vector(st_ref[]))
        end

        st = try
            solve_gap(
                model,
                T_fm,
                μ;
                solver_backend=:models,
                solver=solver_primary,
                initial_guess=initial_guess,
                residual_norm_max=residual_norm_max,
                xi=xi,
                p_num=p_num,
                t_num=t_num,
            )
        catch
            if solver_secondary === nothing
                nothing
            else
                try
                    solve_gap(
                        model,
                        T_fm,
                        μ;
                        solver_backend=:models,
                        solver=solver_secondary,
                        initial_guess=initial_guess,
                        residual_norm_max=residual_norm_max,
                        xi=xi,
                        p_num=p_num,
                        t_num=t_num,
                    )
                catch
                    nothing
                end
            end
        end

        if st === nothing
            F[1] = 1e6
            return nothing
        end

        x_state = SVector{5}(Tuple(state_vector(st)))
        pressure = -omega(model, x_state, T_fm, μ_vec; p_num=p_num, t_num=t_num, xi=xi)

        pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
        rho_vec = ForwardDiff.gradient(pressure_mu, μ_vec)
        rho_norm = sum(rho_vec) / (3.0 * rho0)

        pressure_T = τ -> -omega(model, x_state, τ, μ_vec; p_num=p_num, t_num=t_num, xi=xi)
        entropy = ForwardDiff.derivative(pressure_T, T_fm)

        energy = -pressure + sum(μ_vec .* rho_vec) + T_fm * entropy
        masses = calculate_mass_vec(model, SVector{3}(x_state[1], x_state[2], x_state[3]))

        if mass_positive_constraint && any(m -> !isfinite(m) || m <= 0.0, masses)
            F[1] = 1e6
            return nothing
        end

        st_ref[] = st
        x_state_ref[] = x_state
        mu_vec_ref[] = SVector{3}(Tuple(μ_vec))
        pressure_ref[] = pressure
        rho_norm_ref[] = rho_norm
        entropy_ref[] = entropy
        energy_ref[] = energy
        masses_ref[] = masses

        F[1] = entropy - s_target
        return nothing
    end

    res = nlsolve(
        residual_fn!,
        [Float64(mu0)];
        autodiff=:finite,
        method=:trust_region,
        xtol=1e-9,
        ftol=1e-9,
        nlsolve_kwargs...,
    )

    if st_ref[] === nothing
        tmpF = zeros(1)
        residual_fn!(tmpF, res.zero)
    end

    gap_vec = gap_residual(model, st_ref[], T_fm, mu_vec_ref[]; xi=xi, p_num=p_num, t_num=t_num)
    gap_norm = sqrt(sum(abs2, gap_vec))
    entropy_residual = abs(entropy_ref[] - s_target)
    residual_norm = max(gap_norm, entropy_residual)

    omega_val = -pressure_ref[]
    thermo_finite = isfinite(omega_val) && isfinite(pressure_ref[]) && isfinite(rho_norm_ref[]) && isfinite(entropy_ref[]) && isfinite(energy_ref[])
    phys = physicality_check(x_state_ref[], masses_ref[]) && thermo_finite

    Φ = x_state_ref[][4]
    Φbar = x_state_ref[][5]
    soft_phys = thermo_finite &&
        isfinite(Φ) && isfinite(Φbar) &&
        (-1e-3 <= Φ <= 1 + 1e-3) && (-1e-3 <= Φbar <= 1 + 1e-3) &&
        all(isfinite, masses_ref[]) && all(m -> m >= -1e-8, masses_ref[])

    accept_tol = max(Float64(residual_norm_max), 1e-3)
    near_accept_tol = min(1e-4, 10 * accept_tol)
    converged = isfinite(residual_norm) && (
        (res.f_converged && phys && residual_norm <= accept_tol) ||
        ((phys || soft_phys) && residual_norm <= near_accept_tol)
    )

    return (
        converged=converged,
        solution=Float64[
            x_state_ref[][1], x_state_ref[][2], x_state_ref[][3], x_state_ref[][4], x_state_ref[][5],
            mu_vec_ref[][1], mu_vec_ref[][2], mu_vec_ref[][3],
        ],
        x_state=x_state_ref[],
        mu_vec=mu_vec_ref[],
        omega=omega_val,
        pressure=pressure_ref[],
        rho_norm=rho_norm_ref[],
        entropy=entropy_ref[],
        energy=energy_ref[],
        masses=masses_ref[],
        iterations=res.iterations,
        residual_norm=residual_norm,
    )
end

function solve_fixedsigma_constraint(
    model::AbstractQCDModel,
    T_fm::Real,
    sigma_target::Real;
    seed_guess::AbstractVector,
    mu0::Real=0.2,
    solver_primary::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:finite, xtol=1e-9, ftol=1e-9),
    solver_secondary::Union{Nothing, AbstractGapSolver}=nothing,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    rho0::Real,
    physicality_check::Function=((_, _) -> true),
    nlsolve_kwargs...,
)
    st_ref = Ref{Any}(nothing)
    x_state_ref = Ref(SVector{5, Float64}(fill(NaN, 5)))
    mu_vec_ref = Ref(SVector{3, Float64}(fill(NaN, 3)))
    pressure_ref = Ref(NaN)
    rho_norm_ref = Ref(NaN)
    entropy_ref = Ref(NaN)
    energy_ref = Ref(NaN)
    masses_ref = Ref(SVector{3, Float64}(fill(NaN, 3)))

    residual_fn! = (F, x) -> begin
        μ = Float64(x[1])
        μ_vec = normalize_mu_vec(μ)

        initial_guess = if st_ref[] === nothing
            Float64.(seed_guess)
        else
            Vector{Float64}(state_vector(st_ref[]))
        end

        st = try
            solve_gap(
                model,
                T_fm,
                μ;
                solver_backend=:models,
                solver=solver_primary,
                initial_guess=initial_guess,
                residual_norm_max=residual_norm_max,
                xi=xi,
                p_num=p_num,
                t_num=t_num,
            )
        catch
            if solver_secondary === nothing
                nothing
            else
                try
                    solve_gap(
                        model,
                        T_fm,
                        μ;
                        solver_backend=:models,
                        solver=solver_secondary,
                        initial_guess=initial_guess,
                        residual_norm_max=residual_norm_max,
                        xi=xi,
                        p_num=p_num,
                        t_num=t_num,
                    )
                catch
                    nothing
                end
            end
        end

        if st === nothing
            F[1] = 1e6
            return nothing
        end

        x_state = SVector{5}(Tuple(state_vector(st)))
        pressure = -omega(model, x_state, T_fm, μ_vec; p_num=p_num, t_num=t_num, xi=xi)

        pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
        rho_vec = ForwardDiff.gradient(pressure_mu, μ_vec)
        rho_norm = sum(rho_vec) / (3.0 * rho0)

        pressure_T = τ -> -omega(model, x_state, τ, μ_vec; p_num=p_num, t_num=t_num, xi=xi)
        entropy = ForwardDiff.derivative(pressure_T, T_fm)

        energy = -pressure + sum(μ_vec .* rho_vec) + T_fm * entropy
        masses = calculate_mass_vec(model, SVector{3}(x_state[1], x_state[2], x_state[3]))

        st_ref[] = st
        x_state_ref[] = x_state
        mu_vec_ref[] = SVector{3}(Tuple(μ_vec))
        pressure_ref[] = pressure
        rho_norm_ref[] = rho_norm
        entropy_ref[] = entropy
        energy_ref[] = energy
        masses_ref[] = masses

        n_B = rho_norm * rho0
        if !isfinite(n_B) || abs(n_B) <= 1e-12
            F[1] = 1e6
            return nothing
        end

        F[1] = entropy / n_B - sigma_target
        return nothing
    end

    res = nlsolve(
        residual_fn!,
        [Float64(mu0)];
        autodiff=:finite,
        method=:trust_region,
        xtol=1e-9,
        ftol=1e-9,
        nlsolve_kwargs...,
    )

    if st_ref[] === nothing
        tmpF = zeros(1)
        residual_fn!(tmpF, res.zero)
    end

    gap_vec = gap_residual(model, st_ref[], T_fm, mu_vec_ref[]; xi=xi, p_num=p_num, t_num=t_num)
    gap_norm = sqrt(sum(abs2, gap_vec))
    n_B_ref = rho_norm_ref[] * rho0
    sigma_residual = if isfinite(n_B_ref) && abs(n_B_ref) > 1e-12
        abs(entropy_ref[] / n_B_ref - sigma_target)
    else
        Inf
    end
    residual_norm = max(gap_norm, sigma_residual)

    omega_val = -pressure_ref[]
    thermo_finite = isfinite(omega_val) && isfinite(pressure_ref[]) && isfinite(rho_norm_ref[]) && isfinite(entropy_ref[]) && isfinite(energy_ref[])
    phys = physicality_check(x_state_ref[], masses_ref[]) && thermo_finite
    converged = res.f_converged && phys && isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-3)

    return (
        converged=converged,
        solution=Float64[
            x_state_ref[][1], x_state_ref[][2], x_state_ref[][3], x_state_ref[][4], x_state_ref[][5],
            mu_vec_ref[][1], mu_vec_ref[][2], mu_vec_ref[][3],
        ],
        x_state=x_state_ref[],
        mu_vec=mu_vec_ref[],
        omega=omega_val,
        pressure=pressure_ref[],
        rho_norm=rho_norm_ref[],
        entropy=entropy_ref[],
        energy=energy_ref[],
        masses=masses_ref[],
        iterations=res.iterations,
        residual_norm=residual_norm,
    )
end

function solve_fixedasymrho_constraint(
    model::AbstractQCDModel,
    T_fm::Real,
    rho_target::Real,
    ud_ratio_target::Real,
    s_target::Real;
    seed_guess::AbstractVector,
    mu0::AbstractVector=Float64[0.2, 0.2, 0.2],
    solver_primary::AbstractGapSolver=NLsolveGapSolver(method=:newton, jacobian=:finite, xtol=1e-9, ftol=1e-9),
    solver_secondary::Union{Nothing, AbstractGapSolver}=nothing,
    xi::Real=0.0,
    p_num::Int=24,
    t_num::Int=8,
    residual_norm_max::Real=1e-6,
    rho0::Real,
    enforce_physicality::Bool=false,
    physicality_check::Function=((_, _) -> true),
    nlsolve_kwargs...,
)
    st_ref = Ref{Any}(nothing)
    x_state_ref = Ref(SVector{5, Float64}(fill(NaN, 5)))
    mu_vec_ref = Ref(SVector{3, Float64}(fill(NaN, 3)))
    pressure_ref = Ref(NaN)
    rho_norm_ref = Ref(NaN)
    entropy_ref = Ref(NaN)
    energy_ref = Ref(NaN)
    masses_ref = Ref(SVector{3, Float64}(fill(NaN, 3)))
    rho_vec_ref = Ref(SVector{3, Float64}(fill(NaN, 3)))

    residual_fn! = (F, x) -> begin
        μ_vec = SVector{3}(Float64(x[1]), Float64(x[2]), Float64(x[3]))

        initial_guess = if st_ref[] === nothing
            Float64.(seed_guess)
        else
            Vector{Float64}(state_vector(st_ref[]))
        end

        st = try
            solve_gap(
                model,
                T_fm,
                μ_vec;
                solver_backend=:models,
                solver=solver_primary,
                initial_guess=initial_guess,
                residual_norm_max=residual_norm_max,
                xi=xi,
                p_num=p_num,
                t_num=t_num,
            )
        catch
            if solver_secondary === nothing
                nothing
            else
                try
                    solve_gap(
                        model,
                        T_fm,
                        μ_vec;
                        solver_backend=:models,
                        solver=solver_secondary,
                        initial_guess=initial_guess,
                        residual_norm_max=residual_norm_max,
                        xi=xi,
                        p_num=p_num,
                        t_num=t_num,
                    )
                catch
                    nothing
                end
            end
        end

        if st === nothing
            F[1] = 1e6
            F[2] = 1e6
            F[3] = 1e6
            return nothing
        end

        x_state = SVector{5}(Tuple(state_vector(st)))
        pressure = -omega(model, x_state, T_fm, μ_vec; p_num=p_num, t_num=t_num, xi=xi)

        pressure_mu = μtrial -> -omega(model, x_state, T_fm, μtrial; p_num=p_num, t_num=t_num, xi=xi)
        rho_vec = ForwardDiff.gradient(pressure_mu, μ_vec)
        rho_norm = sum(rho_vec) / (3.0 * rho0)

        pressure_T = τ -> -omega(model, x_state, τ, μ_vec; p_num=p_num, t_num=t_num, xi=xi)
        entropy = ForwardDiff.derivative(pressure_T, T_fm)

        energy = -pressure + sum(μ_vec .* rho_vec) + T_fm * entropy
        masses = calculate_mass_vec(model, SVector{3}(x_state[1], x_state[2], x_state[3]))

        rho_u, rho_d, rho_s = rho_vec[1], rho_vec[2], rho_vec[3]
        nB = sum(rho_vec) / (3.0 * rho0)
        ud_ratio = if abs(rho_d) > 1e-12
            rho_u / rho_d
        else
            rho_u / (rho_d >= 0 ? 1e-12 : -1e-12)
        end

        st_ref[] = st
        x_state_ref[] = x_state
        mu_vec_ref[] = μ_vec
        pressure_ref[] = pressure
        rho_norm_ref[] = rho_norm
        entropy_ref[] = entropy
        energy_ref[] = energy
        masses_ref[] = masses
        rho_vec_ref[] = SVector{3}(Tuple(rho_vec))

        F[1] = nB - rho_target
        F[2] = ud_ratio - ud_ratio_target
        F[3] = rho_s - s_target
        return nothing
    end

    res = nlsolve(
        residual_fn!,
        Float64.(mu0);
        autodiff=:finite,
        method=:trust_region,
        xtol=1e-9,
        ftol=1e-9,
        nlsolve_kwargs...,
    )

    if st_ref[] === nothing
        tmpF = zeros(3)
        residual_fn!(tmpF, res.zero)
    end

    gap_vec = gap_residual(model, st_ref[], T_fm, mu_vec_ref[]; xi=xi, p_num=p_num, t_num=t_num)
    gap_norm = sqrt(sum(abs2, gap_vec))
    rho_u, rho_d, rho_s = rho_vec_ref[][1], rho_vec_ref[][2], rho_vec_ref[][3]
    nB = sum(rho_vec_ref[]) / (3.0 * rho0)
    ud_ratio = if abs(rho_d) > 1e-12
        rho_u / rho_d
    else
        rho_u / (rho_d >= 0 ? 1e-12 : -1e-12)
    end
    nB_residual = abs(nB - rho_target)
    ud_residual = abs(ud_ratio - ud_ratio_target)
    s_residual = abs(rho_s - s_target)
    residual_norm = max(gap_norm, nB_residual, ud_residual, s_residual)

    omega_val = -pressure_ref[]
    thermo_finite = isfinite(omega_val) && isfinite(pressure_ref[]) && isfinite(rho_norm_ref[]) && isfinite(entropy_ref[]) && isfinite(energy_ref[])
    phys = physicality_check(x_state_ref[], masses_ref[]) && thermo_finite
    converged = if enforce_physicality
        res.f_converged && phys && isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-3)
    else
        isfinite(residual_norm) && residual_norm <= max(Float64(residual_norm_max), 1e-3)
    end

    return (
        converged=converged,
        solution=Float64[
            x_state_ref[][1], x_state_ref[][2], x_state_ref[][3], x_state_ref[][4], x_state_ref[][5],
            mu_vec_ref[][1], mu_vec_ref[][2], mu_vec_ref[][3],
        ],
        x_state=x_state_ref[],
        mu_vec=mu_vec_ref[],
        omega=omega_val,
        pressure=pressure_ref[],
        rho_norm=rho_norm_ref[],
        entropy=entropy_ref[],
        energy=energy_ref[],
        masses=masses_ref[],
        iterations=res.iterations,
        residual_norm=residual_norm,
    )
end
