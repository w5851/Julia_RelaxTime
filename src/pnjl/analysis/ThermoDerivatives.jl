module ThermoDerivatives

using ForwardDiff
using StaticArrays
using NLsolve

using ..Constants_PNJL: ħc_MeV_fm
using ..AnisoGapSolver:
    cached_nodes,
    calculate_thermo,
    calculate_rho,
    calculate_mass_vec,
    calculate_energy_isotropic,
    residual_mu!,
    DEFAULT_MU_GUESS,
    DEFAULT_MOMENTUM_COUNT,
    DEFAULT_THETA_COUNT

export solve_equilibrium_mu, thermo_derivatives, dP_dT, dP_dmu, dE_dT, dE_dmu, dEpsilon_dT, dEpsilon_dmu, dn_dT, dn_dmu, bulk_derivative_coeffs, quasiparticle_energy, mass_derivatives

"""将种子向量转换到当前参数的数值类型，便于在自动微分中复用"""
function _convert_seed(seed_state, T_mev, mu_mev)
    seed_type = promote_type(typeof(T_mev), typeof(mu_mev), Float64)
    return collect(seed_type.(seed_state))
end

"""在指定 (T, μ, ξ) 下求解能隙方程，并返回热力学量。

自动微分调用时会把 T、μ 包装成 Dual，这里不强制转换为 Float64，
以便导数能穿透求解过程。
"""
function solve_equilibrium_mu(T_mev::Real, mu_mev::Real; xi::Real=0.0, seed_state=DEFAULT_MU_GUESS, p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT, nlsolve_kwargs...)
    T_fm = T_mev / ħc_MeV_fm
    mu_fm = mu_mev / ħc_MeV_fm
    thermal_nodes = cached_nodes(p_num, t_num)
    seed_vec = _convert_seed(seed_state, T_mev, mu_mev)
    mu_vec = SVector{3}(mu_fm, mu_fm, mu_fm)

    f! = (F, x) -> residual_mu!(F, x, mu_vec, T_fm, thermal_nodes, xi)
    extra = (; nlsolve_kwargs...)
    res = nlsolve(f!, seed_vec; autodiff=:forward, method=:newton, xtol=1e-9, ftol=1e-9, extra...)

    x_state = SVector{5}(Tuple(res.zero))
    pressure, rho_norm, entropy, energy = calculate_thermo(x_state, mu_vec, T_fm, thermal_nodes, xi)
    rho_vec = calculate_rho(x_state, mu_vec, T_fm, thermal_nodes, xi)
    rho_total = sum(rho_vec) / 3

    return (
        pressure=pressure,
        energy=energy,
        rho=rho_total,
        rho_norm=rho_norm,
        entropy=entropy,
        x_state=x_state,
        mu_vec=mu_vec,
        converged=res.f_converged,
        iterations=res.iterations,
        residual_norm=res.residual_norm,
    )
end

"""递归计算 n 阶导数，使用 ForwardDiff 嵌套调用"""
function _nth_derivative(f, x, n::Int)
    n >= 1 || error("order must be ≥ 1")
    n == 1 && return ForwardDiff.derivative(f, x)
    inner = y -> _nth_derivative(f, y, n - 1)
    return ForwardDiff.derivative(inner, x)
end

"""总压强对温度的 n 阶导数，内部会重新求解能隙以保持平衡"""
function dP_dT(T_mev::Real, mu_mev::Real; order::Int=1, xi::Real=0.0, seed_state=DEFAULT_MU_GUESS, kwargs...)
    g = t -> solve_equilibrium_mu(t, mu_mev; xi=xi, seed_state=seed_state, kwargs...).pressure
    return _nth_derivative(g, T_mev, order)
end

"""总压强对化学势的 n 阶导数"""
function dP_dmu(T_mev::Real, mu_mev::Real; order::Int=1, xi::Real=0.0, seed_state=DEFAULT_MU_GUESS, kwargs...)
    g = μ -> solve_equilibrium_mu(T_mev, μ; xi=xi, seed_state=seed_state, kwargs...).pressure
    return _nth_derivative(g, mu_mev, order)
end

"""能量密度 ε 对温度的 n 阶导数"""
function dEpsilon_dT(T_mev::Real, mu_mev::Real; order::Int=1, xi::Real=0.0, seed_state=DEFAULT_MU_GUESS, kwargs...)
    g = t -> solve_equilibrium_mu(t, mu_mev; xi=xi, seed_state=seed_state, kwargs...).energy
    return _nth_derivative(g, T_mev, order)
end

"""能量密度 ε 对化学势的 n 阶导数"""
function dEpsilon_dmu(T_mev::Real, mu_mev::Real; order::Int=1, xi::Real=0.0, seed_state=DEFAULT_MU_GUESS, kwargs...)
    g = μ -> solve_equilibrium_mu(T_mev, μ; xi=xi, seed_state=seed_state, kwargs...).energy
    return _nth_derivative(g, mu_mev, order)
end

"""给定动量 p_fm 和味道，返回各向同性色散关系 E = sqrt(p^2 + m^2)（m 来自平衡解）。"""
function quasiparticle_energy(T_mev::Real, mu_mev::Real, p_fm::Real; flavor::Int=1, xi::Real=0.0, seed_state=DEFAULT_MU_GUESS, p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT, kwargs...)
    1 <= flavor <= 3 || error("flavor must be 1,2,3")
    base = solve_equilibrium_mu(T_mev, mu_mev; xi=xi, seed_state=seed_state, p_num=p_num, t_num=t_num, kwargs...)
    masses = calculate_mass_vec(base.x_state)
    return calculate_energy_isotropic(masses[flavor], p_fm)
end

"""返回质量及其对 T/μ 的导数（三味）。"""
function mass_derivatives(T_mev::Real, mu_mev::Real; xi::Real=0.0, seed_state=DEFAULT_MU_GUESS, p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT, kwargs...)
    mass_fn_T = t -> begin
        res = solve_equilibrium_mu(t, mu_mev; xi=xi, seed_state=seed_state, p_num=p_num, t_num=t_num, kwargs...)
        calculate_mass_vec(SVector{3}(res.x_state[1:3]))
    end
    mass_fn_mu = μ -> begin
        res = solve_equilibrium_mu(T_mev, μ; xi=xi, seed_state=seed_state, p_num=p_num, t_num=t_num, kwargs...)
        calculate_mass_vec(SVector{3}(res.x_state[1:3]))
    end
    masses = mass_fn_T(T_mev)
    dM_dT = SVector{3}(ForwardDiff.derivative(t -> mass_fn_T(t)[i], T_mev) for i in 1:3)
    dM_dmu = SVector{3}(ForwardDiff.derivative(μ -> mass_fn_mu(μ)[i], mu_mev) for i in 1:3)
    return (masses=masses, dM_dT=dM_dT, dM_dmu=dM_dmu)
end

"""单粒子能量对温度的导数（链式法则，传入已知 m 与 dM/dT）。"""
function dE_dT(T_mev::Real, mu_mev::Real, p_fm::Real; m::Real, dM_dT::Real)
    E = calculate_energy_isotropic(m, p_fm)
    return (m / E) * dM_dT
end

"""单粒子能量对化学势的导数（链式法则，传入已知 m 与 dM/dμ）。"""
function dE_dmu(T_mev::Real, mu_mev::Real, p_fm::Real; m::Real, dM_dmu::Real)
    E = calculate_energy_isotropic(m, p_fm)
    return (m / E) * dM_dmu
end

"""总粒子数密度 (未归一化) 对温度的一阶导数"""
function dn_dT(T_mev::Real, mu_mev::Real; xi::Real=0.0, seed_state=DEFAULT_MU_GUESS, kwargs...)
    g = t -> solve_equilibrium_mu(t, mu_mev; xi=xi, seed_state=seed_state, kwargs...).rho
    return ForwardDiff.derivative(g, T_mev)
end

"""总粒子数密度 (未归一化) 对化学势的一阶导数"""
function dn_dmu(T_mev::Real, mu_mev::Real; xi::Real=0.0, seed_state=DEFAULT_MU_GUESS, kwargs...)
    g = μ -> solve_equilibrium_mu(T_mev, μ; xi=xi, seed_state=seed_state, kwargs...).rho
    return ForwardDiff.derivative(g, mu_mev)
end

"""同时返回基础热力学量及其一阶导数组合，便于体粘滞系数公式直接调用"""
function thermo_derivatives(T_mev::Real, mu_mev::Real; xi::Real=0.0, seed_state=DEFAULT_MU_GUESS, p_num::Int=DEFAULT_MOMENTUM_COUNT, t_num::Int=DEFAULT_THETA_COUNT, kwargs...)
    base = solve_equilibrium_mu(T_mev, mu_mev; xi=xi, seed_state=seed_state, p_num=p_num, t_num=t_num, kwargs...)
    seed_for_diff = base.x_state

    P_T = dP_dT(T_mev, mu_mev; xi=xi, seed_state=seed_for_diff, p_num=p_num, t_num=t_num, kwargs...)
    P_mu = dP_dmu(T_mev, mu_mev; xi=xi, seed_state=seed_for_diff, p_num=p_num, t_num=t_num, kwargs...)
    E_T = dEpsilon_dT(T_mev, mu_mev; xi=xi, seed_state=seed_for_diff, p_num=p_num, t_num=t_num, kwargs...)
    E_mu = dEpsilon_dmu(T_mev, mu_mev; xi=xi, seed_state=seed_for_diff, p_num=p_num, t_num=t_num, kwargs...)
    n_T = dn_dT(T_mev, mu_mev; xi=xi, seed_state=seed_for_diff, p_num=p_num, t_num=t_num, kwargs...)
    n_mu = dn_dmu(T_mev, mu_mev; xi=xi, seed_state=seed_for_diff, p_num=p_num, t_num=t_num, kwargs...)

    denom_eps = E_T * n_mu - E_mu * n_T
    denom_n = n_T * E_mu - n_mu * E_T
    dP_depsilon_n = denom_eps == 0 ? NaN : (P_T * n_mu - P_mu * n_T) / denom_eps
    dP_dn_epsilon = denom_n == 0 ? NaN : (P_T * E_mu - P_mu * E_T) / denom_n

    md = mass_derivatives(T_mev, mu_mev; xi=xi, seed_state=seed_for_diff, p_num=p_num, t_num=t_num, kwargs...)

    return (
        pressure=base.pressure,
        energy=base.energy,
        rho=base.rho,
        rho_norm=base.rho_norm,
        entropy=base.entropy,
        dP_dT=P_T,
        dP_dmu=P_mu,
        dEpsilon_dT=E_T,
        dEpsilon_dmu=E_mu,
        dn_dT=n_T,
        dn_dmu=n_mu,
        dP_depsilon_n=dP_depsilon_n,
        dP_dn_epsilon=dP_dn_epsilon,
        masses=md.masses,
        dM_dT=md.dM_dT,
        dM_dmu=md.dM_dmu,
        converged=base.converged,
        iterations=base.iterations,
        residual_norm=base.residual_norm,
    )
end

"""返回体粘滞系数公式中常用的两个组合导数"""
function bulk_derivative_coeffs(T_mev::Real, mu_mev::Real; xi::Real=0.0, seed_state=DEFAULT_MU_GUESS, kwargs...)
    derivs = thermo_derivatives(T_mev, mu_mev; xi=xi, seed_state=seed_state, kwargs...)
    return (
        dP_depsilon_n=derivs.dP_depsilon_n,
        dP_dn_epsilon=derivs.dP_dn_epsilon,
        dM_dT=derivs.dM_dT,
        dM_dmu=derivs.dM_dmu,
    )
end

end # module
