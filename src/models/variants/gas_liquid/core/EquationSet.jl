module GasLiquidEquationSet

using NLsolve
using StaticArrays

export GasLiquidCoreParams
export GasLiquidState
export state_vector
export solve_equilibrium

Base.@kwdef struct GasLiquidCoreParams
    m_nucleon_inv_fm::Float64 = 939.0 / 197.3269804
    g_sigma::Float64 = 2.2
    g_delta::Float64 = 0.6
    sigma_target_scale::Float64 = 0.22
end

struct GasLiquidState{T}
    sigma::T
    delta::T
    mu_p::T
    mu_n::T
end

@inline state_vector(st::GasLiquidState{T}) where {T} = SVector{4, T}(st.sigma, st.delta, st.mu_p, st.mu_n)

@inline function _mu_baryon(mu_vec)
    if mu_vec isa Real
        return float(mu_vec)
    end
    length(mu_vec) == 3 || throw(ArgumentError("mu_vec must be Real or length-3 vector"))
    # Models 主线通常用三味夸克化学势；这里按 baryon-like 口径先取平均值。
    return float((mu_vec[1] + mu_vec[2] + mu_vec[3]) / 3)
end

@inline function _sigma_target(T::Real, muB::Real, p::GasLiquidCoreParams)
    T_safe = max(float(T), 1e-8)
    return p.sigma_target_scale * tanh((muB - T_safe) / (T_safe + 0.3))
end

function solve_equilibrium(
    T::Real,
    mu_vec;
    params::GasLiquidCoreParams=GasLiquidCoreParams(),
    initial_guess::Union{Nothing,AbstractVector}=nothing,
)
    muB = _mu_baryon(mu_vec)
    x0 = initial_guess === nothing ? [0.05, 0.0, muB, muB] : collect(float.(initial_guess))
    length(x0) == 4 || throw(ArgumentError("initial_guess must have length 4"))

    function residual!(F, x)
        sigma_t = _sigma_target(T, muB, params)
        F[1] = x[1] - sigma_t
        F[2] = x[2]
        F[3] = x[3] - muB
        F[4] = x[4] - muB
        return nothing
    end

    res = nlsolve(residual!, x0; autodiff=:forward, method=:newton, xtol=1e-10, ftol=1e-10)
    if !(res.f_converged && isfinite(res.residual_norm))
        error("gas-liquid equilibrium solve failed: f_converged=$(res.f_converged), residual_norm=$(res.residual_norm)")
    end

    x = res.zero
    return GasLiquidState(float(x[1]), float(x[2]), float(x[3]), float(x[4]))
end

end # module
