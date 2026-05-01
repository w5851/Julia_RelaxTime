raw"""
MesonDensity

介子数密度最小实现入口。当前阶段仅覆盖稳定粒子极限：

- 玻色分布函数
- `π/K` 聚合通道默认简并因子
- 稳定粒子极限数密度
- `K/π` 比值与最小温度扫描

后续 BU / BW / 各向异性扩展在此基础上演进。
"""
module MesonDensity

using ..GaussLegendre: gauleg

export DEFAULT_MESON_DENSITY_Q_NODES
export bose_distribution, meson_degeneracy
export stable_meson_number_density, stable_kpi_ratio, stable_kpi_scan

const DEFAULT_MESON_DENSITY_Q_NODES = 256

@inline function _require_nonnegative(name::AbstractString, value::Float64)
    value >= 0.0 && return
    throw(ArgumentError("$(name) must be nonnegative, got $(value)"))
end

@inline function _default_qmax(mass::Float64, T::Float64, μ::Float64)::Float64
    gap = max(mass - μ, 0.0)
    return max(8.0, 20.0 * T + 10.0 * gap, 8.0 * mass + 10.0 * T)
end

"""
    meson_degeneracy(meson; charge_resolved=false) -> Int

返回当前主线下的 `π/K` 简并因子。

- 聚合通道（默认）：`d_π = 3`、`d_K = 4`
- 电荷分辨通道：`d = 1`
"""
@inline function meson_degeneracy(meson::Symbol; charge_resolved::Bool=false)::Int
    if meson === :pi
        return charge_resolved ? 1 : 3
    elseif meson === :K
        return charge_resolved ? 1 : 4
    else
        throw(ArgumentError("Unsupported meson $(meson). Use :pi or :K."))
    end
end

raw"""
    bose_distribution(E, μ, T) -> Float64

玻色分布函数：

```math
g(E) = 1 / (\exp((E-\mu)/T) - 1)
```
"""
function bose_distribution(E::Float64, μ::Float64, T::Float64)::Float64
    _require_nonnegative("temperature T", T)
    T == 0.0 && return 0.0
    E > μ || throw(ArgumentError("Bose distribution requires E > μ to avoid pole, got E=$(E), μ=$(μ)"))

    exponent = (E - μ) / T
    exponent > 700.0 && return 0.0
    return 1.0 / expm1(exponent)
end

raw"""
    stable_meson_number_density(mass, T; μ=0.0, degeneracy=1,
                                qmax=nothing, num_q_nodes=DEFAULT_MESON_DENSITY_Q_NODES) -> Float64

稳定粒子极限介子数密度：

```math
n_M = d_M \\int_0^\\infty \\frac{dq\\,q^2}{2\\pi^2}
      \\frac{1}{\\exp((E_M-\\mu_M)/T)-1},
\\qquad
E_M = \\sqrt{q^2 + m_M^2}.
```
"""
function stable_meson_number_density(mass::Float64, T::Float64;
                                     μ::Float64=0.0,
                                     degeneracy::Integer=1,
                                     qmax::Union{Nothing,Float64}=nothing,
                                     num_q_nodes::Int=DEFAULT_MESON_DENSITY_Q_NODES)::Float64
    _require_nonnegative("mass", mass)
    _require_nonnegative("temperature T", T)
    degeneracy > 0 || throw(ArgumentError("degeneracy must be positive, got $(degeneracy)"))
    num_q_nodes > 1 || throw(ArgumentError("num_q_nodes must be > 1, got $(num_q_nodes)"))
    mass > μ || throw(ArgumentError("Stable boson density requires mass > μ, got mass=$(mass), μ=$(μ)"))
    T == 0.0 && return 0.0

    q_upper = qmax === nothing ? _default_qmax(mass, T, μ) : Float64(qmax)
    q_upper > 0.0 || throw(ArgumentError("qmax must be positive, got $(q_upper)"))

    nodes, weights = gauleg(0.0, q_upper, num_q_nodes)
    integral = 0.0
    @inbounds for i in eachindex(nodes, weights)
        q = nodes[i]
        E = hypot(q, mass)
        integral += weights[i] * q^2 * bose_distribution(E, μ, T)
    end
    return degeneracy * integral / (2.0 * π^2)
end

"""
    stable_kpi_ratio(m_pi, m_K, T; μ_pi=0.0, μ_K=0.0, d_pi=3, d_K=4, kwargs...) -> Float64

返回稳定粒子极限的 `K/π` 数密度比值。
"""
function stable_kpi_ratio(m_pi::Float64, m_K::Float64, T::Float64;
                          μ_pi::Float64=0.0, μ_K::Float64=0.0,
                          d_pi::Integer=3, d_K::Integer=4,
                          kwargs...)::Float64
    n_pi = stable_meson_number_density(m_pi, T; μ=μ_pi, degeneracy=d_pi, kwargs...)
    n_K = stable_meson_number_density(m_K, T; μ=μ_K, degeneracy=d_K, kwargs...)
    return iszero(n_pi) ? NaN : n_K / n_pi
end

"""
    stable_kpi_scan(temperatures; m_pi, m_K, μ_pi=0.0, μ_K=0.0, d_pi=3, d_K=4, kwargs...)
        -> NamedTuple

对一组温度执行稳定粒子极限 `π/K` 数密度扫描，返回：

- `temperatures`
- `n_pi`
- `n_K`
- `kpi_ratio`
"""
function stable_kpi_scan(temperatures::AbstractVector{<:Real};
                         m_pi::Float64, m_K::Float64,
                         μ_pi::Float64=0.0, μ_K::Float64=0.0,
                         d_pi::Integer=3, d_K::Integer=4,
                         kwargs...)
    Ts = Float64[Float64(T) for T in temperatures]
    n_pi = Vector{Float64}(undef, length(Ts))
    n_K = Vector{Float64}(undef, length(Ts))
    ratios = Vector{Float64}(undef, length(Ts))

    for i in eachindex(Ts)
        Ti = Ts[i]
        n_pi[i] = stable_meson_number_density(m_pi, Ti; μ=μ_pi, degeneracy=d_pi, kwargs...)
        n_K[i] = stable_meson_number_density(m_K, Ti; μ=μ_K, degeneracy=d_K, kwargs...)
        ratios[i] = iszero(n_pi[i]) ? NaN : n_K[i] / n_pi[i]
    end

    return (temperatures=Ts, n_pi=n_pi, n_K=n_K, kpi_ratio=ratios)
end

end # module MesonDensity
