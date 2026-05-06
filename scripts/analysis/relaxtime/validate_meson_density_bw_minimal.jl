"""
Phase D3: 最小 BW（Breit-Wigner）过渡验证脚本。

目标：
- 复用 meson workflow 的 `M_pi/M_K/Gamma_pi/Gamma_K`
- 对比稳定粒子极限与最小 BW 质量展宽近似数密度
- 生成 `stable vs BW` 的最小 CSV

当前固定口径：
- muB = 0
- xi = 0
- pi/K 两道
- 对脚本级最小实现，先采用 BW 质量分布 smearing 代理
"""

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))

using .Constants_PNJL: ħc_MeV_fm
using .Models: solve_gap_and_meson_point
using Main.MesonDensity: meson_degeneracy, stable_meson_number_density
using QuadGK

const DEFAULT_OUTPUT = joinpath(
    PROJECT_ROOT,
    "data",
    "outputs",
    "results",
    "relaxtime",
    "analysis",
    "meson_density_bw_minimal.csv",
)

const DEFAULT_TS_MEV = collect(170.0:5.0:210.0)
const DEFAULT_Q_NODES = 160
const DEFAULT_RTOL = 1e-5

@inline function _fmt(x)
    x isa Bool && return x ? "true" : "false"
    return string(x)
end

@inline function _default_qmax(mass::Float64, T::Float64, μ::Float64)::Float64
    gap = max(mass - μ, 0.0)
    return max(8.0, 20.0 * T + 10.0 * gap, 8.0 * mass + 10.0 * T)
end

@inline function _default_mmax(mass::Float64, gamma::Float64, T::Float64)::Float64
    return max(8.0, mass + 20.0 * abs(gamma) + 20.0 * T, 8.0 * mass + 20.0 * T)
end

@inline function _bw_mass_kernel(m::Float64, mass0::Float64, gamma::Float64)::Float64
    return (gamma / 2.0) / ((m - mass0)^2 + gamma^2 / 4.0)
end

function bw_meson_number_density(
    mass::Float64,
    gamma::Float64,
    T::Float64;
    degeneracy::Int,
    qmax::Union{Nothing,Float64}=nothing,
    mmax::Union{Nothing,Float64}=nothing,
    mass_rtol::Float64=DEFAULT_RTOL,
)
    mass >= 0.0 || throw(ArgumentError("mass must be nonnegative, got $(mass)"))
    T >= 0.0 || throw(ArgumentError("temperature T must be nonnegative, got $(T)"))
    degeneracy > 0 || throw(ArgumentError("degeneracy must be positive, got $(degeneracy)"))
    T == 0.0 && return 0.0

    gamma_eff = abs(Float64(gamma))
    gamma_eff <= 1e-8 && return stable_meson_number_density(mass, T; degeneracy=degeneracy, num_q_nodes=DEFAULT_Q_NODES)

    q_upper = qmax === nothing ? _default_qmax(mass, T, 0.0) : Float64(qmax)
    m_upper = mmax === nothing ? _default_mmax(mass, gamma_eff, T) : Float64(mmax)
    q_upper > 0.0 || throw(ArgumentError("qmax must be positive, got $(q_upper)"))
    m_upper > 0.0 || throw(ArgumentError("mmax must be positive, got $(m_upper)"))

    norm_val, = quadgk(m -> _bw_mass_kernel(m, mass, gamma_eff), 0.0, m_upper; rtol=mass_rtol)
    norm_val > 0.0 || throw(ArgumentError("BW mass kernel normalization must be positive, got $(norm_val)"))

    mass_integrand = function (m)
        rho = _bw_mass_kernel(m, mass, gamma_eff) / norm_val
        n_m = stable_meson_number_density(m, T; degeneracy=1, qmax=q_upper, num_q_nodes=DEFAULT_Q_NODES)
        return rho * n_m
    end

    mass_val, = quadgk(mass_integrand, 0.0, m_upper; rtol=mass_rtol)
    return degeneracy * mass_val
end

function _parse_t_values(raw::AbstractString)
    vals = Float64[]
    for seg in split(raw, ',')
        s = strip(seg)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    isempty(vals) && throw(ArgumentError("temperature list cannot be empty"))
    return vals
end

function _selected_temperatures()
    raw = strip(get(ENV, "MESON_D3_T_VALUES", ""))
    isempty(raw) && return DEFAULT_TS_MEV
    return _parse_t_values(raw)
end

function main()
    output = isempty(ARGS) ? DEFAULT_OUTPUT : abspath(ARGS[1])
    mkpath(dirname(output))

    Ts_MeV = _selected_temperatures()
    continuation_state = nothing
    rows = NamedTuple[]

    for T_MeV in Ts_MeV
        T_fm = T_MeV / ħc_MeV_fm
        res = solve_gap_and_meson_point(
            T_fm,
            0.0;
            xi=0.0,
            mesons=(:pi, :K),
            continuation_state=continuation_state,
            mixed_branch_align=:strict_sign_binding,
            p_num=8,
            t_num=4,
            solver_kwargs=(; iterations=20),
            mass_kwargs=(; iterations=20),
        )

        mpi = res.meson_results[:pi]
        mK = res.meson_results[:K]

        d_pi = meson_degeneracy(:pi)
        d_K = meson_degeneracy(:K)

        n_pi_stable = stable_meson_number_density(Float64(mpi.mass), T_fm; degeneracy=d_pi, num_q_nodes=DEFAULT_Q_NODES)
        n_K_stable = stable_meson_number_density(Float64(mK.mass), T_fm; degeneracy=d_K, num_q_nodes=DEFAULT_Q_NODES)
        n_pi_bw = bw_meson_number_density(Float64(mpi.mass), Float64(mpi.gamma), T_fm; degeneracy=d_pi)
        n_K_bw = bw_meson_number_density(Float64(mK.mass), Float64(mK.gamma), T_fm; degeneracy=d_K)

        ratio_stable = iszero(n_pi_stable) ? NaN : n_K_stable / n_pi_stable
        ratio_bw = iszero(n_pi_bw) ? NaN : n_K_bw / n_pi_bw

        push!(rows, (
            T_MeV=T_MeV,
            m_pi=Float64(mpi.mass),
            gamma_pi=Float64(mpi.gamma),
            n_pi_stable=n_pi_stable,
            n_pi_bw=n_pi_bw,
            delta_n_pi=n_pi_bw - n_pi_stable,
            bw_over_stable_pi=iszero(n_pi_stable) ? NaN : n_pi_bw / n_pi_stable,
            m_K=Float64(mK.mass),
            gamma_K=Float64(mK.gamma),
            n_K_stable=n_K_stable,
            n_K_bw=n_K_bw,
            delta_n_K=n_K_bw - n_K_stable,
            bw_over_stable_K=iszero(n_K_stable) ? NaN : n_K_bw / n_K_stable,
            kpi_stable=ratio_stable,
            kpi_bw=ratio_bw,
            delta_kpi=ratio_bw - ratio_stable,
        ))

        continuation_state = res.continuation_state
    end

    open(output, "w") do io
        println(io, "T_MeV,m_pi,gamma_pi,n_pi_stable,n_pi_bw,delta_n_pi,bw_over_stable_pi,m_K,gamma_K,n_K_stable,n_K_bw,delta_n_K,bw_over_stable_K,kpi_stable,kpi_bw,delta_kpi")
        for r in rows
            println(io, join((_fmt(getproperty(r, k)) for k in (
                :T_MeV, :m_pi, :gamma_pi, :n_pi_stable, :n_pi_bw, :delta_n_pi, :bw_over_stable_pi,
                :m_K, :gamma_K, :n_K_stable, :n_K_bw, :delta_n_K, :bw_over_stable_K,
                :kpi_stable, :kpi_bw, :delta_kpi
            )), ','))
        end
    end

    println("Wrote Phase D3 BW validation CSV: ", output)
    for r in rows
        println(
            "T=$(r.T_MeV) MeV | ",
            "K/pi stable=$(r.kpi_stable), ",
            "BW=$(r.kpi_bw), ",
            "delta=$(r.delta_kpi)"
        )
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
