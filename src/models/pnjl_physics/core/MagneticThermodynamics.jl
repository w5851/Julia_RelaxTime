"""
    MagneticThermodynamics

外磁场下 PNJL 热力学主接口。

核心公式对应：
- Ω_f^0, Ω_f^T (Landau 求和)
- E_{f,n} = sqrt(2n|q_f|eB + pz^2 + M_f^2)
- ρ_f = ∂P/∂μ_f（提供数值导数与低温近似两种路径）
- G(B) IMC 参数化
"""
module MagneticThermodynamics

using StaticArrays

const _CONSTANTS_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "constants", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PATH)
end
using Main.Constants_PNJL:
    G_fm2,
    K_fm5,
    m_ud0_inv_fm,
    m_s0_inv_fm,
    ħc_MeV_fm,
    ρ0_inv_fm3

const _PNJL_CORE_PATH = normpath(joinpath(@__DIR__, "..", "PNJLCore.jl"))
if !isdefined(@__MODULE__, :PNJLCore)
    include(_PNJL_CORE_PATH)
end

include("MagneticIntegrals.jl")
using .MagneticIntegrals:
    QUARK_CHARGE_ABS,
    MAGNETIC_EB_MIN_MEV2,
    MAGNETIC_EB_MIN_FM2,
    validate_magnetic_eB,
    resolve_nmax_from_cutoff,
    omega0_flavor_landau,
    omegat_flavor_landau,
    density_flavor_landau

export MagneticIMCParams, default_imc_params, coupling_GB
export MagneticConfig, default_magnetic_config
export MAGNETIC_EB_MIN_MEV2, MAGNETIC_EB_MIN_FM2, validate_magnetic_eB
export calculate_magnetic_omega_components, calculate_magnetic_omega
export calculate_magnetic_pressure, calculate_magnetic_rho
export calculate_magnetic_number_densities
export magnetic_nmax_convergence_report

const ρ0 = ρ0_inv_fm3
const _PNJL_PARAMS_REF = Ref{Any}(nothing)

@inline function _pnjl_params()
    p = _PNJL_PARAMS_REF[]
    if p === nothing
        p = PNJLCore.pnjl_params()
        _PNJL_PARAMS_REF[] = p
    end
    return p
end

Base.@kwdef struct MagneticIMCParams
    a::Float64 = 0.108805
    b::Float64 = -1.0133e-4
    c::Float64 = 0.02228
    d::Float64 = 1.84558e-4
    Λ_QCD_MeV::Float64 = 300.0
end

default_imc_params() = MagneticIMCParams()

struct MagneticConfig
    eB_fm2::Float64
    n_max::Union{Int, Nothing}
    p_num::Int
    pz_max::Float64
    cutoff_N::Int
    imc::MagneticIMCParams

    function MagneticConfig(
        eB_fm2::Real,
        n_max::Union{Nothing, Int}=nothing,
        p_num::Int=96,
        pz_max::Real=0.0,
        cutoff_N::Int=10,
        imc::MagneticIMCParams=default_imc_params(),
    )
        eB_value = validate_magnetic_eB(eB_fm2)
        new(eB_value, n_max, p_num, Float64(pz_max), cutoff_N, imc)
    end
end

function MagneticConfig(
    ;
    eB_fm2::Real=0.0,
    n_max::Union{Nothing, Int}=nothing,
    p_num::Int=96,
    pz_max::Real=0.0,
    cutoff_N::Int=10,
    imc::MagneticIMCParams=default_imc_params(),
)
    return MagneticConfig(eB_fm2, n_max, p_num, pz_max, cutoff_N, imc)
end

default_magnetic_config(; eB_fm2::Real=0.0) = MagneticConfig(eB_fm2=float(eB_fm2))

@inline function _Λ_QCD_inv_fm(params::MagneticIMCParams)
    return params.Λ_QCD_MeV / ħc_MeV_fm
end

@inline function coupling_GB(eB_fm2::Real; G0::Real=G_fm2, imc::MagneticIMCParams=default_imc_params())
    eB_value = validate_magnetic_eB(eB_fm2)
    ζ = eB_value / (_Λ_QCD_inv_fm(imc)^2)
    num = 1 + imc.a * ζ^2 + imc.b * ζ^3
    den = 1 + imc.c * ζ^2 + imc.d * ζ^4
    return float(G0) * num / den
end

@inline function _calculate_mass_vec_with_GB(φ::SVector{3, T}, G_B::Real) where {T}
    φ_u, φ_d, φ_s = φ
    Gv = convert(T, G_B)
    Kv = convert(T, K_fm5)
    return SVector{3, T}(
        m_ud0_inv_fm - 4 * Gv * φ_u + 2 * Kv * φ_d * φ_s,
        m_ud0_inv_fm - 4 * Gv * φ_d + 2 * Kv * φ_u * φ_s,
        m_s0_inv_fm - 4 * Gv * φ_s + 2 * Kv * φ_u * φ_d,
    )
end

@inline function _chiral_with_GB(φ::SVector{3, T}, G_B::Real) where {T}
    Gv = convert(T, G_B)
    Kv = convert(T, K_fm5)
    return 2 * Gv * sum(φ .^ 2) - 4 * Kv * prod(φ)
end

@inline function _resolve_nmax(masses::SVector{3, <:Real}, mu_vec::SVector{3, <:Real}, eB::Real, conf::MagneticConfig)
    conf.n_max !== nothing && return conf.n_max::Int
    eB_value = validate_magnetic_eB(eB)
    nvals = ntuple(i -> resolve_nmax_from_cutoff(masses[i], mu_vec[i], QUARK_CHARGE_ABS[i], eB_value), 3)
    return max(maximum(nvals), 3)
end

@inline function _validate_magnetic_controls(
    T_fm::Real,
    xi::Real,
    p_num::Int,
    t_num::Int,
)
    T_value = Float64(T_fm)
    isfinite(T_value) && T_value > 0.0 || throw(ArgumentError(
        "magnetic thermodynamics requires finite T_fm > 0, got $(T_fm)",
    ))
    isfinite(Float64(xi)) || throw(ArgumentError("magnetic route requires finite xi, got $(xi)"))
    abs(Float64(xi)) <= 1e-14 || throw(ArgumentError(
        "PNJLMagneticModel currently supports only xi=0; got xi=$(xi)",
    ))
    p_num >= 4 || throw(ArgumentError("magnetic p_num must be >= 4, got $(p_num)"))
    t_num >= 1 || throw(ArgumentError("magnetic t_num must be >= 1, got $(t_num)"))
    return nothing
end

@inline function _effective_magnetic_config(
    magnetic::MagneticConfig;
    p_num::Union{Nothing, Int}=nothing,
    pz_max::Union{Nothing, Real}=nothing,
    n_max::Union{Nothing, Int}=nothing,
    cutoff_N::Union{Nothing, Int}=nothing,
)
    p_num_value = p_num === nothing ? magnetic.p_num : Int(p_num)
    pz_max_value = pz_max === nothing ? magnetic.pz_max : Float64(pz_max)
    n_max_value = n_max === nothing ? magnetic.n_max : Int(n_max)
    cutoff_value = cutoff_N === nothing ? magnetic.cutoff_N : Int(cutoff_N)
    p_num_value >= 4 || throw(ArgumentError("magnetic p_num must be >= 4, got $(p_num_value)"))
    pz_max_value >= 0.0 || throw(ArgumentError("magnetic pz_max must be nonnegative, got $(pz_max_value)"))
    n_max_value === nothing || n_max_value >= 0 || throw(ArgumentError("magnetic n_max must be >= 0, got $(n_max_value)"))
    cutoff_value >= 1 || throw(ArgumentError("magnetic cutoff_N must be >= 1, got $(cutoff_value)"))
    return MagneticConfig(
        eB_fm2=magnetic.eB_fm2,
        n_max=n_max_value,
        p_num=p_num_value,
        pz_max=pz_max_value,
        cutoff_N=cutoff_value,
        imc=magnetic.imc,
    )
end

function calculate_magnetic_omega_components(
    x_state::SVector{5, T},
    mu_vec::SVector{3, M},
    T_fm::Real,
    magnetic::MagneticConfig=default_magnetic_config();
    xi::Real=0.0,
    p_num::Union{Nothing, Int}=nothing,
    t_num::Int=8,
    pz_max::Union{Nothing, Real}=nothing,
    n_max::Union{Nothing, Int}=nothing,
    cutoff_N::Union{Nothing, Int}=nothing,
) where {T, M}
    _validate_magnetic_controls(T_fm, xi, p_num === nothing ? magnetic.p_num : p_num, t_num)
    conf = _effective_magnetic_config(
        magnetic;
        p_num=p_num,
        pz_max=pz_max,
        n_max=n_max,
        cutoff_N=cutoff_N,
    )
    φ = SVector{3, T}(x_state[1], x_state[2], x_state[3])
    Φ, Φbar = x_state[4], x_state[5]

    G_B = coupling_GB(conf.eB_fm2; imc=conf.imc)
    masses = _calculate_mass_vec_with_GB(φ, G_B)
    n_max = _resolve_nmax(masses, mu_vec, conf.eB_fm2, conf)

    pz_max_val = conf.pz_max > 0 ? conf.pz_max : max(8 * Main.Constants_PNJL.Λ_inv_fm, 25.0)

    vac = 0.0
    therm = 0.0
    @inbounds for i in 1:3
        vac += omega0_flavor_landau(
            masses[i],
            QUARK_CHARGE_ABS[i],
            conf.eB_fm2;
            n_max=n_max,
            p_num=conf.p_num,
            pz_max=pz_max_val,
            cutoff_N=conf.cutoff_N,
        )
        therm += omegat_flavor_landau(
            masses[i],
            mu_vec[i],
            T_fm,
            Φ,
            Φbar,
            QUARK_CHARGE_ABS[i],
            conf.eB_fm2;
            n_max=n_max,
            p_num=conf.p_num,
            pz_max=pz_max_val,
        )
    end

    chi = _chiral_with_GB(φ, G_B)
    poly = PNJLCore.polyakov_potential(_pnjl_params(), Φ, Φbar, T_fm)
    ω = chi + poly + vac + therm
    return (chi=chi, poly=poly, vac=vac, therm=therm, masses=masses, omega=ω, n_max=n_max, G_B=G_B)
end

@inline function calculate_magnetic_omega(x_state, mu_vec, T_fm, magnetic::MagneticConfig=default_magnetic_config(); kwargs...)
    st = SVector{5, Float64}(Tuple(x_state))
    μ = SVector{3, Float64}(Tuple(mu_vec))
    return calculate_magnetic_omega_components(st, μ, T_fm, magnetic; kwargs...).omega
end

@inline function calculate_magnetic_pressure(x_state, mu_vec, T_fm, magnetic::MagneticConfig=default_magnetic_config(); kwargs...)
    return -calculate_magnetic_omega(x_state, mu_vec, T_fm, magnetic; kwargs...)
end

function _calculate_magnetic_net_density_vector(
    x_state,
    mu_vec,
    T_fm,
    magnetic::MagneticConfig=default_magnetic_config();
    xi::Real=0.0,
    p_num::Union{Nothing, Int}=nothing,
    t_num::Int=8,
    pz_max::Union{Nothing, Real}=nothing,
    n_max::Union{Nothing, Int}=nothing,
    cutoff_N::Union{Nothing, Int}=nothing,
)
    st = SVector{5, Float64}(Tuple(x_state))
    μ = SVector{3, Float64}(Tuple(mu_vec))
    comp = calculate_magnetic_omega_components(st, μ, T_fm, magnetic;
        xi=xi, p_num=p_num, t_num=t_num, pz_max=pz_max, n_max=n_max, cutoff_N=cutoff_N)
    masses = comp.masses
    n_max = comp.n_max
    conf = _effective_magnetic_config(magnetic;
        p_num=p_num, pz_max=pz_max, n_max=n_max, cutoff_N=cutoff_N)
    pz_max_val = conf.pz_max > 0 ? conf.pz_max : max(8 * Main.Constants_PNJL.Λ_inv_fm, 25.0)

    q = MVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for i in 1:3
        q[i] = density_flavor_landau(
            masses[i],
            μ[i],
            T_fm,
            st[4],
            st[5],
            QUARK_CHARGE_ABS[i],
            conf.eB_fm2;
            n_max=n_max,
            p_num=conf.p_num,
            pz_max=pz_max_val,
        )
    end
    return SVector{3}(q)
end

function calculate_magnetic_rho(x_state, mu_vec, T_fm, magnetic::MagneticConfig=default_magnetic_config(); kwargs...)
    return _calculate_magnetic_net_density_vector(x_state, mu_vec, T_fm, magnetic; kwargs...)
end

function calculate_magnetic_number_densities(x_state, mu_vec, T_fm, magnetic::MagneticConfig=default_magnetic_config(); kwargs...)
    q = _calculate_magnetic_net_density_vector(x_state, mu_vec, T_fm, magnetic; kwargs...)
    # `quark` is retained as the historical field name, but its formal
    # magnetic meaning is the net quark density q - qbar.
    return (quark=q, antiquark=nothing, net=q, baryon=sum(q) / 3.0 / ρ0)
end

"""magnetic_nmax_convergence_report(x_state, mu_vec, T_fm, magnetic; delta_n=6, rtol=3e-2)

给出 Landau 截断收敛报告：比较 `n_base` 与 `n_base + delta_n` 下的 Ω 相对差异。
"""
function magnetic_nmax_convergence_report(
    x_state,
    mu_vec,
    T_fm,
    magnetic::MagneticConfig=default_magnetic_config();
    delta_n::Int=6,
    rtol::Real=3e-2,
)
    delta_n >= 1 || throw(ArgumentError("delta_n must be >= 1, got $delta_n"))
    st = SVector{5, Float64}(Tuple(x_state))
    μ = SVector{3, Float64}(Tuple(mu_vec))

    comp0 = calculate_magnetic_omega_components(st, μ, T_fm, magnetic)
    n_base = comp0.n_max
    conf_probe = MagneticConfig(
        eB_fm2=magnetic.eB_fm2,
        n_max=n_base + delta_n,
        p_num=magnetic.p_num,
        pz_max=magnetic.pz_max,
        cutoff_N=magnetic.cutoff_N,
        imc=magnetic.imc,
    )
    comp1 = calculate_magnetic_omega_components(st, μ, T_fm, conf_probe)

    denom = max(1.0, abs(comp1.omega))
    rel_diff = abs(comp1.omega - comp0.omega) / denom
    return (
        converged=rel_diff <= rtol,
        rtol=Float64(rtol),
        rel_diff=rel_diff,
        n_base=n_base,
        n_probe=n_base + delta_n,
        omega_base=comp0.omega,
        omega_probe=comp1.omega,
    )
end

end # module MagneticThermodynamics
