"""
    MagneticThermodynamics

外磁场下 PNJL 热力学主接口。

核心公式对应：
- MFIR：零场三动量截断真空项 + Hurwitz-zeta 磁场真空修正
- legacy：完整 Landau 真空项（仅显式 diagnostic 路由）
- Ω_f^T (Landau 求和)
- E_{f,n} = sqrt(2n|q_f|eB + pz^2 + M_f^2)
- ρ_f = ∂P/∂μ_f（提供数值导数与低温近似两种路径）
- G(B) IMC 参数化
"""
module MagneticThermodynamics

using StaticArrays
using TOML

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
    parent = parentmodule(@__MODULE__)
    if isdefined(parent, :PNJLCore)
        const PNJLCore = getproperty(parent, :PNJLCore)
    else
        include(_PNJL_CORE_PATH)
    end
end

include("MagneticIntegrals.jl")
using .MagneticIntegrals:
    QUARK_CHARGE_ABS,
    MAGNETIC_EB_MIN_MEV2,
    MAGNETIC_EB_MIN_FM2,
    validate_magnetic_eB,
    resolve_nmax_from_cutoff,
    omega_magnetic_mfir,
    omega0_flavor_landau,
    omegat_flavor_landau,
    density_flavor_landau

export MagneticIMCParams, default_imc_params, coupling_GB
export MagneticConfig, default_magnetic_config
export MAGNETIC_ROUTE_DEFAULT, MAGNETIC_ROUTES
export MAGNETIC_EB_MIN_MEV2, MAGNETIC_EB_MIN_FM2, validate_magnetic_eB
export calculate_magnetic_omega_components, calculate_magnetic_omega
export calculate_magnetic_pressure, calculate_magnetic_rho
export calculate_magnetic_number_densities
export magnetic_nmax_convergence_report
export resolve_magnetic_nmax

const ρ0 = ρ0_inv_fm3
const _PNJL_PARAMS_REF = Ref{Any}(nothing)
const MAGNETIC_ROUTE_DEFAULT = :mfir
const MAGNETIC_ROUTES = (:mfir, :landau_legacy)

@inline function _pnjl_params()
    p = _PNJL_PARAMS_REF[]
    if p === nothing
        p = PNJLCore.pnjl_params()
        _PNJL_PARAMS_REF[] = p
    end
    return p
end

Base.@kwdef struct MagneticIMCParams
    a::Float64 = 0.0108805
    b::Float64 = -1.0133e-4
    c::Float64 = 0.02228
    d::Float64 = 1.84558e-4
    Λ_QCD_MeV::Float64 = 300.0
end

const _MAGNETIC_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "config", "models", "pnjl"))
const _MAGNETIC_CONFIG_CACHE = Dict{String, Dict{String, Any}}()
const _MAGNETIC_DEFAULTS = Dict{String, Any}(
    "magnetic" => Dict{String, Any}(
        "route" => "mfir",
        "n_max" => -1,
        "n_max_policy" => "thermal_tail",
        "thermal_tail_factor" => 30.0,
        "n_max_floor" => 3,
        "n_max_cap" => 10000,
        "p_num" => 128,
        "pz_max" => 40.0,
        "cutoff_N" => 10,
        "zeta_num" => 256,
        "eB_fm2" => 0.0,
        "imc" => Dict{String, Any}(
            "a" => 0.0108805,
            "b" => -1.0133e-4,
            "c" => 0.02228,
            "d" => 1.84558e-4,
            "Lambda_QCD_MeV" => 300.0,
        ),
    ),
)

function _magnetic_profile_config(profile::String)
    return get!(_MAGNETIC_CONFIG_CACHE, profile) do
        cfg = deepcopy(_MAGNETIC_DEFAULTS)
        path = joinpath(_MAGNETIC_CONFIG_DIR, string(profile, ".toml"))
        if isfile(path)
            parsed = TOML.parsefile(path)
            raw = get(parsed, "magnetic", Dict{String, Any}())
            raw isa Dict || throw(ArgumentError("[magnetic] must be a TOML table in $(path)"))
            merged = deepcopy(cfg["magnetic"])
            for (key, value) in raw
                if key == "imc" && value isa Dict
                    merged["imc"] = merge(merged["imc"], value)
                else
                    merged[key] = value
                end
            end
            cfg["magnetic"] = merged
        elseif profile != "magnetic_default"
            @warn "magnetic profile not found; using built-in defaults" profile path
        end
        cfg
    end
end

function default_imc_params(; profile::String=get(ENV, "PNJL_MAGNETIC_PROFILE", "magnetic_default"))
    raw = get(_magnetic_profile_config(profile), "magnetic", Dict{String, Any}())
    imc = get(raw, "imc", Dict{String, Any}())
    return MagneticIMCParams(
        a=Float64(get(imc, "a", 0.0108805)),
        b=Float64(get(imc, "b", -1.0133e-4)),
        c=Float64(get(imc, "c", 0.02228)),
        d=Float64(get(imc, "d", 1.84558e-4)),
        Λ_QCD_MeV=Float64(get(imc, "Lambda_QCD_MeV", 300.0)),
    )
end

function _magnetic_route(route::Symbol)
    route in MAGNETIC_ROUTES || throw(ArgumentError(
        "unsupported magnetic route=$(route); expected one of $(MAGNETIC_ROUTES)",
    ))
    return route
end

struct MagneticConfig
    eB_fm2::Float64
    n_max::Union{Int, Nothing}
    p_num::Int
    pz_max::Float64
    cutoff_N::Int
    imc::MagneticIMCParams
    route::Symbol
    zeta_num::Int
    params::PNJLCore.PNJLParams
    n_max_policy::Symbol
    thermal_tail_factor::Float64
    n_max_floor::Int
    n_max_cap::Int

    function MagneticConfig(
        eB_fm2::Real,
        n_max::Union{Nothing, Int}=nothing,
        p_num::Int=128,
        pz_max::Real=0.0,
        cutoff_N::Int=10,
        imc::MagneticIMCParams=default_imc_params(),
        route::Symbol=MAGNETIC_ROUTE_DEFAULT,
        zeta_num::Int=256,
        params::PNJLCore.PNJLParams=_pnjl_params(),
        n_max_policy::Symbol=:thermal_tail,
        thermal_tail_factor::Real=30.0,
        n_max_floor::Int=3,
        n_max_cap::Int=10000,
    )
        eB_value = validate_magnetic_eB(eB_fm2)
        n_max === nothing || n_max >= 0 || throw(ArgumentError("magnetic n_max must be >= 0, got $(n_max)"))
        p_num >= 4 || throw(ArgumentError("magnetic p_num must be >= 4, got $(p_num)"))
        pz_max >= 0 || throw(ArgumentError("magnetic pz_max must be nonnegative, got $(pz_max)"))
        cutoff_N >= 1 || throw(ArgumentError("magnetic cutoff_N must be >= 1, got $(cutoff_N)"))
        zeta_num >= 8 || throw(ArgumentError("magnetic zeta_num must be >= 8, got $(zeta_num)"))
        n_max_policy in (:thermal_tail, :vacuum_cutoff) || throw(ArgumentError(
            "unsupported magnetic n_max_policy=$(n_max_policy); expected :thermal_tail or :vacuum_cutoff",
        ))
        tail_value = Float64(thermal_tail_factor)
        isfinite(tail_value) && tail_value > 0.0 || throw(ArgumentError(
            "magnetic thermal_tail_factor must be finite and > 0, got $(thermal_tail_factor)",
        ))
        n_max_floor >= 0 || throw(ArgumentError("magnetic n_max_floor must be >= 0, got $(n_max_floor)"))
        n_max_cap >= n_max_floor || throw(ArgumentError(
            "magnetic n_max_cap must be >= n_max_floor, got $(n_max_cap) < $(n_max_floor)",
        ))
        new(
            eB_value, n_max, p_num, Float64(pz_max), cutoff_N, imc,
            _magnetic_route(route), zeta_num, params, n_max_policy,
            tail_value, n_max_floor, n_max_cap,
        )
    end
end

function MagneticConfig(
    ;
    eB_fm2::Real=0.0,
    n_max::Union{Nothing, Int}=nothing,
    p_num::Int=128,
    pz_max::Real=0.0,
    cutoff_N::Int=10,
    imc::MagneticIMCParams=default_imc_params(),
    route::Symbol=MAGNETIC_ROUTE_DEFAULT,
    zeta_num::Int=256,
    n_max_policy::Symbol=:thermal_tail,
    thermal_tail_factor::Real=30.0,
    n_max_floor::Int=3,
    n_max_cap::Int=10000,
    params::PNJLCore.PNJLParams=_pnjl_params(),
)
    return MagneticConfig(
        eB_fm2, n_max, p_num, pz_max, cutoff_N, imc, route, zeta_num, params,
        n_max_policy, thermal_tail_factor, n_max_floor, n_max_cap,
    )
end

function default_magnetic_config(
    ;
    eB_fm2::Real=0.0,
    profile::String=get(ENV, "PNJL_MAGNETIC_PROFILE", "magnetic_default"),
    n_max::Union{Nothing, Int}=nothing,
    p_num::Union{Nothing, Int}=nothing,
    pz_max::Union{Nothing, Real}=nothing,
    cutoff_N::Union{Nothing, Int}=nothing,
    route::Union{Nothing, Symbol}=nothing,
    zeta_num::Union{Nothing, Int}=nothing,
    imc::Union{Nothing, MagneticIMCParams}=nothing,
    n_max_policy::Union{Nothing, Symbol}=nothing,
    thermal_tail_factor::Union{Nothing, Real}=nothing,
    n_max_floor::Union{Nothing, Int}=nothing,
    n_max_cap::Union{Nothing, Int}=nothing,
    params::PNJLCore.PNJLParams=_pnjl_params(),
)
    raw = get(_magnetic_profile_config(profile), "magnetic", Dict{String, Any}())
    configured_n = Int(get(raw, "n_max", -1))
    configured_n >= -1 || throw(ArgumentError("magnetic config n_max must be >= -1, got $(configured_n)"))
    n_value = n_max === nothing ? (configured_n < 0 ? nothing : configured_n) : n_max
    return MagneticConfig(
        eB_fm2,
        n_value,
        p_num === nothing ? Int(get(raw, "p_num", 128)) : p_num,
        pz_max === nothing ? Float64(get(raw, "pz_max", 40.0)) : pz_max,
        cutoff_N === nothing ? Int(get(raw, "cutoff_N", 10)) : cutoff_N,
        imc === nothing ? default_imc_params(profile=profile) : imc,
        route === nothing ? Symbol(lowercase(String(get(raw, "route", "mfir")))) : route,
        zeta_num === nothing ? Int(get(raw, "zeta_num", 256)) : zeta_num,
        params,
        n_max_policy === nothing ? Symbol(lowercase(String(get(raw, "n_max_policy", "thermal_tail")))) : n_max_policy,
        thermal_tail_factor === nothing ? Float64(get(raw, "thermal_tail_factor", 30.0)) : thermal_tail_factor,
        n_max_floor === nothing ? Int(get(raw, "n_max_floor", 3)) : n_max_floor,
        n_max_cap === nothing ? Int(get(raw, "n_max_cap", 10000)) : n_max_cap,
    )
end

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

@inline function _calculate_mass_vec_with_GB(
    φ::SVector{3, T},
    G_B::Real;
    params::PNJLCore.PNJLParams=_pnjl_params(),
) where {T}
    φ_u, φ_d, φ_s = φ
    Gv = convert(T, G_B)
    Kv = convert(T, params.K_fm5)
    m_ud0 = convert(T, params.m_ud0_inv_fm)
    m_s0 = convert(T, params.m_s0_inv_fm)
    return SVector{3, T}(
        m_ud0 - 4 * Gv * φ_u + 2 * Kv * φ_d * φ_s,
        m_ud0 - 4 * Gv * φ_d + 2 * Kv * φ_u * φ_s,
        m_s0 - 4 * Gv * φ_s + 2 * Kv * φ_u * φ_d,
    )
end

@inline function _chiral_with_GB(
    φ::SVector{3, T},
    G_B::Real;
    params::PNJLCore.PNJLParams=_pnjl_params(),
) where {T}
    Gv = convert(T, G_B)
    Kv = convert(T, params.K_fm5)
    return 2 * Gv * sum(φ .^ 2) - 4 * Kv * prod(φ)
end

@inline function _validate_nmax_budget(candidate::Int, conf::MagneticConfig, T_fm, eB)
    candidate <= conf.n_max_cap || throw(ArgumentError(
        "magnetic n_max budget exceeded: candidate=$(candidate), cap=$(conf.n_max_cap), " *
        "T_fm=$(T_fm), eB_fm2=$(eB); increase n_max_cap or use a stronger field",
    ))
    return max(candidate, conf.n_max_floor)
end

"""resolve_magnetic_nmax(T_fm, mu_vec, magnetic; masses=nothing) -> Int

Resolve the finite Landau-layer budget used by the magnetic thermal term.
An explicit `magnetic.n_max` wins. Otherwise `:thermal_tail` uses a shared
temperature/chemical-potential tail estimate and applies the configured floor
and cap; `:vacuum_cutoff` is an explicit legacy policy and requires masses.
"""
function resolve_magnetic_nmax(
    T_fm::Real,
    mu_vec::AbstractVector{<:Real},
    magnetic::MagneticConfig;
    masses::Union{Nothing, SVector{3, <:Real}}=nothing,
)
    T_value = Float64(T_fm)
    isfinite(T_value) && T_value > 0.0 || throw(ArgumentError(
        "magnetic n_max resolution requires finite T_fm > 0, got $(T_fm)",
    ))
    all(isfinite, mu_vec) || throw(ArgumentError("magnetic n_max resolution requires finite mu_vec"))
    eB_value = validate_magnetic_eB(magnetic.eB_fm2)
    magnetic.n_max !== nothing && return magnetic.n_max::Int

    candidate = if magnetic.n_max_policy === :thermal_tail
        μ_abs = maximum(abs, mu_vec)
        E_tail = μ_abs + magnetic.thermal_tail_factor * T_value
        q_min = minimum(QUARK_CHARGE_ABS)
        ceil(Int, max(E_tail, 0.0)^2 / (2 * q_min * eB_value))
    else
        masses === nothing && throw(ArgumentError(
            "magnetic vacuum_cutoff n_max policy requires masses or explicit n_max",
        ))
        nvals = ntuple(i -> resolve_nmax_from_cutoff(
            masses[i], mu_vec[i], QUARK_CHARGE_ABS[i], eB_value,
        ), 3)
        maximum(nvals)
    end
    return _validate_nmax_budget(candidate, magnetic, T_fm, eB_value)
end

@inline function _resolve_nmax(
    masses::SVector{3, <:Real},
    mu_vec::SVector{3, <:Real},
    T_fm::Real,
    eB::Real,
    conf::MagneticConfig,
)
    conf.n_max !== nothing && return conf.n_max::Int
    conf.eB_fm2 == Float64(eB) || throw(ArgumentError("magnetic n_max resolution received mismatched eB"))
    return resolve_magnetic_nmax(T_fm, mu_vec, conf; masses=masses)
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
        route=magnetic.route,
        zeta_num=magnetic.zeta_num,
        n_max_policy=magnetic.n_max_policy,
        thermal_tail_factor=magnetic.thermal_tail_factor,
        n_max_floor=magnetic.n_max_floor,
        n_max_cap=magnetic.n_max_cap,
        params=magnetic.params,
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

    params = conf.params
    G_B = coupling_GB(conf.eB_fm2; G0=params.G_fm2, imc=conf.imc)
    masses = _calculate_mass_vec_with_GB(φ, G_B; params=params)
    n_max = _resolve_nmax(masses, mu_vec, T_fm, conf.eB_fm2, conf)

    pz_max_val = conf.pz_max > 0 ? conf.pz_max : max(8 * Main.Constants_PNJL.Λ_inv_fm, 25.0)

    vac = 0.0
    therm = 0.0
    @inbounds for i in 1:3
        if conf.route === :mfir
            Λ = convert(typeof(masses[i]), params.Λ_inv_fm)
            vac += -2 * params.N_color * PNJLCore.vacuum_integral_with_cutoff(masses[i], Λ)
            vac += omega_magnetic_mfir(
                masses[i],
                QUARK_CHARGE_ABS[i],
                conf.eB_fm2;
                zeta_num=conf.zeta_num,
                N_c=params.N_color,
            )
        else
            vac += omega0_flavor_landau(
                masses[i],
                QUARK_CHARGE_ABS[i],
                conf.eB_fm2;
                n_max=n_max,
                p_num=conf.p_num,
                pz_max=pz_max_val,
                cutoff_N=conf.cutoff_N,
            )
        end
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

    chi = _chiral_with_GB(φ, G_B; params=params)
    poly = PNJLCore.polyakov_potential(params, Φ, Φbar, T_fm)
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
        route=magnetic.route,
        zeta_num=magnetic.zeta_num,
        n_max_policy=magnetic.n_max_policy,
        thermal_tail_factor=magnetic.thermal_tail_factor,
        n_max_floor=magnetic.n_max_floor,
        n_max_cap=magnetic.n_max_cap,
        params=magnetic.params,
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
