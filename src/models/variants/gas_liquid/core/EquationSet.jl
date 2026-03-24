module GasLiquidEquationSet

using NLsolve
using StaticArrays

const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "config", "ConfigLoader.jl"))
if !isdefined(@__MODULE__, :ConfigLoader)
    include(_CONFIG_LOADER_PATH)
end
using .ConfigLoader: load_config

export GasLiquidCoreParams
export GasLiquidState
export state_vector
export solve_equilibrium
export field_residuals
export mu_baryon
export effective_masses
export density_bundle

const _DEFAULT_PROFILE = "default"
const _GAS_CONFIG_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "..", "config", "models", "gas_liquid"))
const _DEFAULT_GAS_CONFIG = Dict{String, Any}(
    "model" => Dict{String, Any}(
        "m_nucleon_MeV" => 939.0,
        "m_sigma_MeV" => 500.0,
        "m_omega_MeV" => 783.0,
        "m_rho_MeV" => 770.0,
        "m_delta_MeV" => 983.0,
        "f_sigma" => 10.329,
        "f_omega" => 5.423,
        "f_rho" => 3.15,
        "g_delta" => 0.0,
        "b" => 0.00692,
        "c" => -0.0048,
        "hbarc_MeV_fm" => 197.3269804,
    ),
)

@inline _mev_to_fminv(x_MeV::Real, hbarc_MeV_fm::Real) = float(x_MeV) / float(hbarc_MeV_fm)

@inline function _load_gasliquid_config(; profile::String=get(ENV, "GAS_LIQUID_PARAM_PROFILE", _DEFAULT_PROFILE))
    return load_config(_GAS_CONFIG_DIR, _DEFAULT_GAS_CONFIG; profile=profile).config
end

Base.@kwdef struct GasLiquidCoreParams
    hbarc_MeV_fm::Float64 = 197.3269804

    m_nucleon_inv_fm::Float64 = 939.0 / 197.3269804
    m_sigma_inv_fm::Float64 = 500.0 / 197.3269804
    m_omega_inv_fm::Float64 = 783.0 / 197.3269804
    m_rho_inv_fm::Float64 = 770.0 / 197.3269804
    m_delta_inv_fm::Float64 = 983.0 / 197.3269804

    f_sigma::Float64 = 10.329
    f_omega::Float64 = 5.423
    f_rho::Float64 = 3.15

    g_sigma::Float64 = f_sigma * m_sigma_inv_fm
    g_omega::Float64 = f_omega * m_omega_inv_fm
    g_rho::Float64 = f_rho * m_rho_inv_fm
    g_delta::Float64 = 0.0

    b::Float64 = 0.00692
    c::Float64 = -0.0048
end

function GasLiquidCoreParams(; profile::Union{Nothing, String}=nothing, kwargs...)
    profile_eff = something(profile, get(ENV, "GAS_LIQUID_PARAM_PROFILE", _DEFAULT_PROFILE))
    cfg_model = get(_load_gasliquid_config(profile=profile_eff), "model", Dict{String, Any}())

    hbarc = Float64(get(cfg_model, "hbarc_MeV_fm", 197.3269804))
    mN = _mev_to_fminv(Float64(get(cfg_model, "m_nucleon_MeV", 939.0)), hbarc)
    mσ = _mev_to_fminv(Float64(get(cfg_model, "m_sigma_MeV", 500.0)), hbarc)
    mω = _mev_to_fminv(Float64(get(cfg_model, "m_omega_MeV", 783.0)), hbarc)
    mρ = _mev_to_fminv(Float64(get(cfg_model, "m_rho_MeV", 770.0)), hbarc)
    mδ = _mev_to_fminv(Float64(get(cfg_model, "m_delta_MeV", 983.0)), hbarc)

    fσ = Float64(get(cfg_model, "f_sigma", 10.329))
    fω = Float64(get(cfg_model, "f_omega", 5.423))
    fρ = Float64(get(cfg_model, "f_rho", 3.15))

    base = (
        hbarc_MeV_fm=hbarc,
        m_nucleon_inv_fm=mN,
        m_sigma_inv_fm=mσ,
        m_omega_inv_fm=mω,
        m_rho_inv_fm=mρ,
        m_delta_inv_fm=mδ,
        f_sigma=fσ,
        f_omega=fω,
        f_rho=fρ,
        g_sigma=fσ * mσ,
        g_omega=fω * mω,
        g_rho=fρ * mρ,
        g_delta=Float64(get(cfg_model, "g_delta", 0.0)),
        b=Float64(get(cfg_model, "b", 0.00692)),
        c=Float64(get(cfg_model, "c", -0.0048)),
    )

    merged = merge(base, kwargs)
    return GasLiquidCoreParams(
        Float64(merged[:hbarc_MeV_fm]),
        Float64(merged[:m_nucleon_inv_fm]),
        Float64(merged[:m_sigma_inv_fm]),
        Float64(merged[:m_omega_inv_fm]),
        Float64(merged[:m_rho_inv_fm]),
        Float64(merged[:m_delta_inv_fm]),
        Float64(merged[:f_sigma]),
        Float64(merged[:f_omega]),
        Float64(merged[:f_rho]),
        Float64(merged[:g_sigma]),
        Float64(merged[:g_omega]),
        Float64(merged[:g_rho]),
        Float64(merged[:g_delta]),
        Float64(merged[:b]),
        Float64(merged[:c]),
    )
end

struct GasLiquidState{T}
    sigma::T
    delta::T
    mu_p::T
    mu_n::T
end

@inline state_vector(st::GasLiquidState{T}) where {T} = SVector{4, T}(st.sigma, st.delta, st.mu_p, st.mu_n)

@inline function mu_baryon(mu_vec)
    if mu_vec isa Real
        return float(mu_vec)
    end
    length(mu_vec) == 3 || throw(ArgumentError("mu_vec must be Real or length-3 vector"))
    return float((mu_vec[1] + mu_vec[2] + mu_vec[3]) / 3)
end

@inline function effective_masses(st::GasLiquidState{T}, p::GasLiquidCoreParams) where {T}
    mp = p.m_nucleon_inv_fm - p.g_sigma * st.sigma - p.g_delta * st.delta
    mn = p.m_nucleon_inv_fm - p.g_sigma * st.sigma + p.g_delta * st.delta
    return (mp=mp, mn=mn)
end

@inline function _fermi_pair(E::TE, μ::TM, Tval::TT) where {TE, TM, TT}
    Tp = promote_type(TE, TM, TT)
    Eτ = convert(Tp, E)
    μτ = convert(Tp, μ)
    Tτ = convert(Tp, Tval)
    f = inv(one(Tp) + exp((Eτ - μτ) / Tτ))
    fb = inv(one(Tp) + exp((Eτ + μτ) / Tτ))
    return f, fb
end

function density_bundle(st::GasLiquidState{Tst}, T::Real, p::GasLiquidCoreParams; p_num::Int=96) where {Tst}
    Tm = promote_type(Tst, Float64)
    Tv = T
    pmax = max(convert(Tm, 7.5), abs(convert(Tm, st.mu_p)), abs(convert(Tm, st.mu_n)), convert(Tm, p.m_nucleon_inv_fm + 4.0))
    dp = pmax / convert(Tm, p_num)

    m = effective_masses(st, p)

    rho_p = zero(Tm)
    rho_n = zero(Tm)
    rho_s_p = zero(Tm)
    rho_s_n = zero(Tm)
    kin = zero(Tm)
    logsum = zero(Tm)

    @inbounds for k in 0:p_num
        mom = convert(Tm, k) * dp
        w = (k == 0 || k == p_num) ? convert(Tm, 0.5) : one(Tm)
        measure = mom^2 * dp / convert(Tm, 2 * pi^2)

        Ep = sqrt(mom^2 + m.mp^2)
        En = sqrt(mom^2 + m.mn^2)

        fp, fbp = _fermi_pair(Ep, st.mu_p, Tv)
        fn, fbn = _fermi_pair(En, st.mu_n, Tv)

        rho_p += w * convert(Tm, 2.0) * measure * (fp - fbp)
        rho_n += w * convert(Tm, 2.0) * measure * (fn - fbn)
        rho_s_p += w * convert(Tm, 2.0) * measure * (m.mp / Ep) * (fp + fbp)
        rho_s_n += w * convert(Tm, 2.0) * measure * (m.mn / En) * (fn + fbn)
        kin += w * convert(Tm, 2.0) * measure * (Ep * (fp + fbp) + En * (fn + fbn))

        logsum += w * convert(Tm, 2.0) * measure * (
            log1p(exp(-(Ep - st.mu_p) / Tv)) +
            log1p(exp(-(Ep + st.mu_p) / Tv)) +
            log1p(exp(-(En - st.mu_n) / Tv)) +
            log1p(exp(-(En + st.mu_n) / Tv))
        )
    end

    return (
        masses=m,
        rho_p=rho_p,
        rho_n=rho_n,
        rho_B=rho_p + rho_n,
        rho_s_p=rho_s_p,
        rho_s_n=rho_s_n,
        kinetic_energy=kin,
        logsum=logsum,
    )
end

function field_residuals(st::GasLiquidState{Tst}, T::Real, p::GasLiquidCoreParams; p_num::Int=96, mu_target::Union{Nothing, Real}=nothing) where {Tst}
    dens = density_bundle(st, T, p; p_num=p_num)
    Tm = promote_type(Tst, typeof(T), typeof(dens.rho_B))
    μtar = isnothing(mu_target) ? convert(Tm, 0.5) * (st.mu_p + st.mu_n) : convert(Tm, mu_target)

    σ = st.sigma
    δ = st.delta

    r_sigma = p.m_sigma_inv_fm^2 * σ + p.b * p.m_nucleon_inv_fm * (p.g_sigma * σ)^2 + p.c * (p.g_sigma * σ)^3 - p.g_sigma * (dens.rho_s_p + dens.rho_s_n)
    r_delta = p.m_delta_inv_fm^2 * δ - p.g_delta * (dens.rho_s_p - dens.rho_s_n)
    r_mu_p = st.mu_p - μtar
    r_mu_n = st.mu_n - μtar

    return SVector{4, Tm}(Tm(r_sigma), Tm(r_delta), Tm(r_mu_p), Tm(r_mu_n))
end

function solve_equilibrium(
    T::Real,
    mu_vec;
    params::GasLiquidCoreParams=GasLiquidCoreParams(),
    initial_guess::Union{Nothing, AbstractVector}=nothing,
    p_num::Int=96,
)
    muB = mu_baryon(mu_vec)
    x0 = initial_guess === nothing ? [0.02, 0.0, muB, muB] : collect(float.(initial_guess))
    length(x0) == 4 || throw(ArgumentError("initial_guess must have length 4"))

    function residual!(F, x)
        st = GasLiquidState(x[1], x[2], x[3], x[4])
        r = field_residuals(st, T, params; p_num=p_num, mu_target=muB)
        F[1] = r[1]
        F[2] = r[2]
        F[3] = r[3]
        F[4] = r[4]
        return nothing
    end

    res = nlsolve(residual!, x0; autodiff=:forward, method=:trust_region, xtol=1e-10, ftol=1e-10)
    if !(res.f_converged && isfinite(res.residual_norm))
        error("gas-liquid equilibrium solve failed: f_converged=$(res.f_converged), residual_norm=$(res.residual_norm)")
    end

    x = res.zero
    return GasLiquidState(float(x[1]), float(x[2]), float(x[3]), float(x[4]))
end

end # module
