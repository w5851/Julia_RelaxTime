"""Finite-temperature RMF/RMFT equation set for the gas-liquid model.

The core uses natural units (fm and fm^-1). Coupling ratios are defined as
`f_i = g_i^2 / m_i^2` in fm^2. The state stores the field contributions
`S = g_sigma*sigma`, `D = g_delta*delta` and the effective nucleon chemical
potentials. Vector fields are generated from the densities.
"""
module GasLiquidEquationSet

using NLsolve
using StaticArrays

const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "config", "ConfigLoader.jl"))
if !isdefined(@__MODULE__, :ConfigLoader)
    include(_CONFIG_LOADER_PATH)
end
using .ConfigLoader: load_config

const _GAUSS_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "..", "integration", "GaussLegendre.jl"))
if !isdefined(@__MODULE__, :GaussLegendre)
    include(_GAUSS_PATH)
end
using .GaussLegendre: gauleg

export GasLiquidCoreParams, GasLiquidState, RMFSolverResult
export state_vector, coupling_ratio, parameter_summary
export solve_equilibrium, solve_rmf
export field_residuals, effective_masses, effective_chemical_potentials
export density_bundle, field_contributions, physical_mu_pair
export mu_baryon

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
        "b" => 0.00692,
        "c" => -0.0048,
        "hbarc_MeV_fm" => 197.3269804,
    ),
)

@inline _mev_to_fminv(x_MeV::Real, hbarc_MeV_fm::Real) = float(x_MeV) / float(hbarc_MeV_fm)

@inline function _load_gasliquid_config(; profile::String=get(ENV, "GAS_LIQUID_PARAM_PROFILE", _DEFAULT_PROFILE))
    return load_config(_GAS_CONFIG_DIR, _DEFAULT_GAS_CONFIG; profile=profile).config
end

"""RMF parameters in natural units.

`f_sigma`, `f_omega`, `f_rho` and `f_delta` are all fm^2 and obey
`f_i = g_i^2/m_i^2`. `g_delta` is retained as a derived compatibility field;
new code should use `f_delta`.
"""
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
    f_delta::Float64 = 0.0
    g_sigma::Float64 = sqrt(max(f_sigma, 0.0)) * m_sigma_inv_fm
    g_omega::Float64 = sqrt(max(f_omega, 0.0)) * m_omega_inv_fm
    g_rho::Float64 = sqrt(max(f_rho, 0.0)) * m_rho_inv_fm
    g_delta::Float64 = sqrt(max(f_delta, 0.0)) * m_delta_inv_fm
    b::Float64 = 0.00692
    c::Float64 = -0.0048
end

@inline coupling_ratio(g::Real, m_inv_fm::Real) = (float(g) / float(m_inv_fm))^2

function _validate_params(p::GasLiquidCoreParams)
    p.hbarc_MeV_fm > 0 || throw(ArgumentError("hbarc_MeV_fm must be positive"))
    all(>(0.0), (p.m_nucleon_inv_fm, p.m_sigma_inv_fm, p.m_omega_inv_fm, p.m_rho_inv_fm, p.m_delta_inv_fm)) ||
        throw(ArgumentError("all meson and nucleon masses must be positive"))
    all(>=(0.0), (p.f_sigma, p.f_omega, p.f_rho, p.f_delta)) ||
        throw(ArgumentError("coupling ratios f_i must be non-negative"))
    return p
end

function GasLiquidCoreParams(; profile::Union{Nothing, String}=nothing, kwargs...)
    profile_eff = something(profile, get(ENV, "GAS_LIQUID_PARAM_PROFILE", _DEFAULT_PROFILE))
    cfg_model = get(_load_gasliquid_config(profile=profile_eff), "model", Dict{String, Any}())

    hbarc = Float64(get(cfg_model, "hbarc_MeV_fm", 197.3269804))
    mN = _mev_to_fminv(Float64(get(cfg_model, "m_nucleon_MeV", 939.0)), hbarc)
    msigma = _mev_to_fminv(Float64(get(cfg_model, "m_sigma_MeV", 500.0)), hbarc)
    momega = _mev_to_fminv(Float64(get(cfg_model, "m_omega_MeV", 783.0)), hbarc)
    mrho = _mev_to_fminv(Float64(get(cfg_model, "m_rho_MeV", 770.0)), hbarc)
    mdelta = _mev_to_fminv(Float64(get(cfg_model, "m_delta_MeV", 983.0)), hbarc)

    fsigma = Float64(get(cfg_model, "f_sigma", 10.329))
    fomega = Float64(get(cfg_model, "f_omega", 5.423))
    frho = Float64(get(cfg_model, "f_rho", 3.15))
    # Legacy profiles may specify a bare g_delta. Convert it only at the
    # configuration boundary; the core never uses the bare coupling.
    fdelta = if haskey(cfg_model, "f_delta")
        Float64(cfg_model["f_delta"])
    else
        gd = Float64(get(cfg_model, "g_delta", 0.0))
        (gd / mdelta)^2
    end

    base = (
        hbarc_MeV_fm=hbarc,
        m_nucleon_inv_fm=mN,
        m_sigma_inv_fm=msigma,
        m_omega_inv_fm=momega,
        m_rho_inv_fm=mrho,
        m_delta_inv_fm=mdelta,
        f_sigma=fsigma,
        f_omega=fomega,
        f_rho=frho,
        f_delta=fdelta,
        g_sigma=sqrt(max(fsigma, 0.0)) * msigma,
        g_omega=sqrt(max(fomega, 0.0)) * momega,
        g_rho=sqrt(max(frho, 0.0)) * mrho,
        g_delta=sqrt(max(fdelta, 0.0)) * mdelta,
        b=Float64(get(cfg_model, "b", 0.00692)),
        c=Float64(get(cfg_model, "c", -0.0048)),
    )

    merged = merge(base, kwargs)
    # Accept a bare g_delta only for legacy callers. If both are supplied,
    # f_delta is authoritative.
    fdelta_eff = haskey(kwargs, :f_delta) ? Float64(merged[:f_delta]) :
        (haskey(kwargs, :g_delta) ? coupling_ratio(merged[:g_delta], merged[:m_delta_inv_fm]) : Float64(merged[:f_delta]))
    gdelta_eff = sqrt(max(fdelta_eff, 0.0)) * Float64(merged[:m_delta_inv_fm])
    gsigma_eff = haskey(kwargs, :g_sigma) ? Float64(merged[:g_sigma]) : sqrt(max(Float64(merged[:f_sigma]), 0.0)) * Float64(merged[:m_sigma_inv_fm])
    gomega_eff = haskey(kwargs, :g_omega) ? Float64(merged[:g_omega]) : sqrt(max(Float64(merged[:f_omega]), 0.0)) * Float64(merged[:m_omega_inv_fm])
    grho_eff = haskey(kwargs, :g_rho) ? Float64(merged[:g_rho]) : sqrt(max(Float64(merged[:f_rho]), 0.0)) * Float64(merged[:m_rho_inv_fm])
    p = GasLiquidCoreParams(
        Float64(merged[:hbarc_MeV_fm]),
        Float64(merged[:m_nucleon_inv_fm]),
        Float64(merged[:m_sigma_inv_fm]),
        Float64(merged[:m_omega_inv_fm]),
        Float64(merged[:m_rho_inv_fm]),
        Float64(merged[:m_delta_inv_fm]),
        Float64(merged[:f_sigma]),
        Float64(merged[:f_omega]),
        Float64(merged[:f_rho]),
        fdelta_eff,
        gsigma_eff,
        gomega_eff,
        grho_eff,
        gdelta_eff,
        Float64(merged[:b]),
        Float64(merged[:c]),
    )
    return _validate_params(p)
end

"""Compact, serializable parameter semantics used by manifests."""
@inline function parameter_summary(p::GasLiquidCoreParams)
    return (
        coupling_definition="f_i=g_i^2/m_i^2",
        units="natural_units_fm",
        hbarc_MeV_fm=p.hbarc_MeV_fm,
        m_nucleon_inv_fm=p.m_nucleon_inv_fm,
        m_sigma_inv_fm=p.m_sigma_inv_fm,
        m_omega_inv_fm=p.m_omega_inv_fm,
        m_rho_inv_fm=p.m_rho_inv_fm,
        m_delta_inv_fm=p.m_delta_inv_fm,
        f_sigma_fm2=p.f_sigma,
        f_omega_fm2=p.f_omega,
        f_rho_fm2=p.f_rho,
        f_delta_fm2=p.f_delta,
        b=p.b,
        c=p.c,
    )
end

"""Core state: `(S,D,mu_p_star,mu_n_star)` in fm^-1."""
struct GasLiquidState{T}
    S::T
    D::T
    mu_p::T
    mu_n::T
end

function Base.getproperty(st::GasLiquidState, name::Symbol)
    name === :sigma && return getfield(st, :S)
    name === :delta && return getfield(st, :D)
    name === :mu_p_star && return getfield(st, :mu_p)
    name === :mu_n_star && return getfield(st, :mu_n)
    return getfield(st, name)
end

Base.propertynames(::GasLiquidState, private::Bool=false) = (:S, :D, :mu_p, :mu_n, :sigma, :delta, :mu_p_star, :mu_n_star)

@inline state_vector(st::GasLiquidState{T}) where {T} = SVector{4, T}(st.S, st.D, st.mu_p, st.mu_n)

@inline function mu_baryon(mu_vec)
    mu_vec isa Real && return float(mu_vec)
    length(mu_vec) == 3 || throw(ArgumentError("mu_vec must be Real or length-3 vector"))
    return float((mu_vec[1] + mu_vec[2] + mu_vec[3]) / 3)
end

@inline function _resolve_mu_pair(mu_vec; mu_B=nothing, mu_3=nothing)
    if mu_B !== nothing
        mb = float(mu_B)
        m3 = mu_3 === nothing ? 0.0 : float(mu_3)
        return mb + m3, mb - m3, mb, m3
    elseif mu_vec isa Real
        mb = float(mu_vec)
        m3 = mu_3 === nothing ? 0.0 : float(mu_3)
        return mb + m3, mb - m3, mb, m3
    elseif length(mu_vec) == 2
        mb = float(mu_vec[1])
        m3 = float(mu_vec[2])
        return mb + m3, mb - m3, mb, m3
    else
        mb = mu_baryon(mu_vec)
        m3 = mu_3 === nothing ? 0.0 : float(mu_3)
        return mb + m3, mb - m3, mb, m3
    end
end

@inline function effective_masses(st::GasLiquidState, p::GasLiquidCoreParams)
    return (mp=p.m_nucleon_inv_fm - st.S - st.D, mn=p.m_nucleon_inv_fm - st.S + st.D)
end

@inline _fermi(E, mu, T) = T <= 0.0 ? (E < mu ? one(E) : zero(E)) : begin
    z = (E - mu) / T
    z > 40.0 ? 0.0 : (z < -40.0 ? 1.0 : inv(1.0 + exp(z)))
end

@inline _antifermi(E, mu, T) = T <= 0.0 ? zero(E) : begin
    z = (E + mu) / T
    z > 40.0 ? 0.0 : (z < -40.0 ? 1.0 : inv(1.0 + exp(z)))
end

@inline _log1pexp(x) = x > 40.0 ? x + log1p(exp(-x)) : log1p(exp(x))

function _density_integrals(st::GasLiquidState, T::Real, p::GasLiquidCoreParams; p_num::Int=96, p_max_inv_fm::Real=10.0)
    # Preserve ForwardDiff dual values for the generic Models thermodynamics
    # facade; only normalize integer/real inputs, do not force Float64 here.
    Tval = float(T)
    Tval >= 0.0 || throw(ArgumentError("temperature must be non-negative in fm^-1"))
    p_num > 0 || throw(ArgumentError("p_num must be positive"))
    pmax = Float64(p_max_inv_fm)
    pmax > 0.0 || throw(ArgumentError("p_max_inv_fm must be positive"))
    nodes, weights = gauleg(0.0, pmax, p_num)
    masses = effective_masses(st, p)

    rho_p = 0.0; rho_n = 0.0
    rho_s_p = 0.0; rho_s_n = 0.0
    kinetic = 0.0; log_integral = 0.0; entropy = 0.0
    @inbounds for i in eachindex(nodes)
        k = nodes[i]
        measure = k^2 * weights[i] / (2.0 * pi^2)
        Ep = sqrt(k^2 + masses.mp^2)
        En = sqrt(k^2 + masses.mn^2)
        fp = _fermi(Ep, st.mu_p, Tval); fbp = _antifermi(Ep, st.mu_p, Tval)
        fn = _fermi(En, st.mu_n, Tval); fbn = _antifermi(En, st.mu_n, Tval)

        rho_p += 2.0 * measure * (fp - fbp)
        rho_n += 2.0 * measure * (fn - fbn)
        rho_s_p += 2.0 * measure * (masses.mp / Ep) * (fp + fbp)
        rho_s_n += 2.0 * measure * (masses.mn / En) * (fn + fbn)
        kinetic += 2.0 * measure * (Ep * (fp + fbp) + En * (fn + fbn))
        if Tval > 0.0
            xp = (Ep - st.mu_p) / Tval; xbp = (Ep + st.mu_p) / Tval
            xn = (En - st.mu_n) / Tval; xbn = (En + st.mu_n) / Tval
            log_integral += measure * (_log1pexp(-xp) + _log1pexp(-xbp) + _log1pexp(-xn) + _log1pexp(-xbn))
            # The x*n terms are the beta(E +/- mu) n pieces. The log terms
            # are required by -dOmega/dT and keep the entropy non-negative.
            entropy += 2.0 * measure * (_log1pexp(-xp) + xp * fp + _log1pexp(-xbp) + xbp * fbp +
                                         _log1pexp(-xn) + xn * fn + _log1pexp(-xbn) + xbn * fbn)
        end
    end

    return (
        masses=masses,
        rho_p=rho_p,
        rho_n=rho_n,
        rho_B=rho_p + rho_n,
        rho_3=rho_p - rho_n,
        rho_s_p=rho_s_p,
        rho_s_n=rho_s_n,
        rho_s=rho_s_p + rho_s_n,
        rho_s3=rho_s_p - rho_s_n,
        kinetic_energy=kinetic,
        thermal_log_integral=log_integral,
        entropy=entropy,
        p_num=p_num,
        p_max_inv_fm=pmax,
        mu_p_star=st.mu_p,
        mu_n_star=st.mu_n,
    )
end

"""Densities for a state whose chemical potentials are effective mu* values."""
@inline density_bundle(st::GasLiquidState, T::Real, p::GasLiquidCoreParams; p_num::Int=96, p_max_inv_fm::Real=10.0) =
    _density_integrals(st, T, p; p_num=p_num, p_max_inv_fm=p_max_inv_fm)

@inline function field_contributions(dens, p::GasLiquidCoreParams)
    return (W=p.f_omega * dens.rho_B, R=p.f_rho * dens.rho_3)
end

@inline function physical_mu_pair(st::GasLiquidState, dens, p::GasLiquidCoreParams)
    fields = field_contributions(dens, p)
    return (mu_p=st.mu_p + fields.W + fields.R, mu_n=st.mu_n + fields.W - fields.R)
end

@inline function effective_chemical_potentials(mu_p::Real, mu_n::Real, dens, p::GasLiquidCoreParams)
    fields = field_contributions(dens, p)
    return (mu_p_star=mu_p - fields.W - fields.R, mu_n_star=mu_n - fields.W + fields.R)
end

@inline effective_chemical_potentials(st::GasLiquidState, dens, p::GasLiquidCoreParams) = (mu_p_star=st.mu_p, mu_n_star=st.mu_n)

"""Solve the vector-field fixed point for legacy states supplied with physical mu."""
function _density_bundle_from_physical_mu(S, D, mu_p::Real, mu_n::Real, T::Real, p::GasLiquidCoreParams;
    p_num::Int=96, p_max_inv_fm::Real=10.0, max_iter::Int=100, tol::Float64=1e-11)
    W = 0.0; R = 0.0
    dens = nothing
    st_eff = GasLiquidState(Float64(S), Float64(D), Float64(mu_p), Float64(mu_n))
    converged = false
    for _ in 1:max_iter
        st_eff = GasLiquidState(Float64(S), Float64(D), Float64(mu_p) - W - R, Float64(mu_n) - W + R)
        dens = _density_integrals(st_eff, T, p; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
        fields = field_contributions(dens, p)
        err = max(abs(fields.W - W), abs(fields.R - R))
        W = 0.5 * W + 0.5 * fields.W
        R = 0.5 * R + 0.5 * fields.R
        if err < tol
            converged = true
            break
        end
    end
    return (densities=dens, state=st_eff, fields=(W=W, R=R), converged=converged)
end

@inline function _field_residual_pair(st::GasLiquidState, dens, p::GasLiquidCoreParams)
    return (st.S - p.f_sigma * (dens.rho_s - p.b * p.m_nucleon_inv_fm * st.S^2 - p.c * st.S^3),
            st.D - p.f_delta * dens.rho_s3)
end

"""Return scalar/delta residuals plus either mu or density constraints.

With no target keywords the historical diagnostic meaning is preserved:
`st.mu_p/st.mu_n` are interpreted as physical mu values and the vector fixed
point is solved internally. New solver code always supplies explicit targets.
"""
function field_residuals(st::GasLiquidState, T::Real, p::GasLiquidCoreParams;
    p_num::Int=96, p_max_inv_fm::Real=10.0,
    mu_target::Union{Nothing, Real}=nothing,
    mu_p_target::Union{Nothing, Real}=nothing,
    mu_n_target::Union{Nothing, Real}=nothing,
    rho_B_target::Union{Nothing, Real}=nothing,
    rho_3_target::Union{Nothing, Real}=nothing)
    if mu_target !== nothing
        mu_p_target = mu_target
        mu_n_target = mu_target
    end

    if mu_p_target === nothing && mu_n_target === nothing && rho_B_target === nothing && rho_3_target === nothing
        legacy = _density_bundle_from_physical_mu(st.S, st.D, st.mu_p, st.mu_n, T, p; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
        rS, rD = _field_residual_pair(GasLiquidState(st.S, st.D, legacy.state.mu_p, legacy.state.mu_n), legacy.densities, p)
        return SVector{4, Float64}(rS, rD, 0.0, 0.0)
    end

    dens = density_bundle(st, T, p; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
    rS, rD = _field_residual_pair(st, dens, p)
    fields = field_contributions(dens, p)
    r3, r4 = if mu_p_target !== nothing || mu_n_target !== nothing
        mu_p_target === nothing && throw(ArgumentError("mu_p_target and mu_n_target must be supplied together"))
        mu_n_target === nothing && throw(ArgumentError("mu_p_target and mu_n_target must be supplied together"))
        (st.mu_p + fields.W + fields.R - Float64(mu_p_target), st.mu_n + fields.W - fields.R - Float64(mu_n_target))
    elseif rho_B_target !== nothing || rho_3_target !== nothing
        rho_B_target === nothing && throw(ArgumentError("rho_B_target and rho_3_target must be supplied together"))
        rho_3_target === nothing && throw(ArgumentError("rho_B_target and rho_3_target must be supplied together"))
        (dens.rho_B - Float64(rho_B_target), dens.rho_3 - Float64(rho_3_target))
    else
        (0.0, 0.0)
    end
    return SVector{4, Float64}(rS, rD, r3, r4)
end

struct RMFSolverResult{S, D, F}
    state::S
    densities::D
    fields::F
    mu_p::Float64
    mu_n::Float64
    mode::Symbol
    converged::Bool
    solver_status::Symbol
    iterations::Int
    residual_norm::Float64
    attempts::Int
    fallback_used::Bool
    failure_reason::String
end

@inline function _result_from_state(st, T, p; p_num, p_max_inv_fm, mu_p, mu_n, mode, converged, status, iterations, residual_norm, attempts, fallback_used, failure_reason="")
    dens = density_bundle(st, T, p; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
    fields = field_contributions(dens, p)
    return RMFSolverResult(st, dens, fields, Float64(mu_p), Float64(mu_n), mode, converged, status, iterations, residual_norm, attempts, fallback_used, String(failure_reason))
end

function _initial_guess_fixed_rho(p, rho_B_target, rho_3_target)
    rho_p = max((rho_B_target + rho_3_target) / 2, 0.0)
    rho_n = max((rho_B_target - rho_3_target) / 2, 0.0)
    kFp = (3.0 * pi^2 * rho_p)^(1 / 3)
    kFn = (3.0 * pi^2 * rho_n)^(1 / 3)
    m = p.m_nucleon_inv_fm
    return [0.20, 0.0, sqrt(kFp^2 + m^2), sqrt(kFn^2 + m^2)]
end

"""Solve either a fixed-mu or fixed-density RMF point."""
function solve_equilibrium(T::Real, mu_vec=0.0;
    params::GasLiquidCoreParams=GasLiquidCoreParams(),
    mode::Symbol=:fixed_mu,
    mu_B=nothing,
    mu_3=nothing,
    rho_B_target=nothing,
    rho_3_target=nothing,
    alpha=nothing,
    initial_guess::Union{Nothing, AbstractVector}=nothing,
    p_num::Int=96,
    p_max_inv_fm::Real=10.0,
    max_iter::Int=200,
    xtol::Float64=1e-10,
    ftol::Float64=1e-10,
    return_result::Bool=false)
    mode in (:fixed_mu, :fixed_rho) || throw(ArgumentError("mode must be :fixed_mu or :fixed_rho"))
    T >= 0 || throw(ArgumentError("temperature must be non-negative in fm^-1"))

    target_mu_p, target_mu_n, _, _ = _resolve_mu_pair(mu_vec; mu_B=mu_B, mu_3=mu_3)
    target_rho_B = nothing
    target_rho_3 = nothing
    if mode === :fixed_rho
        target_rho_B = rho_B_target === nothing ? throw(ArgumentError("rho_B_target is required for :fixed_rho")) : Float64(rho_B_target)
        target_rho_B >= 0 || throw(ArgumentError("rho_B_target must be non-negative"))
        if rho_3_target !== nothing
            target_rho_3 = Float64(rho_3_target)
        elseif alpha !== nothing
            a = Float64(alpha)
            -1.0 <= a <= 1.0 || throw(ArgumentError("alpha must lie in [-1,1]"))
            target_rho_3 = -a * target_rho_B
        else
            target_rho_3 = 0.0
        end
        abs(target_rho_3) <= target_rho_B + 1e-12 || throw(ArgumentError("|rho_3_target| must not exceed rho_B_target"))
    end

    x0 = if initial_guess === nothing
        mode === :fixed_mu ? [0.20, 0.0, target_mu_p, target_mu_n] : _initial_guess_fixed_rho(params, target_rho_B, target_rho_3)
    else
        collect(Float64.(initial_guess))
    end
    length(x0) == 4 || throw(ArgumentError("initial_guess must have length 4"))

    guesses = Vector{Vector{Float64}}()
    push!(guesses, x0)
    push!(guesses, [0.0, 0.0, x0[3], x0[4]])
    push!(guesses, [min(max(x0[1], 0.0), 1.2), 0.0, x0[3], x0[4]])

    last_state = GasLiquidState(x0[1], x0[2], x0[3], x0[4])
    last_status = :not_started
    last_norm = Inf
    last_iterations = 0
    last_failure_reason = ""
    for (attempt, guess) in enumerate(guesses)
        function residual!(F, x)
            st = GasLiquidState(x[1], x[2], x[3], x[4])
            r = if mode === :fixed_mu
                field_residuals(st, T, params; p_num=p_num, p_max_inv_fm=p_max_inv_fm, mu_p_target=target_mu_p, mu_n_target=target_mu_n)
            else
                field_residuals(st, T, params; p_num=p_num, p_max_inv_fm=p_max_inv_fm, rho_B_target=target_rho_B, rho_3_target=target_rho_3)
            end
            F[1] = r[1]; F[2] = r[2]; F[3] = r[3]; F[4] = r[4]
            return nothing
        end

        try
            res = nlsolve(residual!, guess; autodiff=:finite, method=:trust_region, xtol=xtol, ftol=ftol, iterations=max_iter)
            x = Float64.(res.zero)
            last_state = GasLiquidState(x[1], x[2], x[3], x[4])
            last_norm = Float64(res.residual_norm)
            last_iterations = hasproperty(res, :iterations) ? Int(res.iterations) : 0
            converged = Bool(getproperty(res, :f_converged)) && isfinite(last_norm) && last_norm <= max(ftol * 10, 1e-8)
            last_status = converged ? :converged : :max_iterations
            if converged
                dens = density_bundle(last_state, T, params; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
                masses = effective_masses(last_state, params)
                if !(isfinite(masses.mp) && isfinite(masses.mn) && masses.mp > 0.0 && masses.mn > 0.0)
                    last_status = :unphysical_mass
                    last_failure_reason = "effective nucleon mass is non-positive or non-finite"
                    converged = false
                else
                    solved_mu = physical_mu_pair(last_state, dens, params)
                    result = _result_from_state(last_state, T, params; p_num=p_num, p_max_inv_fm=p_max_inv_fm,
                        mu_p=mode === :fixed_mu ? target_mu_p : solved_mu.mu_p,
                        mu_n=mode === :fixed_mu ? target_mu_n : solved_mu.mu_n,
                        mode=mode, converged=true, status=:converged, iterations=last_iterations,
                        residual_norm=last_norm, attempts=attempt, fallback_used=attempt > 1, failure_reason="")
                    return return_result ? result : result.state
                end
            end
        catch err
            last_status = :exception
            last_norm = Inf
            last_iterations = 0
            last_failure_reason = sprint(showerror, err)
        end
    end

    failed_dens = density_bundle(last_state, T, params; p_num=p_num, p_max_inv_fm=p_max_inv_fm)
    failed_phys = physical_mu_pair(last_state, failed_dens, params)
    isempty(last_failure_reason) && (last_failure_reason = "solver did not meet residual tolerance")
    result = _result_from_state(last_state, T, params; p_num=p_num, p_max_inv_fm=p_max_inv_fm,
        mu_p=mode === :fixed_mu ? target_mu_p : failed_phys.mu_p,
        mu_n=mode === :fixed_mu ? target_mu_n : failed_phys.mu_n,
        mode=mode, converged=false, status=last_status, iterations=last_iterations,
        residual_norm=last_norm, attempts=length(guesses), fallback_used=length(guesses) > 1,
        failure_reason=last_failure_reason)
    return return_result ? result : throw(ErrorException("gas-liquid RMF solve failed: status=$(result.solver_status), residual_norm=$(result.residual_norm)"))
end

@inline solve_rmf(args...; kwargs...) = solve_equilibrium(args...; return_result=true, kwargs...)

end # module
