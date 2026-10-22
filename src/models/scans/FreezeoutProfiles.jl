"""
    FreezeoutProfiles

Chemical freeze-out baseline parameterizations used by charged/freeze-out
meson-density workflows.
"""
module FreezeoutProfiles

using Main.Constants_PNJL: ħc_MeV_fm
using TOML

const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "config", "ConfigLoader.jl"))
if !isdefined(@__MODULE__, :ConfigLoader)
    Base.include(@__MODULE__, _CONFIG_LOADER_PATH)
end
using .ConfigLoader: deep_merge

export AbstractFreezeoutParameterProfile, FreezeoutParameterProfile, AndronicFreezeoutProfile
export load_freezeout_profile, freezeout_profile_dir
export freezeout_muB_GeV, freezeout_temperature_GeV
export freezeout_point_GeV, freezeout_point_MeV
export build_freezeout_scan_points

const _FREEZEOUT_PROFILE_DIR = normpath(joinpath(
    @__DIR__, "..", "..", "..",
    "config", "physics", "freezeout",
))

abstract type AbstractFreezeoutParameterProfile end

"""Cleymans polynomial in muB, with the historical constructor and fields."""
struct FreezeoutParameterProfile <: AbstractFreezeoutParameterProfile
    profile_name::String
    family::String
    source_tag::String
    a_GeV::Float64
    b_GeV_inv1::Float64
    c_GeV_inv3::Float64
    d_GeV::Float64
    e_GeV_inv1::Float64
end

"""Andronic logistic in log(sqrt(sNN)/GeV); not a crossover parameterization."""
struct AndronicFreezeoutProfile <: AbstractFreezeoutParameterProfile
    profile_name::String
    family::String
    source_tag::String
    T_lim_GeV::Float64
    log_offset::Float64
    log_width::Float64
    d_GeV::Float64
    e_GeV_inv1::Float64
end

@inline freezeout_profile_dir() = _FREEZEOUT_PROFILE_DIR

function _default_freezeout_config()
    return Dict{String, Any}(
        "profile_name" => "default",
        "family" => "cleymans_polynomial_rational_v1",
        "source_tag" => "cleymans2006_like",
        "coefficients" => Dict{String, Any}(
            "a_GeV" => 0.166,
            "b_GeV_inv1" => 0.139,
            "c_GeV_inv3" => 0.053,
            "d_GeV" => 1.308,
            "e_GeV_inv1" => 0.273,
        ),
    )
end

@inline function _require_string(cfg::Dict{String, Any}, key::String)
    haskey(cfg, key) || throw(ArgumentError("freezeout profile missing key: $(key)"))
    value = cfg[key]
    value isa AbstractString || throw(ArgumentError("freezeout profile key $(key) must be string, got $(typeof(value))"))
    return String(value)
end

@inline function _require_coeff(coeffs::Dict{String, Any}, key::String)
    haskey(coeffs, key) || throw(ArgumentError("freezeout profile missing coefficient: $(key)"))
    coeffs[key] isa Real && !(coeffs[key] isa Bool) ||
        throw(ArgumentError("freezeout coefficient $(key) must be a real number"))
    value = Float64(coeffs[key])
    isfinite(value) || throw(ArgumentError("freezeout coefficient $(key) must be finite, got $(coeffs[key])"))
    return value
end

function _coerce_profile(cfg::Dict{String, Any})
    coeffs = get(cfg, "coefficients", nothing)
    coeffs isa Dict{String, Any} || throw(ArgumentError("freezeout profile requires [coefficients] table"))
    name = _require_string(cfg, "profile_name")
    family = _require_string(cfg, "family")
    source = _require_string(cfg, "source_tag")
    d = _require_coeff(coeffs, "d_GeV")
    e = _require_coeff(coeffs, "e_GeV_inv1")
    d > 0 && e > 0 || throw(ArgumentError("freezeout d_GeV and e_GeV_inv1 must be positive"))
    if family == "cleymans_polynomial_rational_v1"
        a = _require_coeff(coeffs, "a_GeV")
        b = _require_coeff(coeffs, "b_GeV_inv1")
        c = _require_coeff(coeffs, "c_GeV_inv3")
        a > 0 && b >= 0 && c >= 0 ||
            throw(ArgumentError("freezeout requires a_GeV>0, b_GeV_inv1>=0 and c_GeV_inv3>=0"))
        return FreezeoutParameterProfile(name, family, source, a, b, c, d, e)
    elseif family == "andronic_logistic_rational_v1"
        Tlim = _require_coeff(coeffs, "T_lim_GeV")
        offset = _require_coeff(coeffs, "log_offset")
        width = _require_coeff(coeffs, "log_width")
        Tlim > 0 && width > 0 ||
            throw(ArgumentError("freezeout T_lim_GeV and log_width must be positive"))
        return AndronicFreezeoutProfile(name, family, source, Tlim, offset, width, d, e)
    end
    throw(ArgumentError("unsupported chemical freezeout family: $(family)"))
end

function load_freezeout_profile(; profile::String="default", dir::AbstractString=_FREEZEOUT_PROFILE_DIR)
    path = joinpath(dir, profile * ".toml")
    if !isfile(path)
        profile == "default" && return _coerce_profile(_default_freezeout_config())
        throw(ArgumentError("freezeout profile not found: $(path)"))
    end
    cfg = try
        TOML.parsefile(path)
    catch err
        throw(ArgumentError("invalid freezeout profile $(path): $(sprint(showerror, err))"))
    end
    # Preserve legacy Cleymans overrides. A different family must supply all of
    # its own coefficients, so a partial Andronic file cannot inherit d/e from
    # Cleymans and silently create a hybrid chemical freezeout line.
    if get(cfg, "family", "cleymans_polynomial_rational_v1") == "cleymans_polynomial_rational_v1"
        cfg = deep_merge(_default_freezeout_config(), cfg)
    end
    return _coerce_profile(cfg)
end

@inline function freezeout_muB_GeV(profile::AbstractFreezeoutParameterProfile, sqrt_s_NN_GeV::Real)
    s = Float64(sqrt_s_NN_GeV)
    isfinite(s) && s >= 0.0 || throw(ArgumentError("sqrt_s_NN_GeV must be finite and nonnegative, got $(sqrt_s_NN_GeV)"))
    denom = 1.0 + profile.e_GeV_inv1 * s
    denom > 0.0 || throw(ArgumentError("freezeout muB denominator must stay positive, got $(denom)"))
    return profile.d_GeV / denom
end

@inline function freezeout_temperature_GeV(profile::FreezeoutParameterProfile, muB_GeV::Real)
    μ = Float64(muB_GeV)
    isfinite(μ) && μ >= 0.0 || throw(ArgumentError("muB_GeV must be finite and nonnegative, got $(muB_GeV)"))
    return profile.a_GeV - profile.b_GeV_inv1 * μ^2 - profile.c_GeV_inv3 * μ^4
end

@inline function _temperature_at_energy(profile::AndronicFreezeoutProfile, s, muB)
    s == 0 && return 0.0
    # The division acts on log(s) only: 2.60 - log(s)/0.45.
    x = profile.log_offset - log(s) / profile.log_width
    return x > 0 ? profile.T_lim_GeV * exp(-x) / (1 + exp(-x)) :
        profile.T_lim_GeV / (1 + exp(x))
end

@inline _temperature_at_energy(profile::FreezeoutParameterProfile, s, muB) =
    freezeout_temperature_GeV(profile, muB)

@inline function freezeout_temperature_GeV(profile::AndronicFreezeoutProfile, muB_GeV::Real)
    μ = Float64(muB_GeV)
    isfinite(μ) && 0 <= μ <= profile.d_GeV ||
        throw(ArgumentError("Andronic muB_GeV must be finite and in [0, d_GeV], got $(muB_GeV)"))
    μ == 0 && return profile.T_lim_GeV
    s = (profile.d_GeV / μ - 1) / profile.e_GeV_inv1
    return _temperature_at_energy(profile, s, μ)
end

@inline function freezeout_point_GeV(profile::AbstractFreezeoutParameterProfile, sqrt_s_NN_GeV::Real)
    μB_GeV = freezeout_muB_GeV(profile, sqrt_s_NN_GeV)
    T_GeV = _temperature_at_energy(profile, Float64(sqrt_s_NN_GeV), μB_GeV)
    return (
        sqrt_s_NN_GeV=Float64(sqrt_s_NN_GeV),
        muB_GeV=μB_GeV,
        T_GeV=T_GeV,
    )
end

@inline function freezeout_point_MeV(profile::AbstractFreezeoutParameterProfile, sqrt_s_NN_GeV::Real)
    pt = freezeout_point_GeV(profile, sqrt_s_NN_GeV)
    return (
        sqrt_s_NN_GeV=pt.sqrt_s_NN_GeV,
        muB_GeV=pt.muB_GeV,
        T_GeV=pt.T_GeV,
        muB_MeV=pt.muB_GeV * 1000.0,
        T_MeV=pt.T_GeV * 1000.0,
        muq_MeV=pt.muB_GeV * 1000.0 / 3.0,
        muB_fm=pt.muB_GeV * 1000.0 / ħc_MeV_fm,
        T_fm=pt.T_GeV * 1000.0 / ħc_MeV_fm,
    )
end

function _sort_points(points::Vector{<:NamedTuple}, traversal::Symbol)
    if traversal === :as_given
        return points
    elseif traversal === :sqrts_ascending
        return sort(points; by=p -> p.sqrt_s_NN_GeV)
    elseif traversal === :sqrts_descending
        return sort(points; by=p -> p.sqrt_s_NN_GeV, rev=true)
    elseif traversal === :muB_descending
        return sort(points; by=p -> p.muB_GeV, rev=true)
    elseif traversal === :muB_ascending
        return sort(points; by=p -> p.muB_GeV)
    end
    throw(ArgumentError("unsupported freezeout traversal: $(traversal)"))
end

function build_freezeout_scan_points(
    sqrt_s_NN_values;
    profile::AbstractFreezeoutParameterProfile=load_freezeout_profile(),
    traversal::Symbol=:sqrts_ascending,
)
    values = Float64.(collect(sqrt_s_NN_values))
    isempty(values) && throw(ArgumentError("sqrt_s_NN_values must not be empty"))
    points = [freezeout_point_MeV(profile, s) for s in values]
    return _sort_points(points, traversal)
end

end # module FreezeoutProfiles
