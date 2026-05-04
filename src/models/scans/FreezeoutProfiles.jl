"""
    FreezeoutProfiles

Chemical freeze-out baseline parameterizations used by charged/freeze-out
meson-density workflows.
"""
module FreezeoutProfiles

using Main.Constants_PNJL: ħc_MeV_fm

const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "config", "ConfigLoader.jl"))
if !isdefined(@__MODULE__, :ConfigLoader)
    Base.include(@__MODULE__, _CONFIG_LOADER_PATH)
end
using .ConfigLoader: load_config

export FreezeoutParameterProfile
export load_freezeout_profile, freezeout_profile_dir
export freezeout_muB_GeV, freezeout_temperature_GeV
export freezeout_point_GeV, freezeout_point_MeV
export build_freezeout_scan_points

const _FREEZEOUT_PROFILE_DIR = normpath(joinpath(
    @__DIR__, "..", "..", "..",
    "config", "physics", "freezeout",
))

struct FreezeoutParameterProfile
    profile_name::String
    family::String
    source_tag::String
    a_GeV::Float64
    b_GeV_inv1::Float64
    c_GeV_inv3::Float64
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
    value = Float64(coeffs[key])
    isfinite(value) || throw(ArgumentError("freezeout coefficient $(key) must be finite, got $(coeffs[key])"))
    return value
end

function _coerce_profile(cfg::Dict{String, Any})::FreezeoutParameterProfile
    coeffs = get(cfg, "coefficients", nothing)
    coeffs isa Dict{String, Any} || throw(ArgumentError("freezeout profile requires [coefficients] table"))
    return FreezeoutParameterProfile(
        _require_string(cfg, "profile_name"),
        _require_string(cfg, "family"),
        _require_string(cfg, "source_tag"),
        _require_coeff(coeffs, "a_GeV"),
        _require_coeff(coeffs, "b_GeV_inv1"),
        _require_coeff(coeffs, "c_GeV_inv3"),
        _require_coeff(coeffs, "d_GeV"),
        _require_coeff(coeffs, "e_GeV_inv1"),
    )
end

function load_freezeout_profile(; profile::String="default", dir::AbstractString=_FREEZEOUT_PROFILE_DIR)
    loaded = load_config(dir, _default_freezeout_config(); profile=profile)
    return _coerce_profile(loaded.config)
end

@inline function freezeout_muB_GeV(profile::FreezeoutParameterProfile, sqrt_s_NN_GeV::Real)
    s = Float64(sqrt_s_NN_GeV)
    s >= 0.0 || throw(ArgumentError("sqrt_s_NN_GeV must be nonnegative, got $(sqrt_s_NN_GeV)"))
    denom = 1.0 + profile.e_GeV_inv1 * s
    denom > 0.0 || throw(ArgumentError("freezeout muB denominator must stay positive, got $(denom)"))
    return profile.d_GeV / denom
end

@inline function freezeout_temperature_GeV(profile::FreezeoutParameterProfile, muB_GeV::Real)
    μ = Float64(muB_GeV)
    μ >= 0.0 || throw(ArgumentError("muB_GeV must be nonnegative, got $(muB_GeV)"))
    return profile.a_GeV - profile.b_GeV_inv1 * μ^2 - profile.c_GeV_inv3 * μ^4
end

@inline function freezeout_point_GeV(profile::FreezeoutParameterProfile, sqrt_s_NN_GeV::Real)
    μB_GeV = freezeout_muB_GeV(profile, sqrt_s_NN_GeV)
    T_GeV = freezeout_temperature_GeV(profile, μB_GeV)
    return (
        sqrt_s_NN_GeV=Float64(sqrt_s_NN_GeV),
        muB_GeV=μB_GeV,
        T_GeV=T_GeV,
    )
end

@inline function freezeout_point_MeV(profile::FreezeoutParameterProfile, sqrt_s_NN_GeV::Real)
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
    profile::FreezeoutParameterProfile=load_freezeout_profile(),
    traversal::Symbol=:sqrts_ascending,
)
    values = Float64.(collect(sqrt_s_NN_values))
    isempty(values) && throw(ArgumentError("sqrt_s_NN_values must not be empty"))
    points = [freezeout_point_MeV(profile, s) for s in values]
    return _sort_points(points, traversal)
end

end # module FreezeoutProfiles
