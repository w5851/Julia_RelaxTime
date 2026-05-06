"""
    FlavorChemicalProfiles

Flavor-level chemical-potential profiles used by charged/freeze-out
meson-density workflows.
"""
module FlavorChemicalProfiles

using Main.Constants_PNJL: ħc_MeV_fm

const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "config", "ConfigLoader.jl"))
if !isdefined(@__MODULE__, :ConfigLoader)
    Base.include(@__MODULE__, _CONFIG_LOADER_PATH)
end
using .ConfigLoader: load_config

export FlavorChemicalProfile
export flavor_chemical_profile_dir, load_flavor_chemical_profile
export flavor_mu_profile_MeV, flavor_mu_profile_fm

const _FLAVOR_CHEMICAL_PROFILE_DIR = normpath(joinpath(
    @__DIR__, "..", "..", "..",
    "config", "physics", "flavor_chemical",
))

struct FlavorChemicalProfile
    profile_name::String
    source_tag::String
    family::String
    apply_to_equilibrium::Bool
    mu_u_over_muq::Float64
    mu_d_over_muq::Float64
    mu_s_over_muq::Float64
end

@inline flavor_chemical_profile_dir() = _FLAVOR_CHEMICAL_PROFILE_DIR

function _default_flavor_chemical_config()
    return Dict{String, Any}(
        "profile_name" => "default",
        "source_tag" => "equal_flavor_fixedmu",
        "family" => "flavor_ratio_to_muq_v1",
        "apply_to_equilibrium" => true,
        "coefficients" => Dict{String, Any}(
            "mu_u_over_muq" => 1.0,
            "mu_d_over_muq" => 1.0,
            "mu_s_over_muq" => 1.0,
        ),
    )
end

@inline function _require_string(cfg::Dict{String, Any}, key::String)
    haskey(cfg, key) || throw(ArgumentError("flavor chemical profile missing key: $(key)"))
    value = cfg[key]
    value isa AbstractString || throw(ArgumentError("flavor chemical profile key $(key) must be string, got $(typeof(value))"))
    return String(value)
end

@inline function _require_bool(cfg::Dict{String, Any}, key::String)
    haskey(cfg, key) || throw(ArgumentError("flavor chemical profile missing key: $(key)"))
    value = cfg[key]
    value isa Bool || throw(ArgumentError("flavor chemical profile key $(key) must be bool, got $(typeof(value))"))
    return Bool(value)
end

@inline function _require_coeff(coeffs::Dict{String, Any}, key::String)
    haskey(coeffs, key) || throw(ArgumentError("flavor chemical profile missing coefficient: $(key)"))
    value = Float64(coeffs[key])
    isfinite(value) || throw(ArgumentError("flavor chemical coefficient $(key) must be finite, got $(coeffs[key])"))
    return value
end

function _coerce_profile(cfg::Dict{String, Any})::FlavorChemicalProfile
    coeffs = get(cfg, "coefficients", nothing)
    coeffs isa Dict{String, Any} || throw(ArgumentError("flavor chemical profile requires [coefficients] table"))

    family = _require_string(cfg, "family")
    family == "flavor_ratio_to_muq_v1" || throw(ArgumentError(
        "unsupported flavor chemical profile family: $(family)"
    ))

    return FlavorChemicalProfile(
        _require_string(cfg, "profile_name"),
        _require_string(cfg, "source_tag"),
        family,
        _require_bool(cfg, "apply_to_equilibrium"),
        _require_coeff(coeffs, "mu_u_over_muq"),
        _require_coeff(coeffs, "mu_d_over_muq"),
        _require_coeff(coeffs, "mu_s_over_muq"),
    )
end

function load_flavor_chemical_profile(; profile::String="default", dir::AbstractString=_FLAVOR_CHEMICAL_PROFILE_DIR)
    loaded = load_config(dir, _default_flavor_chemical_config(); profile=profile)
    return _coerce_profile(loaded.config)
end

function flavor_mu_profile_MeV(profile::FlavorChemicalProfile, muq_MeV::Real)
    μq = Float64(muq_MeV)
    isfinite(μq) || throw(ArgumentError("muq_MeV must be finite, got $(muq_MeV)"))

    return (
        profile_name=profile.profile_name,
        source_tag=profile.source_tag,
        family=profile.family,
        apply_to_equilibrium=profile.apply_to_equilibrium,
        muq_MeV=μq,
        mu_u_MeV=profile.mu_u_over_muq * μq,
        mu_d_MeV=profile.mu_d_over_muq * μq,
        mu_s_MeV=profile.mu_s_over_muq * μq,
    )
end

function flavor_mu_profile_fm(profile::FlavorChemicalProfile, muq_MeV::Real)
    mev = flavor_mu_profile_MeV(profile, muq_MeV)
    return merge(mev, (
        mu_u_fm=mev.mu_u_MeV / ħc_MeV_fm,
        mu_d_fm=mev.mu_d_MeV / ħc_MeV_fm,
        mu_s_fm=mev.mu_s_MeV / ħc_MeV_fm,
    ))
end

end # module FlavorChemicalProfiles
