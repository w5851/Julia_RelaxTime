"""
    MesonChemicalProfiles

Meson-level chemical-potential / degeneracy profiles used by charged and
freeze-out meson-density scans.
"""
module MesonChemicalProfiles

using Main.Constants_PNJL: ħc_MeV_fm

const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "config", "ConfigLoader.jl"))
if !isdefined(@__MODULE__, :ConfigLoader)
    Base.include(@__MODULE__, _CONFIG_LOADER_PATH)
end
using .ConfigLoader: load_config

export MesonChemicalProfile
export meson_chemical_profile_dir, load_meson_chemical_profile
export meson_chemical_profile_MeV, meson_chemical_profile_fm

const _MESON_CHEMICAL_PROFILE_DIR = normpath(joinpath(
    @__DIR__, "..", "..", "..",
    "config", "physics", "meson_chemical",
))

struct MesonChemicalProfile
    profile_name::String
    source_tag::String
    pi_label::String
    k_label::String
    charge_resolved::Bool
    mu_pi_MeV::Float64
    mu_K_MeV::Float64
    d_pi::Int
    d_K::Int
end

@inline meson_chemical_profile_dir() = _MESON_CHEMICAL_PROFILE_DIR

function _default_meson_chemical_config()
    return Dict{String, Any}(
        "profile_name" => "default",
        "source_tag" => "neutral_aggregate_default",
        "pi_label" => "pi",
        "k_label" => "K",
        "charge_resolved" => false,
        "chemical_potentials" => Dict{String, Any}(
            "mu_pi_MeV" => 0.0,
            "mu_K_MeV" => 0.0,
        ),
        "degeneracies" => Dict{String, Any}(
            "d_pi" => 3,
            "d_K" => 4,
        ),
    )
end

@inline function _require_string(cfg::Dict{String, Any}, key::String)
    haskey(cfg, key) || throw(ArgumentError("meson chemical profile missing key: $(key)"))
    value = cfg[key]
    value isa AbstractString || throw(ArgumentError("meson chemical profile key $(key) must be string, got $(typeof(value))"))
    return String(value)
end

@inline function _require_bool(cfg::Dict{String, Any}, key::String)
    haskey(cfg, key) || throw(ArgumentError("meson chemical profile missing key: $(key)"))
    value = cfg[key]
    value isa Bool || throw(ArgumentError("meson chemical profile key $(key) must be bool, got $(typeof(value))"))
    return Bool(value)
end

@inline function _require_float(cfg::Dict{String, Any}, key::String)
    haskey(cfg, key) || throw(ArgumentError("meson chemical profile missing key: $(key)"))
    value = Float64(cfg[key])
    isfinite(value) || throw(ArgumentError("meson chemical profile key $(key) must be finite, got $(cfg[key])"))
    return value
end

@inline function _require_positive_int(cfg::Dict{String, Any}, key::String)
    haskey(cfg, key) || throw(ArgumentError("meson chemical profile missing key: $(key)"))
    value = Int(cfg[key])
    value > 0 || throw(ArgumentError("meson chemical profile key $(key) must be positive, got $(cfg[key])"))
    return value
end

function _coerce_profile(cfg::Dict{String, Any})::MesonChemicalProfile
    chemical = get(cfg, "chemical_potentials", nothing)
    chemical isa Dict{String, Any} || throw(ArgumentError("meson chemical profile requires [chemical_potentials] table"))
    degeneracies = get(cfg, "degeneracies", nothing)
    degeneracies isa Dict{String, Any} || throw(ArgumentError("meson chemical profile requires [degeneracies] table"))

    return MesonChemicalProfile(
        _require_string(cfg, "profile_name"),
        _require_string(cfg, "source_tag"),
        _require_string(cfg, "pi_label"),
        _require_string(cfg, "k_label"),
        _require_bool(cfg, "charge_resolved"),
        _require_float(chemical, "mu_pi_MeV"),
        _require_float(chemical, "mu_K_MeV"),
        _require_positive_int(degeneracies, "d_pi"),
        _require_positive_int(degeneracies, "d_K"),
    )
end

function load_meson_chemical_profile(; profile::String="default", dir::AbstractString=_MESON_CHEMICAL_PROFILE_DIR)
    loaded = load_config(dir, _default_meson_chemical_config(); profile=profile)
    return _coerce_profile(loaded.config)
end

@inline function meson_chemical_profile_MeV(profile::MesonChemicalProfile)
    return (
        profile_name=profile.profile_name,
        source_tag=profile.source_tag,
        pi_label=profile.pi_label,
        k_label=profile.k_label,
        charge_resolved=profile.charge_resolved,
        mu_pi_MeV=profile.mu_pi_MeV,
        mu_K_MeV=profile.mu_K_MeV,
        d_pi=profile.d_pi,
        d_K=profile.d_K,
    )
end

@inline function meson_chemical_profile_fm(profile::MesonChemicalProfile)
    return (
        profile_name=profile.profile_name,
        source_tag=profile.source_tag,
        pi_label=profile.pi_label,
        k_label=profile.k_label,
        charge_resolved=profile.charge_resolved,
        mu_pi_fm=profile.mu_pi_MeV / ħc_MeV_fm,
        mu_K_fm=profile.mu_K_MeV / ħc_MeV_fm,
        d_pi=profile.d_pi,
        d_K=profile.d_K,
    )
end

end # module MesonChemicalProfiles
