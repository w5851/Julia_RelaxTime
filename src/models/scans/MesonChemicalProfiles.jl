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
    pi_channel::Symbol
    k_channel::Symbol
    charge_resolved::Bool
    mu_pi_rule::Symbol
    mu_K_rule::Symbol
    mu_pi_MeV::Float64
    mu_K_MeV::Float64
    d_pi::Int
    d_K::Int
end

function MesonChemicalProfile(
    profile_name::String,
    source_tag::String,
    pi_label::String,
    k_label::String,
    pi_channel::Symbol,
    k_channel::Symbol,
    charge_resolved::Bool,
    mu_K_rule::Symbol,
    mu_pi_MeV::Float64,
    mu_K_MeV::Float64,
    d_pi::Int,
    d_K::Int,
)
    return MesonChemicalProfile(
        profile_name,
        source_tag,
        pi_label,
        k_label,
        pi_channel,
        k_channel,
        charge_resolved,
        :constant,
        mu_K_rule,
        mu_pi_MeV,
        mu_K_MeV,
        d_pi,
        d_K,
    )
end

@inline meson_chemical_profile_dir() = _MESON_CHEMICAL_PROFILE_DIR

function _default_meson_chemical_config()
    return Dict{String, Any}(
        "profile_name" => "default",
        "source_tag" => "neutral_aggregate_default",
        "pi_label" => "pi",
        "k_label" => "K",
        "charge_resolved" => false,
        "chemical_potential_rules" => Dict{String, Any}(
            "mu_pi_rule" => "constant",
            "mu_K_rule" => "constant",
        ),
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

@inline function _resolve_meson_channel(label::AbstractString, family::Symbol)
    value = lowercase(strip(String(label)))
    if family === :pi
        value in ("pi", "π") && return :pi
        value in ("pi_plus", "pi+", "π+") && return :pi_plus
        value in ("pi_minus", "pi-", "π-") && return :pi_minus
        throw(ArgumentError("unsupported pion label $(label); use pi, pi_plus, or pi_minus"))
    elseif family === :K
        value in ("k",) && return :K
        value in ("k_plus", "k+") && return :K_plus
        value in ("k_minus", "k-") && return :K_minus
        throw(ArgumentError("unsupported kaon label $(label); use K, K_plus, or K_minus"))
    end
    throw(ArgumentError("unsupported meson family $(family)"))
end

@inline function _resolve_mu_pi_rule(value::AbstractString)
    key = lowercase(strip(String(value)))
    key == "constant" && return :constant
    key == "mu_u_minus_mu_d_signed" && return :mu_u_minus_mu_d_signed
    throw(ArgumentError("unsupported mu_pi_rule $(value); use constant or mu_u_minus_mu_d_signed"))
end

@inline function _resolve_mu_K_rule(value::AbstractString)
    key = lowercase(strip(String(value)))
    key == "constant" && return :constant
    key == "mu_u_minus_mu_s_signed" && return :mu_u_minus_mu_s_signed
    throw(ArgumentError("unsupported mu_K_rule $(value); use constant or mu_u_minus_mu_s_signed"))
end

@inline function _signed_pion_mu_from_flavor(pi_channel::Symbol, flavor_mev)::Float64
    Δ = Float64(flavor_mev.mu_u_MeV) - Float64(flavor_mev.mu_d_MeV)
    if pi_channel === :pi_plus
        return Δ
    elseif pi_channel === :pi_minus
        return -Δ
    elseif pi_channel === :pi
        return Δ
    end
    throw(ArgumentError("mu_pi signed flavor rule requires pion channel, got $(pi_channel)"))
end

@inline function _signed_kaon_mu_from_flavor(k_channel::Symbol, flavor_mev)::Float64
    Δ = Float64(flavor_mev.mu_u_MeV) - Float64(flavor_mev.mu_s_MeV)
    if k_channel === :K_plus
        return Δ
    elseif k_channel === :K_minus
        return -Δ
    elseif k_channel === :K
        return Δ
    end
    throw(ArgumentError("mu_K signed flavor rule requires kaon channel, got $(k_channel)"))
end

function _coerce_profile(cfg::Dict{String, Any})::MesonChemicalProfile
    chemical = get(cfg, "chemical_potentials", nothing)
    chemical isa Dict{String, Any} || throw(ArgumentError("meson chemical profile requires [chemical_potentials] table"))
    degeneracies = get(cfg, "degeneracies", nothing)
    degeneracies isa Dict{String, Any} || throw(ArgumentError("meson chemical profile requires [degeneracies] table"))
    rules = get(cfg, "chemical_potential_rules", Dict{String, Any}())
    rules isa Dict{String, Any} || throw(ArgumentError("meson chemical profile [chemical_potential_rules] must be table"))

    pi_label = _require_string(cfg, "pi_label")
    k_label = _require_string(cfg, "k_label")

    return MesonChemicalProfile(
        _require_string(cfg, "profile_name"),
        _require_string(cfg, "source_tag"),
        pi_label,
        k_label,
        _resolve_meson_channel(pi_label, :pi),
        _resolve_meson_channel(k_label, :K),
        _require_bool(cfg, "charge_resolved"),
        _resolve_mu_pi_rule(String(get(rules, "mu_pi_rule", "constant"))),
        _resolve_mu_K_rule(String(get(rules, "mu_K_rule", "constant"))),
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
        pi_channel=profile.pi_channel,
        k_channel=profile.k_channel,
        charge_resolved=profile.charge_resolved,
        mu_pi_rule=profile.mu_pi_rule,
        mu_K_rule=profile.mu_K_rule,
        mu_pi_MeV=profile.mu_pi_MeV,
        mu_K_MeV=profile.mu_K_MeV,
        d_pi=profile.d_pi,
        d_K=profile.d_K,
    )
end

function meson_chemical_profile_fm(profile::MesonChemicalProfile; flavor_mev=nothing)
    mu_pi_MeV = if profile.mu_pi_rule === :constant
        profile.mu_pi_MeV
    else
        flavor_mev === nothing && throw(ArgumentError("meson chemical profile $(profile.profile_name) requires flavor_mev to resolve mu_pi_rule=$(profile.mu_pi_rule)"))
        _signed_pion_mu_from_flavor(profile.pi_channel, flavor_mev)
    end
    mu_K_MeV = if profile.mu_K_rule === :constant
        profile.mu_K_MeV
    else
        flavor_mev === nothing && throw(ArgumentError("meson chemical profile $(profile.profile_name) requires flavor_mev to resolve mu_K_rule=$(profile.mu_K_rule)"))
        _signed_kaon_mu_from_flavor(profile.k_channel, flavor_mev)
    end
    return (
        profile_name=profile.profile_name,
        source_tag=profile.source_tag,
        pi_label=profile.pi_label,
        k_label=profile.k_label,
        pi_channel=profile.pi_channel,
        k_channel=profile.k_channel,
        charge_resolved=profile.charge_resolved,
        mu_pi_rule=profile.mu_pi_rule,
        mu_K_rule=profile.mu_K_rule,
        mu_pi_fm=mu_pi_MeV / ħc_MeV_fm,
        mu_K_fm=mu_K_MeV / ħc_MeV_fm,
        d_pi=profile.d_pi,
        d_K=profile.d_K,
    )
end

end # module MesonChemicalProfiles
