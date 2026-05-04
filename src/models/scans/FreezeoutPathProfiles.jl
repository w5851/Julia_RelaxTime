"""
    FreezeoutPathProfiles

Path-level profiles layered on top of freeze-out baseline parameterizations.
They describe how `(sqrt(s_NN) -> (T, mu_B))` should be interpreted for workflow
scans when the actual literature scan path is not identical to the raw chemical
freeze-out baseline.
"""
module FreezeoutPathProfiles

using Main.Constants_PNJL: ħc_MeV_fm

const _CONFIG_LOADER_PATH = normpath(joinpath(@__DIR__, "..", "..", "config", "ConfigLoader.jl"))
if !isdefined(@__MODULE__, :ConfigLoader)
    Base.include(@__MODULE__, _CONFIG_LOADER_PATH)
end
using .ConfigLoader: load_config

using ..FreezeoutProfiles

export FreezeoutPathProfile
export freezeout_path_profile_dir, load_freezeout_path_profile
export build_freezeout_path_points

const _FREEZEOUT_PATH_PROFILE_DIR = normpath(joinpath(
    @__DIR__, "..", "..", "..",
    "config", "physics", "freezeout_path",
))

struct FreezeoutPathProfile
    profile_name::String
    family::String
    source_tag::String
    switch_sqrt_s_NN_GeV::Float64
    constant_T_MeV::Float64
end

@inline freezeout_path_profile_dir() = _FREEZEOUT_PATH_PROFILE_DIR

function _default_path_profile_config()
    return Dict{String, Any}(
        "profile_name" => "baseline_freezeout",
        "family" => "baseline_freezeout_v1",
        "source_tag" => "raw_chemical_freezeout_baseline",
        "parameters" => Dict{String, Any}(
            "switch_sqrt_s_NN_GeV" => -1.0,
            "constant_T_MeV" => -1.0,
        ),
    )
end

@inline function _require_string(cfg::Dict{String, Any}, key::String)
    haskey(cfg, key) || throw(ArgumentError("freezeout path profile missing key: $(key)"))
    value = cfg[key]
    value isa AbstractString || throw(ArgumentError("freezeout path profile key $(key) must be string, got $(typeof(value))"))
    return String(value)
end

@inline function _optional_param(parameters::Dict{String, Any}, key::String, fallback::Float64)
    haskey(parameters, key) || return fallback
    value = Float64(parameters[key])
    isfinite(value) || throw(ArgumentError("freezeout path parameter $(key) must be finite, got $(parameters[key])"))
    return value
end

function _coerce_profile(cfg::Dict{String, Any})::FreezeoutPathProfile
    params = get(cfg, "parameters", Dict{String, Any}())
    params isa Dict{String, Any} || throw(ArgumentError("freezeout path profile [parameters] must be a table"))
    return FreezeoutPathProfile(
        _require_string(cfg, "profile_name"),
        _require_string(cfg, "family"),
        _require_string(cfg, "source_tag"),
        _optional_param(params, "switch_sqrt_s_NN_GeV", -1.0),
        _optional_param(params, "constant_T_MeV", -1.0),
    )
end

function load_freezeout_path_profile(; profile::String="baseline_freezeout", dir::AbstractString=_FREEZEOUT_PATH_PROFILE_DIR)
    loaded = load_config(dir, _default_path_profile_config(); profile=profile)
    return _coerce_profile(loaded.config)
end

@inline function _segment_label(profile::FreezeoutPathProfile, sqrt_s_NN_GeV::Float64)
    if profile.family == "baseline_freezeout_v1"
        return "baseline_freezeout"
    elseif profile.family == "freezeout_plus_constT_v1"
        if profile.switch_sqrt_s_NN_GeV >= 0.0 && sqrt_s_NN_GeV > profile.switch_sqrt_s_NN_GeV
            return "constant_T"
        end
        return "baseline_freezeout"
    end
    throw(ArgumentError("unsupported freezeout path family: $(profile.family)"))
end

function _apply_profile_to_point(profile::FreezeoutPathProfile, pt)
    seg = _segment_label(profile, Float64(pt.sqrt_s_NN_GeV))
    if seg == "constant_T"
        T_MeV = profile.constant_T_MeV
        T_MeV > 0.0 || throw(ArgumentError("constant_T_MeV must be positive for path profile $(profile.profile_name)"))
        return merge(pt, (
            T_MeV=T_MeV,
            T_GeV=T_MeV / 1000.0,
            T_fm=T_MeV / ħc_MeV_fm,
            path_profile=profile.profile_name,
            path_segment=seg,
        ))
    end
    return merge(pt, (
        path_profile=profile.profile_name,
        path_segment=seg,
    ))
end

function build_freezeout_path_points(
    sqrt_s_NN_values;
    freezeout_profile::FreezeoutProfiles.FreezeoutParameterProfile=FreezeoutProfiles.load_freezeout_profile(),
    path_profile::FreezeoutPathProfile=load_freezeout_path_profile(),
    traversal::Symbol=:sqrts_ascending,
)
    base_points = FreezeoutProfiles.build_freezeout_scan_points(
        sqrt_s_NN_values;
        profile=freezeout_profile,
        traversal=traversal,
    )
    return [_apply_profile_to_point(path_profile, pt) for pt in base_points]
end

end # module FreezeoutPathProfiles
