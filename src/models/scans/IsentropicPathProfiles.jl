"""
    IsentropicPathProfiles

最小 fixed-σ（等熵线）路径 profile 与点列生成器。
第一版只固化：

- `sigma_target`
- 温度点列
- 路径排序与路径元数据

不在此层提前抽象为更通用的外部 path family。
"""
module IsentropicPathProfiles

using Printf

export IsentropicPathProfile
export load_isentropic_path_profile
export build_isentropic_path_points

struct IsentropicPathProfile
    profile_name::String
    family::String
    source_tag::String
    sigma_target::Float64
end

@inline function _require_sigma_target(sigma_target::Real)
    σ = Float64(sigma_target)
    isfinite(σ) || throw(ArgumentError("sigma_target must be finite, got $(sigma_target)"))
    σ > 0.0 || throw(ArgumentError("sigma_target must be positive, got $(sigma_target)"))
    return σ
end

function load_isentropic_path_profile(;
    sigma_target::Real,
    profile::String="default",
)
    profile == "default" || throw(ArgumentError("unsupported isentropic path profile: $(profile)"))
    return IsentropicPathProfile(
        profile,
        "fixed_sigma_v1",
        "models.FixedSigma",
        _require_sigma_target(sigma_target),
    )
end

function _sort_temperature_values(values::Vector{Float64}, traversal::Symbol)
    if traversal === :as_given
        return values
    elseif traversal === :T_ascending
        return sort(values)
    elseif traversal === :T_descending
        return sort(values; rev=true)
    end
    throw(ArgumentError("unsupported isentropic traversal: $(traversal)"))
end

function build_isentropic_path_points(
    T_MeV_values;
    sigma_target::Real,
    profile::IsentropicPathProfile=load_isentropic_path_profile(; sigma_target=sigma_target),
    traversal::Symbol=:T_ascending,
)
    values = Float64.(collect(T_MeV_values))
    isempty(values) && throw(ArgumentError("T_MeV_values must not be empty"))
    for (idx, value) in pairs(values)
        isfinite(value) || throw(ArgumentError("T_MeV_values[$(idx)] must be finite, got $(value)"))
        value > 0.0 || throw(ArgumentError("T_MeV_values[$(idx)] must be positive, got $(value)"))
    end

    sorted = _sort_temperature_values(values, traversal)
    label = @sprintf("sigma=%.6f", profile.sigma_target)
    return map(enumerate(sorted)) do (idx, T_MeV)
        (
            path_family="isentropic",
            path_profile=profile.profile_name,
            path_segment="fixed_sigma",
            path_point_index=idx - 1,
            path_order_key=T_MeV,
            path_label=label,
            sigma_target=profile.sigma_target,
            T_MeV=T_MeV,
        )
    end
end

end # module IsentropicPathProfiles
