const _MANIFEST_EXTENSION_ALLOWED_KEYS = Set{String}([
    "pipeline_family",
    "baseline_suite",
    "physics_profile",
    "adapter_version",
])

function _validate_manifest_extension_keys(keys_iter)
    invalid = String[]
    for key in keys_iter
        key_str = String(key)
        key_str in _MANIFEST_EXTENSION_ALLOWED_KEYS || push!(invalid, key_str)
    end
    isempty(invalid) || throw(ArgumentError("unsupported manifest extension key(s): $(join(sort!(unique!(invalid)), ", "))"))
    return nothing
end

function normalize_manifest_extensions(meta)::Dict{String, Any}
    _validate_manifest_extension_keys(keys(meta))
    normalized = Dict{String, Any}()
    for (k, v) in pairs(meta)
        normalized[String(k)] = v
    end
    return normalized
end

function build_manifest_extensions(meta::NamedTuple)
    defaults = Dict{String, Any}(
        "pipeline_family" => "generic",
        "baseline_suite" => "none",
        "physics_profile" => "default",
        "adapter_version" => "v1",
    )

    merged = Dict{String, Any}()
    merge!(merged, defaults)

    merge!(merged, normalize_manifest_extensions(meta))

    return merged
end
