"""ConfigLoader

一个小而通用的 TOML 配置加载器。

目标：
- 让不同模型（NJL/PNJL/rPNJL/…）复用同一套“默认配置 + profile 覆盖”的加载逻辑
- 保持 include 驱动项目结构下的可用性（不要求 package 化）

约定：
- `dir` 下的配置文件命名为 `<profile>.toml`（默认 `default`）
- 若文件不存在：使用 `default_config`
- 若文件存在但解析失败：警告并回退到 `default_config`

返回：`(config, profile, path)`
"""
module ConfigLoader

using Printf: @sprintf
using SHA: sha256
using TOML

export deep_merge, load_config, reset_config_loader_cache!, config_loader_cache_stats

const _CONFIG_CACHE_LOCK = ReentrantLock()
const _CONFIG_CACHE = Dict{String, NamedTuple{(:config, :profile, :path), Tuple{Dict{String, Any}, String, Union{Nothing, String}}}}()
const _CONFIG_CACHE_STATS = Dict{Symbol, Int}(:hit => 0, :miss => 0)

function _normalize_dir(dir::AbstractString)
    candidate = abspath(normpath(dir))
    return isdir(candidate) ? realpath(candidate) : candidate
end

function deep_merge(base::Dict{String, Any}, override::Dict{String, Any})
    result = deepcopy(base)
    for (k, v) in override
        if v isa Dict{String, Any} && haskey(result, k) && result[k] isa Dict{String, Any}
            result[k] = deep_merge(result[k]::Dict{String, Any}, v)
        else
            result[k] = v
        end
    end
    return result
end

function _json_escape(s::AbstractString)
    io = IOBuffer()
    print(io, '"')
    print(io, Base.escape_string(String(s)))
    print(io, '"')
    return String(take!(io))
end

function _canonical_json(io::IO, value)
    if value isa Dict
        pairs = collect(value)
        normalized = map(pairs) do (k, v)
            (String(k), v)
        end
        sort!(normalized; by=first)
        print(io, '{')
        for (idx, (k, v)) in enumerate(normalized)
            idx > 1 && print(io, ',')
            print(io, _json_escape(k), ':')
            _canonical_json(io, v)
        end
        print(io, '}')
        return
    end

    if value isa NamedTuple
        _canonical_json(io, Dict{String, Any}(string(k) => v for (k, v) in pairs(value)))
        return
    end

    if value isa AbstractVector || value isa Tuple
        print(io, '[')
        for (idx, element) in enumerate(value)
            idx > 1 && print(io, ',')
            _canonical_json(io, element)
        end
        print(io, ']')
        return
    end

    if value isa AbstractString
        print(io, _json_escape(value))
        return
    end

    if value isa Bool
        print(io, value ? "true" : "false")
        return
    end

    if value isa Integer
        print(io, value)
        return
    end

    if value isa AbstractFloat
        if !isfinite(value)
            throw(ArgumentError("Non-finite float value found in config fingerprint input"))
        end
        normalized = iszero(value) ? zero(value) : value
        print(io, @sprintf("%.17g", normalized))
        return
    end

    if value === nothing
        print(io, "null")
        return
    end

    throw(ArgumentError("Unsupported fingerprint value type: $(typeof(value))"))
end

function _canonical_json(value)
    io = IOBuffer()
    _canonical_json(io, value)
    return String(take!(io))
end

function _hex_digest(bytes::Vector{UInt8})
    io = IOBuffer()
    for b in bytes
        print(io, lowercase(string(b; base=16, pad=2)))
    end
    return String(take!(io))
end

function _profile_content_hash(path::AbstractString)
    if !isfile(path)
        return "__MISSING__"
    end
    return _hex_digest(collect(sha256(read(path, String))))
end

function _config_fingerprint(
    dir::AbstractString,
    default_config::Dict{String, Any},
    profile::String,
    inherited_configs::AbstractVector,
)
    normalized_dir = _normalize_dir(dir)
    path = joinpath(normalized_dir, string(profile, ".toml"))
    inherited_payload = Any[
        inherited isa Dict{String, Any} ? inherited : Dict{String, Any}(inherited)
        for inherited in inherited_configs
    ]
    payload = Dict{String, Any}(
        "dir" => normalized_dir,
        "profile" => profile,
        "default_config" => default_config,
        "inherited_configs" => inherited_payload,
        "profile_content_hash" => _profile_content_hash(path),
    )
    canonical = _canonical_json(payload)
    return _hex_digest(collect(sha256(canonical)))
end

function reset_config_loader_cache!()
    lock(_CONFIG_CACHE_LOCK) do
        empty!(_CONFIG_CACHE)
        _CONFIG_CACHE_STATS[:hit] = 0
        _CONFIG_CACHE_STATS[:miss] = 0
    end
    return nothing
end

function config_loader_cache_stats()
    lock(_CONFIG_CACHE_LOCK) do
        return (
            entries=length(_CONFIG_CACHE),
            hit=_CONFIG_CACHE_STATS[:hit],
            miss=_CONFIG_CACHE_STATS[:miss],
        )
    end
end

"""load_config(dir, default_config; profile="default", inherited_configs=Dict[])

从 `dir/<profile>.toml` 读取配置，并与 `default_config` deep-merge。

可通过 `inherited_configs` 传入一组“父配置”Dict，按顺序先合并父配置再合并 profile 文件，
用于实现跨目录/跨层级的参数继承（例如 physics 层共享参数注入模型层）。

返回 `(config, profile, path)`，其中 `path` 若文件不存在则为 `nothing`。
"""
function load_config(
    dir::AbstractString,
    default_config::Dict{String, Any};
    profile::String="default",
    inherited_configs::AbstractVector=Any[],
)
    normalized_dir = _normalize_dir(dir)
    path = joinpath(normalized_dir, string(profile, ".toml"))
    fingerprint = _config_fingerprint(normalized_dir, default_config, profile, inherited_configs)

    cached = lock(_CONFIG_CACHE_LOCK) do
        entry = get(_CONFIG_CACHE, fingerprint, nothing)
        if entry === nothing
            _CONFIG_CACHE_STATS[:miss] += 1
            return nothing
        end
        _CONFIG_CACHE_STATS[:hit] += 1
        return deepcopy(entry)
    end

    if cached !== nothing
        return cached
    end

    cfg = deepcopy(default_config)
    for inherited in inherited_configs
        inherited_dict = inherited isa Dict{String, Any} ? inherited : Dict{String, Any}(inherited)
        cfg = deep_merge(cfg, inherited_dict)
    end

    if isfile(path)
        try
            parsed = TOML.parsefile(path)
            cfg = deep_merge(cfg, parsed)
        catch err
            @warn "解析配置失败，使用内置默认值" profile path exception=(err, catch_backtrace())
        end
    else
        if profile != "default"
            @warn "未找到指定配置文件，退回默认参数" profile path
        end
    end

    result = (config=cfg, profile=profile, path=isfile(path) ? path : nothing)

    lock(_CONFIG_CACHE_LOCK) do
        _CONFIG_CACHE[fingerprint] = deepcopy(result)
    end

    return deepcopy(result)
end

end # module ConfigLoader
