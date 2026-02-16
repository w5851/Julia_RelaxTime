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

using TOML

export deep_merge, load_config

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

"""load_config(dir, default_config; profile="default")

从 `dir/<profile>.toml` 读取配置，并与 `default_config` deep-merge。

返回 `(config, profile, path)`，其中 `path` 若文件不存在则为 `nothing`。
"""
function load_config(dir::AbstractString, default_config::Dict{String, Any}; profile::String="default")
    path = joinpath(dir, string(profile, ".toml"))

    cfg = deepcopy(default_config)
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

    return (config=cfg, profile=profile, path=isfile(path) ? path : nothing)
end

end # module ConfigLoader
