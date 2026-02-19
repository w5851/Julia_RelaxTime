const MODULE_REGISTRY = [
    Dict(
        "id" => "pnjl-gap",
        "name" => "PNJL Gap Single Point",
        "description" => "Run PNJL single-point solve via exported PNJL.solve interface",
        "params_schema" => Dict(
            "T_mev" => "Float64 (required)",
            "mu_mev" => "Float64 (required for FixedMu)",
            "rho_target" => "Float64 (optional, when provided uses FixedRho)",
            "xi" => "Float64 (optional, default 0.0)",
            "p_num" => "Int (optional)",
            "t_num" => "Int (optional)",
        ),
    ),
    Dict(
        "id" => "pnjl-scan",
        "name" => "PNJL Scan Job",
        "description" => "Run PNJL T-mu / T-rho scans as background jobs with status/result endpoints",
        "params_schema" => Dict(
            "kind" => "String (required): tmu | trho",
            "params" => "Dict (optional): scan kwargs, e.g. T_values/mu_values/rho_values/output_path",
            "params.xi" => "Float64 (optional): single xi",
            "params.xi_values" => "Vector{Float64} (optional): explicit xi list",
            "params.xi_grid" => "Dict(start, stop, step) (optional): generate xi list",
            "params.max_retries" => "Int (optional, default 0, max 3)",
        ),
    ),
]

@inline function _to_float64(x)
    return x isa Number ? Float64(x) : parse(Float64, String(x))
end

@inline function _to_int(x, default::Int)
    x === nothing && return default
    return x isa Integer ? Int(x) : parse(Int, String(x))
end

function _to_symbol_dict(obj)
    data = Dict{Symbol,Any}()
    for (k, v) in pairs(obj)
        key = k isa Symbol ? k : Symbol(string(k))
        data[key] = v
    end
    return data
end
