const MODULE_REGISTRY = [
    Dict(
        "id" => "pnjl-gap",
        "name" => "PNJL Gap Single Point",
        "description" => "Run PNJL-compatible single-point solve via Models unified interfaces",
        "invocation_style" => "sync",
        "service_surface" => "point",
        "default_client_surface" => "service",
        "stable_entrypoint" => "Models.solve_pnjl_point",
        "http" => Dict(
            "method" => "POST",
            "path" => "/api/modules/pnjl-gap/run",
        ),
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
        "invocation_style" => "async",
        "service_surface" => "job",
        "default_client_surface" => "service",
        "stable_entrypoint" => "Models.run_scan_pipeline",
        "http" => Dict(
            "create" => Dict("method" => "POST", "path" => "/api/modules/pnjl-scan/jobs"),
            "status" => Dict("method" => "GET", "path" => "/api/modules/pnjl-scan/jobs/{job_id}"),
            "result" => Dict("method" => "GET", "path" => "/api/modules/pnjl-scan/jobs/{job_id}/result"),
            "cancel" => Dict("method" => "POST", "path" => "/api/modules/pnjl-scan/jobs/{job_id}/cancel"),
        ),
        "params_schema" => Dict(
            "kind" => "String (required): tmu | trho",
            "params" => "Dict (optional): scan kwargs, e.g. T_values/mu_values/rho_values/output_path",
            "params.xi" => "Float64 (optional): single xi",
            "params.xi_values" => "Vector{Float64} (optional): explicit xi list",
            "params.xi_grid" => "Dict(start, stop, step) (optional): generate xi list",
            "params.max_retries" => "Int (optional, default 0, max 3)",
        ),
    ),
    Dict(
        "id" => "transport-point",
        "name" => "Transport Point",
        "description" => "Run transport post-processing for a single PNJL equilibrium point via Models stable workflow APIs",
        "invocation_style" => "sync",
        "service_surface" => "point",
        "default_client_surface" => "service",
        "stable_entrypoint" => "Models.solve_transport_from_equilibrium",
        "http" => Dict(
            "method" => "POST",
            "path" => "/api/modules/transport-point/run",
        ),
        "params_schema" => Dict(
            "T_mev" => "Float64 (required)",
            "mu_mev" => "Float64 (required for FixedMu; alias mu)",
            "xi" => "Float64 (optional, default 0.0)",
            "tau" => "Float64 or Dict(u,d,s,ubar,dbar,sbar) (required, finite, non-negative)",
            "compute_bulk" => "Bool (optional, default false)",
            "p_num" => "Int (optional): equilibrium momentum nodes",
            "t_num" => "Int (optional): equilibrium angular nodes",
            "transport.p_nodes" => "Int (optional): transport integration momentum nodes",
            "transport.p_max" => "Float64 (optional): transport momentum cutoff in fm^-1",
            "transport.cos_nodes" => "Int (optional): transport angular nodes",
        ),
    ),
]

@inline _new_message_id() = string(uuid4())

function _error_payload(code::String, message::String; message_id::String=_new_message_id(), extras::AbstractDict=Dict{String, Any}())
    payload = Dict{String, Any}(
        "status" => "error",
        "error_code" => code,
        "error" => message,
        "message_id" => message_id,
    )
    for (k, v) in pairs(extras)
        payload[string(k)] = v
    end
    return payload
end

@inline function _to_float64(x)
    return x isa Number ? Float64(x) : parse(Float64, String(x))
end

@inline function _to_int(x, default::Int)
    x === nothing && return default
    return x isa Integer ? Int(x) : parse(Int, String(x))
end

@inline function _to_bool(x, default::Bool)
    x === nothing && return default
    if x isa Bool
        return x
    end
    if x isa Number
        return Int(x) != 0
    end
    s = lowercase(String(x))
    if s in ("1", "true", "yes", "y", "on")
        return true
    elseif s in ("0", "false", "no", "n", "off")
        return false
    end
    error("Invalid boolean value: $(x)")
end

function _to_symbol_dict(obj)
    data = Dict{Symbol,Any}()
    for (k, v) in pairs(obj)
        key = k isa Symbol ? k : Symbol(string(k))
        data[key] = v
    end
    return data
end
