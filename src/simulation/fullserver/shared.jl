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
            "T_mev" => "Float64 (required; alias t_mev)",
            "mu_mev" => "Float64 (required; alias mu)",
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
    Dict(
        "id" => "script-tasks",
        "name" => "Stable Script Task Catalog",
        "description" => "Discover and run curated long-running run_* script workflows as background jobs with safe presets, log tails, and artifact discovery",
        "invocation_style" => "async",
        "service_surface" => "job",
        "default_client_surface" => "service",
        "stable_entrypoint" => "scripts/dev/run_with_sysimage + curated run_* scripts",
        "http" => Dict(
            "catalog" => Dict("method" => "GET", "path" => "/api/modules/script-tasks"),
            "create" => Dict("method" => "POST", "path" => "/api/modules/script-tasks/jobs"),
            "status" => Dict("method" => "GET", "path" => "/api/modules/script-tasks/jobs/{job_id}"),
            "result" => Dict("method" => "GET", "path" => "/api/modules/script-tasks/jobs/{job_id}/result"),
            "cancel" => Dict("method" => "POST", "path" => "/api/modules/script-tasks/jobs/{job_id}/cancel"),
        ),
        "params_schema" => Dict(
            "task_id" => "String (required): id from GET /api/modules/script-tasks",
            "preset" => "String (optional): smoke | canonical | custom; defaults to task default_preset",
            "confirm_heavy" => "Bool (required for heavy presets)",
            "custom_args" => "Vector{String} (required for custom preset): one argv item per element",
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
    x isa Real && !(x isa Bool) && return Float64(x)
    x isa AbstractString && return parse(Float64, x)
    throw(ArgumentError("Expected a number or numeric string"))
end

@inline function _to_int(x, default::Int)
    x === nothing && return default
    x isa Integer && !(x isa Bool) && return Int(x)
    x isa AbstractString && return parse(Int, x)
    throw(ArgumentError("Expected an integer or integer string"))
end

@inline function _to_bool(x, default::Bool)
    x === nothing && return default
    if x isa Bool
        return x
    end
    if x isa Real && isfinite(x) && isinteger(x)
        return !iszero(x)
    end
    x isa AbstractString || throw(ArgumentError("Invalid boolean value"))
    s = lowercase(strip(x))
    if s in ("1", "true", "yes", "y", "on")
        return true
    elseif s in ("0", "false", "no", "n", "off")
        return false
    end
    throw(ArgumentError("Invalid boolean value: $(x)"))
end

function _to_symbol_dict(obj)
    obj isa AbstractDict || obj isa NamedTuple || throw(ArgumentError("Expected a JSON object"))
    data = Dict{Symbol,Any}()
    for (k, v) in pairs(obj)
        k isa Symbol || k isa AbstractString || throw(ArgumentError("Object keys must be strings"))
        key = k isa Symbol ? k : Symbol(k)
        data[key] = v
    end
    return data
end

function _point_request_params(req::HTTP.Request)
    body = isempty(req.body) ? Dict{Symbol,Any}() : _to_symbol_dict(JSON3.read(req.body))
    return haskey(body, :params) ? _to_symbol_dict(body[:params]) : body
end
