function handle_modules_list()
    headers = [
        "Content-Type" => "application/json; charset=utf-8",
        "Access-Control-Allow-Origin" => "*",
    ]
    return HTTP.Response(200, headers, JSON3.write(MODULE_REGISTRY))
end

@inline function _pnjl_json_response(status::Int, payload::AbstractDict)
    headers = [
        "Content-Type" => "application/json; charset=utf-8",
        "Access-Control-Allow-Origin" => "*",
    ]
    return HTTP.Response(status, headers, JSON3.write(payload))
end

@inline function _pnjl_error_response(status::Int, code::String, message::String)
    return _pnjl_json_response(status, _error_payload(code, message))
end

function handle_pnjl_single_point(req::HTTP.Request)
    if req.method != "POST"
        return HTTP.Response(405, ["Content-Type" => "text/plain"], "Method Not Allowed")
    end
    body = isempty(req.body) ? Dict{Symbol,Any}() : JSON3.read(String(req.body))
    params_obj = haskey(body, :params) ? body[:params] : body
    params_dict = params_obj isa Dict ? params_obj : _to_symbol_dict(params_obj)

    try
        scan_opts = Models.default_scan_numeric_options()
        solver_result = Models.solve_pnjl_point(
            T_mev=get(params_dict, :T_mev, nothing),
            t_mev=get(params_dict, :t_mev, nothing),
            mu_mev=get(params_dict, :mu_mev, nothing),
            mu=get(params_dict, :mu, nothing),
            rho_target=get(params_dict, :rho_target, nothing),
            xi=get(params_dict, :xi, 0.0),
            p_num=get(params_dict, :p_num, scan_opts.p_num),
            t_num=get(params_dict, :t_num, scan_opts.t_num),
            allow_seed_fallback=get(params_dict, :allow_seed_fallback, true),
        )

        result = Dict(
            "converged" => solver_result.converged,
            "omega" => solver_result.omega,
            "pressure" => solver_result.pressure,
            "rho_norm" => solver_result.rho_norm,
            "entropy" => solver_result.entropy,
            "energy" => solver_result.energy,
            "iterations" => solver_result.iterations,
            "residual_norm" => solver_result.residual_norm,
            "xi" => solver_result.xi,
            "seed_fallback_used" => solver_result.seed_fallback_used,
            "x_state" => solver_result.x_state,
            "mu_vec" => solver_result.mu_vec,
            "masses" => solver_result.masses,
        )

        return _pnjl_json_response(200, Dict("status" => "ok", "result" => result))
    catch e
        @error "PNJL single-point request failed" exception=(e, catch_backtrace())
        if e isa ArgumentError || e isa DomainError
            return _pnjl_error_response(400, "INVALID_REQUEST", "Invalid PNJL request parameters")
        end
        return _pnjl_error_response(500, "PNJL_SINGLE_POINT_FAILED", "PNJL single-point solve failed")
    end
end
