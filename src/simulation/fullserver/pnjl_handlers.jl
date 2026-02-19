function handle_modules_list()
    headers = [
        "Content-Type" => "application/json; charset=utf-8",
        "Access-Control-Allow-Origin" => "*",
    ]
    return HTTP.Response(200, headers, JSON3.write(MODULE_REGISTRY))
end

function handle_pnjl_single_point(req::HTTP.Request)
    if req.method != "POST"
        return HTTP.Response(405, ["Content-Type" => "text/plain"], "Method Not Allowed")
    end
    body = isempty(req.body) ? Dict{Symbol,Any}() : JSON3.read(String(req.body))
    params_obj = haskey(body, :params) ? body[:params] : body
    params_dict = params_obj isa Dict ? params_obj : _to_symbol_dict(params_obj)

    try
        t_mev_raw = get(params_dict, :T_mev, get(params_dict, :t_mev, nothing))
        t_mev_raw === nothing && error("Missing required parameter: T_mev")

        xi = _to_float64(get(params_dict, :xi, 0.0))
        t_mev = _to_float64(t_mev_raw)
        p_num = _to_int(get(params_dict, :p_num, nothing), 24)
        t_num = _to_int(get(params_dict, :t_num, nothing), 12)

        t_fm = t_mev / ħc_MeV_fm

        solver_result = if haskey(params_dict, :rho_target)
            rho_target = _to_float64(params_dict[:rho_target])
            PNJL.solve(PNJL.FixedRho(rho_target), t_fm; xi=xi, p_num=p_num, t_num=t_num)
        else
            mu_mev_raw = get(params_dict, :mu_mev, get(params_dict, :mu, nothing))
            mu_mev_raw === nothing && error("Missing required parameter: mu_mev (for FixedMu mode)")
            mu_fm = _to_float64(mu_mev_raw) / ħc_MeV_fm
            PNJL.solve(PNJL.FixedMu(), t_fm, mu_fm; xi=xi, p_num=p_num, t_num=t_num)
        end

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
            "x_state" => collect(solver_result.x_state),
            "mu_vec" => collect(solver_result.mu_vec),
            "masses" => collect(solver_result.masses),
        )

        headers = [
            "Content-Type" => "application/json; charset=utf-8",
            "Access-Control-Allow-Origin" => "*",
        ]
        return HTTP.Response(200, headers, JSON3.write(Dict("status" => "ok", "result" => result)))
    catch e
        error_msg = sprint(showerror, e, catch_backtrace())
        headers = [
            "Content-Type" => "application/json; charset=utf-8",
            "Access-Control-Allow-Origin" => "*",
        ]
        payload = Dict("status" => "error", "error" => error_msg)
        return HTTP.Response(500, headers, JSON3.write(payload))
    end
end
