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

        model = Models.create_model(:PNJL)

        solver_result = if haskey(params_dict, :rho_target)
            rho_target = _to_float64(params_dict[:rho_target])
            mu_seed_mev = _to_float64(get(params_dict, :mu_mev, get(params_dict, :mu, 0.0)))
            mu_seed_fm = mu_seed_mev / ħc_MeV_fm
            seed_state = try
                st_seed = Models.solve_gap(model, t_fm, mu_seed_fm; xi=xi, p_num=p_num, t_num=t_num)
                Models.state_vector(st_seed)
            catch
                [0.02, 0.02, 0.03, 0.5, 0.5]
            end

            solver_primary = Models.NLsolveGapSolver(method=:newton, jacobian=:finite, xtol=1e-9, ftol=1e-9)
            solver_secondary = Models.NLsolveGapSolver(method=:trust_region, jacobian=:finite, xtol=1e-9, ftol=1e-9)
            r = Models.solve_fixedrho_constraint(
                model,
                t_fm,
                rho_target;
                seed_guess=seed_state,
                mu0=mu_seed_fm,
                solver_primary=solver_primary,
                solver_secondary=solver_secondary,
                xi=xi,
                p_num=p_num,
                t_num=t_num,
                residual_norm_max=1e-6,
            )
            (
                converged=r.converged,
                omega=r.omega,
                pressure=r.pressure,
                rho_norm=r.rho_norm,
                entropy=r.entropy,
                energy=r.energy,
                iterations=r.iterations,
                residual_norm=r.residual_norm,
                xi=xi,
                x_state=r.x_state,
                mu_vec=r.mu_vec,
                masses=r.masses,
            )
        else
            mu_mev_raw = get(params_dict, :mu_mev, get(params_dict, :mu, nothing))
            mu_mev_raw === nothing && error("Missing required parameter: mu_mev (for FixedMu mode)")
            mu_fm = _to_float64(mu_mev_raw) / ħc_MeV_fm
            st = Models.solve_gap(model, t_fm, mu_fm; xi=xi, p_num=p_num, t_num=t_num)
            x_state = Models.state_vector(st)
            mu_vec = Models.normalize_mu_vec(mu_fm)
            omega_val = Models.omega(model, x_state, t_fm, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
            pressure = -omega_val
            rho = Models.number_densities(model, x_state, t_fm, mu_vec; p_num=p_num, t_num=t_num, xi=xi)
            rho_norm = sum(rho.quark) / 3.0
            masses = Models.calculate_mass_vec(model, x_state)
            (
                converged=true,
                omega=omega_val,
                pressure=pressure,
                rho_norm=rho_norm,
                entropy=NaN,
                energy=NaN,
                iterations=0,
                residual_norm=NaN,
                xi=xi,
                x_state=x_state,
                mu_vec=mu_vec,
                masses=masses,
            )
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
