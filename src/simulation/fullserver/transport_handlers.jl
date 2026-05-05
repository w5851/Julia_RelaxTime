@inline function _transport_error_response(status::Int, code::String, message::String)
    return _pnjl_json_response(status, _error_payload(code, message))
end

@inline function _transport_real_option(value, name::String)
    value === nothing && throw(ArgumentError("Missing required parameter: $(name)"))
    val = _to_float64(value)
    isfinite(val) || throw(ArgumentError("$(name) must be finite"))
    return Float64(val)
end

@inline function _transport_int_option(value, fallback::Union{Nothing,Int}, name::String)
    value === nothing && fallback !== nothing && return fallback
    value === nothing && throw(ArgumentError("Missing required parameter: $(name)"))
    parsed = _to_int(value, fallback === nothing ? 0 : fallback)
    parsed > 0 || throw(ArgumentError("$(name) must be positive"))
    return parsed
end

@inline function _json_safe_value(value)
    if value isa AbstractFloat
        return isfinite(value) ? value : nothing
    elseif value isa Integer || value isa Bool || value isa Nothing || value isa AbstractString
        return value
    elseif value isa AbstractDict
        return Dict(string(k) => _json_safe_value(v) for (k, v) in pairs(value))
    elseif value isa NamedTuple
        return Dict(string(k) => _json_safe_value(v) for (k, v) in pairs(value))
    elseif value isa AbstractVector
        return [_json_safe_value(v) for v in value]
    elseif value isa Tuple
        return [_json_safe_value(v) for v in value]
    end
    return value
end

function _transport_tau_payload(tau_value)
    tau_value === nothing && throw(ArgumentError("Missing required parameter: tau"))
    if tau_value isa Number
        tau = _transport_real_option(tau_value, "tau")
        tau >= 0.0 || throw(ArgumentError("tau must be non-negative"))
        return (u=tau, d=tau, s=tau, ubar=tau, dbar=tau, sbar=tau)
    end

    tau_dict = tau_value isa Dict ? tau_value : _to_symbol_dict(tau_value)
    required = (:u, :d, :s, :ubar, :dbar, :sbar)
    values = Dict{Symbol, Float64}()
    for key in required
        haskey(tau_dict, key) || throw(ArgumentError("tau is missing :$(key)"))
        τ = _transport_real_option(tau_dict[key], "tau.$(key)")
        τ >= 0.0 || throw(ArgumentError("tau.$(key) must be non-negative"))
        values[key] = τ
    end
    return (
        u=values[:u],
        d=values[:d],
        s=values[:s],
        ubar=values[:ubar],
        dbar=values[:dbar],
        sbar=values[:sbar],
    )
end

function _transport_integration_config(params_dict::Dict{Symbol,Any})
    haskey(params_dict, :transport) || return nothing
    transport_obj = params_dict[:transport]
    transport_dict = transport_obj isa Dict ? transport_obj : _to_symbol_dict(transport_obj)
    isempty(transport_dict) && return nothing

    workflow = Models.transport_workflow_module()
    config_kwargs = (
        p_nodes=haskey(transport_dict, :p_nodes) ? _transport_int_option(transport_dict[:p_nodes], nothing, "transport.p_nodes") : nothing,
        p_max=haskey(transport_dict, :p_max) ? _transport_real_option(transport_dict[:p_max], "transport.p_max") : nothing,
        cos_nodes=haskey(transport_dict, :cos_nodes) ? _transport_int_option(transport_dict[:cos_nodes], nothing, "transport.cos_nodes") : nothing,
    )
    filtered = (; (k => v for (k, v) in pairs(config_kwargs) if v !== nothing)...)
    isempty(keys(filtered)) && return nothing
    return workflow.TransportIntegrationConfig(; filtered...)
end

@inline function _transport_result_payload(result, T_mev::Float64, mu_mev::Float64)
    eq = result.equilibrium
    tr = result.transport
    return Dict(
        "inputs" => Dict(
            "T_mev" => T_mev,
            "mu_mev" => mu_mev,
            "xi" => result.thermo_params.ξ,
            "tau" => Dict(
                "u" => result.tau.u,
                "d" => result.tau.d,
                "s" => result.tau.s,
                "ubar" => result.tau.ubar,
                "dbar" => result.tau.dbar,
                "sbar" => result.tau.sbar,
            ),
        ),
        "equilibrium" => Dict(
            "converged" => eq.converged,
            "iterations" => eq.iterations,
            "residual_norm" => _json_safe_value(eq.residual_norm),
            "x_state" => _json_safe_value(eq.x_state),
            "mu_vec" => _json_safe_value(eq.mu_vec),
            "masses" => _json_safe_value(eq.masses),
        ),
        "thermo_background" => Dict(
            "pressure" => _json_safe_value(result.thermo_background.pressure),
            "entropy" => _json_safe_value(result.thermo_background.entropy),
            "energy" => _json_safe_value(result.thermo_background.energy),
            "rho_mass" => _json_safe_value(result.thermo_background.rho_mass),
            "c_p" => _json_safe_value(result.thermo_background.c_p),
        ),
        "densities" => Dict(
            "u" => result.densities.u,
            "d" => result.densities.d,
            "s" => result.densities.s,
            "ubar" => result.densities.ubar,
            "dbar" => result.densities.dbar,
            "sbar" => result.densities.sbar,
        ),
        "transport" => Dict(
            "eta" => _json_safe_value(tr.eta),
            "zeta" => _json_safe_value(tr.zeta),
            "sigma" => _json_safe_value(tr.sigma),
            "kappa_BB" => _json_safe_value(tr.kappa_BB),
            "kappa_BQ" => _json_safe_value(tr.kappa_BQ),
            "kappa_BS" => _json_safe_value(tr.kappa_BS),
            "kappa_QQ" => _json_safe_value(tr.kappa_QQ),
            "kappa_QS" => _json_safe_value(tr.kappa_QS),
            "kappa_SS" => _json_safe_value(tr.kappa_SS),
            "lambda" => _json_safe_value(tr.lambda),
            "lorenz_number" => _json_safe_value(tr.lorenz_number),
            "lorentz_legacy" => _json_safe_value(tr.lorentz_legacy),
            "viscous_conductive_coupling_ratio" => _json_safe_value(tr.viscous_conductive_coupling_ratio),
            "prandtl_number" => _json_safe_value(tr.prandtl_number),
            "bulk_to_shear_viscosity_ratio" => _json_safe_value(tr.bulk_to_shear_viscosity_ratio),
        ),
        "reproducibility" => Dict(
            "physics_profile" => result.reproducibility.physics_profile,
            "physics_config_path" => result.reproducibility.physics_config_path,
        ),
    )
end

function handle_transport_point(req::HTTP.Request)
    if req.method != "POST"
        return _transport_error_response(405, "METHOD_NOT_ALLOWED", "Method Not Allowed")
    end

    body = isempty(req.body) ? Dict{Symbol,Any}() : JSON3.read(req.body)
    params_obj = haskey(body, :params) ? body[:params] : body
    params_dict = params_obj isa Dict ? params_obj : _to_symbol_dict(params_obj)

    try
        numeric = Models.default_scan_numeric_options()
        T_mev = _transport_real_option(get(params_dict, :T_mev, get(params_dict, :t_mev, nothing)), "T_mev")
        mu_mev = _transport_real_option(get(params_dict, :mu_mev, get(params_dict, :mu, nothing)), "mu_mev")
        xi = haskey(params_dict, :xi) ? _transport_real_option(params_dict[:xi], "xi") : 0.0
        tau = _transport_tau_payload(get(params_dict, :tau, nothing))
        compute_bulk = _to_bool(get(params_dict, :compute_bulk, false), false)
        p_num = _transport_int_option(get(params_dict, :p_num, nothing), numeric.p_num, "p_num")
        t_num = _transport_int_option(get(params_dict, :t_num, nothing), numeric.t_num, "t_num")
        transport_config = _transport_integration_config(params_dict)

        workflow = Models.transport_workflow_module()
        T_fm = T_mev / ħc_MeV_fm
        mu_fm = mu_mev / ħc_MeV_fm
        equilibrium = workflow.EquilibriumFacade.solve_equilibrium_backend(
            T_fm,
            mu_fm;
            xi=xi,
            solver_backend=:auto,
            p_num=p_num,
            t_num=t_num,
            solver_kwargs=(iterations=30,),
        )
        result = Models.solve_transport_from_equilibrium(
            equilibrium,
            T_fm,
            mu_fm;
            xi=xi,
            tau=tau,
            compute_tau=false,
            compute_bulk=compute_bulk,
            p_num=p_num,
            t_num=t_num,
            transport_config=transport_config,
        )

        return _pnjl_json_response(200, Dict(
            "status" => "ok",
            "result" => _transport_result_payload(result, T_mev, mu_mev),
        ))
    catch e
        if e isa ArgumentError || e isa DomainError
            return _transport_error_response(400, "INVALID_REQUEST", "Invalid transport-point request parameters")
        end
        @error "Transport point request failed" exception=(e, catch_backtrace())
        return _transport_error_response(500, "TRANSPORT_POINT_FAILED", "Transport point solve failed")
    end
end
