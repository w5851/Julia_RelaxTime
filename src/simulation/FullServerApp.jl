module FullServerApp

using HTTP
using JSON3
using LinearAlgebra

include(joinpath(@__DIR__, "MomentumMapping.jl"))
using .MomentumMapping

include(joinpath(@__DIR__, "..", "pnjl", "PNJL.jl"))
using .PNJL

const _CONSTANTS_PATH = normpath(joinpath(@__DIR__, "..", "Constants_PNJL.jl"))
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, _CONSTANTS_PATH)
end
using Main.Constants_PNJL: ħc_MeV_fm

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
]

@inline function _to_float64(x)
    return x isa Number ? Float64(x) : parse(Float64, String(x))
end

@inline function _to_int(x, default::Int)
    x === nothing && return default
    return x isa Integer ? Int(x) : parse(Int, String(x))
end

function handle_compute(req::HTTP.Request)
    try
        body = JSON3.read(String(req.body))

        p1 = [Float64(body.p1x), Float64(body.p1y), Float64(body.p1z)]
        p2 = [Float64(body.p2x), Float64(body.p2y), Float64(body.p2z)]
        m1 = Float64(body.m1)
        m2 = Float64(body.m2)
        m3 = Float64(body.m3)
        m4 = Float64(body.m4)

        theta_star = haskey(body, :theta_star) ? Float64(body.theta_star) : π / 4
        phi_star = haskey(body, :phi_star) ? Float64(body.phi_star) : π / 6

        if any(isnan.([p1; p2; m1; m2; m3; m4; theta_star; phi_star]))
            return HTTP.Response(
                400,
                ["Content-Type" => "application/json"],
                JSON3.write(Dict("success" => false, "error" => "Invalid input: NaN detected")),
            )
        end

        result = calculate_outgoing_momenta(p1, p2, m1, m2, m3, m4, theta_star, phi_star)
        is_valid, checks = validate_kinematics(result, m1, m2, m3, m4, tol = 1e-9)

        response_data = Dict(
            "success" => true,
            "data" => Dict(
                "ellipsoid" => Dict(
                    "center" => result.ellipsoid.center,
                    "axes_directions" => [result.ellipsoid.axes_directions[:, i] for i in 1:3],
                    "half_lengths" => result.ellipsoid.half_lengths,
                ),
                "momenta" => Dict(
                    "p1" => result.p1_lab,
                    "p2" => result.p2_lab,
                    "p3" => result.p3_lab,
                    "p4" => result.p4_lab,
                    "E1" => result.E1,
                    "E2" => result.E2,
                    "E3" => result.E3,
                    "E4" => result.E4,
                ),
                "physics" => Dict(
                    "s" => result.s,
                    "sqrt_s" => sqrt(result.s),
                    "p_star" => result.p_star,
                    "beta" => norm(result.beta),
                    "beta_vector" => result.beta,
                    "gamma" => result.gamma,
                    "theta_star" => result.theta_star,
                    "phi_star" => result.phi_star,
                ),
                "validation" => Dict(
                    "is_valid" => is_valid,
                    "energy_conservation" => checks["energy_conservation"][1],
                    "momentum_conservation" => checks["momentum_conservation"][1],
                ),
            ),
            "error" => nothing,
        )

        headers = [
            "Content-Type" => "application/json; charset=utf-8",
            "Access-Control-Allow-Origin" => "*",
        ]
        return HTTP.Response(200, headers, JSON3.write(response_data))
    catch e
        error_msg = sprint(showerror, e, catch_backtrace())
        @error "Computation error" exception = (e, catch_backtrace())

        response_data = Dict("success" => false, "data" => nothing, "error" => error_msg)
        headers = [
            "Content-Type" => "application/json; charset=utf-8",
            "Access-Control-Allow-Origin" => "*",
        ]
        return HTTP.Response(400, headers, JSON3.write(response_data))
    end
end

function cors_middleware(handler)
    return function (req::HTTP.Request)
        if req.method == "OPTIONS"
            return HTTP.Response(
                200,
                [
                    "Access-Control-Allow-Origin" => "*",
                    "Access-Control-Allow-Methods" => "POST, GET, OPTIONS",
                    "Access-Control-Allow-Headers" => "Content-Type",
                    "Access-Control-Max-Age" => "86400",
                ],
            )
        end

        resp = handler(req)
        HTTP.setheader(resp, "Access-Control-Allow-Origin" => "*")
        HTTP.setheader(resp, "Access-Control-Allow-Methods", "POST, GET, OPTIONS")
        HTTP.setheader(resp, "Access-Control-Allow-Headers", "Content-Type")
        return resp
    end
end

function serve_static_file(path::String, repo_root::String)
    if contains(path, "..") || contains(path, "\\")
        return HTTP.Response(403, "Forbidden")
    end

    clean_path = startswith(path, "/") ? path[2:end] : path
    file_path = if path == "/" || path == ""
        joinpath(repo_root, "web", "index.html")
    else
        joinpath(repo_root, "web", clean_path)
    end

    if !isfile(file_path)
        return HTTP.Response(404, "Not Found: $path")
    end

    content_type = if endswith(file_path, ".html")
        "text/html; charset=utf-8"
    elseif endswith(file_path, ".css")
        "text/css; charset=utf-8"
    elseif endswith(file_path, ".js")
        "application/javascript; charset=utf-8"
    elseif endswith(file_path, ".json")
        "application/json; charset=utf-8"
    else
        "application/octet-stream"
    end

    try
        content = read(file_path)
        return HTTP.Response(200, ["Content-Type" => content_type], body = content)
    catch e
        @error "Error serving file" file = file_path exception = e
        return HTTP.Response(500, "Internal Server Error")
    end
end

function handle_modules_list()
    headers = [
        "Content-Type" => "application/json; charset=utf-8",
        "Access-Control-Allow-Origin" => "*",
    ]
    return HTTP.Response(200, headers, JSON3.write(MODULE_REGISTRY))
end

function _to_symbol_dict(obj)
    data = Dict{Symbol,Any}()
    for (k, v) in pairs(obj)
        key = k isa Symbol ? k : Symbol(string(k))
        data[key] = v
    end
    return data
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

function route_request(req::HTTP.Request, repo_root::String)
    path = String(HTTP.URIs.unescapeuri(req.target))

    if path == "/health"
        return HTTP.Response(200, ["Content-Type" => "text/plain"], "OK")
    elseif path == "/compute"
        return handle_compute(req)
    elseif path == "/api/modules" && req.method == "GET"
        return handle_modules_list()
    elseif path == "/api/modules/pnjl-gap/run"
        return handle_pnjl_single_point(req)
    else
        static_path = String(split(path, '?')[1])
        return serve_static_file(static_path, repo_root)
    end
end

function build_app(repo_root::String)
    return cors_middleware(req -> route_request(req, repo_root))
end

end
