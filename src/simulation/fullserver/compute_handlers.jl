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
