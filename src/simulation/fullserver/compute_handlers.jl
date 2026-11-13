function handle_compute(req::HTTP.Request)
    try
        (; p1, p2, m1, m2, m3, m4, theta_star, phi_star) = _parse_compute_request(req)

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
        ]
        return HTTP.Response(200, headers, JSON3.write(response_data))
    catch e
        @error "Computation error" exception = (e, catch_backtrace())

        message_id = _new_message_id()
        response_data = if e isa ArgumentError || e isa DomainError
            merge(Dict(
                "success" => false,
                "data" => nothing,
            ), _error_payload("INVALID_INPUT", "Invalid request parameters"; message_id=message_id))
        else
            merge(Dict(
                "success" => false,
                "data" => nothing,
            ), _error_payload("COMPUTATION_ERROR", "Computation failed"; message_id=message_id))
        end

        status = e isa ArgumentError || e isa DomainError ? 400 : 500
        headers = [
            "Content-Type" => "application/json; charset=utf-8",
        ]
        return HTTP.Response(status, headers, JSON3.write(response_data))
    end
end
