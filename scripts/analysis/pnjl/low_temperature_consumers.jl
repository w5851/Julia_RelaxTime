#!/usr/bin/env julia
"""Small, consumer-oriented low-temperature checks. Run numerical work on quote-remote."""
module LowTemperatureConsumers

using Dates
using ForwardDiff
using JSON3
using LinearAlgebra
using StaticArrays
using TaylorDiff

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    Base.include(Main, joinpath(ROOT, "src", "models", "Models.jl"))
end
const M = Main.Models
const PI = M.PNJLIntegrals
const TD = M.PNJLChiBTaylorDiff
const HC = 197.327

function fixed_rule(f, a, b, nodes, weights)
    half = (b - a) / 2
    center = (a + b) / 2
    value = weights[1] * f(center + half * nodes[1])
    for i in 2:length(nodes)
        value += weights[i] * f(center + half * nodes[i])
    end
    return half * value
end

function big_gauss_legendre(n)
    nodes, weights = zeros(BigFloat, n), zeros(BigFloat, n)
    for i in 1:cld(n, 2)
        z = cos(big(pi) * (BigFloat(i) - big"0.25") / (BigFloat(n) + big"0.5"))
        converged = false
        for _ in 1:40
            p, previous = one(z), zero(z)
            for k in 1:n
                p, previous = ((2k - 1) * z * p - (k - 1) * previous) / k, p
            end
            derivative = n * (z * p - previous) / (z^2 - 1)
            change = p / derivative
            z -= change
            if abs(change) <= 8eps(z)
                weight = 2 / ((1 - z^2) * derivative^2)
                nodes[i], nodes[n + 1 - i] = -z, z
                weights[i] = weights[n + 1 - i] = weight
                converged = true
                break
            end
        end
        converged || error("BigFloat Gauss-Legendre root $(i)/$(n) did not converge")
    end
    return nodes, weights
end

function repeated_derivative(f, x, order)
    order == 0 && return f(x)
    return ForwardDiff.derivative(y -> repeated_derivative(f, y, order - 1), x)
end

function high_precision_window_pressure(masses, phi, phibar, mu_vec, T, xi;
    order=4, nodes=96, tail=56)
    # Original q-integral, independent of the production comoving-energy identity.
    # BigFloat removes cancellation of its large positive/negative derivatives.
    return setprecision(BigFloat, 192) do
        masses_big = BigFloat.(masses)
        xs, weights = big_gauss_legendre(nodes)
        t = BigFloat(T)
        intervals = Vector{Vector{BigFloat}}()
        for i in eachindex(masses_big)
            m = masses_big[i]
            a = abs(BigFloat(mu_vec[i]))
            energies = [m, max(m, a) + tail * t]
            for width in (-tail, -20, -8, -3, 0, 3, 8, 20, tail)
                push!(energies, max(m, a + width * t))
            end
            bounds = sort!(unique!(sqrt.(max.(energies .^ 2 .- m^2, big"0"))))
            push!(intervals, bounds)
        end
        evaluate = delta -> begin
            total = zero(delta)
            for i in eachindex(masses_big)
                m = masses_big[i]
                f = q -> q^2 * 2t * PI.calculate_log_term(sqrt(q^2 + m^2),
                    BigFloat(mu_vec[i]) + delta, t, BigFloat(phi), BigFloat(phibar))
                bounds = intervals[i]
                for j in 1:(length(bounds) - 1)
                    total += fixed_rule(f, bounds[j], bounds[j + 1], xs, weights)
                end
            end
            total / (2big(pi)^2 * sqrt(1 + BigFloat(xi)))
        end
        [Float64(repeated_derivative(evaluate, BigFloat(0), k)) for k in 0:order]
    end
end

"""Independent fixed-node oracle; windows use E-|mu|, including mu approximately M."""
function window_pressure(masses, phi, phibar, mu_vec, T, xi; nodes=96, tail=44.0)
    T0 = PI._primal_float(T)
    T0 > 0 || throw(ArgumentError("window oracle requires T > 0"))
    xs, ws = Main.GaussLegendre.standard_nodes_weights(nodes)
    total = zero(T + mu_vec[1] + masses[1])
    for i in eachindex(masses)
        m0 = abs(PI._primal_float(masses[i]))
        mu0 = abs(PI._primal_float(mu_vec[i]))
        energies = [m0, max(m0, mu0) + tail * T0]
        for width in (-tail, -16.0, -6.0, -2.0, 0.0, 2.0, 6.0, 16.0, tail)
            push!(energies, max(m0, mu0 + width * T0))
        end
        bounds = sort!(unique!(sqrt.(max.(energies .^ 2 .- m0^2, 0.0))))
        f = q -> q^2 * 2T * PI.calculate_log_term(
            sqrt(q^2 + masses[i]^2), mu_vec[i], T, phi, phibar)
        for j in 1:(length(bounds) - 1)
            total += fixed_rule(f, bounds[j], bounds[j + 1], xs, ws)
        end
    end
    return total / (2pi^2 * sqrt(1 + xi))
end

"""Zero-temperature ideal Fermi-sea pressure and mu derivatives (fm units)."""
function zero_reference(masses, mu_vec, xi)
    return setprecision(BigFloat, 192) do
        values = zeros(BigFloat, 5)
        for (mass, mu) in zip(masses, mu_vec)
            m, a = abs(BigFloat(mass)), abs(BigFloat(mu))
            a > m || continue
            q = sqrt(a^2 - m^2)
            pref = BigFloat(3) / (big(pi)^2 * sqrt(1 + BigFloat(xi)))
            radial = m == 0 ? a^4 / 12 :
                (a * q * (2a^2 - 5m^2) + 3m^4 * asinh(q / m)) / 24
            values[1] += pref * radial
            values[2] += sign(mu) * pref * q^3 / 3
            values[3] += pref * a * q
            values[4] += sign(mu) * pref * (2a^2 - m^2) / q
            values[5] += pref * a * (2a^2 - 3m^2) / q^3
        end
        result = Union{Nothing, Float64}[Float64(value) for value in values]
        if any(abs(mu_vec[i]) == abs(masses[i]) && masses[i] != 0 for i in eachindex(masses))
            result[4] = nothing
            result[5] = nothing
        end
        result
    end
end

function pressure_derivatives(evaluate, mu0, T0)
    delta = TD.taylor_variable(0.0, 4)
    mu_series = mu0 .+ delta
    pressure = evaluate(mu_series, T0)
    chemical = [TD.nth_derivative_from_series(pressure, k) for k in 0:4]
    thermal = if T0 > 0
        series = evaluate(mu0, T0 + TD.taylor_variable(0.0, 2))
        [TD.nth_derivative_from_series(series, k) for k in 1:2]
    else
        Float64[]
    end
    return (; chemical, thermal)
end

function oracle_consumer_check(model, state, T, mu, xi, quadrature; order=4)
    phi = SVector(state[1], state[2], state[3])
    masses = M.calculate_mass_vec(model, phi)
    pressure = (trial, temperature) -> -PI.calculate_log_sum_rs_reduced_adaptive(
        masses, state[4], state[5], trial, temperature, xi;
        rtol=quadrature.thermo_quadrature_rtol,
        atol=quadrature.thermo_quadrature_atol, maxevals=200000)
    adaptive = pressure_derivatives(pressure, mu, T)
    oracle = pressure_derivatives((trial, temperature) ->
        window_pressure(masses, state[4], state[5], trial, temperature, xi), mu, T)
    high = T * HC <= 0.1 && any(mu .> masses) ?
        high_precision_window_pressure(masses, state[4], state[5], mu, T, xi; order) : nothing
    return (; adaptive, oracle, high)
end

function residual_reference(model, state, T, mu, xi; nodes=96)
    pressure = trial -> begin
        phi = SVector(trial[1], trial[2], trial[3])
        masses = M.calculate_mass_vec(model, phi)
        -M.calculate_chiral(model, phi) - M.vacuum_contribution(model, masses) -
            M.polyakov_potential(model, trial[4], trial[5], T) +
            window_pressure(masses, trial[4], trial[5], mu, T, xi; nodes)
    end
    return ForwardDiff.gradient(pressure, SVector{5}(state))
end

function capture(f)
    start = time_ns()
    try
        value = f()
        return (; status="ok", elapsed_seconds=(time_ns() - start) / 1e9, value)
    catch err
        return (; status="error", elapsed_seconds=(time_ns() - start) / 1e9,
            error=sprint(showerror, err), exception=string(typeof(err)))
    end
end

json_safe(x) = x
json_safe(x::AbstractFloat) = isfinite(x) ? x : string(x)
json_safe(x::NamedTuple) = map(json_safe, x)
json_safe(x::AbstractArray) = map(json_safe, x)
json_safe(x::Tuple) = map(json_safe, x)

function write_report(output, rows; stage, label)
    mkpath(dirname(abspath(output)))
    payload = json_safe((; schema="pnjl_low_temperature_consumers_v1",
            created_utc=string(now(UTC)), julia=string(VERSION), host=get(ENV, "COMPUTERNAME", ""),
            stage, label, units="T and mu in fm^-1 internally; pressure derivatives d^nP/dmu_q^n",
            nonfinite_encoding="NaN and +/-Inf are strings, preserving undefined consumer outputs", rows))
    open(output * ".tmp", "w") do io
        JSON3.pretty(io, payload)
    end
    mv(output * ".tmp", output; force=true)
end

function run_fixed(output; label="initial", quick=false)
    rows = Any[]
    masses = SVector(1.5, 1.5, 2.5)
    phi, phibar = 0.1, 0.1
    temperatures = quick ? (0.0, 0.01, 0.1, 1.0, 20.0) : (0.0, 0.01, 0.1, 1.0, 5.0, 20.0, 60.0)
    cases = [(mu=1.25, xi=0.0), (mu=1.5, xi=0.0), (mu=1.8, xi=0.0)]
    quick || append!(cases, [(mu=1.8, xi=-0.3), (mu=1.8, xi=0.3)])
    for point in cases, T_MeV in temperatures
        T0 = T_MeV / HC
        mu0 = SVector(point.mu, point.mu, point.mu)
        zero_limit = zero_reference(masses, mu0, point.xi)
        methods = Any[]
        for rtol in (1e-6, 1e-8, 1e-10)
            result = capture() do
                evaluate = (mu, T) -> -PI.calculate_log_sum_rs_reduced_adaptive(
                    masses, phi, phibar, mu, T, point.xi; rtol, atol=rtol * 1e-3, maxevals=200000)
                if T0 == 0 && point.mu in masses
                    (; chemical=[evaluate(mu0, T0)], thermal=Float64[])
                else
                    pressure_derivatives(evaluate, mu0, T0)
                end
            end
            push!(methods, (; method="adaptive", rtol, result))
        end
        if T0 > 0
            for nodes in (64, 96)
                result = capture() do
                    evaluate = (mu, T) -> window_pressure(masses, phi, phibar, mu, T, point.xi; nodes)
                    pressure_derivatives(evaluate, mu0, T0)
                end
                push!(methods, (; method="energy_window_oracle", nodes, result))
            end
            if T_MeV <= 0.1 && point.mu > minimum(masses)
                result = capture(() -> high_precision_window_pressure(masses, phi, phibar,
                    mu0, T0, point.xi))
                push!(methods, (; method="bigfloat_window_oracle", bits=192, nodes=96, result))
            end
        end
        push!(rows, (; T_MeV, mu_q_inv_fm=point.mu, xi=point.xi, masses_inv_fm=collect(masses),
            phi, phibar, at_fermi_threshold=point.mu in masses, zero_limit, methods))
        write_report(output, rows; stage="fixed", label)
        println("fixed T=", T_MeV, " mu=", point.mu, " xi=", point.xi,
            " errors=", count(x -> x.result.status != "ok", methods))
        flush(stdout)
    end
    return rows
end

function quadrature_options(rtol)
    return (; thermo_quadrature_policy=:rs_reduced_adaptive,
        thermo_quadrature_rtol=rtol, thermo_quadrature_atol=rtol * 1e-3,
        thermo_quadrature_maxevals=200000)
end

function consumer_values(model, state, T, mu, xi, quadrature; residual_tol=1e-9)
    context = TD.build_taylor_base_context(T, mu; model, xi, p_num=24, t_num=8,
        base_state=state, base_residual_tol=residual_tol, polish_max_iterations=8,
        thermo_quadrature_kwargs=quadrature)
    common = (; model, xi, p_num=24, t_num=8, base_context=context,
        series_residual_tol=residual_tol, thermo_quadrature_kwargs=quadrature)
    mass_result = capture(() -> M.mass_derivatives(T, mu[1]; order=1, common...))
    thermo_result = capture(() -> M.thermo_derivatives(T, mu[1]; common...))
    pressure_mu = capture() do
        series = TD.pressure_series_direction(T, mu, SVector(1.0, 1.0, 1.0);
            order=4, common...)
        [TD.nth_derivative_from_series(series, k) for k in 0:4]
    end
    chi_result = capture(() -> M.chi2_B(T, 3mu[1]; common...))
    chi4_result = capture(() -> M.chi4_B(T, 3mu[1]; common...))
    bulk_result = capture(() -> M.bulk_viscosity_coefficients(T, mu[1]; model, xi,
        p_num=24, t_num=8, base_state=context.x0, base_residual_tol=residual_tol,
        series_residual_tol=residual_tol, thermo_quadrature_kwargs=quadrature))
    density = M.number_densities(model, context.x0, T, mu; xi, quadrature...)
    masses = M.calculate_mass_vec(model, SVector(context.x0[1], context.x0[2], context.x0[3]))
    zero_pressure = M.model_pressure(model, context.x0, mu, 0.0; xi, quadrature...)
    zero_density = M.number_densities(model, context.x0, 0.0, mu; xi, quadrature...)
    return (; state=collect(context.x0), masses=collect(masses),
        fermi_offsets=collect((mu .- masses) ./ T), residual_norm=context.base_residual_norm,
        jacobian_condition=cond(context.J0), density,
        zero_fixed_state=(; pressure=zero_pressure, density=zero_density),
        mass_result, thermo_result, pressure_mu, chi_result, chi4_result, bulk_result)
end

function run_self_consistent(output; label="consumers", quick=false)
    rows = Any[]
    model = M.create_model(:PNJL)
    points = quick ? [(mu=1.25, xi=0.0), (mu=2.2, xi=0.0)] :
        [(mu=1.25, xi=0.0), (mu=1.8, xi=0.0), (mu=2.2, xi=0.0),
         (mu=2.2, xi=-0.3), (mu=2.2, xi=0.3)]
    temperatures = quick ? (60.0, 5.0, 0.1, 0.01) : (60.0, 20.0, 5.0, 1.0, 0.1, 0.01)
    for point in points
        seed = nothing
        for T_MeV in temperatures
            T = T_MeV / HC
            mu = SVector(point.mu, point.mu, point.mu)
            solved = capture() do
                q = quadrature_options(1e-8)
                state = M.solve_gap(model, T, mu; xi=point.xi,
                    initial_guess=seed, residual_norm_max=1e-9,
                    solver=M.NLsolveGapSolver(ftol=1e-11, xtol=1e-12), q...)
                collect(M.state_vector(state))
            end
            methods = Any[]
            if solved.status == "ok"
                seed = solved.value
                for rtol in (1e-6, 1e-8, 1e-10)
                    result = capture(() -> consumer_values(model, seed, T, mu, point.xi,
                        quadrature_options(rtol)))
                    push!(methods, (; rtol, result))
                end
            end
            oracle = if solved.status == "ok" && point.mu >= 1.8 && T_MeV in (0.01, 5.0, 60.0)
                capture() do
                    state = methods[end].result.status == "ok" ? methods[end].result.value.state : seed
                    (; fixed=oracle_consumer_check(model, state, T, mu, point.xi, quadrature_options(1e-10)),
                       residual=collect(residual_reference(model, state, T, mu, point.xi)))
                end
            else
                nothing
            end
            push!(rows, (; T_MeV, mu_q_inv_fm=point.mu, xi=point.xi, solved, methods, oracle))
            write_report(output, rows; stage="self_consistent", label)
            println("consumer T=", T_MeV, " mu=", point.mu, " xi=", point.xi,
                " solve=", solved.status, " checked=", length(methods))
            flush(stdout)
        end
    end
    return rows
end

function timed_best(f; repeats=3)
    f() # compile before measuring
    seconds = Inf
    value = nothing
    for _ in 1:repeats
        start = time_ns()
        value = f()
        seconds = min(seconds, (time_ns() - start) / 1e9)
    end
    return (; seconds, value)
end

function tensor_check(model, row)
    state = SVector{5, Float64}(row.methods[end].result.value.state)
    T = Float64(row.T_MeV) / HC
    mu = Float64(row.mu_q_inv_fm)
    xi = Float64(row.xi)
    masses = M.calculate_mass_vec(model, SVector(state[1], state[2], state[3]))
    q = quadrature_options(1e-10)
    adaptive = timed_best() do
        pressure_derivatives((trial, temperature) -> M.thermal_contribution(
            model, masses, state[4], state[5], trial, temperature; xi, q...) * -1,
            SVector(mu, mu, mu), T)
    end
    fixed = map((24, 96, 384)) do nodes
        result = timed_best() do
            pressure_derivatives((trial, temperature) -> M.thermal_contribution(
                model, masses, state[4], state[5], trial, temperature;
                xi, p_num=nodes, t_num=16) * -1, SVector(mu, mu, mu), T)
        end
        (; nodes, angular_nodes=16, result)
    end
    return (; T_MeV=Float64(row.T_MeV), mu_q_inv_fm=mu, xi,
        state=collect(state), masses=collect(masses), adaptive, fixed)
end

function point_observables(model, state, T, mu, xi, quadrature)
    evaluate = (trial_mu, trial_T) -> M.model_pressure(
        model, state, trial_mu, trial_T; xi, quadrature...)
    p = evaluate(mu, T)
    # At a stationary mean field these are total first derivatives.
    density = M.number_densities(model, state, T, mu; xi, quadrature...)
    entropy = ForwardDiff.derivative(t -> evaluate(mu, t), T)
    masses = M.calculate_mass_vec(model, SVector(state[1], state[2], state[3]))
    return vcat(p, sum(density.quark - density.antiquark), entropy, collect(masses))
end

function finite_difference_consumers(model, row)
    T = Float64(row.T_MeV) / HC
    mu = Float64(row.mu_q_inv_fm)
    xi = Float64(row.xi)
    mu_vec = SVector(mu, mu, mu)
    state = SVector{5, Float64}(row.methods[end].result.value.state)
    q = quadrature_options(1e-10)
    strict = consumer_values(model, state, T, mu_vec, xi, q; residual_tol=1e-13)
    center = SVector{5, Float64}(strict.state)
    fd = Any[]
    for direction in (:mu, :T), fraction in (0.01, 0.003, 0.001)
        # Chemical steps follow the branch's physical scale; thermal steps stay
        # strictly positive even at T=0.01 MeV. No global phase reselection.
        h = direction === :mu ? 0.1 * fraction : T * fraction
        values = Dict{Int, Vector{Float64}}()
        residuals = Float64[]
        drift = 0.0
        offsets = direction === :mu ? (-3, -2, -1, 1, 2, 3) : (-2, -1, 1, 2)
        for k in offsets
            trial_mu = direction === :mu ? mu_vec .+ k * h : mu_vec
            trial_T = direction === :T ? T + k * h : T
            solved = M.solve_gap(model, trial_T, trial_mu; xi, initial_guess=center,
                residual_norm_max=1e-12,
                solver=M.NLsolveGapSolver(ftol=1e-14, xtol=1e-14), q...)
            trial_state = SVector{5, Float64}(M.state_vector(solved))
            TD._assert_same_branch(model, center, trial_state)
            drift = max(drift, maximum(abs, trial_state - center))
            _, residual = TD._primal_residual(model, trial_state, trial_T, trial_mu,
                xi, 24, 8; thermo_quadrature_kwargs=q)
            push!(residuals, residual)
            values[k] = point_observables(model, trial_state, trial_T, trial_mu, xi, q)
        end
        derivative = (values[-2] - 8values[-1] + 8values[1] - values[2]) / (12h)
        # Third derivative of the independently integrated density = P_mu^4.
        fourth_pressure = direction === :mu ?
            (values[2][2] - 2values[1][2] + 2values[-1][2] - values[-2][2]) / (2h^3) : nothing
        fourth_pressure_order4 = direction === :mu ?
            (values[-3][2] - 8values[-2][2] + 13values[-1][2] - 13values[1][2] +
             8values[2][2] - values[3][2]) / (8h^3) : nothing
        push!(fd, (; direction=string(direction), h, first_derivatives=derivative,
            fourth_pressure, fourth_pressure_order4,
            max_state_drift=drift, max_residual=maximum(residuals)))
    end
    return (; T_MeV=Float64(row.T_MeV), mu_q_inv_fm=mu, xi,
        original=row.methods[end].result.value, strict,
        observable_order=("pressure", "net_quark_density", "entropy", "M_u", "M_d", "M_s"), fd)
end

function run_audit(output, consumer_input; label="audit")
    input = JSON3.read(read(consumer_input, String))
    model = M.create_model(:PNJL)
    rows = Any[]
    select(T, mu, xi=0.0) = only(filter(r -> r.T_MeV == T &&
        r.mu_q_inv_fm == mu && r.xi == xi, input.rows))
    for (T, mu, xi) in ((60.0, 2.2, 0.0), (5.0, 2.2, 0.0),
        (0.01, 2.2, 0.0), (0.01, 2.2, 0.3))
        push!(rows, (; check="tensor", result=capture(() -> tensor_check(model, select(T, mu, xi)))))
        write_report(output, rows; stage="audit", label)
        println("audit tensor T=", T, " mu=", mu, " xi=", xi)
        flush(stdout)
    end
    for (T, mu, xi) in ((60.0, 2.2, 0.0), (5.0, 2.2, 0.0),
        (0.01, 2.2, 0.0), (0.01, 1.8, 0.0))
        push!(rows, (; check="total_derivatives",
            result=capture(() -> finite_difference_consumers(model, select(T, mu, xi)))))
        write_report(output, rows; stage="audit", label)
        println("audit total derivatives T=", T, " mu=", mu)
        flush(stdout)
    end
    for (T, ratio) in ((0.01 / HC, 0.3 / (0.01 / HC)), (0.3 / 63.9, 63.9),
        (0.3 / 64.1, 64.1))
        masses = SVector(1.5, 1.5, 2.5)
        mu = SVector(1.8, 1.8, 1.8)
        result = capture() do
            adaptive = pressure_derivatives((trial, temperature) ->
                -PI.calculate_log_sum_rs_reduced_adaptive(masses, 0.12, 0.28,
                    trial, temperature, 0.3; rtol=1e-10, atol=1e-13, maxevals=200000), mu, T)
            medium = high_precision_window_pressure(masses, 0.12, 0.28, mu, T, 0.3)
            strict = high_precision_window_pressure(masses, 0.12, 0.28, mu, T, 0.3;
                nodes=128, tail=72)
            (; T_MeV=T * HC, mu_q_inv_fm=1.8, xi=0.3, ratio, adaptive,
                oracle_96_56T=medium, oracle_128_72T=strict)
        end
        push!(rows, (; check="oracle_self_convergence", result))
        write_report(output, rows; stage="audit", label)
        println("audit oracle ratio=", ratio)
        flush(stdout)
    end
    return rows
end

function main(args=ARGS)
    length(args) >= 1 || throw(ArgumentError("usage: low_temperature_consumers.jl OUTPUT.json [LABEL] [--quick]"))
    if "--audit" in args
        length(args) == 4 || throw(ArgumentError("usage: low_temperature_consumers.jl OUTPUT.json LABEL --audit CONSUMERS.json"))
        return run_audit(args[1], args[4]; label=args[2])
    end
    run = "--consumers" in args ? run_self_consistent : run_fixed
    run(args[1]; label=length(args) >= 2 ? args[2] : "initial", quick="--quick" in args)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    LowTemperatureConsumers.main()
end
