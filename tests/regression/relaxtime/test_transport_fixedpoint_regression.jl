using Test
using Printf

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const TRANSPORT_BASELINE_PATH = joinpath(PROJECT_ROOT, "tests", "baselines", "relaxtime", "baseline_transport_fixedpoints_v1.csv")

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const TransportWorkflow = Models.transport_workflow_module()
using .TransportWorkflow

function _point_key(T::Float64, mu::Float64, xi::Float64)
    return @sprintf("%.6f|%.6f|%.6f", T, mu, xi)
end

function _load_baseline(path::String)
    isfile(path) || error("baseline CSV not found: $path")
    rows = NamedTuple[]
    lines = readlines(path)
    isempty(lines) && error("baseline CSV is empty: $path")

    for line in lines[2:end]
        s = strip(line)
        isempty(s) && continue
        startswith(s, "#") && continue
        cols = split(s, ',')
        length(cols) == 6 || error("invalid baseline row: $line")
        T = parse(Float64, strip(cols[1]))
        mu = parse(Float64, strip(cols[2]))
        xi = parse(Float64, strip(cols[3]))
        eta = parse(Float64, strip(cols[4]))
        sigma = parse(Float64, strip(cols[5]))
        zeta = parse(Float64, strip(cols[6]))
        push!(rows, (T=T, mu=mu, xi=xi, eta=eta, sigma=sigma, zeta=zeta))
    end

    return rows
end

function _run_models_case(T, mu, xi, tau, models_solver)
    return TransportWorkflow.solve_gap_and_transport(
        T,
        mu;
        xi=xi,
        tau=tau,
        compute_tau=false,
        compute_bulk=true,
        p_num=8,
        t_num=4,
        transport_config=TransportIntegrationConfig(p_nodes=8, p_max=3.5),
        solver_backend=:models,
        models_solver=models_solver,
        models_residual_norm_max=1e-4,
        seed_state=nothing,
    )
end

@testset "Transport fixedpoint regression" begin
    rows = _load_baseline(TRANSPORT_BASELINE_PATH)
    baseline = Dict(_point_key(row.T, row.mu, row.xi) => row for row in rows)

    tau = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)
    models_solver = Models.NLsolveGapSolver(
        method=:trust_region,
        jacobian=:finite,
        xtol=1e-10,
        ftol=1e-10,
    )

    rtol = 8e-2
    atol = 1e-6

    for row in rows
        point_label = "T=$(row.T), mu=$(row.mu), xi=$(row.xi)"
        @testset "fixedpoint $(point_label)" begin
            res_models = _run_models_case(row.T, row.mu, row.xi, tau, models_solver)

            @test res_models.equilibrium.converged isa Bool
            @test all(isfinite, res_models.masses)
            @test isfinite(res_models.transport.eta)
            @test isfinite(res_models.transport.sigma)
            @test isfinite(res_models.transport.zeta)

            expected_key = _point_key(row.T, row.mu, row.xi)
            @test haskey(baseline, expected_key)
            expected = baseline[expected_key]

            @test isapprox(res_models.transport.eta, expected.eta; rtol=rtol, atol=atol)
            @test isapprox(res_models.transport.sigma, expected.sigma; rtol=rtol, atol=atol)
            @test isapprox(res_models.transport.zeta, expected.zeta; rtol=rtol, atol=atol)
        end
    end
end
