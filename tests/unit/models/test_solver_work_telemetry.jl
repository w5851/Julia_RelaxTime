using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Constants_PNJL)
    include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
end
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "request-scoped solver work telemetry" begin
    telemetry = Models.SolverWorkTelemetry()
    @test telemetry.equilibrium_requests == 0
    @test telemetry.fixedrho_requests == 0

    Models.record_solver_request!(telemetry; fixedrho=true)
    Models.record_governed_attempt!(telemetry; origin=:primary)
    Models.record_governed_attempt!(telemetry; origin=:method_rescue)
    Models.record_postprocess_residual!(telemetry)
    Models.record_attempt_outcome!(telemetry; converged=false)
    Models.record_solver_exception!(telemetry)
    Models.record_scan_retry!(telemetry)
    Models.record_nlsolve_work!(telemetry,
        (f_calls=4, g_calls=2, iterations=3),
        :newton,
    )

    snapshot = Models.solver_work_snapshot(telemetry)
    @test snapshot.equilibrium_requests == 1
    @test snapshot.fixedrho_requests == 1
    @test snapshot.governed_attempts == 2
    @test snapshot.root_fallbacks == 1
    @test snapshot.governed_rescues == 1
    @test snapshot.nlsolve_f_calls == 4
    @test snapshot.nlsolve_g_calls == 2
    @test snapshot.newton_iterations == 3
    @test snapshot.postprocess_residual_calls == 1
    @test snapshot.nonconverged_attempts == 1
    @test snapshot.exceptions == 1
    @test snapshot.scan_retries == 1

    Models.reset_solver_work!(telemetry)
    @test all(value == 0 for value in values(Models.solver_work_snapshot(telemetry)))
    @test_throws ArgumentError Models.record_solver_request!(Ref(0); fixedrho=false)
end

@testset "FixedAsymmetricRho inner work and residual instrumentation" begin
    telemetry = Models.SolverWorkTelemetry()
    Models.record_solver_request!(telemetry; fixedasymrho=true)
    Models.record_fixedasymrho_inner_solve!(telemetry)
    Models.record_fixedasymrho_inner_solve!(telemetry)
    residual! = (F, x) -> (F .= x .^ 2 .- 1)
    @test Models._instrument_fixedasymrho_residual(residual!, nothing) === residual!
    instrumented! = Models._instrument_fixedasymrho_residual(residual!, telemetry)
    output = zeros(2)
    instrumented!(output, [2.0, 3.0])
    @test output == [3.0, 8.0]
    instrumented!(output, [1.0, 2.0])
    @test output == [0.0, 3.0]
    snapshot = Models.solver_work_snapshot(telemetry)
    @test snapshot.equilibrium_requests == snapshot.fixedasymrho_requests == 1
    @test snapshot.fixedrho_requests == 0
    @test snapshot.fixedasymrho_inner_solves == 2
    @test snapshot.fixedasymrho_residual_evaluations == 2
    @test_throws ArgumentError Models._instrument_fixedasymrho_residual(residual!, Ref(0))
    Models.reset_solver_work!(telemetry)
    @test all(iszero, values(Models.solver_work_snapshot(telemetry)))
end

@testset "SolverResult v1 remains unchanged with telemetry" begin
    @test Models.SOLVER_CONTRACT_VERSION_V1 == :v1
    @test :contract_version in Models.SOLVER_RESULT_REQUIRED_FIELDS
    @test !(:work_telemetry in Models.SOLVER_RESULT_REQUIRED_FIELDS)
end
