# Small numerical parity check, executed by the Actions pilot.
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(REPO_ROOT, "src", "models", "Models.jl"))
using Test

@testset "FixedAsymmetricRho telemetry preserves results" begin
    model = Models.create_model(:PNJL)
    mode = Models.FixedAsymmetricRho(0.35, 0.876, 0.0)
    T = 120.0 / 197.327
    options = (; xi=0.0, p_num=8, t_num=4, iterations=300,
        residual_norm_max=1e-6, semantic_mode=:ground_state, evaluate_all_attempts=true)
    plain = Models.solve_multi(model, mode, T; options...)
    telemetry = Models.SolverWorkTelemetry()
    measured = Models.solve_multi(model, mode, T; work_telemetry=telemetry, options...)
    @test measured.converged == plain.converged
    @test measured.solution ≈ plain.solution atol=1e-10 rtol=1e-10
    @test measured.pressure ≈ plain.pressure atol=1e-10 rtol=1e-10
    @test measured.residual_norm ≈ plain.residual_norm atol=1e-10
    work = Models.solver_work_snapshot(telemetry)
    @test work.equilibrium_requests == work.fixedasymrho_requests == 1
    @test work.fixedrho_requests == 0
    @test work.governed_attempts > 0
    @test work.fixedasymrho_inner_solves >= work.governed_attempts
    @test work.fixedasymrho_residual_evaluations > work.nlsolve_f_calls > 0
    @test work.nlsolve_g_calls > 0
    @test work.trust_region_iterations > 0
    @test work.newton_iterations == 0 # the current inner implementation uses trust region
    @test work.postprocess_residual_calls > 0
    println(work)
end
