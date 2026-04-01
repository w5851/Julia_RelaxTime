using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Wave-B compat routing smoke" begin
    @test isdefined(Main.Models, :solver_migration_status)

    status = Main.Models.solver_migration_status("Models.solve_fixedmu_constraint")
    @test status.status == :active
    @test status.route == :compat_shim
    @test status.unified_entry == "Models.solve_constraint(model, FixedMu(), T; μ_fm=...)"

    m = Main.Models.create_model(:NJL)
    T = 0.5
    seed = Float64.(Main.Models.gap_initial_guess(m, T, SVector{3}(0.0, 0.0, 0.0)))

    via_unified = Main.Models.solve_constraint(m, Main.Models.FixedMu(), T; μ_fm=0.0, seed_guess=seed, p_num=24, t_num=6)
    via_legacy = Main.Models.solve_fixedmu_constraint(m, T, 0.0; seed_guess=seed, p_num=24, t_num=6)

    @test via_unified.converged == via_legacy.converged
    @test isapprox(via_unified.pressure, via_legacy.pressure; rtol=1e-10, atol=1e-12)
    @test isapprox(via_unified.rho_norm, via_legacy.rho_norm; rtol=1e-10, atol=1e-12)
end
