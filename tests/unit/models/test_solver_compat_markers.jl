using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "solver compatibility markers" begin
    @test isdefined(Models, :solver_migration_map)

    mapping = Models.solver_migration_map()
    @test mapping isa AbstractVector
    @test any(e -> e.old_entry == "Models.solve_fixedmu_constraint" && e.new_entry == "Models.solve_constraint(model, FixedMu(), T; μ_fm=...)", mapping)
    @test any(e -> e.old_entry == "Models.solve_fixedrho_constraint" && e.new_entry == "Models.solve_constraint(model, FixedRho(...), T)", mapping)
    @test all(e -> e.status == :hard_deprecated, mapping)

    @test isdefined(Models, :solve_fixedmu_constraint)
    @test isdefined(Models, :solve_fixedrho_constraint)
    @test isdefined(Models, :solve_fixedentropy_constraint)
    @test isdefined(Models, :solve_fixedsigma_constraint)
    @test isdefined(Models, :solve_fixedasymrho_constraint)

    @test isdefined(Models, :solve_constraint)

    m = Models.create_model(:NJL)
    T = 0.5
    seed = Float64.(Models.gap_initial_guess(m, T, SVector{3}(0.0, 0.0, 0.0)))

    via_unified = Models.solve_constraint(m, Models.FixedMu(), T; μ_fm=0.0, seed_guess=seed, p_num=24, t_num=6)
    legacy_err = try
        Models.solve_fixedmu_constraint(m, T, 0.0; seed_guess=seed, p_num=24, t_num=6)
        nothing
    catch exc
        exc
    end

    @test via_unified.converged
    @test isfinite(via_unified.pressure)
    @test isfinite(via_unified.rho_norm)
    @test legacy_err isa ArgumentError
    @test occursin("hard-deprecated", sprint(showerror, legacy_err))
end
