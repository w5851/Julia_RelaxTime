using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Wave-D compat cleanup smoke" begin
    @test isdefined(Main.Models, :solver_migration_status)

    status = Main.Models.solver_migration_status("Models.solve_fixedmu_constraint")
    @test status.status == :hard_deprecated
    @test status.route == :compat_shim
    @test occursin("Models.solve_constraint", status.unified_entry)

    m = Main.Models.create_model(:NJL)
    T = 0.5
    seed = Float64.(Main.Models.gap_initial_guess(m, T, SVector{3}(0.0, 0.0, 0.0)))

    err = try
        Main.Models.solve_fixedmu_constraint(m, T, 0.0; seed_guess=seed, p_num=24, t_num=6)
        nothing
    catch exc
        exc
    end

    @test err isa ArgumentError
    @test occursin("hard-deprecated", sprint(showerror, err))
    @test occursin("Models.solve_constraint", sprint(showerror, err))
end
