using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "residual spine contract" begin
    model = Models.create_model(:PNJL)
    x_state = SVector{5}(-1.84, -1.84, -2.23, 0.5, 0.5)
    mu_vec = SVector{3}(0.0, 0.0, 0.0)
    T_fm = 0.5
    xi = 0.0
    p_num = 8
    t_num = 4

    gap_norm = Models._gap_norm_from_state(model, x_state, mu_vec, T_fm; xi=xi, p_num=p_num, t_num=t_num)

    thermal_nodes = Models.cached_nodes(p_num, t_num)
    params = Models.GapParams(T_fm, thermal_nodes, xi; p_num=p_num, t_num=t_num, model_kind=:PNJL)
    F_core = zeros(5)
    Models.gap_core_residual!(F_core, model, x_state, mu_vec, params)
    expected_norm = sqrt(sum(abs2, F_core))

    @test isapprox(gap_norm, expected_norm; rtol=1e-9, atol=1e-9)

    solver_path = joinpath(PROJECT_ROOT, "src", "models", "solver", "ConstraintSolverCommon.jl")
    source = read(solver_path, String)
    @test occursin("_gap_norm_from_state", source)
    @test !occursin("gap_residual(model", source)
end
