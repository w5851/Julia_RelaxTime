using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

_models_entry = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !(isdefined(Main, :Models) && isdefined(Main.Models, :omega) && isdefined(Main.Models, :solve_gap))
    Base.include(Main, _models_entry)
end

@testset "Models PNJL params injection (minimal slice)" begin
    # Default profile (from config/models/pnjl/default.toml): N_color=3
    m_default = Main.Models.create_model(:PNJL; profile="default", physics_profile="default")

    # unittest profile overrides N_color in config/models/pnjl/unittest.toml
    m_unittest = Main.Models.create_model(:PNJL; profile="unittest", physics_profile="default")

    @test m_default.consts.N_color == 3
    @test m_unittest.consts.N_color == 7

    # number_densities should scale linearly with N_color since it uses pref=2*N_color.
    st = Main.Models.MeanFieldState(SVector{3, Float64}(0.10, 0.10, 0.10); Phi=0.6, PhiBar=0.6)
    T_fm = 0.20
    mu = 0.0

    d_default = Main.Models.number_densities(m_default, st, T_fm, mu; p_num=10, t_num=2, xi=0.0)
    d_unittest = Main.Models.number_densities(m_unittest, st, T_fm, mu; p_num=10, t_num=2, xi=0.0)

    s_default = sum(d_default.quark) + sum(d_default.antiquark)
    s_unittest = sum(d_unittest.quark) + sum(d_unittest.antiquark)

    @test isfinite(s_default)
    @test isfinite(s_unittest)
    @test s_default != 0.0

    expected_ratio = Float64(m_unittest.consts.N_color) / Float64(m_default.consts.N_color)
    @test isapprox(s_unittest / s_default, expected_ratio; rtol=1e-10, atol=0.0)
end
