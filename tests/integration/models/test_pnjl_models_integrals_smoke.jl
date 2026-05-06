using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

_models_entry = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !(isdefined(Main, :Models) && isdefined(Main.Models, :omega) && isdefined(Main.Models, :solve_gap))
    Base.include(Main, _models_entry)
end

@testset "Models PNJL integrals (no legacy Thermodynamics include)" begin
    m = Models.create_model(:PNJL; profile="default", physics_profile="default")

    φ = SVector{3, Float64}(0.010, 0.020, 0.030)
    masses = Models.calculate_mass_vec(m, φ)

    Φ = 0.12
    Φbar = 0.34
    T_fm = 0.20
    mu_vec = SVector{3, Float64}(0.0, 0.0, 0.0)

    therm1 = Models.thermal_contribution(m, masses, Φ, Φbar, mu_vec, T_fm; p_num=10, t_num=2, xi=0.3)
    therm2 = Models.thermal_contribution(m, masses, Φ, Φbar, mu_vec, T_fm; p_num=10, t_num=2, xi=0.3)

    @test isfinite(therm1)
    @test therm1 == therm2

    # Node cache should be populated with the current (p_num, t_num, p_max_inv_fm) key shape.
    key = (10, 2, Models.PNJLIntegrals.DEFAULT_THERMAL_P_MAX_INV_FM)
    @test haskey(Models.PNJLIntegrals.NODE_CACHE, key)
end
