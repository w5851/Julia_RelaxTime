using Test

const _MODELS_ENTRY = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_MODELS_ENTRY)
end
using .Models

using StaticArrays

@testset "Models omega: NJL2" begin
    m = Models.create_model(:NJL2; profile="default")

    p_num = 24
    t_num = 6

    T = 0.8
    φ = @SVector [0.01, 0.02, 0.0]

    mu_vec_a = @SVector [0.0, 0.0, 0.0]
    mu_vec_b = @SVector [0.0, 0.0, 1.23]

    comp_a = Models.omega_components(m, φ, T, mu_vec_a; p_num=p_num, t_num=t_num, xi=0.0)
    comp_b = Models.omega_components(m, φ, T, mu_vec_b; p_num=p_num, t_num=t_num, xi=0.0)

    @test all(isfinite, (comp_a.chi, comp_a.poly, comp_a.vac, comp_a.therm, comp_a.omega))
    @test comp_a.omega ≈ (comp_a.chi + comp_a.poly + comp_a.vac + comp_a.therm)
    @test comp_a.poly == 0.0

    # NJL2 ignores μ_s (third flavor placeholder)
    @test comp_a.omega == comp_b.omega

    st = Models.MeanFieldState(φ; Phi=0.123, PhiBar=0.456)
    comp_s = Models.omega_components(m, st, T, mu_vec_a; p_num=p_num, t_num=t_num, xi=0.0)
    @test comp_s.omega == comp_a.omega

    x5 = SVector{5, Float64}(φ[1], φ[2], φ[3], 0.3, 0.4)
    comp_5 = Models.omega_components(m, x5, T, mu_vec_a; p_num=p_num, t_num=t_num, xi=0.0)
    @test comp_5.omega == comp_a.omega

    dens = Models.number_densities(m, st, T, mu_vec_a; p_num=16, t_num=4, xi=0.0)
    @test dens.quark[3] == 0.0
    @test dens.antiquark[3] == 0.0
end
