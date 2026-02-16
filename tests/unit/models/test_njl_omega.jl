using Test

const _MODELS_ENTRY = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_MODELS_ENTRY)
end
using .Models

using StaticArrays

@testset "Models omega: NJL" begin
    m = Models.create_model(:NJL; profile="default")

    # Use small node counts to keep unit test fast
    p_num = 24
    t_num = 6

    T = 0.8
    mu_vec = @SVector [0.0, 0.0, 0.0]

    φ = @SVector [0.01, 0.02, 0.03]

    # x_state length=3 path (implicit Φ=Φbar=1)
    comp3 = Models.omega_components(m, φ, T, mu_vec; p_num=p_num, t_num=t_num, xi=0.0)
    @test all(isfinite, (comp3.chi, comp3.poly, comp3.vac, comp3.therm, comp3.omega))
    @test comp3.omega ≈ (comp3.chi + comp3.poly + comp3.vac + comp3.therm)
    @test comp3.poly == 0.0

    # MeanFieldState path
    st = Models.MeanFieldState(φ; Phi=0.123, PhiBar=0.456)
    compS = Models.omega_components(m, st, T, mu_vec; p_num=p_num, t_num=t_num, xi=0.0)
    @test compS.omega == comp3.omega

    # x_state length=5 path (Φ/Φbar ignored by NJLModel, but accepted)
    x5 = SVector{5, Float64}(φ[1], φ[2], φ[3], 0.123, 0.456)
    comp5 = Models.omega_components(m, x5, T, mu_vec; p_num=p_num, t_num=t_num, xi=0.0)
    @test comp5.omega ≈ (comp5.chi + comp5.poly + comp5.vac + comp5.therm)
    @test comp5.poly == 0.0

    # NJL ignores Φ/Φbar, so changing them should not affect the result
    x5b = SVector{5, Float64}(φ[1], φ[2], φ[3], 0.9, 0.1)
    comp5b = Models.omega_components(m, x5b, T, mu_vec; p_num=p_num, t_num=t_num, xi=0.0)
    @test comp5.omega == comp5b.omega
end
