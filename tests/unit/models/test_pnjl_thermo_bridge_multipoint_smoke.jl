using Test

# New models entry
const _MODELS_ENTRY = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_MODELS_ENTRY)
end
using .Models

# Legacy PNJL thermodynamics
include(joinpath(@__DIR__, "..", "..", "..", "src", "pnjl", "core", "Thermodynamics.jl"))
using .Thermodynamics

# Model-based thermo bridge
include(joinpath(@__DIR__, "..", "..", "..", "src", "pnjl", "core", "ModelThermodynamics.jl"))
using .ModelThermodynamics

using StaticArrays

@testset "PNJL thermo bridge multipoint (legacy vs Models)" begin
    model = Models.create_model(:PNJL)

    # Keep this cheap: same representative state for all points.
    x_state = @SVector [0.02, 0.02, 0.03, 0.5, 0.5]
    xi = 0.0

    p_num = 24
    t_num = 6
    thermal_nodes = Thermodynamics.Integrals.cached_nodes(p_num, t_num)

    points = [
        (T=0.15, μ=0.0),
        (T=0.80, μ=0.10),
        (T=1.10, μ=0.20),
    ]

    for pt in points
        T = pt.T
        mu_vec = @SVector [pt.μ, pt.μ, pt.μ]

        ω_legacy = Thermodynamics.calculate_omega(x_state, mu_vec, T, thermal_nodes, xi)
        ω_models = Models.omega(model, x_state, T, mu_vec; p_num=p_num, t_num=t_num, xi=xi)

        comp_legacy = Thermodynamics.calculate_omega_components(x_state, mu_vec, T, thermal_nodes, xi)
        comp_models = Models.omega_components(model, x_state, T, mu_vec; p_num=p_num, t_num=t_num, xi=xi)

        @test isfinite(ω_legacy)
        @test isfinite(ω_models)
        @test isapprox(ω_models, ω_legacy; rtol=1e-9, atol=1e-11)

        @test isfinite(comp_legacy.omega)
        @test isfinite(comp_models.omega)
        @test isapprox(comp_models.omega, ω_models; rtol=1e-12, atol=1e-12)
        @test isapprox(comp_legacy.omega, ω_legacy; rtol=1e-12, atol=1e-12)

        for k in (:chi, :poly, :vac, :therm)
            @test isfinite(getfield(comp_legacy, k))
            @test isfinite(getfield(comp_models, k))
            @test isapprox(getfield(comp_models, k), getfield(comp_legacy, k); rtol=1e-9, atol=1e-11)
        end

        P_legacy, rho_norm_legacy, s_legacy, e_legacy = Thermodynamics.calculate_thermo(x_state, mu_vec, T, thermal_nodes, xi)
        P_m, rho_norm_m, s_m, e_m = ModelThermodynamics.thermo_model(model, x_state, mu_vec, T; p_num=p_num, t_num=t_num, xi=xi)

        @test isfinite(P_legacy)
        @test isfinite(P_m)
        @test isapprox(P_m, P_legacy; rtol=1e-9, atol=1e-11)

        @test isfinite(rho_norm_legacy)
        @test isfinite(rho_norm_m)
        @test isapprox(rho_norm_m, rho_norm_legacy; rtol=1e-9, atol=1e-11)

        @test isfinite(s_legacy)
        @test isfinite(s_m)
        @test isapprox(s_m, s_legacy; rtol=5e-8, atol=5e-10)

        @test isfinite(e_legacy)
        @test isfinite(e_m)
        @test isapprox(e_m, e_legacy; rtol=5e-8, atol=5e-10)
    end
end
