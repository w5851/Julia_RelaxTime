using Test

# New models entry
_models_entry = joinpath(@__DIR__, "..", "..", "..", "src", "models", "Models.jl")
if !isdefined(Main, :Models)
    include(_models_entry)
end
using .Models

# Model-based thermo facade
include(joinpath(@__DIR__, "..", "..", "..", "src", "models", "core", "ThermoFacade.jl"))
using .ThermoFacade

using StaticArrays

@testset "PNJL thermo bridge multipoint (unified vs Models)" begin
    model = Models.create_model(:PNJL)

    # Keep this cheap: same representative state for all points.
    x_state = @SVector [0.02, 0.02, 0.03, 0.5, 0.5]
    xi = 0.0

    p_num = 24
    t_num = 6
    thermal_nodes = ThermoFacade.cached_nodes_legacy(p_num, t_num)

    points = [
        (T=0.15, μ=0.0),
        (T=0.80, μ=0.10),
        (T=1.10, μ=0.20),
    ]

    for pt in points
        T = pt.T
        mu_vec = @SVector [pt.μ, pt.μ, pt.μ]

        ω_bridge = ThermoFacade.calculate_omega_backend(
            x_state,
            mu_vec,
            T;
            model=model,
            model_kind=:PNJL,
            p_num=p_num,
            t_num=t_num,
            xi=xi,
        )
        ω_models = Models.omega(model, x_state, T, mu_vec; p_num=p_num, t_num=t_num, xi=xi)

        comp_bridge = ThermoFacade.calculate_omega_components_backend(
            x_state,
            mu_vec,
            T;
            model=model,
            model_kind=:PNJL,
            p_num=p_num,
            t_num=t_num,
            xi=xi,
        )
        comp_models = Models.omega_components(model, x_state, T, mu_vec; p_num=p_num, t_num=t_num, xi=xi)

        @test isfinite(ω_bridge)
        @test isfinite(ω_models)
        @test isapprox(ω_models, ω_bridge; rtol=1e-9, atol=1e-11)

        @test isfinite(comp_bridge.omega)
        @test isfinite(comp_models.omega)
        @test isapprox(comp_models.omega, ω_models; rtol=1e-12, atol=1e-12)
        @test isapprox(comp_bridge.omega, ω_bridge; rtol=1e-12, atol=1e-12)

        for k in (:chi, :poly, :vac, :therm)
            @test isfinite(getfield(comp_bridge, k))
            @test isfinite(getfield(comp_models, k))
            @test isapprox(getfield(comp_models, k), getfield(comp_bridge, k); rtol=1e-9, atol=1e-11)
        end

        P_bridge, rho_norm_bridge, s_bridge, e_bridge = ThermoFacade.calculate_thermo_backend(
            x_state,
            mu_vec,
            T;
            model=model,
            model_kind=:PNJL,
            p_num=p_num,
            t_num=t_num,
            xi=xi,
        )

        @test isfinite(P_bridge)
        @test isapprox(P_bridge, -ω_bridge; rtol=1e-9, atol=1e-11)

        @test isfinite(rho_norm_bridge)
        @test rho_norm_bridge >= 0

        @test isfinite(s_bridge)
        @test isfinite(s_bridge)

        @test isfinite(e_bridge)
        @test isfinite(e_bridge)
    end
end
