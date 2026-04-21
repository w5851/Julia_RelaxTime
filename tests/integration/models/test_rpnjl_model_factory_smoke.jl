using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
_models_entry = joinpath(PROJECT_ROOT, "src", "models", "Models.jl")
if !(isdefined(Main, :Models) && isdefined(Main.Models, :create_model) && isdefined(Main.Models, :omega))
    Base.include(Main, _models_entry)
end

@testset "Models RPNJL factory smoke" begin
    model = Models.create_model(:RPNJL)
    model_deg = Models.create_model(:RPNJL; use_rpnjl_extensions=false)
    @test model isa Models.RPNJLModel

    x_state = Models.MeanFieldState([0.001, 0.001, 0.06, 0.2, 0.2])
    T_fm = 0.15
    mu = 0.0

    ω = Models.omega(model, x_state, T_fm, mu; p_num=8, t_num=4, xi=0.0)
    @test isfinite(ω)

    comps = Models.omega_components(model, x_state, T_fm, mu; p_num=8, t_num=4, xi=0.0)
    @test all(isfinite, (comps.chi, comps.poly, comps.vac, comps.therm, comps.omega))

    @test hasproperty(model.ext, :g1_fm8)
    @test hasproperty(model.ext, :g2_fm8)
    @test hasproperty(model.ext, :kappa)
    @test isapprox(model.ext.kappa, 0.1; rtol=0.0, atol=1e-12)

    Φ = 0.26
    Φbar = 0.18
    poly_deg = Models.polyakov_potential(model_deg, Φ, Φbar, T_fm)
    poly_ext = Models.polyakov_potential(model, Φ, Φbar, T_fm)
    @test isfinite(poly_deg)
    @test isfinite(poly_ext)

    model_k0 = Models.RPNJLModel(model.base, merge(model.ext, (kappa=0.0,)), true)
    poly_k0 = Models.polyakov_potential(model_k0, Φ, Φbar, T_fm)
    @test isfinite(poly_k0)
    @test !isapprox(poly_ext, poly_k0; rtol=1e-8, atol=1e-12)
end
