using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !(isdefined(Main, :Models) && isdefined(Main.Models, :create_model))
    Base.include(Main, joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
if !isdefined(Main, :Constants_PNJL)
    Base.include(Main, joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
end

@testset "PNJL/RPNJL critical params validation smoke" begin
    cfg = Main.Constants_PNJL.pnjl_constants(profile="default", physics_profile="default", log_config=true)

    @test cfg.hbarc_MeV_fm > 0
    @test cfg.alpha_em > 0
    @test cfg.N_color >= 1
    @test cfg.N_flavor >= 1
    @test cfg.ρ0_fm3 > 0
    @test cfg.Λ_inv_fm > 0
    @test cfg.G_fm2 > 0
    @test cfg.K_fm5 > 0
    @test cfg.T0_inv_fm > 0

    model = Main.Models.create_model(:RPNJL; log_config=true)
    @test model.ext.kappa >= 0
    @test isfinite(model.ext.g1_fm8)
    @test isfinite(model.ext.g2_fm8)
end
