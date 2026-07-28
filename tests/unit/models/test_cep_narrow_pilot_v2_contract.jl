using Test

const V2_PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Models)
    include(joinpath(V2_PROJECT_ROOT, "src", "models", "Models.jl"))
end
if !isdefined(Main, :PilotV2Config)
    include(joinpath(V2_PROJECT_ROOT, "scripts", "analysis", "pnjl_cep_narrow_pilot_v2.jl"))
end

@testset "CEP narrow pilot v2 contract" begin
    cfg = Main._v2_parse_config([
        "--xi", "0.5",
        "--method", "rho_support_cascade",
        "--stage", "discovery",
        "--calculation-sha", repeat("a", 40),
    ])
    @test cfg.method == :rho_support_cascade
    @test cfg.stage == :discovery
    @test cfg.cascade_coarse_step == 0.05
    @test cfg.cascade_fine_step == 0.025
    @test cfg.targeted_max_points == 12
    @test cfg.temperature_resolution_target_MeV == 0.125

    @test Main._v2_target_config(12).target_point_count == 9
    @test Main._v2_target_config(8).target_point_count == 7
    @test Main._v2_target_config(4) === nothing

    canonical = (T=100.0, muq=300.0)
    local_cfg = Main.PilotV2Config(xi=0.5, calculation_sha=repeat("a", 40))
    window = Main._v2_window(local_cfg, canonical)
    @test window.T_min == 92.0
    @test window.T_max == 132.0
    @test isempty(window.required_anchors)
    @test Main._v2_in_window(100.0, window)
    @test !Main._v2_in_window(91.999, window)
    @test !Main._v2_in_window(132.001, window)
end
