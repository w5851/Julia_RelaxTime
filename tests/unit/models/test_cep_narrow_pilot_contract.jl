using Test

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
if !isdefined(Main, :Constants_PNJL)
    include(joinpath(PROJECT_ROOT, "src", "constants", "Constants_PNJL.jl"))
end
if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end
const PILOT_SCRIPT = joinpath(PROJECT_ROOT, "scripts", "analysis", "pnjl_cep_narrow_pilot.jl")
if !isdefined(Main, :PilotConfig)
    include(PILOT_SCRIPT)
end

@testset "CEP narrow pilot configuration and curve contract" begin
    cfg = Main.parse_config([
        "--xi", "0.0",
        "--method", "rho_support_cascade",
        "--calculation-sha", repeat("a", 40),
    ])
    @test cfg.xi == 0.0
    @test cfg.method == :rho_support_cascade
    @test cfg.p_num == 24
    @test cfg.t_num == 8
    @test cfg.targeted_max_points == 12
    @test cfg.oracle_refine_step == 0.003125

    rho = collect(0.0:0.05:2.0)
    x = rho .- 1.0
    mu = x .^ 3 .- 0.1 .* x
    curve = Main._curve_status(rho, mu)
    @test curve.status == :resolved_s_shape
    @test isfinite(curve.mu_transition)

    support = Main.Criticality.RhoSupportConfig(max_extra_points=12)
    @test support.max_extra_points <= cfg.targeted_max_points
end
