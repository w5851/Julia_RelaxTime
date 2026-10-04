using Test

isdefined(Main, :ChargedGBUQ0WindowAudit) || include(joinpath(@__DIR__, "..", "..", "..",
    "scripts", "analysis", "relaxtime", "audit_charged_gbu_q0_window.jl"))
const Q0_WINDOW_AUDIT_TEST = Main.ChargedGBUQ0WindowAudit

@testset "Finite-window sensitivity compares explicit cutoffs and separate resolutions" begin
    a = Q0_WINDOW_AUDIT_TEST
    vs = a.variants()
    @test [v.settings.omega_lower_inv_fm for v in vs[1:3]] == [1e-5,1e-4,1e-6]
    @test vs[4].settings.omega_nodes == 128
    @test vs[5].settings.mesh == 128
    @test vs[6].settings.q_nodes == 16
    @test last(vs).settings.endpoint_policy == "strict_zero_limit"
    @test a.compare(.1,.10001).passed
    @test !a.compare(.1,.11).passed
    @test !a.compare(NaN,.1).passed
    @test a.compare(0.,1e-13).passed
    @test a.compare(-.1,-.10001).passed
    @test length(a.CASES) == 4
end
