using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

include(joinpath(PROJECT_ROOT, "src", "Constants_PNJL.jl"))
include(joinpath(PROJECT_ROOT, "src", "pnjl", "PNJL.jl"))

using .PNJL

@testset "magnetic nmax convergence criterion" begin
    x_state = SVector{5, Float64}(-0.03, -0.03, -0.04, 0.2, 0.2)
    T_fm = 150.0 / Constants_PNJL.ħc_MeV_fm
    mu_fm = 80.0 / Constants_PNJL.ħc_MeV_fm
    eB_fm2 = 2.0e4 / (Constants_PNJL.ħc_MeV_fm^2)
    mu_vec = SVector{3, Float64}(mu_fm, mu_fm, mu_fm)

    conf = Main.PNJL.default_magnetic_config(eB_fm2=eB_fm2)
    rep = Main.PNJL.magnetic_nmax_convergence_report(x_state, mu_vec, T_fm, conf; delta_n=6, rtol=3e-2)

    @test rep.n_probe > rep.n_base
    @test rep.rel_diff >= 0.0
    @test rep.rtol == 3e-2
    @test rep.converged isa Bool
end
