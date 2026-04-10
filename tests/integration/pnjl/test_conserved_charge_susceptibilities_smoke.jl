using Test
using StaticArrays

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

const PNJL = Models.pnjl_module()

@testset "PNJL conserved-charge susceptibilities smoke" begin
    T_fm = 0.55
    muB_fm = 0.24
    kwargs = (; xi=0.0, p_num=6, t_num=3)

    vals = (
        PNJL.chi2_B(T_fm, muB_fm; kwargs...),
        PNJL.chi3_B(T_fm, muB_fm; kwargs...),
        PNJL.chi4_B(T_fm, muB_fm; kwargs...),
        PNJL.baryon_Ssigma(T_fm, muB_fm; kwargs...),
        PNJL.baryon_kappa_sigma2(T_fm, muB_fm; kwargs...),
    )
    @test all(isfinite, vals)

    chi2 = PNJL.chi2_B(T_fm, muB_fm; kwargs...)
    chi3 = PNJL.chi3_B(T_fm, muB_fm; kwargs...)
    chi4 = PNJL.chi4_B(T_fm, muB_fm; kwargs...)

    @test isapprox(PNJL.baryon_Ssigma(T_fm, muB_fm; kwargs...), chi3 / chi2; rtol=1e-8, atol=1e-10)
    @test isapprox(PNJL.baryon_kappa_sigma2(T_fm, muB_fm; kwargs...), chi4 / chi2; rtol=1e-8, atol=1e-10)
    @test abs(PNJL.chi1_B(0.5, 0.0; kwargs...)) <= 1e-8
    @test abs(PNJL.chi3_B(0.5, 0.0; kwargs...)) <= 1e-8
end
