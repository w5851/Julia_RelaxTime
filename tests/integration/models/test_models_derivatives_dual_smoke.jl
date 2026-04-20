using Test
using ForwardDiff

const PROJECT_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const HBARC_MEV_FM = 197.327

if !isdefined(Main, :Models)
    include(joinpath(PROJECT_ROOT, "src", "models", "Models.jl"))
end

@testset "Models derivatives dual compatibility smoke" begin
    T0 = 100.0 / HBARC_MEV_FM
    μ0 = 20.0 / HBARC_MEV_FM

    md = Models.mass_derivatives(T0, μ0; order=1, p_num=8, t_num=4)
    @test length(md.masses) == 3
    @test all(isfinite, md.masses)
    @test length(md.dM_dT) == 3
    @test length(md.dM_dmu) == 3

    td = Models.thermo_derivatives(T0, μ0; p_num=8, t_num=4)
    @test isfinite(td.pressure)
    @test isfinite(td.energy)

    coeff = Models.bulk_viscosity_coefficients(T0, μ0; p_num=8, t_num=4)
    @test haskey(coeff, :v_n_sq)
    @test haskey(coeff, :dμB_dT_sigma)
    @test haskey(coeff, :dM_dT)
    @test haskey(coeff, :dM_dμB)

    f_mass(T) = first(Models.mass_derivatives(T, μ0; order=1, p_num=8, t_num=4).masses)
    f_pressure(μ) = Models.thermo_derivatives(T0, μ; p_num=8, t_num=4).pressure

    @test isfinite(ForwardDiff.derivative(f_mass, T0))
    @test isfinite(ForwardDiff.derivative(f_pressure, μ0))

    chi2B = Models.chi2_B(T0, μ0; p_num=8, t_num=4)
    chi4B = Models.chi4_B(T0, μ0; p_num=8, t_num=4)
    kappa_sigma2 = Models.baryon_kappa_sigma2(T0, μ0; p_num=8, t_num=4)
    @test isfinite(chi2B)
    @test isfinite(chi4B)
    @test isfinite(kappa_sigma2)
    @test kappa_sigma2 ≈ chi4B / chi2B

    h_mu = 1e-3
    chi3_p = Models.chi3_B(T0, μ0 + h_mu; p_num=8, t_num=4)
    chi3_m = Models.chi3_B(T0, μ0 - h_mu; p_num=8, t_num=4)
    chi4_fd = T0 * (chi3_p - chi3_m) / (2h_mu)
    @test isfinite(chi4_fd)
    @test isapprox(chi4B, chi4_fd; rtol=5e-2, atol=5e-3)

    @test !isdefined(Models.ConservedChargeSusceptibilities, :_nth_derivative)
    @test !isdefined(Models.ConservedChargeSusceptibilities, :_susceptibility_scale)
end
