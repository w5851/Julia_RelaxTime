using Test

const _RELAXTIME_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAXTIME_PATH)
end
using Main.TransportCoefficients
const TCV = Main.RelaxTime.TransportCoefficientsValidation

const QUARK_PARAMS = (m=(u=0.3,d=0.3,s=0.5), μ=(u=0.2,d=0.2,s=0.2))
const THERMO_PARAMS = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
const QUARK_PARAMS_STRUCT = QuarkParams(QUARK_PARAMS)
const THERMO_PARAMS_STRUCT = ThermoParams(THERMO_PARAMS)

const TAU_ZERO = (u=0.0,d=0.0,s=0.0,ubar=0.0,dbar=0.0,sbar=0.0)
const TAU_ONE = (u=1.0,d=1.0,s=1.0,ubar=1.0,dbar=1.0,sbar=1.0)
const DENSITIES_ONE = (u=0.12,d=0.09,s=0.03,ubar=0.02,dbar=0.01,sbar=0.005)
const PRESSURE_ONE = 0.08
const ENERGY_ONE = 0.60

@testset "TransportCoefficients: eta/sigma basic" begin
    eta0 = shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ZERO, p_nodes=16, p_max=10.0)
    sigma0 = electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ZERO, p_nodes=16, p_max=10.0)
    @test eta0 == 0.0
    @test sigma0 == 0.0

    eta1 = shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, p_nodes=16, p_max=10.0)
    sigma1 = electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, p_nodes=16, p_max=10.0)
    @test isfinite(eta1)
    @test isfinite(sigma1)
    @test eta1 > 0.0
    @test sigma1 > 0.0

    cfg = TransportIntegrationConfig(p_nodes=16, p_max=10.0)
    eta1_cfg = shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, config=cfg)
    sigma1_cfg = electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, config=cfg)
    @test isapprox(eta1_cfg, eta1; rtol=1e-12, atol=0.0)
    @test isapprox(sigma1_cfg, sigma1; rtol=1e-12, atol=0.0)

    eta1_struct = shear_viscosity(QUARK_PARAMS_STRUCT, THERMO_PARAMS_STRUCT; tau=TAU_ONE, p_nodes=16, p_max=10.0)
    sigma1_struct = electric_conductivity(QUARK_PARAMS_STRUCT, THERMO_PARAMS_STRUCT; tau=TAU_ONE, p_nodes=16, p_max=10.0)
    @test isapprox(eta1_struct, eta1; rtol=1e-12, atol=0.0)
    @test isapprox(sigma1_struct, sigma1; rtol=1e-12, atol=0.0)
end

@testset "TransportCoefficients: conserved-charge helpers" begin
    @test conserved_charge_value(:u, :B) == 1 / 3
    @test conserved_charge_value(:ubar, :B) == -1 / 3
    @test conserved_charge_value(:u, :Q) == 2 / 3
    @test conserved_charge_value(:dbar, :Q) == 1 / 3
    @test conserved_charge_value(:s, :S) == -1
    @test conserved_charge_value(:sbar, :S) == 1

    ch = conserved_charge_densities(DENSITIES_ONE)
    @test isapprox(ch.B, (0.10 + 0.08 + 0.025) / 3; rtol=1e-12, atol=0.0)
    @test isapprox(ch.Q, (2 * 0.10 - 0.08 - 0.025) / 3; rtol=1e-12, atol=0.0)
    @test isapprox(ch.S, -0.025; rtol=1e-12, atol=0.0)
    @test enthalpy_density(PRESSURE_ONE, ENERGY_ONE) == PRESSURE_ONE + ENERGY_ONE
    @test isapprox(rho_mass_from_densities(QUARK_PARAMS.m, DENSITIES_ONE), (0.12 + 0.02) * 0.3 + (0.09 + 0.01) * 0.3 + (0.03 + 0.005) * 0.5; rtol=1e-12, atol=0.0)
end

@testset "TransportCoefficients: kappa/lambda basic" begin
    cfg = TransportIntegrationConfig(p_nodes=16, p_max=10.0)

    κBB0 = kappa_BB(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ZERO, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)
    κBQ0 = kappa_BQ(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ZERO, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)
    κBS0 = kappa_BS(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ZERO, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)
    κQQ0 = kappa_QQ(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ZERO, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)
    κQS0 = kappa_QS(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ZERO, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)
    κSS0 = kappa_SS(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ZERO, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)
    @test κBB0 == 0.0
    @test κBQ0 == 0.0
    @test κBS0 == 0.0
    @test κQQ0 == 0.0
    @test κQS0 == 0.0
    @test κSS0 == 0.0

    κBB = kappa_BB(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)
    κBQ = kappa_BQ(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)
    κBS = kappa_BS(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)
    κQQ = kappa_QQ(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)
    κQS = kappa_QS(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)
    κSS = kappa_SS(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)
    κmat = diffusion_matrix(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)
    λ = lambda_from_kappa_BB(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, densities=DENSITIES_ONE, pressure=PRESSURE_ONE, energy=ENERGY_ONE, config=cfg)

    @test isfinite(κBB)
    @test isfinite(κBQ)
    @test isfinite(κBS)
    @test isfinite(κQQ)
    @test isfinite(κQS)
    @test isfinite(κSS)
    @test κBB > 0.0
    @test κQQ > 0.0
    @test κSS > 0.0
    @test κmat.BQ == κmat.QB
    @test κmat.BS == κmat.SB
    @test κmat.QS == κmat.SQ
    @test isapprox(κmat.BB, κBB; rtol=1e-12, atol=0.0)
    @test isapprox(κmat.BQ, κBQ; rtol=1e-12, atol=0.0)
    @test isapprox(κmat.BS, κBS; rtol=1e-12, atol=0.0)
    @test isapprox(κmat.QQ, κQQ; rtol=1e-12, atol=0.0)
    @test isapprox(κmat.QS, κQS; rtol=1e-12, atol=0.0)
    @test isapprox(κmat.SS, κSS; rtol=1e-12, atol=0.0)
    @test size(κmat.matrix) == (3, 3)
    @test isfinite(λ)
    @test λ > 0.0
end

@testset "TransportCoefficients: diffusion follows legacy occupancy factor" begin
    cfg = TransportIntegrationConfig(p_nodes=12, p_max=8.0)
    base_provider = default_transport_provider()

    provider_f02 = merge(base_provider, (
        quark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.2,
        antiquark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.2,
        quark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.2,
        antiquark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.2,
    ))
    provider_f08 = merge(base_provider, (
        quark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.8,
        antiquark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.8,
        quark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.8,
        antiquark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.8,
    ))

    κ02 = kappa_BB(
        QUARK_PARAMS,
        THERMO_PARAMS;
        tau=TAU_ONE,
        densities=DENSITIES_ONE,
        pressure=PRESSURE_ONE,
        energy=ENERGY_ONE,
        provider=provider_f02,
        config=cfg,
    )

    κ08 = kappa_BB(
        QUARK_PARAMS,
        THERMO_PARAMS;
        tau=TAU_ONE,
        densities=DENSITIES_ONE,
        pressure=PRESSURE_ONE,
        energy=ENERGY_ONE,
        provider=provider_f08,
        config=cfg,
    )

    @test isfinite(κ02)
    @test isfinite(κ08)
    # Legacy/doc formula uses occupancy kernel F=f; for fixed f=0.2 and f=0.8,
    # diffusion scales approximately linearly with f.
    @test isapprox(κ08 / κ02, 4.0; rtol=1e-2, atol=0.0)
end

@testset "TransportCoefficients: derived ratio helpers" begin
    λ = 2.4
    σ = 1.5
    T = 0.2
    η = 0.9
    s = 3.0
    ζ = 0.18
    c_p = 4.0
    rho_mass = 2.0

    @test isapprox(lorenz_number(λ, σ, T), λ / (σ * T); rtol=1e-12, atol=0.0)
    @test isapprox(lorentz_legacy(λ, σ, T), λ / (σ / T); rtol=1e-12, atol=0.0)
    @test isapprox(viscous_conductive_coupling_ratio(η, s, σ, T), (η / s) / (σ / T); rtol=1e-12, atol=0.0)
    @test isapprox(prandtl_number(η, c_p, λ, rho_mass), η * c_p / (λ * rho_mass); rtol=1e-12, atol=0.0)
    @test isapprox(bulk_to_shear_viscosity_ratio(ζ, η), ζ / η; rtol=1e-12, atol=0.0)

    @test isnan(lorenz_number(λ, 0.0, T))
    @test isnan(lorentz_legacy(λ, 0.0, T))
    @test isnan(viscous_conductive_coupling_ratio(η, 0.0, σ, T))
    @test isnan(prandtl_number(η, c_p, 0.0, rho_mass))
    @test isnan(bulk_to_shear_viscosity_ratio(ζ, 0.0))
end

@testset "TransportCoefficients: input validation guards" begin
    cfg = TransportIntegrationConfig(p_nodes=8, p_max=4.0)

    bad_T = merge(THERMO_PARAMS, (T=0.0,))
    @test_throws ArgumentError shear_viscosity(QUARK_PARAMS, bad_T; tau=TAU_ONE, config=cfg)

    bad_mass = (m=(u=-0.1, d=0.3, s=0.5), μ=QUARK_PARAMS.μ)
    @test_throws ArgumentError shear_viscosity(bad_mass, THERMO_PARAMS; tau=TAU_ONE, config=cfg)

    tau_missing = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0)
    @test_throws ErrorException shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=tau_missing, config=cfg)

    tau_bad_nonfinite = (u=1.0, d=NaN, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)
    @test_throws ArgumentError shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=tau_bad_nonfinite, config=cfg)

    tau_bad_negative = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=-1.0, sbar=1.0)
    @test_throws ArgumentError shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=tau_bad_negative, config=cfg)

    bad_charges = (u=NaN, d=-1 / 3, s=-1 / 3)
    @test_throws ArgumentError electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, config=cfg, charges=bad_charges)

    bad_cfg = TransportIntegrationConfig(p_nodes=8, p_max=4.0)
    bad_bulk = (
        v_n_sq=0.3,
        dμB_dT_sigma=0.1,
        masses=[0.3, 0.3],
        dM_dT=[0.0, 0.0, 0.0],
        dM_dμB=[0.0, 0.0, 0.0],
    )
    @test_throws ArgumentError bulk_viscosity_isentropic(
        QUARK_PARAMS,
        THERMO_PARAMS;
        tau=TAU_ONE,
        config=bad_cfg,
        bulk_coeffs_isentropic=bad_bulk,
    )

    bad_bulk_nonfinite = (
        v_n_sq=0.3,
        dμB_dT_sigma=Inf,
        masses=[0.3, 0.3, 0.5],
        dM_dT=[0.0, 0.0, NaN],
        dM_dμB=[0.0, 0.0, 0.0],
    )
    @test_throws ArgumentError bulk_viscosity_isentropic(
        QUARK_PARAMS,
        THERMO_PARAMS;
        tau=TAU_ONE,
        config=bad_cfg,
        bulk_coeffs_isentropic=bad_bulk_nonfinite,
    )
end

@testset "TransportCoefficients: numerical safeguards" begin
    thermo_params_aniso = merge(THERMO_PARAMS, (ξ=0.2,))

    # 数值保护只覆盖物理占据数范围内的极端输入；这里用接近 1 的边界值
    # 验证极小能量下系数仍保持有限且非负，而不改变核心物理表达式。
    provider_edge_f = (
        energy_from_p=(p::Float64, m::Float64) -> 1e-30,
        energy_from_p_aniso=(p::Float64, m::Float64, ξ::Float64, c::Float64) -> 1e-30,
        quark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 1.0 - 1e-12,
        antiquark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 1.0 - 1e-12,
        quark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 1.0 - 1e-12,
        antiquark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 1.0 - 1e-12,
        prefer_energy_aniso=true,
    )

    η = shear_viscosity(
        QUARK_PARAMS,
        thermo_params_aniso;
        tau=TAU_ONE,
        provider=provider_edge_f,
        p_nodes=6,
        cos_nodes=6,
        p_max=4.0,
    )
    @test isfinite(η)
    @test η >= 0.0
end

@testset "TransportCoefficients: bulk_viscosity unified entry" begin
    coeffs = (
        v_n_sq=0.3,
        dμB_dT_sigma=0.1,
        masses=[0.3, 0.3, 0.5],
        dM_dT=[0.0, 0.0, 0.0],
        dM_dμB=[0.0, 0.0, 0.0],
    )

    zeta_iso = bulk_viscosity_isentropic(
        QUARK_PARAMS,
        THERMO_PARAMS;
        tau=TAU_ONE,
        p_nodes=8,
        p_max=4.0,
        bulk_coeffs_isentropic=coeffs,
    )

    zeta_unified = bulk_viscosity(
        QUARK_PARAMS,
        THERMO_PARAMS;
        formula=:isentropic,
        tau=TAU_ONE,
        p_nodes=8,
        p_max=4.0,
        bulk_coeffs_isentropic=coeffs,
    )

    zeta_legacy_alias = bulk_viscosity(
        QUARK_PARAMS,
        THERMO_PARAMS;
        formula=:isentropic,
        tau=TAU_ONE,
        p_nodes=8,
        p_max=4.0,
        bulk_coeffs=coeffs,
    )

    @test isapprox(zeta_unified, zeta_iso; rtol=1e-12, atol=0.0)
    @test isapprox(zeta_legacy_alias, zeta_iso; rtol=1e-12, atol=0.0)
    @test_throws ErrorException bulk_viscosity(QUARK_PARAMS, THERMO_PARAMS; formula=:unknown, tau=TAU_ONE, bulk_coeffs_isentropic=coeffs)
end

@testset "TransportCoefficients: sigma scales with q^2" begin
    charges1 = (u=2/3, d=-1/3, s=-1/3)
    charges2 = (u=4/3, d=-2/3, s=-2/3)

    σ1 = electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, charges=charges1, p_nodes=16, p_max=10.0)
    σ2 = electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, charges=charges2, p_nodes=16, p_max=10.0)

    # charges2 = 2 * charges1 => q^2 scales by 4
    @test isapprox(σ2 / σ1, 4.0; rtol=1e-3)
end

@testset "TransportCoefficients: config precedence and validation" begin
    cfg_p8 = TransportIntegrationConfig(p_nodes=8, p_max=10.0)
    eta_p16 = shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, p_nodes=16, p_max=10.0)
    eta_override = shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, config=cfg_p8, p_nodes=16)
    @test isapprox(eta_override, eta_p16; rtol=1e-12, atol=0.0)

    @test_throws ArgumentError TransportIntegrationConfig(p_grid=[0.0, 1.0])
    @test_throws ArgumentError TransportIntegrationConfig(cos_grid=[-1.0, 1.0])
    @test_throws ArgumentError TransportIntegrationConfig(p_grid=[0.0, 1.0], p_w=[1.0])
    @test_throws ArgumentError TransportIntegrationConfig(cos_grid=[-1.0, 1.0], cos_w=[1.0])
    @test_throws ArgumentError TransportIntegrationConfig(p_max=0.0)

    # positional-config overloads
    eta_pos = shear_viscosity(QUARK_PARAMS, THERMO_PARAMS, cfg_p8; tau=TAU_ONE, p_nodes=16)
    eta_kw = shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, config=cfg_p8, p_nodes=16)
    @test isapprox(eta_pos, eta_kw; rtol=1e-12, atol=0.0)

    sigma_pos = electric_conductivity(QUARK_PARAMS, THERMO_PARAMS, cfg_p8; tau=TAU_ONE)
    sigma_kw = electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, config=cfg_p8)
    @test isapprox(sigma_pos, sigma_kw; rtol=1e-12, atol=0.0)

    # request-based recommended entry
    req = TransportRequest(
        QUARK_PARAMS,
        THERMO_PARAMS;
        tau=TAU_ONE,
        integration=cfg_p8,
    )
    eta_req = shear_viscosity(req; p_nodes=16)
    @test isapprox(eta_req, eta_kw; rtol=1e-12, atol=0.0)

    sigma_req = electric_conductivity(req)
    @test isapprox(sigma_req, sigma_kw; rtol=1e-12, atol=0.0)

    req_struct = TransportRequest(
        QUARK_PARAMS_STRUCT,
        THERMO_PARAMS_STRUCT;
        tau=TAU_ONE,
        integration=cfg_p8,
    )
    eta_req_struct = shear_viscosity(req_struct; p_nodes=16)
    sigma_req_struct = electric_conductivity(req_struct)
    @test isapprox(eta_req_struct, eta_kw; rtol=1e-12, atol=0.0)
    @test isapprox(sigma_req_struct, sigma_kw; rtol=1e-12, atol=0.0)
end

@testset "TransportCoefficients: request adapter layering" begin
    req = TransportCoefficients._build_transport_request(
        QUARK_PARAMS,
        THERMO_PARAMS,
        TAU_ONE;
        charges=Main.TransportCoefficients.default_charges(),
        degeneracy=Main.TransportCoefficients.degeneracy_default(),
        config=nothing,
        kwargs=pairs((; p_nodes=12, p_max=5.0)),
    )

    @test req isa TransportRequest
    @test req.integration.p_nodes == 12
    @test req.integration.p_max == 5.0
    @test req.physics.degeneracy == Main.TransportCoefficients.degeneracy_default()
end

@testset "TransportCoefficients: contract layer guards" begin
    @test_throws ErrorException TCV._validate_transport_request_contract(
        DENSITIES_ONE,
        PRESSURE_ONE,
        nothing,
    )
    @test_throws ErrorException TCV._validate_transport_request_contract(
        DENSITIES_ONE,
        nothing,
        ENERGY_ONE,
    )
    @test_throws ErrorException TCV._validate_transport_request_contract(
        nothing,
        PRESSURE_ONE,
        ENERGY_ONE,
    )

    @test TCV._validate_transport_request_contract(
        DENSITIES_ONE,
        PRESSURE_ONE,
        ENERGY_ONE,
    ) == true
    @test TCV._validate_transport_request_contract(
        nothing,
        nothing,
        nothing,
    ) == false
end

@testset "TransportCoefficients: provider injection smoke" begin
    prov = default_transport_provider()

    eta_default = shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, p_nodes=16, p_max=10.0)
    eta_injected = shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, p_nodes=16, p_max=10.0, provider=prov)
    @test isapprox(eta_injected, eta_default; rtol=1e-12, atol=0.0)

    sigma_default = electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, p_nodes=16, p_max=10.0)
    sigma_injected = electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, p_nodes=16, p_max=10.0, provider=prov)
    @test isapprox(sigma_injected, sigma_default; rtol=1e-12, atol=0.0)

    req = TransportRequest(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, integration=TransportIntegrationConfig(p_nodes=16, p_max=10.0))
    eta_req_default = shear_viscosity(req)
    eta_req_injected = shear_viscosity(req; provider=prov)
    @test isapprox(eta_req_injected, eta_req_default; rtol=1e-12, atol=0.0)

    toy_provider = (
        energy_from_p=(p::Float64, m::Float64) -> sqrt(p * p + m * m),
        quark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.1,
        antiquark_distribution=(E::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64) -> 0.1,
        quark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.1,
        antiquark_distribution_aniso=(p::Float64, m::Float64, μ::Float64, T::Float64, Φ::Float64, Φbar::Float64, ξ::Float64, c::Float64) -> 0.1,
    )

    eta_toy = shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, p_nodes=16, p_max=10.0, provider=toy_provider)
    @test isfinite(eta_toy)
    @test eta_toy > 0.0
    @test !isapprox(eta_toy, eta_default; rtol=1e-6, atol=0.0)
end

@testset "TransportCoefficients: provider mass/mu override smoke" begin
    # Make sure transport integrand can be decoupled from hard-coded quark_params fields.
    qp_bad = (m=(u=NaN, d=NaN, s=NaN), μ=(u=999.0, d=999.0, s=999.0))

    @test_throws ArgumentError shear_viscosity(qp_bad, THERMO_PARAMS; tau=TAU_ONE, p_nodes=16, p_max=10.0)

    cached_m = (u=0.3, d=0.3, s=0.5)
    cached_mu = (u=0.2, d=0.2, s=0.2)

    prov_base = default_transport_provider()
    prov = merge(prov_base, (
        mass_for_species=(sp::Symbol, qp, tp) -> begin
            sp in (:u, :ubar) && return cached_m.u
            sp in (:d, :dbar) && return cached_m.d
            sp in (:s, :sbar) && return cached_m.s
            error("unknown species=$sp")
        end,
        mu_for_species=(sp::Symbol, qp, tp) -> begin
            sp in (:u, :ubar) && return cached_mu.u
            sp in (:d, :dbar) && return cached_mu.d
            sp in (:s, :sbar) && return cached_mu.s
            error("unknown species=$sp")
        end,
    ))

    eta_bad = shear_viscosity(qp_bad, THERMO_PARAMS; tau=TAU_ONE, p_nodes=16, p_max=10.0, provider=prov)
    @test isfinite(eta_bad)
    @test eta_bad > 0.0

    @testset "TransportCoefficients anisotropic energy hook" begin
        # Minimal anisotropic smoke: ensure ξ≠0 path is finite and that overriding
        # `energy_from_p_aniso` affects the result when `prefer_energy_aniso=true`.
        thermo_params_aniso = merge(THERMO_PARAMS, (ξ=0.2,))

        η0 = shear_viscosity(
            QUARK_PARAMS,
            thermo_params_aniso;
            tau=TAU_ONE,
            p_nodes=6,
            cos_nodes=6,
            p_max=6.0,
        )
        @test isfinite(η0)
        @test η0 > 0

        provider = default_transport_provider()
        provider2 = merge(
            provider,
            (
                prefer_energy_aniso=true,
                energy_from_p_aniso=(p, m, ξ, cosθ) -> sqrt(p * p + m * m + 10.0 * ξ * (p * cosθ)^2),
            ),
        )

        η2 = shear_viscosity(
            QUARK_PARAMS,
            thermo_params_aniso;
            tau=TAU_ONE,
            provider=provider2,
            p_nodes=6,
            cos_nodes=6,
            p_max=6.0,
        )
        @test isfinite(η2)
        @test η2 > 0
        @test !isapprox(η2, η0; rtol=1e-10, atol=0.0)
    end

    @testset "TransportCoefficients aniso distribution fallback" begin
        # Provider without *_distribution_aniso should still work for ξ≠0 if
        # energy passthrough is available.
        thermo_params_aniso = merge(THERMO_PARAMS, (ξ=0.2,))

        prov_full = default_transport_provider()
        prov_no_aniso = (
            energy_from_p=prov_full.energy_from_p,
            energy_from_p_aniso=prov_full.energy_from_p_aniso,
            quark_distribution=prov_full.quark_distribution,
            antiquark_distribution=prov_full.antiquark_distribution,
            prefer_energy_aniso=false,
        )

        η = shear_viscosity(
            QUARK_PARAMS,
            thermo_params_aniso;
            tau=TAU_ONE,
            provider=prov_no_aniso,
            p_nodes=6,
            cos_nodes=6,
            p_max=6.0,
        )
        @test isfinite(η)
        @test η > 0
    end
end

@testset "TransportCoefficients: transport_coefficients diffusion extension" begin
    cfg = TransportIntegrationConfig(p_nodes=10, p_max=6.0)

    tr_basic = transport_coefficients(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, config=cfg)
    @test isnan(tr_basic.kappa_BB)
    @test isnan(tr_basic.kappa_BQ)
    @test isnan(tr_basic.kappa_BS)
    @test isnan(tr_basic.kappa_QQ)
    @test isnan(tr_basic.kappa_QS)
    @test isnan(tr_basic.kappa_SS)
    @test tr_basic.diffusion_matrix === nothing
    @test isnan(tr_basic.lambda)
    @test isnan(tr_basic.lorenz_number)
    @test isnan(tr_basic.lorentz_legacy)
    @test isnan(tr_basic.viscous_conductive_coupling_ratio)
    @test isnan(tr_basic.prandtl_number)
    @test isnan(tr_basic.bulk_to_shear_viscosity_ratio)

    tr_full = transport_coefficients(
        QUARK_PARAMS,
        THERMO_PARAMS;
        tau=TAU_ONE,
        config=cfg,
        densities=DENSITIES_ONE,
        pressure=PRESSURE_ONE,
        energy=ENERGY_ONE,
        entropy=0.42,
        c_p=1.8,
        rho_mass=0.75,
    )
    @test isfinite(tr_full.kappa_BB)
    @test isfinite(tr_full.kappa_BQ)
    @test isfinite(tr_full.kappa_BS)
    @test isfinite(tr_full.kappa_QQ)
    @test isfinite(tr_full.kappa_QS)
    @test isfinite(tr_full.kappa_SS)
    @test tr_full.diffusion_matrix !== nothing
    @test isapprox(tr_full.diffusion_matrix.BQ, tr_full.kappa_BQ; rtol=1e-12, atol=0.0)
    @test isfinite(tr_full.lambda)
    @test isfinite(tr_full.lorenz_number)
    @test isfinite(tr_full.lorentz_legacy)
    @test isfinite(tr_full.viscous_conductive_coupling_ratio)
    @test isfinite(tr_full.prandtl_number)
    @test isnan(tr_full.bulk_to_shear_viscosity_ratio)
    @test tr_full.kappa_BB > 0.0
    @test isapprox(tr_full.lorenz_number, tr_full.lambda / (tr_full.sigma * THERMO_PARAMS.T); rtol=1e-12, atol=0.0)
    @test isapprox(tr_full.lorentz_legacy, tr_full.lambda / (tr_full.sigma / THERMO_PARAMS.T); rtol=1e-12, atol=0.0)
    @test isapprox(tr_full.viscous_conductive_coupling_ratio, (tr_full.eta / 0.42) / (tr_full.sigma / THERMO_PARAMS.T); rtol=1e-12, atol=0.0)
    @test isapprox(tr_full.prandtl_number, tr_full.eta * 1.8 / (tr_full.lambda * 0.75); rtol=1e-12, atol=0.0)
    @test isnan(tr_full.bulk_to_shear_viscosity_ratio)
end
