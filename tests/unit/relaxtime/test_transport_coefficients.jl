using Test

const _TRANSPORT_COEFFICIENTS_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "TransportCoefficients.jl"))
if !isdefined(Main, :TransportCoefficients)
    Base.include(Main, _TRANSPORT_COEFFICIENTS_PATH)
end
using Main.TransportCoefficients

const QUARK_PARAMS = (m=(u=0.3,d=0.3,s=0.5), μ=(u=0.2,d=0.2,s=0.2))
const THERMO_PARAMS = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
const QUARK_PARAMS_STRUCT = QuarkParams(QUARK_PARAMS)
const THERMO_PARAMS_STRUCT = ThermoParams(THERMO_PARAMS)

const TAU_ZERO = (u=0.0,d=0.0,s=0.0,ubar=0.0,dbar=0.0,sbar=0.0)
const TAU_ONE = (u=1.0,d=1.0,s=1.0,ubar=1.0,dbar=1.0,sbar=1.0)

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

@testset "TransportCoefficients: input validation guards" begin
    cfg = TransportIntegrationConfig(p_nodes=8, p_max=4.0)

    bad_T = merge(THERMO_PARAMS, (T=0.0,))
    @test_throws ErrorException shear_viscosity(QUARK_PARAMS, bad_T; tau=TAU_ONE, config=cfg)

    bad_mass = (m=(u=-0.1, d=0.3, s=0.5), μ=QUARK_PARAMS.μ)
    @test_throws ErrorException shear_viscosity(bad_mass, THERMO_PARAMS; tau=TAU_ONE, config=cfg)

    tau_missing = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0)
    @test_throws ErrorException shear_viscosity(QUARK_PARAMS, THERMO_PARAMS; tau=tau_missing, config=cfg)

    bad_charges = (u=NaN, d=-1 / 3, s=-1 / 3)
    @test_throws ErrorException electric_conductivity(QUARK_PARAMS, THERMO_PARAMS; tau=TAU_ONE, config=cfg, charges=bad_charges)

    bad_cfg = TransportIntegrationConfig(p_nodes=8, p_max=4.0)
    bad_bulk = (
        v_n_sq=0.3,
        dμB_dT_sigma=0.1,
        masses=[0.3, 0.3],
        dM_dT=[0.0, 0.0, 0.0],
        dM_dμB=[0.0, 0.0, 0.0],
    )
    @test_throws ErrorException bulk_viscosity_isentropic(
        QUARK_PARAMS,
        THERMO_PARAMS;
        tau=TAU_ONE,
        config=bad_cfg,
        bulk_coeffs_isentropic=bad_bulk,
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

    @test_throws ErrorException TransportIntegrationConfig(p_grid=[0.0, 1.0])
    @test_throws ErrorException TransportIntegrationConfig(cos_grid=[-1.0, 1.0])
    @test_throws ErrorException TransportIntegrationConfig(p_grid=[0.0, 1.0], p_w=[1.0])
    @test_throws ErrorException TransportIntegrationConfig(cos_grid=[-1.0, 1.0], cos_w=[1.0])

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

    @test_throws ErrorException shear_viscosity(qp_bad, THERMO_PARAMS; tau=TAU_ONE, p_nodes=16, p_max=10.0)

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

