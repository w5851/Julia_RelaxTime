using Test

include("../../src/relaxtime/RelaxationTime.jl")
using .RelaxationTime

const DENSITIES_SAMPLE = (u=1.0, d=1.0, s=2.0, ubar=3.0, dbar=3.0, sbar=4.0)
const RATES_SAMPLE = (
    uu_to_uu=1.0,
    ud_to_ud=2.0,
    us_to_us=3.0,
    usbar_to_usbar=5.0,
    uubar_to_uubar=7.0,
    uubar_to_ddbar=11.0,
    uubar_to_ssbar=13.0,
    udbar_to_udbar=17.0,
    ss_to_ss=19.0,
    ssbar_to_ssbar=23.0,
    ssbar_to_uubar=29.0,
)

const EXPECTED_TAU_INV = (
    u = 173.0,
    d = 173.0,
    s = 398.0,
    ubar = 142.0,
    dbar = 142.0,
    sbar = 410.0,
)

const EXPECTED_TAU = (
    u = 1 / EXPECTED_TAU_INV.u,
    d = 1 / EXPECTED_TAU_INV.d,
    s = 1 / EXPECTED_TAU_INV.s,
    ubar = 1 / EXPECTED_TAU_INV.ubar,
    dbar = 1 / EXPECTED_TAU_INV.dbar,
    sbar = 1 / EXPECTED_TAU_INV.sbar,
)

# Minimal parameter sets (values unused when existing_rates covers all processes)
const QUARK_PARAMS = (m=(u=0.1,d=0.1,s=0.2), μ=(u=0.0,d=0.0,s=0.0))
const THERMO_PARAMS = (T=0.15, Φ=0.5, Φbar=0.5, ξ=0.0)
const K_COEFFS = (K_σπ=1.0, K_σK=1.0, K_σ=1.0, K_δπ=1.0, K_δK=1.0)

@testset "relaxation_rates algebra" begin
    tau_inv = relaxation_rates(DENSITIES_SAMPLE, RATES_SAMPLE)
    @test tau_inv.u   ≈ EXPECTED_TAU_INV.u
    @test tau_inv.d   ≈ EXPECTED_TAU_INV.d
    @test tau_inv.s   ≈ EXPECTED_TAU_INV.s
    @test tau_inv.ubar ≈ EXPECTED_TAU_INV.ubar
    @test tau_inv.dbar ≈ EXPECTED_TAU_INV.dbar
    @test tau_inv.sbar ≈ EXPECTED_TAU_INV.sbar
end

@testset "relaxation_times uses provided rates" begin
    result = relaxation_times(
        QUARK_PARAMS,
        THERMO_PARAMS,
        K_COEFFS;
        densities=DENSITIES_SAMPLE,
        existing_rates=RATES_SAMPLE,
    )

    @test result.tau_inv.u   ≈ EXPECTED_TAU_INV.u
    @test result.tau_inv.d   ≈ EXPECTED_TAU_INV.d
    @test result.tau_inv.s   ≈ EXPECTED_TAU_INV.s
    @test result.tau_inv.ubar ≈ EXPECTED_TAU_INV.ubar
    @test result.tau_inv.dbar ≈ EXPECTED_TAU_INV.dbar
    @test result.tau_inv.sbar ≈ EXPECTED_TAU_INV.sbar

    @test result.tau.u   ≈ EXPECTED_TAU.u
    @test result.tau.d   ≈ EXPECTED_TAU.d
    @test result.tau.s   ≈ EXPECTED_TAU.s
    @test result.tau.ubar ≈ EXPECTED_TAU.ubar
    @test result.tau.dbar ≈ EXPECTED_TAU.dbar
    @test result.tau.sbar ≈ EXPECTED_TAU.sbar
end
