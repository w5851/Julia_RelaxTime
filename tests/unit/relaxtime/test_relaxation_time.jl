using Test

# Load all relaxtime modules via RelaxTime.jl entry point
const _RELAX_TIME_MODULE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAX_TIME_MODULE_PATH)
end
using Main.RelaxationTime
using Main.ParameterTypes: QuarkParams, ThermoParams

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
    ubar = 79.0,
    dbar = 79.0,
    sbar = 266.0,
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

@testset "rate_value alias table" begin
    @test RelaxationTime.rate_value(RATES_SAMPLE, :dubar_to_dubar) == RATES_SAMPLE.udbar_to_udbar
    @test RelaxationTime.rate_value(RATES_SAMPLE, :subar_to_subar) == RATES_SAMPLE.usbar_to_usbar
    @test RelaxationTime.rate_value(RATES_SAMPLE, :ubardbar_to_ubardbar) == RATES_SAMPLE.ud_to_ud

    rates_dict = Dict{Symbol, Float64}(pairs(RATES_SAMPLE))
    @test RelaxationTime.rate_value(rates_dict, :ubarubar_to_ubarubar) == RATES_SAMPLE.uu_to_uu
    @test RelaxationTime.rate_value(rates_dict, :ubarsbar_to_ubarsbar) == RATES_SAMPLE.us_to_us
    @test RelaxationTime.rate_value(rates_dict, :sbarsbar_to_sbarsbar) == RATES_SAMPLE.ss_to_ss
end

@testset "relaxation_rates warns and clamps negative totals" begin
    densities = (u=1.0, d=1.0, s=1.0, ubar=1.0, dbar=1.0, sbar=1.0)
    rates_negative = (; (k => -abs(v) for (k, v) in pairs(RATES_SAMPLE))...)

    @test_logs (:warn, r"negative relaxation rate encountered; clamping to 0") begin
        tau_inv = relaxation_rates(densities, rates_negative)
        @test tau_inv.u >= 0.0
        @test tau_inv.s >= 0.0
        @test tau_inv.ubar >= 0.0
        @test tau_inv.sbar >= 0.0
    end
end

@testset "relaxation adapters normalize Dict inputs" begin
    densities_dict = Dict{Symbol, Float64}(pairs(DENSITIES_SAMPLE))
    rates_dict = Dict{Symbol, Float64}(pairs(RATES_SAMPLE))

    tau_inv = relaxation_rates(densities_dict, rates_dict)
    @test tau_inv == EXPECTED_TAU_INV

    result = relaxation_times(
        QUARK_PARAMS,
        THERMO_PARAMS,
        K_COEFFS;
        densities=densities_dict,
        existing_rates=rates_dict,
    )
    @test result.tau_inv == EXPECTED_TAU_INV
    @test result.tau == EXPECTED_TAU
end

@testset "relaxation_rates accepts integer-compatible inputs" begin
    densities_int = Dict{Symbol, Int}(k => Int(round(v)) for (k, v) in pairs(DENSITIES_SAMPLE))
    rates_int = Dict{Symbol, Int}(k => Int(round(v)) for (k, v) in pairs(RATES_SAMPLE))

    tau_inv = relaxation_rates(densities_int, rates_int)
    @test tau_inv.u == EXPECTED_TAU_INV.u
    @test tau_inv.s == EXPECTED_TAU_INV.s
    @test tau_inv.ubar == EXPECTED_TAU_INV.ubar
    @test tau_inv.sbar == EXPECTED_TAU_INV.sbar
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

@testset "relaxation_times struct/NamedTuple equivalence on existing rates" begin
    q_struct = QuarkParams(QUARK_PARAMS)
    t_struct = ThermoParams(THERMO_PARAMS)

    res_nt = relaxation_times(
        QUARK_PARAMS,
        THERMO_PARAMS,
        K_COEFFS;
        densities=DENSITIES_SAMPLE,
        existing_rates=RATES_SAMPLE,
    )

    res_struct = relaxation_times(
        q_struct,
        t_struct,
        K_COEFFS;
        densities=DENSITIES_SAMPLE,
        existing_rates=RATES_SAMPLE,
    )

    @test res_struct.tau_inv == res_nt.tau_inv
    @test res_struct.tau == res_nt.tau
    @test res_struct.rates == res_nt.rates
end

@testset "relaxation_times forwards advanced ASR kwargs" begin
    result = relaxation_times(
        QUARK_PARAMS,
        THERMO_PARAMS,
        K_COEFFS;
        densities=DENSITIES_SAMPLE,
        existing_rates=RATES_SAMPLE,
        threshold_subtraction=true,
        asym_window=0.25,
        asym_fit_min_points=6,
        asym_extra_points=9,
        interpolation_mode=:linear,
    )

    @test result.tau_inv == EXPECTED_TAU_INV
    @test result.tau == EXPECTED_TAU
end
