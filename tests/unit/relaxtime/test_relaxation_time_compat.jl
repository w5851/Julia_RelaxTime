using Test

const _RELAX_TIME_COMPAT_MODULE_PATH = normpath(joinpath(@__DIR__, "..", "..", "..", "src", "relaxtime", "RelaxTime.jl"))
if !isdefined(Main, :RelaxTime)
    Base.include(Main, _RELAX_TIME_COMPAT_MODULE_PATH)
end
const RelaxationTimeCompatModule = Main.RelaxationTime

const COMPAT_RATES_SAMPLE = (
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

@testset "RelaxationTime compat aliases" begin
    @testset "rate_lookup deprecation alias" begin
        @test_deprecated RelaxationTimeCompatModule.rate_lookup(COMPAT_RATES_SAMPLE, :dubar_to_dubar)
        @test RelaxationTimeCompatModule.rate_lookup(COMPAT_RATES_SAMPLE, :dubar_to_dubar) == RelaxationTimeCompatModule.rate_value(COMPAT_RATES_SAMPLE, :dubar_to_dubar)
    end

    @testset "density_lookup retired" begin
        @test !isdefined(RelaxationTimeCompatModule, :density_lookup)
    end
end