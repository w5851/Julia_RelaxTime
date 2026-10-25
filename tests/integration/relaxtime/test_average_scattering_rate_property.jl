using Test
include(joinpath(@__DIR__,"thermal_sampling_fixture.jl"))

@testset "Struct and NamedTuple cache construction across channels and anisotropy" begin
    as=Main.RelaxTime.AverageScatteringRate
    for xi in (-0.4,0.0,0.4)
        f=thermal_sampling_fixture(xi)
        q=Main.ParameterTypes.QuarkParams(f.qp)
        t=Main.ParameterTypes.ThermoParams(f.tp)
        for process in (:uu_to_uu,:us_to_us,:ssbar_to_uubar)
            cache=as.build_thermal_weighted_cache(process,q,t,f.kc;f.options...)
            rate=as.average_scattering_rate(process,f.qp,f.tp,f.kc;f.options...,cs_cache=cache)
            @test isfinite(rate) && rate>0
            @test rate==cache.rate
            @test as.average_scattering_rate(process,q,t,f.kc;f.options...,cs_cache=cache)==rate
            @test cache.estimated_error<=f.options.sampling_config.atol+f.options.sampling_config.rtol*rate
        end
    end
end
