using Test
include(joinpath(@__DIR__,"thermal_sampling_fixture.jl"))

@testset "Native cache support, local queries and adapter boundaries" begin
    as=Main.RelaxTime.AverageScatteringRate
    f=thermal_sampling_fixture(-0.4)
    process=:ssbar_to_uubar
    cache=as.build_thermal_weighted_cache(process,f.qp,f.tp,f.kc;f.options...)
    @test first(cache.context.bounds)==(2f.qp.m.s)^2
    @test first(cache.table.nodes)>first(cache.context.bounds)
    @test last(cache.table.nodes)<last(cache.context.bounds)
    @test all(isfinite,cache.table.sigma) && all(>=(0),cache.table.sigma)
    @test as.get_sigma(cache,first(cache.context.bounds))==0
    @test as.get_sigma(cache,last(cache.context.bounds))==0
    s=cache.table.nodes[length(cache.table.nodes)÷2]
    direct=Main.RelaxTime.TotalCrossSection.total_cross_section(process,s,f.qp,f.tp,f.kc;n_points=4)
    @test as.get_sigma(cache,s;direct=true)==direct
    @test_throws ArgumentError as.get_sigma(cache,NaN)
    @test_throws ArgumentError as.average_scattering_rate(process,f.qp,merge(f.tp,(ξ=-1.0,)),f.kc;f.options...)
    q=Main.ParameterTypes.QuarkParams(f.qp)
    for (particle,flavor) in ((:u,:u),(:ubar,:u),(:d,:d),(:sbar,:s))
        @test as.get_mass(particle,q)==getproperty(f.qp.m,flavor)
        @test as.get_mu(particle,q)==getproperty(f.qp.μ,flavor)
    end
end
