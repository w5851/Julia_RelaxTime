using Test
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime", "ThermalWeightedSampling.jl"))
if !isdefined(Main,:StreamedThermalValidation)
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime", "StreamedThermalValidation.jl"))
end
using .ThermalWeightedSampling, .ThermalAverageCandidates, .StreamedThermalValidation

@testset "Streamed outer-grid integrals retain actual weights and support" begin
    for process in (:uubar_to_uubar, :ud_to_ud, :us_to_us), xi in (-0.43, 0.36)
        state = (quark_params=(m=(u=0.8,d=0.8,s=1.7), μ=(u=0.2,d=0.2,s=0.7)),
            thermo_params=(T=0.9,Φ=0.3,Φbar=0.6,ξ=xi), K_coeffs=(unused=0.0,))
        plan = prepare_thermal_plan(process, state; p_nodes=11, angle_nodes=5, phi_nodes=9)
        bounds = ThermalAverageCandidates.process_bounds(process, state)
        measure = thermal_measure(plan, bounds)
        f(s) = 0.2 + 0.03s + 2exp(-((s-8.13)/0.17)^2)
        g(s) = 0.19 + 0.031s
        query(s) = (f(s), g(s), abs(f(s)-g(s)))
        result = streamed_integrals(plan, query, bounds, Val(3))
        expected = (integrate_measure(measure,f), integrate_measure(measure,g),
            integrate_measure(measure,s->abs(f(s)-g(s))))
        @test all(isapprox.(result.values, expected; rtol=2e-13))
        unfolded = streamed_integrals(plan, query, bounds, Val(3); fold_equal_mass=false)
        @test all(isapprox.(result.values, unfolded.values; rtol=2e-13))
        support = thermal_support(plan, bounds)
        @test support.lower == first(measure.nodes)
        @test support.upper == last(measure.nodes)
        @test support.events == result.events == measure.events
        @test_throws DomainError streamed_integrals(plan, s->(NaN,), bounds, Val(1))
        comparison=TableComparison(f,g,f)
        @test @inferred(comparison(2.3)) == (f(2.3),g(2.3),abs(f(2.3)-g(2.3)),abs(f(2.3)-g(2.3)))
        compared=streamed_integrals(plan,comparison,bounds,Val(4))
        @test isapprox(compared.values[3],expected[3];rtol=2e-13)
    end
end

@testset "Common table intervals preserve all reference and L1 queries" begin
    a=WeightedSigmaTable([1e-8,0.01,0.2,0.6,1.1,2.0],[1e6,24.0,1.4,0.3,2.0,0.7],0.0,0.6)
    b=WeightedSigmaTable([1e-8,0.003,0.13,0.5,0.9,1.4,2.0],[8e5,31.0,2.3,0.1,1.2,1.4,0.8],0.0,0.0)
    c=WeightedSigmaTable([1e-8,0.04,0.3,0.8,1.7,2.0],[7e5,7.0,1.0,0.9,1.1,0.9],0.0,0.3)
    raw=TableComparison(a,b,c);compiled=compile_comparison(raw)
    queries=sort!(unique(vcat(a.nodes,b.nodes,c.nodes,
        collect(10.0.^range(-8,log10(2.0);length=101)),collect(range(0.01,1.99;length=71)))))
    for s in queries
        @test all(isapprox.(compiled(s),raw(s);rtol=3e-13,atol=2e-13))
    end
    @test @inferred(compiled(0.123)) isa NTuple{4,Float64}
    @test_throws DomainError compiled(0.0)
    @test_throws ArgumentError compile_comparison(TableComparison(identity,identity,identity))
end
