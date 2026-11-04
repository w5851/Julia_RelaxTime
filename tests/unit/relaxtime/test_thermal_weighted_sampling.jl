using Test
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","ThermalWeightedSampling.jl"))
using .ThermalWeightedSampling
using .ThermalAverageCandidates
using .ThresholdUpgradeCommon

@testset "Discrete thermal measure versus unchanged quadrature" begin
    for process in (:uubar_to_uubar,:us_to_us), xi in (-0.43,0.36)
        state=(quark_params=(m=(u=0.8,d=0.8,s=1.7),μ=(u=0.2,d=0.2,s=0.7)),
            thermo_params=(T=0.9,Φ=0.3,Φbar=0.6,ξ=xi),K_coeffs=(unused=0.0,))
        plan=prepare_thermal_plan(process,state;p_nodes=10,angle_nodes=5,phi_nodes=9)
        bounds=ThermalAverageCandidates.process_bounds(process,state)
        measure=thermal_measure(plan,bounds)
        f=s->0.3+0.1s+3exp(-((s-9.13)/0.2)^2)
        oracle=only(integrate_queries(plan,(f,),(bounds,),Val(:legacy)))
        @test isapprox(integrate_measure(measure,f),oracle;rtol=2e-12,atol=1e-14)
        @test length(measure.nodes)<=measure.events
        unfolded=thermal_measure(plan,bounds;fold_equal_mass=false)
        @test isapprox(integrate_measure(measure,f),integrate_measure(unfolded,f);rtol=2e-13)
        @test measure.events<=unfolded.events
        @test_throws ArgumentError thermal_measure(plan,bounds;max_events=1)
    end
end

@testset "Width-aware discovery and thermal error for translated narrow resonances" begin
    for location in (0.013,0.227,1.217), width in (1e-7,2e-4)
        nodes=sort!(unique(vcat(collect(range(1e-7,3.0;length=5001)),
            location .+ width.*range(-24,24;length=193))))
        edges=vcat(first(nodes),(nodes[1:end-1].+nodes[2:end])./2,last(nodes))
        weights=diff(edges).*exp.(-nodes)
        measure=ThermalMeasure(nodes,weights)
        f=s->0.2+width/((s-location)^2+width^2)
        discovery=real_axis_feature_seeds(s->complex(s-location,width),collect(range(1e-7,3.0;length=65)))
        @test any(r->isapprox(r.location,location;atol=1e-13),discovery.features)
        built=build_thermal_weighted_table(f,measure,0.0;rtol=2e-4,extra_seeds=discovery.seeds)
        @test built.converged
        @test built.guard_phase==2
        @test abs(built.rate/integrate_measure(measure,f)-1)<=2e-4
        @test isapprox(built.rate,integrate_measure(measure,validated_table(built));rtol=2e-12)
        @test built.estimated_error<=1e-12+2e-4abs(built.rate)
    end
end

@testset "Shifted audit detects a feature invisible to quarter probes" begin
    nodes=collect(range(0.01,2.0;length=6001))
    measure=ThermalMeasure(nodes,fill(1/6001,6001))
    seeds=thermal_seed_grid(measure,0.0;order=8)
    a,b=seeds[end-4],seeds[end-3]
    first_i=searchsortedlast(nodes,a)+1;last_i=searchsortedlast(nodes,b)
    median=(nodes[(first_i+last_i)÷2]-a)/(b-a)
    zeros_at=(0.0,0.25,0.5,0.75,1.0,median)
    f=s->a<s<b ? 0.2+1e10prod(((s-a)/(b-a)-z)^2 for z in zeros_at) : 0.2
    built=build_thermal_weighted_table(f,measure,0.0;rtol=1e-4,seed_order=8)
    @test built.converged
    @test any(h->h.guard_phase>0 && h.estimated_error>h.tolerance,built.history)
    @test abs(built.rate/integrate_measure(measure,f)-1)<1e-4
end

@testset "Budget exhaustion cannot become an accepted cache" begin
    measure=ThermalMeasure(collect(range(0.01,2.0;length=1001)),ones(1001))
    result=build_thermal_weighted_table(s->1+exp(-((s-0.13)/0.01)^2),measure,0.0;max_evaluations=10)
    @test !result.converged
    @test_throws ArgumentError validated_table(result)
    @test_throws DomainError build_thermal_weighted_table(s->NaN,measure,0.0)
    @test_throws DomainError build_thermal_weighted_table(s->1e308,measure,0.0)
    @test_throws ArgumentError build_thermal_weighted_table(identity,measure,0.0;rtol=-1.0)
    @test_throws ArgumentError ThermalMeasure([1.0,2.0],[1.0,-1.0])
    zero=build_thermal_weighted_table(s->0.0,measure,0.0)
    @test zero.converged && zero.rate==0.0
    @test_throws DomainError validated_table(zero)(2.1)
end

@testset "Known regulator singularities split feature discovery" begin
    grid=[0.0,0.2,0.5,1.0,1.5,2.0]
    visited=Float64[]
    denominator=s->begin
        push!(visited,s)
        s==1 && error("singular endpoint must not be evaluated")
        complex(log(abs(s-1))+8,1e-3)
    end
    found=real_axis_feature_seeds(denominator,grid;breakpoints=[1.0])
    @test all(!=(1.0),visited)
    for root in (1-exp(-8),1+exp(-8))
        @test any(r->isapprox(r.location,root;atol=1e-12),found.features)
    end
    @test all(r->abs(r.location-1)>1e-5,found.features)
    pole=real_axis_feature_seeds(s->complex((s-0.3)/(s-1),1e-4),grid;breakpoints=[1.0])
    @test length(pole.features)==1
    @test isapprox(only(pole.features).location,0.3;atol=1e-12)
    @test all(isfinite,pole.seeds)
end
