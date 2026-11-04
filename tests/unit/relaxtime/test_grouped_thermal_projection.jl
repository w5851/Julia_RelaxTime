using Test
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime", "GroupedThermalProjection.jl"))
using .ThermalAverageCandidates
using .GroupedThermalProjection

@testset "Grouped angular geometry preserves the full discrete integral" begin
    for xi in (-0.4,0.4), process in (:uubar_to_uubar,:us_to_us), custom in (false,true)
        state=(quark_params=(m=(u=1.0,d=1.0,s=2.0),μ=(u=0.3,d=0.3,s=0.3)),
            thermo_params=(T=0.8,Φ=0.45,Φbar=0.6,ξ=xi),K_coeffs=(unused=0.0,))
        extras=custom ? (cos_grid=[-0.71,-0.1,0.0,0.32,0.84],cos_w=[0.2,0.4,0.6,0.3,0.5],
            phi_grid=[0.07,0.7,pi,2pi-0.7],phi_w=[0.4,1.0,2.0,0.8]) : (;)
        plan=prepare_thermal_plan(process,state;p_nodes=12,angle_nodes=5,phi_nodes=9,extras...)
        groups=angular_groups(plan)
        @test sum(length(g.pairs) for g in groups)==length(plan.cw)^2
        @test isapprox(sum(length(g.pairs)*sum(g.weights) for g in groups),
            length(plan.cw)^2*sum(plan.phiw);rtol=1e-14)
        s0=(plan.mi+plan.mj)^2; split=s0+1.1
        bounds=ThermalAverageCandidates.process_bounds(process,state)
        q1=s->s<split ? 0.2+0.6/(s-s0) : 0.2+0.6/(split-s0)+0.01*(s-split)
        q2=s->0.4+0.03*s
        primary=s->0.7+2exp(-((s-s0-0.07)/0.02)^2)
        grid=[s0+0.1,s0+0.5,split,s0+2.0,bounds[2]-0.1]
        full=project_linear_rates(plan,grid,s0,primary,(q1,q2),(split,0.0),bounds)
        grouped=grouped_project_linear_rates(plan,grid,s0,primary,(q1,q2),(split,0.0),bounds;groups=groups)
        @test all(isapprox.(grouped,full;rtol=1e-12,atol=1e-14))
        exact=integrate_queries(plan,(primary,q1,q2),(bounds,bounds,bounds),Val(:equal_mass))
        @test all(isapprox.(grouped,exact;rtol=1e-12,atol=1e-14))
        @test isapprox(first(grouped),compensated_primary_rate(plan,primary,bounds);rtol=1e-13,atol=1e-14)
    end
end
