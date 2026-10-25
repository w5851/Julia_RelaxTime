using Test
if !isdefined(Main,:ThresholdCacheCandidates)
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime", "ThresholdCacheCandidates.jl"))
end
if !isdefined(Main,:ThermalAverageCandidates)
    include(joinpath(@__DIR__, "..", "..", "..", "scripts", "analysis", "relaxtime", "ThermalAverageCandidates.jl"))
end
using .ThresholdUpgradeCommon
using .ThresholdCacheCandidates
using .ThermalAverageCandidates

@testset "Prepared thermal quadrature: same nodes and physical boundaries" begin
    for xi in (-0.4, 0.0, 0.4), mu in (0.0, 0.7)
        state = (quark_params=(m=(u=1.2,d=1.2,s=2.0),μ=(u=mu,d=mu,s=mu)),
            thermo_params=(T=0.8,Φ=0.45,Φbar=0.6,ξ=xi),K_coeffs=(unused=0.0,))
        for process in (:uubar_to_uubar, :uubar_to_ssbar, :us_to_us, :ssbar_to_uubar)
            query(s) = 0.3 + 0.2exp(-((s-12)/0.7)^2)
            old = functional_rate(query, process, state; p_nodes=12, angle_nodes=4, phi_nodes=8)
            for velocity in (:legacy, :equal_mass)
                new = prepared_rate(query, process, state; velocity=velocity,
                    p_nodes=12, angle_nodes=4, phi_nodes=8)
                @test isapprox(new, old; rtol=1e-12, atol=1e-14)
            end
        end
        processes = (:uubar_to_uubar, :uubar_to_ssbar)
        qs = (s->0.4, s->0.6)
        joint = shared_rates(qs, processes, state; p_nodes=12, angle_nodes=4, phi_nodes=8)
        for k in 1:2
            separate = functional_rate(qs[k], processes[k], state; p_nodes=12, angle_nodes=4, phi_nodes=8)
            @test isapprox(joint[k], separate; rtol=1e-12, atol=1e-14)
        end
    end
    state=(quark_params=(m=(u=1.0,d=1.0,s=2.0),μ=(u=0.0,d=0.0,s=0.0)),
        thermo_params=(T=0.8,Φ=0.5,Φbar=0.5,ξ=0.0),K_coeffs=(unused=0.0,))
    @test_throws ArgumentError shared_rates((identity,identity),(:uu_to_uu,:us_to_us),state)
    @test_throws ArgumentError shared_rates((identity,),(:uu_to_uu,:uu_to_uu),state)
    @test_throws ArgumentError prepared_rate(identity,:uu_to_uu,state;velocity=:unknown)
end

@testset "Union-grid projection preserves raw and weighted interpolation" begin
    state=(quark_params=(m=(u=1.0,d=1.0,s=2.0),μ=(u=0.3,d=0.3,s=0.3)),
        thermo_params=(T=0.8,Φ=0.45,Φbar=0.6,ξ=0.3),K_coeffs=(unused=0.0,))
    process=:uubar_to_uubar
    plan=prepare_thermal_plan(process,state;p_nodes=12,angle_nodes=4,phi_nodes=8)
    s0=4.0;split=5.1
    # Distinct knots, threshold fallback and a change of representation.
    q1=s->s<split ? 0.2+0.6/(s-s0) : 0.2+0.6/(split-s0)+0.01*(s-split)
    q2=s->s<7.0 ? 0.4+0.03*s : 0.4+0.03*7+0.01*(s-7)
    primary=s->0.7
    grid=[4.1,4.5,5.1,6.2,7.0,12.0,48.0]
    bounds=ThermalAverageCandidates.process_bounds(process,state)
    projected=project_linear_rates(plan,grid,s0,primary,(q1,q2),(split,0.0),bounds)
    exact=shared_rates((primary,q1,q2),(process,process,process),state;p_nodes=12,angle_nodes=4,phi_nodes=8)
    @test all(isapprox.(projected,exact;rtol=1e-12,atol=1e-14))
    @test_throws ArgumentError project_linear_rates(plan,reverse(grid),s0,primary,(q1,q2),(split,0.0),bounds)
end
