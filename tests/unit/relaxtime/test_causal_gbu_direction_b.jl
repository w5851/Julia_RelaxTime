using Test,ForwardDiff
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","causal_gbu_direction_b.jl"))
const testcausalgbudirectionb_B=CausalGBUDirectionB

@testset "Static thermal extension is a differentiable potential" begin
    for phi in (0.,0.3,1.),bar in (0.,0.4,1.),x in (-20.,-1.,0.,1.,20.)
        @test ForwardDiff.derivative(y->testcausalgbudirectionb_B.log_partition(y,phi,bar),x)≈
            -3testcausalgbudirectionb_B.P.occupation(x,phi,bar) atol=1e-14
    end
    for mu in (-0.2,0.,0.4),tc in (3.,8.),phi in (0.3,1.)
        f(v)=testcausalgbudirectionb_B.static_flavor(v[1],v[2],0.8,phi,0.4,3.,tc;nodes=64).omega_inv_fm4
        x=[1.2,mu]
        d=ForwardDiff.gradient(f,x)
        h=ForwardDiff.hessian(f,x)
        state=testcausalgbudirectionb_B.static_flavor(x...,0.8,phi,0.4,3.,tc;nodes=64)
        @test d[1]≈state.condensate_inv_fm3 atol=1e-12
        @test -d[2]≈state.density_inv_fm3 atol=1e-12
        @test h[1,2]≈h[2,1] atol=1e-12
        @test !state.production_authorized
    end
    @test_throws ArgumentError testcausalgbudirectionb_B.static_flavor(1.,0.,0.,0.3,0.4,3.,8.)
end

@testset "Matched positive tail is not a full RPA stability proof" begin
    for c in (0.05,0.2),S in (4.,20.)
        f(s)=testcausalgbudirectionb_B.matched_linear_tail(s,c,S)
        @test f(0.)==0
        @test ForwardDiff.derivative(f,0.)==0
        @test ForwardDiff.derivative(s->ForwardDiff.derivative(f,s),0.)≈2c/(pi*S)
        us,ws=testcausalgbudirectionb_B.R.gauleg(0.,1.,128)
        for s in (-2.,1.0+0.3im,2.0+0.02im)
            direct=c*s^2/pi*sum(w/(S-s*u) for (u,w) in zip(us,ws))
            @test f(s)≈direct atol=1e-12
        end
        @test imag(f(1.5S))≈1.5c*S atol=1e-12
        @test isreal(f(-S)) && f(-S)>0
        pole=testcausalgbudirectionb_B.toy_uhp_pole(c,S,0.2)
        @test pole.residual<1e-12 && pole.imaginary_frequency_inv_fm>0
        @test !pole.is_project_pole && !pole.production_authorized
    end
    @test_throws ArgumentError testcausalgbudirectionb_B.matched_linear_tail(4.,0.1,4.)
    @test_throws ArgumentError testcausalgbudirectionb_B.toy_uhp_pole(0.,4.,0.2)
end
