using Test
include(joinpath(@__DIR__, "LoopQuadratureReference.jl"))
include(joinpath(LOOP_REFERENCE_ROOT,"scripts","analysis","relaxtime","ThresholdCacheCandidates.jl"))
using .LoopQuadratureReference, .ThresholdCacheCandidates
using QuadGK

@testset "Stable angular PV moment and independent nonsingular integral" begin
    for B in (-2.3, 0.7), ratio in (0.0, 1e-9, 1e-5, 0.12, 0.6)
        A=ratio*B
        oracle,_=quadgk(c->c*c/(A*c+B),-1.0,1.0;rtol=1e-12)
        @test isapprox(angular_pv_second_moment(A,B),oracle;rtol=2e-12)
        zeroth,_=quadgk(c->1/(A*c+B),-1.0,1.0;rtol=1e-12)
        @test isapprox(A==0 ? 2/B : angular_log(A,B)/A,zeroth;rtol=2e-12)
    end
end

@testset "Independent k0 reference matches the paired threshold identity" begin
    tp=(T=0.82,Φ=0.38,Φbar=0.47,ξ=0.31)
    m=0.73;mu=0.29
    plan=PairLoopPlan(m,mu,tp;n_nodes=128)
    for epsilon in (-1e-3,-1e-6,-1e-9,0.0,1e-9,1e-6,1e-3,0.2)
        s=4m^2*(1+epsilon)
        ref=loop_reference(sqrt(s),0.0,m,m,mu,mu,tp)
        @test isfinite(ref.value) && isfinite(ref.error)
        @test abs(ref.value-pair_b0_pv(plan,s))<max(2e-8,4ref.error)
    end
end

@testset "Unequal-mass symmetry and finite-k limit" begin
    for xi in (0.0,0.31)
        tp=(T=0.82,Φ=0.38,Φbar=0.47,ξ=xi)
        a=loop_reference(2.6,0.0,0.73,1.41,0.29,0.57,tp)
        b=loop_reference(-2.6,0.0,1.41,0.73,0.57,0.29,tp)
        @test isapprox(a.value,b.value;rtol=2e-12,atol=2e-12)
        at_zero=loop_reference(1.3,0.0,0.73,0.73,0.29,0.29,tp)
        near_zero=loop_reference(1.3,1e-4,0.73,0.73,0.29,0.29,tp)
        @test isapprox(at_zero.value,near_zero.value;rtol=1e-7,atol=2e-9)
        coarse=loop_reference(0.43,1.21,0.73,1.41,0.29,0.57,tp;rtol=2e-8)
        fine=loop_reference(0.43,1.21,0.73,1.41,0.29,0.57,tp;rtol=2e-10)
        @test isfinite(fine.value) && isfinite(fine.error)
        @test abs(coarse.value-fine.value)<max(2e-8,4coarse.error)
    end
end
