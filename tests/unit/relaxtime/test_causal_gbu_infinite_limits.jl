using Test
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","causal_gbu_infinite_limits.jl"))
const testcausalgbuinfinitelimits_IL=CausalGBUInfiniteLimits

@testset "Threshold phase is an analytic limit, independent of counting" begin
    for f in (-2.,-0.01,0.01,2.),A in (-2.,-0.1,0.1,2.)
        limit=testcausalgbuinfinitelimits_IL.threshold_limit(f,A)
        @test !limit.count_used
        @test limit.phase isa Float64
        @test limit.inverse_margin>0
        previous=Inf
        for u in (1e-8,1e-12,1e-16,1e-20)
            actual=-angle(complex(f,-A*sqrt(u)))
            delta=abs(actual-limit.phase)
            @test delta<=previous
            previous=delta
        end
        @test previous<1e-7
    end
    @test_throws ArgumentError testcausalgbuinfinitelimits_IL.threshold_limit(0.,1.)
    @test_throws ArgumentError testcausalgbuinfinitelimits_IL.threshold_limit(1.,0.)
    @test_throws ArgumentError testcausalgbuinfinitelimits_IL.threshold_limit(1+0.1im,1.)
end

@testset "Gap contour must not cut off a near-Mott real root" begin
    w=testcausalgbuinfinitelimits_IL.gap_window(0.,1.,[1e-8,0.5,1-2e-8])
    @test 0<w.left<1e-8
    @test 1-2e-8<w.right<1
    @test w.left_margin≈1.25e-9
    @test testcausalgbuinfinitelimits_IL.gap_window(0.,1.,Float64[]).left==1e-7
    @test_throws ArgumentError testcausalgbuinfinitelimits_IL.gap_window(0.,1.,[0.])
    @test_throws ArgumentError testcausalgbuinfinitelimits_IL.gap_window(0.,1.,[1.])
end
