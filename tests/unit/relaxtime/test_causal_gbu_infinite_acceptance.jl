using Test
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","causal_gbu_infinite_acceptance.jl"))
const testcausalgbuinfiniteacceptance_IA=CausalGBUInfiniteAcceptance

@testset "Real-axis Nyquist includes the near-axis strip" begin
    for n in (16,32,64)
        stable=testcausalgbuinfiniteacceptance_IA.nyquist(z->(z-0.3)/(z+2im),[-1.,0.,1.],[(0.29,0.31)],4.;nodes=n)
        @test stable.passed && stable.count==0
        for height in (0.5,0.01,0.0001)
            # Place a partition at the narrow instability; no eta floor can hide it.
            f(z)=(z-0.2-height*im)/(z+2im)
            r=testcausalgbuinfiniteacceptance_IA.nyquist(f,[-1.,0.,0.2-height,0.2,0.2+height,1.],Tuple{Float64,Float64}[],4.;nodes=n)
            @test r.passed && r.count==1
        end
    end
    @test_throws ArgumentError testcausalgbuinfiniteacceptance_IA.nyquist(identity,[-1.,1.],[(-5.,-4.)],3.)
end
