using Test
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","causal_gbu_infinite_qgate.jl"))
const testcausalgbuinfiniteqgate_IQ=CausalGBUInfiniteQGate

@testset "Geometric gaps and continuous finite-q spectral nodes" begin
    bg=(m=(u=1.1,d=1.3,s=1.9),mu=(u=0.2,d=0.3,s=0.1),T=0.8,Phi=0.3,PhiBar=0.4,
        vacuum=2.5,coupling=Dict(:K_plus=>0.05))
    for q in (0.,0.6,1.4,3.2),ch in testcausalgbuinfiniteqgate_IQ.R.CHANNELS
        k=testcausalgbuinfiniteqgate_IQ.I.kernel(bg,ch,q)
        for x in (k.landau,(k.landau+k.threshold)/2,k.threshold)
            @test iszero(k.rho(x))
            @test iszero(k.rho(-x))
        end
    end
    k=testcausalgbuinfiniteqgate_IQ.I.kernel(bg,:K_plus,0.6;split_inv_fm=36.)
    p=testcausalgbuinfiniteqgate_IQ.P.profile(k;mesh=128)
    @test all(p.yright[i]==p.yleft[i+1] for i in 1:length(p.left)-1)
    for x in k.edges[2:end-1]
        @test isfinite(testcausalgbuinfiniteqgate_IQ.P.polarization(p,x))
    end
    t=testcausalgbuinfiniteqgate_IQ.topology(p;scan_nodes=16)
    @test t.passed
    @test t.positive_count==t.negative_count==0
end

@testset "Adaptive gap contour retains a synthetic near-endpoint root" begin
    r=1-2e-8
    f(z)=z-r
    old=testcausalgbuinfiniteqgate_IQ.R.contour_count(f,1e-7,1-1e-7,-0.02,0.02;nodes=32,max_nodes=2048)
    @test old.passed
    @test old.count==0 # The fixed window genuinely excludes the known zero.
    window=testcausalgbuinfiniteqgate_IQ.Y.Limits.gap_window(0.,1.,[r])
    adaptive=testcausalgbuinfiniteqgate_IQ.R.contour_count(f,window.left,window.right,-0.02,0.02;nodes=32,max_nodes=2048)
    @test adaptive.passed
    @test adaptive.count==1
    @test adaptive.max_step<pi/2
end
