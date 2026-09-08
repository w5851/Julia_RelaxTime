using Test
include(joinpath(@__DIR__,"..","..","..","scripts","analysis","relaxtime","causal_gbu_continuation_gate.jl"))
const CG=CausalGBUContinuationGate

@testset "Mapped original-cut quadrature retains edges and jumps" begin
    rho(x)=sqrt(max(x,0.))*(1-x)
    @test CG.mapped_cauchy(rho,[0.,1.],0.;nodes=48)≈4/(3pi) atol=1e-12
    @test CG.mapped_cauchy(rho,[0.,1.],1.;nodes=48)≈-2/(3pi) atol=1e-12
    step(x)=x<1 ? 0.4 : -0.2
    for z in (0.3im,2.1,3.,-1.)
        exact=(0.4*(log(complex(1)-z)-log(complex(0)-z))-
            0.2*(log(complex(2)-z)-log(complex(1)-z)))/pi
        @test CG.mapped_cauchy(step,[0.,1.,2.],z;nodes=64)≈exact atol=1e-12
    end
    for x in (0.4,1.5)
        r=CG.mapped_cauchy(step,[0.,1.,2.],x;nodes=64)
        @test imag(r)≈step(x) atol=1e-14
    end
    @test_throws ArgumentError CG.mapped_integral(identity,[0.,0.,1.])
    # Distinct physical panels must survive even when Float64 has no interior.
    a=1.;b=nextfloat(a);c=nextfloat(b)
    @test real(CG.mapped_integral(x->x<BigFloat(b) ? 2. : -1.,[a,b,c]))≈b-a rtol=1e-12
end

@testset "Wide-coordinate cut agrees away from roundoff-scale panels" begin
    bg=(m=(u=1.1,d=1.3,s=1.9),mu=(u=0.2,d=0.3,s=0.1),T=0.8,Phi=0.3,PhiBar=0.4,vacuum=2.5)
    for ch in CG.R.CHANNELS,q in (0.,0.6),x in (-5.,-0.4,0.2,3.5,7.)
        @test CG.wide_cut(bg,ch,q,BigFloat(x),4.;nodes=32)≈CG.A.direct_cut(bg,ch,q,x,4.;nodes=32).imaginary atol=1e-12
    end
end

@testset "A retained source proof cannot silently cross Fermi onset" begin
    b=(m=(u=0.4,d=0.41,s=2.),mu=(u=0.37,d=0.43,s=0.01))
    @test !CG.source_domain(b,:pi_plus).previous_pair_source_proof_applicable
    @test CG.source_domain(b,:K_plus).previous_pair_source_proof_applicable
    @test !CG.source_domain(b,:K_plus).previous_full_source_proof_applicable
end
