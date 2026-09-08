using Test, LinearAlgebra

const _CSB_ROOT = normpath(joinpath(@__DIR__,"..","..",".."))
for (name,path) in ((:Constants_PNJL,"constants/Constants_PNJL.jl"),
                    (:GaussLegendre,"integration/GaussLegendre.jl"),
                    (:PNJLQuarkDistributions,"models/pnjl_physics/QuarkDistribution.jl"),
                    (:OneLoopIntegrals,"relaxtime/OneLoopIntegrals.jl"),
                    (:CausalSpectralBubble,"relaxtime/CausalSpectralBubble.jl"))
    isdefined(Main,name) || Base.include(Main,joinpath(_CSB_ROOT,"src",path))
end
using Main.CausalSpectralBubble: build_spectral_bubble, spectral_bubble, spectral_bubble_cut
using Main.CausalSpectralBubble: PiecewiseSpectralFunction, cauchy_transform, build_bubble_dispersion

@testset "Spectral dispersion: independent loop and kinematic gap" begin
    g = build_spectral_bubble(0.8,1.0,0.2,1.4,-0.1,0.8;Phi=0.3,PhiBar=0.4,
        momentum_nodes=96,angle_nodes=64)
    coarse = build_bubble_dispersion(g;segment_nodes=24)
    fine = build_bubble_dispersion(g;segment_nodes=48)
    for w in (0.0,0.4,2.1)
        target = spectral_bubble(g,w;eta_inv_fm=0.6).value
        ec = abs(cauchy_transform(coarse,w+0.3+0.6im)-target)
        ef = abs(cauchy_transform(fine,w+0.3+0.6im)-target)
        @test ef < ec
        @test ef < 5e-4
    end
    # Nearly equal charged flavors: coordinate roundoff must not smear a cut
    # across the entire analytic gap when a threshold is a mesh knot.
    for q in (0.0,0.5,1.0)
        near = build_spectral_bubble(q,1.6519019722802555,0.3742643948705271,
            1.649895213451757,0.4209954886862639,0.8615142220054972)
        p = build_bubble_dispersion(near;segment_nodes=16)
        lo,hi = hypot(q,near.m1_inv_fm-near.m2_inv_fm),hypot(q,near.m1_inv_fm+near.m2_inv_fm)
        for l in range(lo+1e-6,hi-1e-6;length=9)
            @test imag(cauchy_transform(p,l)) == 0
        end
    end
    for l in (0.3,3.0)
        kwargs = (Φ=0.3,Φbar=0.4)
        f(c) = Main.OneLoopIntegrals.B0_spectral_cut(l,0.8,1.0,0.2,1.4,-0.1,0.8;kwargs...,component=c).imaginary
        @test f(:full) ≈ f(:vacuum)+f(:thermal) atol=1e-12
    end
end

@testset "Cauchy/PV transform of a synthetic triangular spectrum" begin
    p = PiecewiseSpectralFunction([-1,0,1],[0,1,0])
    @test cauchy_transform(p,0.0) ≈ 1im atol=1e-15
    @test cauchy_transform(p,1im) ≈ im*(2atan(1)-log(2))/π atol=1e-15
    @test isfinite(real(cauchy_transform(p,1.0)))
    for w in (-0.7,-0.2,0.0,0.3,0.8)
        @test imag(cauchy_transform(p,w)) ≈ 1-abs(w) atol=1e-15
        @test cauchy_transform(p,w) ≈ -conj(cauchy_transform(p,-w)) atol=1e-14
        @test abs(cauchy_transform(p,w+1e-7im)-cauchy_transform(p,w)) < 2e-6
    end
    # Signed spectra are retained, not made positive by the transform.
    negative = PiecewiseSpectralFunction([-1,0,1],[0,-1,0])
    @test cauchy_transform(negative,0.4+0.2im) ≈ -cauchy_transform(p,0.4+0.2im)
    @test_throws ArgumentError cauchy_transform(p,-0.1im)
    @test_throws ArgumentError PiecewiseSpectralFunction([-1,0,1],[1,1,0])
    @test_throws ArgumentError PiecewiseSpectralFunction([-1,0,0],[0,1,0])
end

@testset "Cauchy cells with nearly coincident knots retain Float64 accuracy" begin
    function oracle(p,z)
        setprecision(256) do
            x,y=BigFloat.(p.energy),BigFloat.(p.imaginary)
            w=Complex{BigFloat}(z)
            j=searchsortedlast(x,real(w))
            rho=imag(w)==0 && 1<=j<length(x) ?
                y[j]+(y[j+1]-y[j])*(real(w)-x[j])/(x[j+1]-x[j]) : BigFloat(0)
            value=zero(w)
            for i in 1:length(x)-1
                a,b=x[i],x[i+1]
                slope=(y[i+1]-y[i])/(b-a)
                value+=slope*(b-a)
                if imag(w)>0
                    value+=(y[i]+slope*(w-a))*log((b-w)/(a-w))
                elseif real(w)!=a && real(w)!=b
                    value+=(y[i]+slope*(w-a)-rho)*log(abs((b-w)/(a-w)))
                end
            end
            if imag(w)==0 && rho!=0
                value+=rho*log(abs((last(x)-w)/(first(x)-w)))
            end
            return ComplexF64(value/big(pi)+rho*im)
        end
    end
    for width in (eps(1.0),1e-12,1e-7)
        p=PiecewiseSpectralFunction([-2.,-1.,-1.0+width,1.,1.0+width,2.],
            [0.,0.,0.7,-0.2,0.,0.])
        for z in (-3.,0.,0.4,3.7810782,1e6,0.4+0.2im,3.7810782+1e-8im,1e6+0.4im)
            @test cauchy_transform(p,z)≈oracle(p,z) atol=1e-12 rtol=1e-12
        end
    end
end

@testset "Causal bubble: Dirac trace and contact algebra" begin
    unit2 = Matrix{ComplexF64}(I,2,2)
    zero2 = zeros(ComplexF64,2,2)
    pauli = ([0 1;1 0],[0 -im;im 0],[1 0;0 -1])
    gamma0 = [unit2 zero2;zero2 -unit2]
    gammas = [[zero2 s;-s zero2] for s in pauli]
    gamma5 = im*gamma0*gammas[1]*gammas[2]*gammas[3]
    slash(E,p) = E*gamma0-sum(p[j]*gammas[j] for j in 1:3)
    m1,m2 = 1.1,1.6
    p,r = [0.3,-0.2,0.7],[-0.1,0.3,0.4]
    E1,E2 = hypot(norm(p),m1),hypot(norm(r),m2)
    q = norm(p-r)
    for channel in (:P,:S), s in (-1,1), t in (-1,1)
        vertex = channel === :P ? im*gamma5 : Matrix{ComplexF64}(I,4,4)
        explicit = -tr((slash(s*E1,p)+m1*I)*vertex*(slash(t*E2,r)+m2*I)*vertex)/2
        mt = channel === :P ? (m1-m2)^2 : (m1+m2)^2
        @test explicit ≈ (s*E1-t*E2)^2-q^2-mt atol=1e-14
    end
    # Arbitrary occupations: polynomial division fixes the contact, independently of PNJL.
    n1,n2 = (0.91,0.13),(0.83,0.07)
    z = 0.4+0.7im
    for mt in ((m1-m2)^2,(m1+m2)^2)
        b0,direct,moment0 = 0.0im,0.0im,0.0
        for (si,s) in enumerate((-1,1)), (ti,t) in enumerate((-1,1))
            u = s*E1-t*E2
            R = s*t*(n2[ti]-n1[si])
            b0 += R/(z-u)
            direct += R*(u^2-q^2-mt)/(z-u)
            moment0 += R
        end
        contact = 2*(E1*(n2[2]-n2[1])+E2*(n1[2]-n1[1]))
        @test abs(moment0) < 1e-15
        @test direct ≈ (z^2-q^2-mt)*b0-contact atol=1e-14
    end
end

@testset "Causal bubble: normalization, symmetry and units" begin
    for channel in (:P,:S), q in (0.0,0.8)
        g = build_spectral_bubble(q,1.0,0.2,1.4,-0.1,0.8;channel=channel,
            Phi=0.3,PhiBar=0.4,momentum_nodes=80,angle_nodes=48)
        rev = build_spectral_bubble(q,1.4,-0.1,1.0,0.2,0.8;channel=channel,
            Phi=0.3,PhiBar=0.4,momentum_nodes=80,angle_nodes=48)
        for w in (0.0,0.3,2.1)
            a = spectral_bubble(g,w;eta_inv_fm=0.4)
            b = spectral_bubble(rev,-w;eta_inv_fm=0.4)
            @test a.value ≈ conj(b.value) atol=2e-10
            @test a.contact_identity_residual < 1e-12
            @test !a.production_authorized
            h = 1e-5
            numerical = (spectral_bubble(g,w+h;eta_inv_fm=0.4).value-
                         spectral_bubble(g,w-h;eta_inv_fm=0.4).value)/(2h)
            @test a.derivative ≈ numerical rtol=1e-8 atol=1e-10
        end
        @test abs(spectral_bubble_cut(g,0.0).imaginary) < 1e-13
        gap = spectral_bubble(g,1.5;eta_inv_fm=0)
        @test imag(gap.value) == 0
        @test gap.analytic_scope === :open_gap
        @test_throws ArgumentError spectral_bubble(g,-0.3;eta_inv_fm=0)
        @test_throws ArgumentError spectral_bubble(g,0.0;eta_inv_fm=-0.1)
    end
    q,m1,mu1,m2,mu2,T = 0.8,1.0,0.2,1.4,-0.1,0.8
    a = build_spectral_bubble(q,m1,mu1,m2,mu2,T;vacuum_cutoff_inv_fm=3.0,thermal_cutoff_inv_fm=6.0)
    b = build_spectral_bubble(2q,2m1,2mu1,2m2,2mu2,2T;vacuum_cutoff_inv_fm=6.0,thermal_cutoff_inv_fm=12.0)
    av,bv = spectral_bubble(a,0.4;eta_inv_fm=0.3),spectral_bubble(b,0.8;eta_inv_fm=0.6)
    @test bv.value ≈ 4av.value rtol=1e-12
    @test bv.B0 ≈ av.B0 rtol=1e-12
    @test bv.derivative ≈ 2av.derivative rtol=1e-12
    @test b.contact_inv_fm2 ≈ 4a.contact_inv_fm2 rtol=1e-12

    # q=0 contact agrees with independently integrated A, including the thermal tail.
    nodes,weights = Main.GaussLegendre.gauleg(0.0,16.0,128)
    for mass in (1.0,1.5)
        g = build_spectral_bubble(0.0,mass,0.2,mass,0.2,0.8;
            Phi=0.3,PhiBar=0.4,thermal_cutoff_inv_fm=16.0,momentum_nodes=128)
        A = Main.OneLoopIntegrals.A(mass,0.2,0.8,0.3,0.4,nodes,weights)
        @test g.contact_inv_fm2 ≈ 2A atol=2e-12
        # Chiral-limit algebra at lambda -> 0, not a new equilibrium solve.
        pi0 = sum(g.pi_residues_inv_fm3[i]/(-g.poles_inv_fm[i]) for i in eachindex(g.poles_inv_fm))
        @test pi0 ≈ -Main.Constants_PNJL.N_color*A/(4π^2) atol=1e-12
        K = -π^2/(Main.Constants_PNJL.N_color*A)
        @test abs(1-4K*pi0) < 1e-12
    end
    for lambda in (0.1,0.5,0.9)
        g = build_spectral_bubble(1.0,1.0,0.0,1.0,0.0,0.001)
        @test spectral_bubble_cut(g,lambda).imaginary == 0
        wide = imag(spectral_bubble(g,lambda;eta_inv_fm=2e-6).value)
        narrow = imag(spectral_bubble(g,lambda;eta_inv_fm=1e-6).value)
        @test wide ≈ 2narrow rtol=1e-9
    end
    for q in (0.0,0.5,1.0), channel in (:P,:S)
        g = build_spectral_bubble(q,1.0,0.0,1.0,0.0,0.001;channel=channel)
        l = 3.0
        mt = channel === :P ? 0.0 : 4.0
        expected = Main.Constants_PNJL.N_color/(8π)*(l^2-q^2-mt)*sqrt(1-4/(l^2-q^2))
        @test spectral_bubble_cut(g,l).imaginary ≈ expected rtol=1e-12
    end
    @test_throws ArgumentError build_spectral_bubble(-1.0,1.0,0.0,1.0,0.0,0.8)
    @test_throws ArgumentError build_spectral_bubble(1.0,1.0,0.0,1.0,0.0,0.8;Phi=-0.1)
    @test_throws ArgumentError build_spectral_bubble(1.0,1.0,0.0,1.0,0.0,0.8;channel=:V)
    @test_throws ArgumentError build_spectral_bubble(1.0,1.0,0.0,1.0,0.0,0.8;thermal_cutoff_inv_fm=1.0)
end
