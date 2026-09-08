"""Independent Pereira et al. (PRC 109, 025206) reference, analysis only.

Eq. (26), (55), (56), (60), and Appendix C18/C19, using RETARDED k0+i0.
The cylindrical lens quadrature does not reuse the project's angular domain,
four-residue sum, cut intervals, or Cauchy reconstruction. All energies use
fm^-1. Lambda is explicit; PNJL occupations are a separately identified
extension, with Phi=PhiBar=1 reducing to the paper's Fermi distributions.
"""
module CausalGBUPereiraReference
const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
isdefined(Main, :GaussLegendre) ||
    Base.include(Main, joinpath(ROOT, "src", "integration", "GaussLegendre.jl"))
using Main.GaussLegendre: gauleg

function validate(q,a,u,b,v,T,phi,phibar,L,n)
    all(isfinite,(q,a,u,b,v,T,phi,phibar,L)) &&
        q>=0 && min(a,b,T,L)>0 && 0<=phi<=1 && 0<=phibar<=1 && n>=4 ||
        throw(ArgumentError("invalid Pereira reference inputs"))
end

"""Independent stable occupation polynomial; no project distribution calls."""
function occupation(x,phi,phibar)
    if x>=0
        y=exp(-x)
        return (phi*y+2phibar*y^2+y^3)/(1+3phi*y+3phibar*y^2+y^3)
    end
    y=exp(x)
    return (phi*y^2+2phibar*y+1)/(y^3+3phi*y^2+3phibar*y+1)
end
occupations(E,mu,T,phi,phibar) =
    (quark=occupation((E-mu)/T,phi,phibar),
     anti=occupation((E+mu)/T,phibar,phi))

function weights(E1,E2,u,v,T,phi,phibar,component)
    component in (:full,:vacuum,:thermal) || throw(ArgumentError("unknown component"))
    n1=component===:vacuum ? (quark=0.,anti=0.) : occupations(E1,u,T,phi,phibar)
    n2=component===:vacuum ? (quark=0.,anti=0.) : occupations(E2,v,T,phi,phibar)
    c=component===:thermal ? 0. : 1.
    return n1,n2,c
end

"""Paper pair/scattering denominators at INTERNAL lambda z in the UHP."""
function pair_scattering(E1,E2,z,n1,n2,c=1.)
    pair=(c-n1.anti-n2.quark)/(E1+E2+z) +
         (c-n1.quark-n2.anti)/(E1+E2-z)
    scattering=-(n1.quark-n2.quark)/(E2-E1+z) -
                (n1.anti-n2.anti)/(E2-E1-z)
    return pair,scattering
end

"""Lens coordinates t=p_z-q/2, y=p_perp^2, p^2 dp dx=dy dt/2.

The two halves t<0,t>0 are integrated separately. A1 and A2 are direct
single-line integrals over this SAME lens, not moments of B0 residues.
"""
function cylindrical_reference(q,a,u,b,v,T,L,zs;Phi=1.,PhiBar=1.,
                               component=:full,nz=64,ny=64)
    validate(q,a,u,b,v,T,Phi,PhiBar,L,min(nz,ny))
    component in (:full,:vacuum,:thermal) || throw(ArgumentError("unknown component"))
    !isempty(zs) && all(z->isfinite(z) && imag(z)>0,zs) ||
        throw(ArgumentError("reference requires upper-half-plane internal lambda"))
    pair=zeros(ComplexF64,length(zs)); scatter=zero(pair)
    A1,A2,volume=0.,0.,0.
    if q<2L
        h=L-q/2
        for (lo,hi) in ((-h,0.),(0.,h))
            ts,tws=gauleg(lo,hi,nz)
            for (t,tw) in zip(ts,tws)
                ymax=(L-abs(t)-q/2)*(L+abs(t)+q/2)
                ys,yws=gauleg(0.,ymax,ny)
                for (y,yw) in zip(ys,yws)
                    E1=sqrt(y+(t+q/2)^2+a^2)
                    E2=sqrt(y+(t-q/2)^2+b^2)
                    n1,n2,c=weights(E1,E2,u,v,T,Phi,PhiBar,component)
                    w=tw*yw/2
                    volume+=2pi*w
                    A1-=2w*(c-n1.quark-n1.anti)/E1
                    A2-=2w*(c-n2.quark-n2.anti)/E2
                    for i in eachindex(zs)
                        p,s=pair_scattering(E1,E2,zs[i],n1,n2,c)
                        pair[i]+=w*p/(E1*E2)
                        scatter[i]+=w*s/(E1*E2)
                    end
                end
            end
        end
    end
    b0=pair.+scatter
    pi_p=[3/(8pi^2)*((z^2-q^2-(a-b)^2)*b0[i]-A1-A2) for (i,z) in enumerate(zs)]
    pi_s=[3/(8pi^2)*((z^2-q^2-(a+b)^2)*b0[i]-A1-A2) for (i,z) in enumerate(zs)]
    return (;pair,scatter,b0,A1,A2,contact=A1+A2,pi_p,pi_s,volume,
             production_authorized=false)
end

"""Closed radial antiderivatives of paper Eq. (60) in vacuum.

Equivalent to integrating the angular width of the lens; no quadrature,
other flavor mass, or spectral reconstruction enters this oracle.
"""
function vacuum_A(q,m,L)
    all(isfinite,(q,m,L)) && q>=0 && min(m,L)>0 || throw(ArgumentError("invalid vacuum A"))
    f2(p)=(p*hypot(p,m)-m^2*asinh(p/m))/2
    q==0 && return -4*f2(L)
    q>=2L && return 0.
    lo=abs(L-q)
    e0,e1=hypot(lo,m),hypot(L,m)
    delta1=e1-e0
    delta3=(e1^3-e0^3)/3-m^2*delta1
    partial=-2*(f2(L)-f2(lo))+(delta3+(q^2-L^2)*delta1)/q
    return partial+(q<L ? -4*f2(L-q) : 0.)
end

"""Allowed epsilon=E2-E1 at fixed E=(E1+E2)/2: Appendix D2-D10."""
function epsilon_interval(E,q,a,b,L)
    s=4E^2-q^2
    s>(a+b)^2 || return nothing
    spread=q*sqrt((s-(a+b)^2)*(s-(a-b)^2))/s
    center=-2E*(a^2-b^2)/s
    lo=max(center-spread,2*(E-hypot(a,L)),2*(b-E))
    hi=min(center+spread,2*(E-a),2*(hypot(b,L)-E))
    return lo<hi ? (lo,hi) : nothing
end

"""Exact quadratic sublevel intervals, clipped to a finite physical interval."""
function quadratic_intervals(A,B,C,lo,hi)
    lo<hi || return Tuple{Float64,Float64}[]
    edges=Float64[lo,hi]
    if A==0
        B!=0 && lo< -C/B <hi && push!(edges,-C/B)
    else
        d=B^2-4A*C
        if d>=0
            Q=-(B+copysign(sqrt(d),B))/2
            roots=Q==0 ? (-B/(2A),) : (Q/A,C/Q)
            for r in roots
                lo<r<hi && push!(edges,r)
            end
        end
    end
    sort!(unique!(edges))
    return [(edges[i],edges[i+1]) for i in 1:length(edges)-1
        if (x=(edges[i]+edges[i+1])/2; (A*x+B)*x+C<0)]
end

function scattering_intervals(epsilon,q,a,b,L)
    lo=max(a+epsilon/2,b-epsilon/2)
    hi=min(hypot(a,L)+epsilon/2,hypot(b,L)-epsilon/2)
    d=a^2-b^2
    return quadratic_intervals(4*(epsilon^2-q^2),4epsilon*d,
        q^4+2q^2*(a^2+b^2)-q^2*epsilon^2+d^2,lo,hi)
end

"""Real-axis B0 cut from Appendix C18/C19 (external k0+i0), not mass-i0.

lambda is INTERNAL; eta in the paper is a +/-1 summation index, not a
finite broadening. q=0 has its own Jacobian. Static degeneracy is flagged.
"""
function reference_cut(lambda,q,a,u,b,v,T,L;Phi=1.,PhiBar=1.,
                       component=:full,nodes=64)
    validate(q,a,u,b,v,T,Phi,PhiBar,L,nodes)
    isfinite(lambda) || throw(ArgumentError("nonfinite lambda"))
    component in (:full,:vacuum,:thermal) || throw(ArgumentError("unknown component"))
    pair,landau=0.,0.
    if q<2L
        for eta in (-1,1)
            if q>0
                E=-eta*lambda/2
                interval=E>0 ? epsilon_interval(E,q,a,b,L) : nothing
                if interval!==nothing
                    es,ws=gauleg(interval...,nodes)
                    for (epsilon,w) in zip(es,ws)
                        E1,E2=E-epsilon/2,E+epsilon/2
                        n1,n2,c=weights(E1,E2,u,v,T,Phi,PhiBar,component)
                        W=eta==1 ? c-n1.anti-n2.quark : c-n1.quark-n2.anti
                        pair-=eta*pi/(2q)*w*W
                    end
                end
                epsilon=-eta*lambda
                for (lo,hi) in scattering_intervals(epsilon,q,a,b,L)
                    es,ws=gauleg(lo,hi,nodes)
                    for (E,w) in zip(es,ws)
                        n1,n2,_=weights(E-epsilon/2,E+epsilon/2,u,v,T,Phi,PhiBar,component)
                        W=eta==1 ? n1.quark-n2.quark : n1.anti-n2.anti
                        landau+=eta*pi/q*w*W
                    end
                end
            elseif lambda!=0
                S=-eta*lambda
                if S>a+b
                    E1=(S^2+a^2-b^2)/(2S); E2=S-E1
                    if a<E1<hypot(a,L) && b<E2<hypot(b,L)
                        p=sqrt((S^2-(a+b)^2)*(S^2-(a-b)^2))/(2S)
                        n1,n2,c=weights(E1,E2,u,v,T,Phi,PhiBar,component)
                        W=eta==1 ? c-n1.anti-n2.quark : c-n1.quark-n2.anti
                        pair-=eta*2pi*p/abs(lambda)*W
                    end
                end
                epsilon=-eta*lambda
                E1=(b^2-a^2-epsilon^2)/(2epsilon); E2=E1+epsilon
                if a<E1<hypot(a,L) && b<E2<hypot(b,L)
                    p=sqrt(E1^2-a^2)
                    n1,n2,_=weights(E1,E2,u,v,T,Phi,PhiBar,component)
                    W=eta==1 ? n1.quark-n2.quark : n1.anti-n2.anti
                    landau+=eta*2pi*p/abs(lambda)*W
                end
            end
        end
    end
    return (;pair,landau,imaginary=pair+landau,
             static_degeneracy_unresolved=q==0 && lambda==0 && a==b,
             production_authorized=false)
end
end
