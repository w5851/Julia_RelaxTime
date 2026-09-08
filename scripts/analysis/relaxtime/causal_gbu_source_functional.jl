"""Fixed-projector source functional: analysis-only independent Hessian oracle.

The signed vacuum + medium - reference determinant is NOT a positive thermal
trace, a selected physical response, or a derivation of a GBU particle yield.
"""
module CausalGBUSourceFunctional
include("causal_gbu_direction_b.jl")
const B=CausalGBUDirectionB
const R=B.R
const P=B.P
using LinearAlgebra,ForwardDiff

const ID2=Matrix{ComplexF64}(I,2,2)
const ID4=Matrix{ComplexF64}(I,4,4)
const ZERO2=zeros(ComplexF64,2,2)
const SIGMA=(ComplexF64[0 1;1 0],ComplexF64[0 -im;im 0],ComplexF64[1 0;0 -1])
const ALPHA=ntuple(i->[ZERO2 SIGMA[i];SIGMA[i] ZERO2],3)
const BETA=[ID2 ZERO2;ZERO2 -ID2]
const GAMMA5=[ZERO2 ID2;ID2 ZERO2]
const VP=im*BETA*GAMMA5

function validate_matrix(d,V,keep)
    n=length(d)
    n>0 && size(V)==(n,n) && length(keep)==n && all(isfinite,d) &&
        all(isfinite,V) && all(>(0),d) && ishermitian(V) &&
        eltype(keep)<:Bool || throw(ArgumentError("invalid projected determinant inputs"))
end

"""-log det(P(D+jV)P on ran(P)), for a positive diagonal synthetic D."""
function projected_action(d,V,keep,j)
    validate_matrix(d,V,keep)
    inds=findall(keep)
    isempty(inds) && return zero(j)
    return -log(det(Matrix(Diagonal(d[inds]))+j*V[inds,inds]))
end

"""Exact second derivative: both internal indices are projected."""
function projected_hessian(d,V,keep)
    validate_matrix(d,V,keep)
    return sum(abs2(V[i,j])/(d[i]*d[j]) for i in eachindex(d),j in eachindex(d)
        if keep[i] && keep[j];init=0.)
end

"""Complexified color eigenvalues; this is not a positive SU(3) heat trace."""
function color_eigenvalues(phi,bar)
    0<=phi<=1 && 0<=bar<=1 || throw(ArgumentError("invalid Polyakov variables"))
    phi==bar==1 && return ones(ComplexF64,3)
    return eigvals(ComplexF64[0 0 1;1 0 -3bar;0 1 3phi])
end

hamiltonian(p,m)=sum(p[i]*ALPHA[i] for i in 1:3)+m*BETA
vertex(channel)=channel===:P ? VP : channel===:S ? BETA :
    throw(ArgumentError("channel must be :P or :S"))

"""Dirac projector trace, independently of the A/B reduction."""
function spin_trace(p,r,a,b,s,t,channel=:P)
    s in (-1,1) && t in (-1,1) || throw(ArgumentError("energy signs must be +/-1"))
    e1=sqrt(sum(abs2,p)+a^2); e2=sqrt(sum(abs2,r)+b^2)
    l1=(ID4+s*hamiltonian(p,a)/e1)/2
    l2=(ID4+t*hamiltonian(r,b)/e2)/2
    V=vertex(channel)
    return real(tr(l1*V*l2*V))
end

function occupation_component(xi,T,phi,bar,component)
    component in (:vacuum,:thermal,:full) || throw(ArgumentError("unknown component"))
    vac=xi<0 ? 1. : 0.
    component===:vacuum && return vac
    full=P.occupation(xi/T,phi,bar)
    return component===:thermal ? full-vac : full
end

"""(n(y)-n(x))/(x-y-z), with the STATIC coincident-state limit explicit."""
function response_quotient(x,y,z,T,phi,bar,component)
    nx=occupation_component(x,T,phi,bar,component)
    ny=occupation_component(y,T,phi,bar,component)
    if iszero(z) && abs(x-y)<1e-10
        # This oracle is restricted below zero-T flavor onset.
        x*y>0 || throw(ArgumentError("zero-T Fermi-surface degeneracy"))
        return component===:vacuum ? 0. :
            -ForwardDiff.derivative(v->P.occupation(v/T,phi,bar),(x+y)/2)
    end
    return (ny-nx)/(x-y-z)
end

"""Thermal/vacuum trace-log for an 8x8 static charged source block.

J multiplies (E12+E21)/sqrt(2); tr_flavor(T^2)=1, hence curvature=-Pi.
Only differences relative to J=0 are used; source-independent constants cancel.
"""
function block_action(p,r,a,u,b,v,T,phi,bar,j,component;channel=:P)
    component in (:vacuum,:thermal,:full) || throw(ArgumentError("unknown component"))
    h1=hamiltonian(p,a)-u*ID4; h2=hamiltonian(r,b)-v*ID4
    V=vertex(channel)/sqrt(2)
    es=eigvals(Hermitian([h1 j*V;j*V h2]))
    vac=3sum(e for e in es if e<0;init=0.)
    component===:vacuum && return vac
    full=-T*sum(B.log_partition(e/T,phi,bar) for e in es)
    return component===:thermal ? full-vac : full
end

"""Independent cylindrical-lens integration; q and all energies in fm^-1.

Pi comes from Dirac projectors, not the B0 polynomial or spectral Cauchy code.
Optional source steps compare an actual trace-log second difference at k0=0.
Finite node atoms and finite differences do not certify cuts or UHP stability.
"""
function source_response(q,a,u,b,v,T,vc,tc,zs;Phi=1.,PhiBar=1.,nodes=32,
                         channel=:P,source_steps=Float64[])
    P.validate(q,a,u,b,v,T,Phi,PhiBar,vc,nodes)
    isfinite(tc) && tc>=vc && abs(u)<a && abs(v)<b ||
        throw(ArgumentError("requires thermal>=vacuum and below zero-T flavor onset"))
    !isempty(zs) && all(z->isfinite(z) && imag(z)>=0 && (imag(z)>0 || iszero(z)),zs) ||
        throw(ArgumentError("use external UHP frequencies or exactly static k0=0"))
    all(h->isfinite(h) && h>0,source_steps) || throw(ArgumentError("invalid source steps"))
    vertex(channel)
    values=zeros(ComplexF64,length(zs))
    curvatures=zeros(length(source_steps))
    for (component,L) in ((:vacuum,vc),(:thermal,tc))
        q<2L || continue
        h=L-q/2
        for (lo,hi) in ((-h,0.),(0.,h))
            ts,tws=R.gauleg(lo,hi,nodes)
            for (t0,tw) in zip(ts,tws)
                ymax=(L-abs(t0)-q/2)*(L+abs(t0)+q/2)
                ys,yws=R.gauleg(0.,ymax,nodes)
                for (y,yw) in zip(ys,yws)
                    p=(sqrt(y),0.,t0+q/2); r=(sqrt(y),0.,t0-q/2)
                    e1=sqrt(sum(abs2,p)+a^2); e2=sqrt(sum(abs2,r)+b^2)
                    w=tw*yw/(8pi^2)
                    for s in (-1,1),t in (-1,1)
                        trace=spin_trace(p,r,a,b,s,t,channel)
                        for k in eachindex(zs)
                            values[k]+=3w*trace*response_quotient(s*e1-u,t*e2-v,zs[k],T,Phi,PhiBar,component)
                        end
                    end
                    if !isempty(source_steps)
                        f(j)=block_action(p,r,a,u,b,v,T,Phi,PhiBar,j,component;channel=channel)
                        f0=f(0.)
                        for (k,j) in enumerate(source_steps)
                            curvatures[k]+=w*(f(j)+f(-j)-2*f0)/j^2
                        end
                    end
                end
            end
        end
    end
    return (pi_inv_fm2=values,source_curvatures_inv_fm2=curvatures,
        source_normalization="real_generator_trace_square_one",
        functional="fixed_projector_vacuum_plus_medium_minus_reference",
        zero_T_reference_same_mu=true,full_stationarity_certified=false,
        positive_physical_trace_certified=false,production_authorized=false)
end
end
