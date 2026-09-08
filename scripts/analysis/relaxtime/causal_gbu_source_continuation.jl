"""Analysis-only extension of the fixed Dirac-sea reference through flavor onset."""
module CausalGBUSourceContinuation
include("causal_gbu_source_functional.jl")
const S=CausalGBUSourceFunctional
const R=S.R
const P=S.P
using LinearAlgebra,ForwardDiff

"""Sufficient band gap, not a condition |mu_i|<M_i on the medium."""
function reference_gap(q,a,u,b,v)
    all(isfinite,(q,a,u,b,v)) && q>=0 && min(a,b)>0 || throw(ArgumentError("invalid reference parameters"))
    return min(2a,2b,hypot(q,a+b)-abs(u-v))
end

function component_occupation(e,mu,s,T,phi,bar,component)
    s in (-1,1) && component in (:vacuum,:thermal,:full) || throw(ArgumentError("invalid band/component"))
    sea=s==-1 ? 1. : 0.
    component===:vacuum && return sea
    full=P.occupation((s*e-mu)/T,phi,bar)
    return component===:thermal ? full-sea : full
end

function quotient(e1,u,s,e2,v,t,z,T,phi,bar,component)
    x,y=s*e1-u,t*e2-v
    nx=component_occupation(e1,u,s,T,phi,bar,component)
    ny=component_occupation(e2,v,t,T,phi,bar,component)
    if iszero(z) && abs(x-y)<1e-10
        s==t || throw(ArgumentError("Dirac sea band gap closed"))
        return component===:vacuum ? 0. : -ForwardDiff.derivative(w->P.occupation(w/T,phi,bar),(x+y)/2)
    end
    return (ny-nx)/(x-y-z)
end

"""Adiabatic rank-four sea: retain the same bands, not all grand-energies <0.

Weyl's bound keeps these bands separated throughout [-abs(j),abs(j)].
This is a local source Hessian oracle, not a full nonuniform finite-J action.
"""
function block_action(p,r,a,u,b,v,T,phi,bar,j,component;channel=:P)
    component in (:vacuum,:thermal,:full) || throw(ArgumentError("unknown component"))
    e1=sqrt(sum(abs2,p)+a^2);e2=sqrt(sum(abs2,r)+b^2)
    gap=min(2e1,2e2,e1+e2-abs(u-v))
    gap>sqrt(2)*abs(j) || throw(ArgumentError("source step can close Dirac reference gap"))
    h1=S.hamiltonian(p,a)-u*S.ID4;h2=S.hamiltonian(r,b)-v*S.ID4
    V=S.vertex(channel)/sqrt(2)
    es=eigvals(Hermitian([h1 j*V;j*V h2]))
    sea=3sum(es[1:4])
    component===:vacuum && return sea
    full=-T*sum(S.B.log_partition(e/T,phi,bar) for e in es)
    return component===:thermal ? full-sea : full
end

function source_response(q,a,u,b,v,T,vc,tc,zs;Phi=1.,PhiBar=1.,nodes=32,
                         channel=:P,source_steps=Float64[])
    P.validate(q,a,u,b,v,T,Phi,PhiBar,vc,nodes)
    gap=reference_gap(q,a,u,b,v)
    isfinite(tc) && tc>=vc && gap>0 || throw(ArgumentError("thermal>=vacuum and open reference band gap required"))
    !isempty(zs) && all(z->isfinite(z) && imag(z)>=0 && (imag(z)>0 || iszero(z)),zs) ||
        throw(ArgumentError("external UHP or static frequencies required"))
    all(h->isfinite(h) && 0<h<gap/sqrt(2),source_steps) || throw(ArgumentError("source steps must preserve the band gap"))
    S.vertex(channel)
    values=zeros(ComplexF64,length(zs));curvatures=zeros(length(source_steps))
    for (component,L) in ((:vacuum,vc),(:thermal,tc))
        q<2L || continue
        h=L-q/2
        for (lo,hi) in ((-h,0.),(0.,h))
            ts,tws=R.gauleg(lo,hi,nodes)
            for (t0,tw) in zip(ts,tws)
                ymax=(L-abs(t0)-q/2)*(L+abs(t0)+q/2)
                ys,yws=R.gauleg(0.,ymax,nodes)
                for (y,yw) in zip(ys,yws)
                    p=(sqrt(y),0.,t0+q/2);r=(sqrt(y),0.,t0-q/2)
                    e1=sqrt(sum(abs2,p)+a^2);e2=sqrt(sum(abs2,r)+b^2)
                    w=tw*yw/(8pi^2)
                    for s in (-1,1),t in (-1,1)
                        tr=S.spin_trace(p,r,a,b,s,t,channel)
                        for k in eachindex(zs)
                            values[k]+=3w*tr*quotient(e1,u,s,e2,v,t,zs[k],T,Phi,PhiBar,component)
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
        reference_gap_lower_bound_inv_fm=gap,reference="adiabatic_rank_four_Dirac_sea",
        below_individual_onset=abs(u)<a && abs(v)<b,
        full_stationarity_certified=false,production_authorized=false)
end
end
