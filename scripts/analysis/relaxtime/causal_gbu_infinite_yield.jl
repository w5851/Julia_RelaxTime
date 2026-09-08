"""Restricted, fail-closed GBU partial-yield candidate on the infinite thermal kernel.

This is not the production provider. Normal-gap counts are independent of
phase extraction. Cases with additional cut winding require further analysis.
"""
module CausalGBUInfiniteYield
include("causal_gbu_infinite_profile.jl")
include("causal_gbu_infinite_limits.jl")
const P=CausalGBUInfiniteProfile
const Limits=CausalGBUInfiniteLimits
const I=P.I
const R=P.R
const O=I.C.E.O

function bracket_roots(f,a,b;nodes=128)
    a<b && nodes>=16 || throw(ArgumentError("ordered root interval and nodes>=16 required"))
    xs=[a+(b-a)*sinpi(i/(2nodes))^2 for i in 0:nodes]
    vals=real.(f.(xs)); roots=Float64[]
    all(isfinite,vals) || throw(ArgumentError("nonfinite root scan"))
    for j in 1:nodes
        l,h=xs[j],xs[j+1]; fl,fh=vals[j],vals[j+1]
        fl==0 && push!(roots,l)
        fl*fh<0 || continue
        for _ in 1:60
            mid=(l+h)/2; fm=real(f(mid))
            if fl*fm<=0; h=mid; else; l=mid;fl=fm;end
        end
        push!(roots,(l+h)/2)
    end
    last(vals)==0 && push!(roots,b)
    return unique(roots)
end

function gap_roots(p;nodes=128)
    k=p.kernel; d,U=k.landau,k.threshold
    gaps=((-U+1e-8,-d-1e-8),(d+1e-8,U-1e-8))
    return [bracket_roots(z->P.inverse(p,z),a,b;nodes=nodes) for (a,b) in gaps]
end

"""Principal cut phase ONLY; a negative-real zero-cut gap is counted separately."""
function cut_phase(p,w;eta=0.)
    f=P.inverse(p,complex(w+p.kernel.shift,eta))
    abs(f)>1e-12 || throw(ArgumentError("unresolved zero in cut integration"))
    return -angle(f)
end

function cut_integral(p,lo,hi,T;nodes=64,eta=0.)
    lo<hi || return 0.
    shift=p.kernel.shift
    edges=sort!(unique!(vcat([lo,hi],filter(x->lo<x<hi,p.kernel.edges.-shift),
        filter(x->lo<x<hi,[-p.threshold.U-p.threshold.H-shift,p.threshold.U+p.threshold.H-shift]))))
    # No Bose clipping: all ordinates are strictly positive and no endpoint is sampled.
    return real(I.C.mapped_integral(edges;nodes=nodes) do w
        g=inv(expm1(w/T))
        g*(1+g)/T*O.gbu_weight(cut_phase(p,w;eta=eta))/pi
    end)
end

"""Root term plus signed continuum with the unitary lower-boundary subtraction."""
function shell(p;nodes=64,upper=nothing,eta=0.)
    k=p.kernel; bg=k.bg; T=bg.T; shift=k.shift
    upper=upper===nothing ? max(24.,k.q+24T) : Float64(upper)
    rootsets=gap_roots(p)
    negative,positive=rootsets
    lower=max(0.,k.landau-shift); threshold=k.threshold-shift
    abs(shift)<k.threshold && threshold<upper<k.split-abs(shift) ||
        throw(ArgumentError("Bose-safe normal-gap geometry required"))
    f0=P.inverse(p,shift)
    real(f0)>0 && abs(imag(f0))<1e-8 || throw(ArgumentError("static instability or unresolved Bose endpoint"))
    roots=[x-shift for x in positive]
    all(>(0),roots) || throw(ArgumentError("nonpositive-energy root"))
    length(negative)<=1 && length(positive)<=1 || throw(ArgumentError("candidate requires reviewed extra-root topology"))
    n=length(roots)
    # Independently compare both gap boundaries with the cut phase limits.
    landau_end=lower>1e-8 ? cut_phase(p,lower-1e-8) : 0.
    # A fixed offset is not the threshold limit near Mott. The sign is
    # certified independently of n, using the original Pi error budget.
    pair_start=Limits.threshold_limit(P.inverse(p,k.threshold),p.threshold.positive;
        inverse_error_budget=4bg.coupling[k.ch]*1e-6).phase
    abs(landau_end)<0.005 && abs(pair_start-n*pi)<0.005 ||
        throw(ArgumentError("independent root count and cut phase limits disagree"))
    highphase=cut_phase(p,upper)
    abs(highphase)<0.005 || throw(ArgumentError("UV phase not near zero"))
    landau=cut_integral(p,0.,lower,T;nodes=nodes,eta=eta)
    pair=cut_integral(p,threshold,upper,T;nodes=nodes,eta=eta)
    bound=sum(x->inv(expm1(x/T)),roots;init=0.)
    boundary=n*inv(expm1(threshold/T))
    measure=k.q^2/(2pi^2)
    return (density=measure*(bound+landau+pair-boundary),
        bound=measure*bound,landau=measure*landau,pair=measure*(pair-boundary),
        root_count=n,negative_root_count=length(negative),roots=roots,
        landau_phase=landau_end,threshold_phase=pair_start,high_phase=highphase,
        omega_tail_conditional_bound=measure*inv(expm1(upper/T)),
        static_inverse=real(f0),production_authorized=false)
end

"""Finite-L numerical comparison, retaining its cut and endpoints (not a new model)."""
function finite_kernel(bg,ch,q,L;cut_nodes=64,split_inv_fm=36.)
    k=I.kernel(bg,ch,q;cut_nodes=cut_nodes,split_inv_fm=max(split_inv_fm,q+28bg.T+2.))
    i,j=R.charged_rpa_spec(ch).pair
    g=R.build_spectral_bubble(q,bg.m[i],bg.mu[i],bg.m[j],bg.mu[j],bg.T;
        Phi=bg.Phi,PhiBar=bg.PhiBar,vacuum_cutoff_inv_fm=bg.vacuum,
        thermal_cutoff_inv_fm=L,momentum_nodes=8,angle_nodes=8)
    edges=sort!(unique!(map(x->x==0 ? 0. : x,vcat(k.edges,I.A.cut_panels(g)))))
    split=max(k.split,last(edges)+1.)
    append!(edges,[-split,split]);sort!(unique!(edges))
    rho=x->k.landau<=abs(x)<=k.threshold ? zero(x) : x isa BigFloat ?
        I.C.wide_cut(bg,ch,q,x,L;nodes=cut_nodes) : I.A.direct_cut(bg,ch,q,x,L;nodes=cut_nodes).imaginary
    return merge(k,(rho=rho,edges=edges,split=split,thermal_limit=L))
end
end
