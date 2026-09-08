"""Analysis-only infinite thermal integral; Lambda remains a vacuum model cutoff.

Semi-infinite quadrature is a coordinate map, not a finite thermal hard edge.
No phase folding, spectral clipping, or production default changes are made.
"""
module CausalGBUInfiniteThermal
include("causal_gbu_continuation_gate.jl")
const C=CausalGBUContinuationGate
const R=C.R
const A=C.A
const P=C.E.P

function integrate_range(f,lo,hi;nodes=48,scale=1.)
    isfinite(lo) && lo<hi && (isfinite(hi) || hi==Inf) && nodes>=8 && isfinite(scale) && scale>0 ||
        throw(ArgumentError("finite lower endpoint, larger upper endpoint, positive scale and nodes>=8 required"))
    ts,ws=R.gauleg(0.,1.,nodes)
    if isfinite(hi)
        return sum(w*(hi-lo)*f(lo+(hi-lo)*t) for (t,w) in zip(ts,ws))
    end
    return sum(w*scale/(1-t)^2*f(lo+scale*t/(1-t)) for (t,w) in zip(ts,ws))
end

"""Energy map resolving thermal decay even for a large finite Landau interval."""
function thermal_range(f,lo,hi,T;nodes=48)
    upper=isfinite(hi) ? (hi-lo)/(T+hi-lo) : one(lo)
    return integrate_range(zero(lo),upper;nodes=nodes) do t
        f(lo+T*t/(1-t))*T/(1-t)^2
    end
end

"""Quadratic sublevel on a semi-infinite domain; handles linear degeneracy."""
function infinite_intervals(A,B,C0,lo)
    all(isfinite,(A,B,C0,lo)) || throw(ArgumentError("finite quadratic parameters required"))
    es=[lo,oftype(lo,Inf)]
    if A==0
        B!=0 && -C0/B>lo && push!(es,-C0/B)
    elseif B^2-4A*C0>=0
        Q=-(B+copysign(sqrt(B^2-4A*C0),B))/2
        for r in (Q==0 ? (-B/(2A),) : (Q/A,C0/Q))
            r>lo && push!(es,r)
        end
    end
    sort!(unique!(es))
    return [(es[i],es[i+1]) for i in 1:length(es)-1 if
        (x=isfinite(es[i+1]) ? (es[i]+es[i+1])/2 : es[i]+max(one(lo),abs(es[i]));
         (A*x+B)*x+C0<0)]
end

"""Pereira on-shell thermal difference on the FULL internal momentum space.

lambda is internal. q0 uses the exact radial Jacobian. Finite-q Landau
energy ranges may be semi-infinite and are never replaced by a finite Lth.
"""
function thermal_cut(bg,ch,q,lambda;nodes=48)
    i,j=R.charged_rpa_spec(ch).pair
    convert_real=lambda isa BigFloat ? BigFloat : Float64
    a,u,b,v,T,phi,bar,Q,z=convert_real.((bg.m[i],bg.mu[i],bg.m[j],bg.mu[j],bg.T,bg.Phi,bg.PhiBar,q,lambda))
    P.validate(Q,a,u,b,v,T,phi,bar,bg.vacuum,nodes)
    isfinite(z) || throw(ArgumentError("finite lambda required"))
    total=zero(z)
    for eta in (-1,1)
        S=-eta*z
        if Q>0
            interval=S>0 ? P.epsilon_interval(S/2,Q,a,b,oftype(z,Inf)) : nothing
            if interval!==nothing
                total-=eta*pi/(2Q)*integrate_range(interval...;nodes=nodes) do epsilon
                    n1,n2,_=P.weights((S-epsilon)/2,(S+epsilon)/2,u,v,T,phi,bar,:thermal)
                    eta==1 ? -n1.anti-n2.quark : -n1.quark-n2.anti
                end
            end
            lo=max(a+S/2,b-S/2);d=a^2-b^2
            for (l,h) in infinite_intervals(4*(S^2-Q^2),4S*d,Q^4+2Q^2*(a^2+b^2)-Q^2*S^2+d^2,lo)
                total+=eta*pi/Q*thermal_range(l,h,T;nodes=nodes) do en
                    n1,n2,_=P.weights(en-S/2,en+S/2,u,v,T,phi,bar,:thermal)
                    eta==1 ? n1.quark-n2.quark : n1.anti-n2.anti
                end
            end
        elseif z!=0
            if S>a+b
                e1=(S^2+a^2-b^2)/(2S);e2=S-e1
                p=sqrt((S^2-(a+b)^2)*(S^2-(a-b)^2))/(2S)
                n1,n2,_=P.weights(e1,e2,u,v,T,phi,bar,:thermal)
                total-=eta*2pi*p/abs(z)*(eta==1 ? -n1.anti-n2.quark : -n1.quark-n2.anti)
            end
            e1=(b^2-a^2-S^2)/(2S);e2=e1+S
            if e1>a && e2>b
                n1,n2,_=P.weights(e1,e2,u,v,T,phi,bar,:thermal)
                total+=eta*2pi*sqrt(e1^2-a^2)/abs(z)*(eta==1 ? n1.quark-n2.quark : n1.anti-n2.anti)
            end
        end
    end
    return 3/(8pi^2)*(z^2-Q^2-(a-b)^2)*total
end

function total_cut(bg,ch,q,lambda;nodes=48)
    i,j=R.charged_rpa_spec(ch).pair
    P.validate(q,bg.m[i],bg.mu[i],bg.m[j],bg.mu[j],bg.T,bg.Phi,bg.PhiBar,bg.vacuum,nodes)
    a0,b0,q0,z0=promote(bg.m[i],bg.m[j],q,lambda)
    # Exact geometric zero, including continuous threshold endpoints. This
    # prevents a rounded squared invariant from manufacturing gap support.
    hypot(q0,a0-b0)<=abs(z0)<=hypot(q0,a0+b0) && return zero(z0)
    # A generic one-sided vacuum cut for roundoff-scale endpoint coordinates.
    if lambda isa BigFloat
        i,j=R.charged_rpa_spec(ch).pair
        # Explicit generic q0/finite-q vacuum has no Landau contribution.
        a,b,Q,z,L=BigFloat.((bg.m[i],bg.m[j],q,lambda,bg.vacuum))
        vacuum=zero(z)
        if Q<2L && z!=0
            if Q==0
                S=abs(z)
                if S>a+b
                    e1=(S^2+a^2-b^2)/(2S);e2=S-e1
                    if a<e1<hypot(a,L) && b<e2<hypot(b,L)
                        vacuum=sign(z)*2pi*sqrt((S^2-(a+b)^2)*(S^2-(a-b)^2))/(2S^2)
                    end
                end
            else
                interval=P.epsilon_interval(abs(z)/2,Q,a,b,L)
                interval===nothing || (vacuum=sign(z)*pi/(2Q)*(interval[2]-interval[1]))
            end
        end
        return 3/(8pi^2)*(z^2-Q^2-(a-b)^2)*vacuum+thermal_cut(bg,ch,q,lambda;nodes=nodes)
    end
    i,j=R.charged_rpa_spec(ch).pair
    v=P.reference_cut(lambda,q,bg.m[i],bg.mu[i],bg.m[j],bg.mu[j],bg.T,bg.vacuum;
        Phi=bg.Phi,PhiBar=bg.PhiBar,component=:vacuum,nodes=nodes)
    return 3/(8pi^2)*(lambda^2-q^2-(bg.m[i]-bg.m[j])^2)*v.imaginary+
        thermal_cut(bg,ch,q,lambda;nodes=nodes)
end

"""Independent radial Dirac-projector integral on the infinite thermal domain.

Only open-UHP frequencies: discrete quadrature poles are not physical modes.
The vacuum is independently integrated in cylindrical lens coordinates.
"""
function radial_polarization(bg,ch,q,z;nodes=96)
    imag(z)>0 || throw(ArgumentError("open UHP required for radial quadrature"))
    i,j=R.charged_rpa_spec(ch).pair
    a,u,b,v,T=bg.m[i],bg.mu[i],bg.m[j],bg.mu[j],bg.T
    vacuum=P.cylindrical_reference(q,a,u,b,v,T,bg.vacuum,[z];Phi=bg.Phi,PhiBar=bg.PhiBar,
        component=:vacuum,nz=nodes,ny=nodes).pi_p[1]
    xs,ws=R.gauleg(-1.,1.,nodes)
    thermal=integrate_range(0.,Inf;nodes=nodes,scale=T) do p
        E1=hypot(p,a)
        n1=P.occupations(E1,u,T,bg.Phi,bg.PhiBar)
        term=zero(complex(z))
        for (x,w) in zip(xs,ws)
            E2=sqrt(p^2+q^2-2p*q*x+b^2)
            n2=P.occupations(E2,v,T,bg.Phi,bg.PhiBar)
            for s in (-1,1),t in (-1,1)
                energy=s*E1-t*E2
                # Explicit multiplication: Julia parses 2E1 as the literal 20.
                trace=-s*t*(energy^2-q^2-(a-b)^2)/(2*E1*E2)
                one=s==1 ? n1.quark : -n1.anti
                two=t==1 ? n2.quark : -n2.anti
                term+=w*trace*(two-one)/(energy-z)
            end
        end
        3p^2/(4pi^2)*term
    end
    return vacuum+thermal
end

function kernel(bg,ch,q;cut_nodes=48,split_inv_fm=0.)
    i,j=R.charged_rpa_spec(ch).pair
    g=R.build_spectral_bubble(q,bg.m[i],bg.mu[i],bg.m[j],bg.mu[j],bg.T;
        Phi=bg.Phi,PhiBar=bg.PhiBar,vacuum_cutoff_inv_fm=bg.vacuum,
        thermal_cutoff_inv_fm=bg.vacuum,momentum_nodes=8,angle_nodes=8)
    threshold=hypot(q,bg.m[i]+bg.m[j]);landau=hypot(q,bg.m[i]-bg.m[j])
    split=max(Float64(split_inv_fm),last(A.cut_panels(g))+4bg.T,threshold+4bg.T)
    # isequal distinguishes signed zeros; a zero-width panel is not a cut.
    edges=sort!(unique!(map(x->x==0 ? 0.0 : x,
        vcat(A.cut_panels(g),[-split,split,-threshold,threshold,-landau,landau,-q,q]))))
    return (bg=bg,ch=ch,q=q,rho=x->total_cut(bg,ch,q,x;nodes=cut_nodes),edges=edges,
        split=split,threshold=threshold,landau=landau,shift=bg.mu[i]-bg.mu[j])
end

"""Original-cut PV/UHP transform with the remaining infinite tails retained."""
function polarization(k,z;nodes=64)
    isfinite(z) && imag(z)>=0 && abs(real(z))<k.split ||
        throw(ArgumentError("closed UHP frequency must lie inside numerical split window"))
    finite=C.mapped_cauchy(k.rho,k.edges,z;nodes=nodes)
    tail=integrate_range(k.split,Inf;nodes=nodes,scale=2k.bg.T) do x
        (k.rho(x)/(x-z)+k.rho(-x)/(-x-z))/pi
    end
    return finite+tail
end

"""Conditional weighted phase-band bound; q measure is included.

Requires caller control of the actual phase branch, not merely small rho.
The derivative representation must include the corresponding boundary terms.
"""
function phase_band_bound(lower,upper,T,q;weight_bound=pi,branch_bound_verified=false)
    all(isfinite,(lower,upper,T,q,weight_bound)) && 0<lower<upper && T>0 && min(q,weight_bound)>=0 ||
        throw(ArgumentError("positive Bose-safe band and nonnegative q/weight bound required"))
    g(w)=inv(expm1(w/T))
    return (shell_absolute_bound_inv_fm2=q^2/(2pi^3)*weight_bound*(g(lower)-g(upper)),
        bound_verified=branch_bound_verified,production_authorized=false)
end
end
