"""Analysis-only regulator oracles; no production provider/default is replaced."""
module CausalGBURegulatorChecks
include("causal_gbu_research_utils.jl")
const R = CausalGBUResearch
using Main.RelaxTime.OneLoopIntegrals: B0_spectral_cut, _spectral_cut_intervals
using Main.PNJLQuarkDistributions: quark_distribution, antiquark_distribution
const NC = Main.Constants_PNJL.N_color

"""Exact center-ball restriction on the angular-delta energy interval.

P^2=(p1^2+p2^2)/2-q^2/4, E2=t*(s*E1-lambda).
Thus P^2<L^2 is a monic quadratic in E1; no cut is removed by a phase rule.
"""
function centered_intervals(l,q,a,b,L,s,t)
    all(isfinite,(l,q,a,b,L)) && q > 0 && L > 0 && a > 0 && b > 0 &&
        s in (-1,1) && t in (-1,1) || throw(ArgumentError("finite positive q,masses,cutoff and signs +/-1 required"))
    disc = -l^2+2(a^2+b^2)+q^2+4L^2
    disc > 0 || return Tuple{Float64,Float64}[]
    lo,hi = (s*l-sqrt(disc))/2,(s*l+sqrt(disc))/2
    intervals = Tuple{Float64,Float64}[]
    # Every center-ball point has both line momenta <= L+q/2.
    for (u,v) in _spectral_cut_intervals(l,q,a,b,L+q/2,s,t)
        left,right = max(u,lo),min(v,hi)
        left < right && push!(intervals,(left,right))
    end
    return intervals
end

function centered_cut(l,q,a,mu1,b,mu2,T,phi,phibar,L;ne=64,component=:thermal)
    all(isfinite,(l,q,a,mu1,b,mu2,T,phi,phibar,L)) && q>=0 && min(a,b,T,L)>0 &&
        0<=phi<=1 && 0<=phibar<=1 && ne>=4 || throw(ArgumentError("invalid centered cut inputs"))
    component in (:vacuum,:thermal) || throw(ArgumentError("unknown component"))
    q == 0 && return B0_spectral_cut(l,q,a,mu1,b,mu2,T;Φ=phi,Φbar=phibar,
        pmax_inv_fm=L,energy_nodes=ne,component=component).imaginary
    occupation(s,E,mu) = component===:vacuum ? (s<0 ? 1.0 : 0.0) :
        (s>0 ? quark_distribution(E,mu,T,phi,phibar) : -antiquark_distribution(E,mu,T,phi,phibar))
    val = 0.0
    support=domain_support(q,a,b,L,:centered)
    for s in (-1,1),t in (-1,1)
        lower,upper=s==t ? support.landau : support.pair
        lower<s*l<upper || continue # exact domain geometry, never a magnitude clip
        s==t && abs(l)>=hypot(q,a-b) && continue
        s!=t && abs(l)<=hypot(q,a+b) && continue
        for (left,right) in centered_intervals(l,q,a,b,L,s,t)
            es,ws=R.gauleg(left,right,ne)
            for (e,w) in zip(es,ws)
                e2=t*(s*e-l)
                val -= pi/q*w*s*t*(occupation(t,e2,mu2)-occupation(s,e,mu1))
            end
        end
    end
    return val
end

"""Independent center-coordinate loop quadrature, including BOTH polynomial moments."""
function centered_atoms(q,a,mu1,b,mu2,T,phi,phibar,vc,tc;np=128,nx=64)
    all(isfinite,(q,a,mu1,b,mu2,T,phi,phibar,vc,tc)) && q>=0 && 0<vc<=tc && min(a,b,T)>0 &&
        0<=phi<=1 && 0<=phibar<=1 && min(np,nx)>=4 || throw(ArgumentError("invalid centered inputs"))
    poles,bweights,weights=Float64[],Float64[],Float64[]
    contact=0.0
    for (vacuum,L) in ((true,vc),(false,tc))
        ps,pws=R.gauleg(0.0,L,np)
        xs,xws=q==0 ? ([0.0],[2.0]) : R.gauleg(-1.0,1.0,nx)
        for (p,pw) in zip(ps,pws),(x,xw) in zip(xs,xws)
            e1=sqrt(p^2+q^2/4+p*q*x+a^2)
            e2=sqrt(p^2+q^2/4-p*q*x+b^2)
            n1=vacuum ? (1.0,0.0) : (-antiquark_distribution(e1,mu1,T,phi,phibar),quark_distribution(e1,mu1,T,phi,phibar))
            n2=vacuum ? (1.0,0.0) : (-antiquark_distribution(e2,mu2,T,phi,phibar),quark_distribution(e2,mu2,T,phi,phibar))
            measure=pw*xw*p^2/(e1*e2)
            contact+=measure*2*(e1*(n2[2]-n2[1])+e2*(n1[2]-n1[1]))
            for (si,s) in enumerate((-1,1)),(ti,t) in enumerate((-1,1))
                u=s*e1-t*e2
                bw=measure*s*t*(n2[ti]-n1[si])
                bw==0 && continue
                push!(poles,u);push!(bweights,bw)
                push!(weights,NC/(8pi^2)*bw*(u^2-q^2-(a-b)^2))
            end
        end
    end
    return (poles=poles,bweights=bweights,weights=weights,contact=contact,q=q,mass_term=(a-b)^2)
end

function atom_value(g,z)
    pi_value=sum(g.weights[i]/(z-g.poles[i]) for i in eachindex(g.poles);init=0.0im)
    b0=sum(g.bweights[i]/(z-g.poles[i]) for i in eachindex(g.poles);init=0.0im)
    reconstructed=NC/(8pi^2)*((z^2-g.q^2-g.mass_term)*b0-g.contact)
    return (value=pi_value,b0=b0,reconstructed=reconstructed,moment0=sum(g.bweights),
        contact_residual=abs(pi_value-reconstructed))
end

"""Kinematic support enclosures from extrema of the SAME compact momentum domain.

These are geometry, not a threshold on tiny occupation weights. All sign branches
are included. At T>0 the thermal domain supplies both pair and scattering cuts.
"""
function domain_support(q,a,b,L,regulator)
    all(isfinite,(q,a,b,L)) && q>=0 && min(a,b,L)>0 || throw(ArgumentError("invalid support geometry"))
    regulator in (:two_line,:centered) || throw(ArgumentError("unknown regulator"))
    regulator===:two_line && q>=2L && return nothing
    lo,hi=regulator===:two_line ? (max(0.0,q-L),L) : (q/2-L,q/2+L)
    pair(p)=hypot(p,a)+hypot(p-q,b)
    pmin=clamp(a*q/(a+b),lo,hi)
    umin=pair(pmin)
    umax=if regulator===:two_line
        hypot(L,a)+hypot(L,b)
    else
        x=L*q==0 ? 0.0 : clamp((b^2-a^2)/(2L*q),-1.0,1.0)
        sqrt(L^2+q^2/4+L*q*x+a^2)+sqrt(L^2+q^2/4-L*q*x+b^2)
    end
    function maxdiff(a,b)
        points=[lo,hi]
        if a!=b
            p=a*q/(a-b)
            lo<=p<=hi && push!(points,p)
        elseif q==0
            push!(points,clamp(0.0,lo,hi))
        end
        return maximum(hypot(p,a)-hypot(p-q,b) for p in points)
    end
    return (pair=(umin,umax),landau=(-maxdiff(b,a),maxdiff(a,b)))
end

function merge_intervals(intervals)
    out=Tuple{Float64,Float64}[]
    for (a,b) in sort(collect(intervals);by=first)
        a<b || continue
        if isempty(out) || a>last(out)[2]
            push!(out,(a,b))
        else
            out[end]=(last(out)[1],max(b,last(out)[2]))
        end
    end
    return out
end

function support_intervals(q,a,b,vc,tc,regulator)
    cuts=Tuple{Float64,Float64}[]
    for (thermal,L) in ((false,vc),(true,tc))
        s=domain_support(q,a,b,L,regulator)
        s===nothing && continue
        l,u=s.pair
        append!(cuts,[(l,u),(-u,-l)])
        if thermal
            l,u=s.landau
            append!(cuts,[(l,u),(-u,-l)])
        end
    end
    return merge_intervals(cuts)
end

function analytic_gaps(cuts,shift;lower=0.0,upper=64.0)
    gaps=Tuple{Float64,Float64}[]
    edge=lower
    for (a,b) in cuts
        a,b=a-shift,b-shift
        b<=lower && continue
        a>=upper && break
        edge<min(a,upper) && push!(gaps,(edge,min(a,upper)))
        edge=max(edge,b)
    end
    edge<upper && push!(gaps,(edge,upper))
    return gaps
end

function centered_profile(bg,channel,q;mesh=256,ne=128)
    a,b=R.charged_rpa_spec(channel).pair
    m1,m2,mu1,mu2=bg.m[a],bg.m[b],bg.mu[a],bg.mu[b]
    vc,tc=bg.vacuum,24.0
    cuts=support_intervals(q,m1,m2,vc,tc,:centered)
    upper=maximum(last,cuts)+1.0
    edges=[-upper,upper,0.0,mu1-mu2]
    for (_,L) in ((false,vc),(true,tc))
        sp=domain_support(q,m1,m2,L,:centered)
        for e in (sp.pair...,sp.landau...,hypot(q,m1+m2),hypot(q,m1-m2))
            append!(edges,[e,-e])
        end
        # Surface tangencies change the cut formula; include all cap boundaries.
        for e in (hypot(L+q/2,m1)+hypot(L-q/2,m2),
                  hypot(L-q/2,m1)+hypot(L+q/2,m2),
                  hypot(L+q/2,m1)-hypot(L-q/2,m2),
                  hypot(L-q/2,m1)-hypot(L+q/2,m2))
            append!(edges,[e,-e])
        end
    end
    sort!(unique!(edges))
    x=Float64[]
    for j in 1:length(edges)-1
        left,right=edges[j],edges[j+1]
        append!(x,[left+(right-left)*(1-cospi(i/mesh))/2 for i in 0:mesh-1])
    end
    push!(x,last(edges));sort!(unique!(x))
    function rho(l)
        b0=sum(centered_cut(l,q,m1,mu1,m2,mu2,bg.T,bg.Phi,bg.PhiBar,L;
            ne=ne,component=c) for (c,L) in ((:vacuum,vc),(:thermal,tc)))
        return NC/(8pi^2)*(l^2-q^2-(m1-m2)^2)*b0
    end
    y=rho.(x)
    y[1]=y[end]=0.0
    p=R.PiecewiseSpectralFunction(x,y)
    shift=mu1-mu2
    return (profile=p,inverse=z->1-4bg.coupling[channel]*R.cauchy_transform(p,z+shift),
        shift=shift,threshold=hypot(q,m1+m2)-shift,landau=hypot(q,m1-m2)-shift,
        grid=(q_inv_fm=q,),support=cuts)
end

"""A bounded full support-gap audit; endpoint margins are explicitly not counted."""
function count_support_gaps(b,s,cuts;upper=64.0,margin=1e-6)
    rows=NamedTuple[]
    for (id,(left,right)) in enumerate(analytic_gaps(cuts,b.shift;upper=upper))
        right-left>2margin || continue
        gap=(left+margin,right-margin)
        # Apply margin ONCE, consistently with the counting rectangle below.
        lo=R.certify_gap_roots((w,_)->b.inverse(w),b.grid.q_inv_fm,[(left,right)];
            physical_sheet=true,real_axis=true,omega_nodes=s.nr,endpoint_margin=margin)
        hi=R.certify_gap_roots((w,_)->b.inverse(w),b.grid.q_inv_fm,[(left,right)];
            physical_sheet=true,real_axis=true,omega_nodes=2s.nr,endpoint_margin=margin)
        f(z)=imag(z)>=0 ? b.inverse(z) : conj(b.inverse(conj(z)))
        height=min(0.02,(gap[2]-gap[1])/20)
        # Begin cheaply in smooth/no-root gaps; retain the SAME adaptive doubling,
        # angle-step threshold and maximum resolution used by the original counter.
        c=R.contour_count(f,gap...,-height,height;nodes=16)
        rootstable=lo.passed && hi.passed && lo.count==hi.count &&
            all(abs(lo.roots[i].omega_inv_fm-hi.roots[i].omega_inv_fm)<1e-7 for i in eachindex(lo.roots))
        push!(rows,(gap_id=id,left_inv_fm=left,right_inv_fm=right,margin_inv_fm=margin,
            root_count=hi.count,contour_count=c.count,passed=rootstable && c.passed && c.count==hi.count,
            status=String(hi.status),coarse_status=String(lo.status),root_nodes_stable=rootstable,
            contour_passed=c.passed,contour_nodes_per_edge=c.nodes_per_edge,
            max_root_shift_inv_fm=lo.count==hi.count ? maximum((abs(lo.roots[i].omega_inv_fm-hi.roots[i].omega_inv_fm) for i in eachindex(lo.roots));init=0.0) : Inf,
            max_contour_step=c.max_step,
            roots_inv_fm=join((r.omega_inv_fm for r in hi.roots),';'),production_authorized=false))
    end
    return rows
end
end
